/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver, Deniz Gunceler

This file is part of JDFTx.

JDFTx is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

JDFTx is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with JDFTx.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#include <electronic/Everything.h>
#include <fluid/NonlinearPCM.h>
#include <fluid/LinearPCM.h>
#include <fluid/PCM_internal.h>
#include <core/ScalarFieldIO.h>
#include <core/SphericalHarmonics.h>
#include <core/Util.h>


//Initialize Kkernel to square-root of the inverse kinetic operator
inline void setPreconditioner(int i, double Gsq, double* preconditioner, double epsBulk, double kappaSq, RadialFunctionG& w1)
{	double epsEff = w1
		? (1 + (epsBulk-1) * std::pow(w1(sqrt(Gsq)), 2))
		: epsBulk;
	preconditioner[i] = (Gsq || kappaSq) ? 1./(epsEff*Gsq + kappaSq) : 0.;
}


NonlinearPCM::NonlinearPCM(const Everything& e, const FluidSolverParams& fsp)
: PCM(e, fsp), pMol(0.), ionNbulk(0.), ionZ(0.), screeningEval(0), dielectricEval(0), isNonlocal(false)
{
	const auto& solvent = fsp.solvents[0];
	pMol = solvent->pMol ? solvent->pMol : solvent->molecule.getDipole().length();
	
	//Initialize dielectric evaluation class:
	dielectricEval = new NonlinearPCMeval::Dielectric(fsp.linearDielectric,
		fsp.T, solvent->Nbulk, pMol, epsBulk, solvent->epsInf);

	//Check and setup ionic screening:
	if(fsp.cations.size() > 1) die("NonlinearPCM currently only supports a single cationic component.\n");
	if(fsp.anions.size() > 1) die("NonlinearPCM currently only supports a single anionic component.\n");
	assert(fsp.anions.size() == fsp.cations.size()); //this should be ensured by charge neutrality check in FluidSolver constructor
	if(fsp.cations.size())
	{	//Ensure charge balanced:
		if(fabs(fsp.cations[0]->molecule.getCharge() + fsp.anions[0]->molecule.getCharge())>1e-12)
			die("NonlinearPCM currently only supports charge-balanced (Z:Z) electrolytes.\n");
		ionNbulk = fsp.cations[0]->Nbulk;
		ionZ = fsp.anions[0]->molecule.getCharge();
		double VhsCation = fsp.cations[0]->molecule.getVhs(); if(!VhsCation) VhsCation = (4*M_PI/3)*pow(fsp.cations[0]->Rvdw,3);
		double VhsAnion = fsp.anions[0]->molecule.getVhs(); if(!VhsAnion) VhsAnion = (4*M_PI/3)*pow(fsp.anions[0]->Rvdw,3);
		screeningEval = new NonlinearPCMeval::Screening(fsp.linearScreening, fsp.T, ionNbulk, ionZ, VhsCation, VhsAnion, epsBulk);
	}
	else
	{	ionNbulk = 0.;
		screeningEval = 0;
	}
	
	//Dielectric lookup table
	double dxMapped = 1./512;
	std::vector<double> samples;
	for(double xMapped=0.; xMapped<=1.; xMapped+=dxMapped)
	{	if(xMapped==0.) //x -> 0
			samples.push_back(0.5 * dielectricEval->NT * (1.0/(3.0 - dielectricEval->alpha) + dielectricEval->X));
		else if(xMapped==1.) //x -> infty
			samples.push_back(0.5 * dielectricEval->NT * dielectricEval->X);
		else
		{	double x = xMapped/(1.-xMapped); //inverse of xMapped = x / (1 + x)
			double eps = dielectricEval->eps_from_x(x), frac, logsinch;
			dielectricEval->calcFunctions(eps, frac, logsinch);
			double energy = dielectricEval->NT * (
				logsinch - 0.5 * dielectricEval->alpha * std::pow(eps * frac, 2)
				+ 0.5 * dielectricEval->X * (x * x)
			);
			samples.push_back(energy/(x*x));
		}
	}
	dielEnergyLookup.init(0, samples, dxMapped);
	
	//Screening lookup table
	if(screeningEval)
	{
		double dVmapped = 1./512;
		std::vector<double> samples;
		for(double Vmapped=-1.; Vmapped<=1.; Vmapped+=dVmapped)
		{	if(fabs(Vmapped)==1)
				samples.push_back(0.);
			else
			{	double V = Vmapped / (1 - Vmapped*Vmapped); //inverse of Vmapped = 2V/(1 + sqrt(1 + 4V^2))
				double x = screeningEval->x_from_V(V);
				
				//Corresponding energy for the nonlinear screening:
				double f_x, f = screeningEval->fHS(x, f_x);
				double energy =
					( exp(+V - f_x * screeningEval->x0minus)
					+ exp(-V - f_x * screeningEval->x0plus)
					+ f_x * x - f - 2.0); //upto factor of NT; zero energy at V=0 and energy -> infty at infinite |V|
				samples.push_back(1./(1. + energy)); //map from [0, infty) -> [1, 0), with zero at Vmapped = +/-1
			}
		}
		ionEnergyLookup.init(1, samples, dVmapped);
	}
	
	//Initialize convolution kernels for nonlocal version, CANON
	isNonlocal = (fsp.pcmVariant == PCM_CANON);
	if(isNonlocal)
	{	const double dG = gInfo.dGradial, Gmax = gInfo.GmaxGrid;
		unsigned nGradial = unsigned(ceil(Gmax/dG))+5;
		std::vector<double> Nw0_samples(nGradial), w1_samples(nGradial);
		for(unsigned i=0; i<nGradial; i++)
		{	double GR = (i*dG) * fsp.Res;
			double j0 = bessel_jl(0, GR);
			double j1_by_x_3 = j0 + bessel_jl(2, GR);  //3 j1(x)/x = j0(x) + j2(x)
			Nw0_samples[i] = solvent->Nbulk * (-fsp.Zcenter) * (1.0 - j0); //Zcenter is in expt charge convention
			w1_samples[i] = j1_by_x_3;
		}
		Nw0.init(0, Nw0_samples, dG);
		w1.init(0, w1_samples, dG);
	}
	
	//Initialize preconditioner:
	preconditioner = std::make_shared<RealKernel>(gInfo);
	applyFuncGsq(gInfo, setPreconditioner, preconditioner->data(), epsBulk, k2factor, w1);
}

NonlinearPCM::~NonlinearPCM()
{	delete dielectricEval;
	dielEnergyLookup.free();
	if(screeningEval)
	{	delete screeningEval;
		ionEnergyLookup.free();
	}
	if(isNonlocal)
	{	Nw0.free();
		w1.free();
	}
}

void NonlinearPCM::loadState(const char* filename)
{	ScalarField Iphi(ScalarFieldData::alloc(gInfo));
	loadRawBinary(Iphi, filename); //saved data is in real space
	phiTot = J(Iphi);
}

void NonlinearPCM::saveState(const char* filename) const
{	if(mpiWorld->isHead()) saveRawBinary(I(phiTot), filename); //saved data is in real space
}

void NonlinearPCM::minimizeFluid()
{	//Info:
	logPrintf("\tNonlinear fluid (bulk dielectric constant: %g) occupying %lf of unit cell\n",
		epsBulk, integral(shape[0])/gInfo.detR);
	if(k2factor)
		logPrintf("\tNonlinear screening (bulk screening length: %g bohrs) occupying %lf of unit cell\n",
		sqrt(epsBulk/k2factor), integral(shape.back())/gInfo.detR);
	logFlush();

	//Minimize:
	minimize(e.fluidMinParams);
	logPrintf("\tNonlinear solve completed after %d iterations at t[s]: %9.2lf\n", iterLast, clock_sec());
}

void NonlinearPCM::step(const ScalarFieldTilde& dir, double alpha)
{	::axpy(alpha, dir, phiTot);
}

double NonlinearPCM::compute(ScalarFieldTilde* grad, ScalarFieldTilde* Kgrad)
{
	//Vacuum term:
	ScalarFieldTilde KinvPhi = L(phiTot) * (-1.0/(4*M_PI*gInfo.detR)); // -nabla^2 phi/(4 pi)
	double Avac = 0.5 * dot(KinvPhi, O(phiTot));
	
	//Dielectric term:
	ScalarField A, A_s_UNUSED; nullToZero(A, gInfo);
	VectorField Dphi = I(gradient(w1 ? w1*phiTot : phiTot)), A_Dphi_null;  //includes CANON nonlocality, if needed
	VectorField& A_Dphi = grad ? Dphi : A_Dphi_null; //Retrieve gradient in place (since Dphi no longer needed)
	(*dielectricEval)(dielEnergyLookup, shape[0], Dphi, A, A_Dphi, A_s_UNUSED);
	
	//Screening term:
	ScalarField A_phi;
	if(screeningEval)
	{	ScalarField phi = I(phiTot);
		if(grad) nullToZero(A_phi, gInfo);
		(*screeningEval)(ionEnergyLookup, shape.back(), phi, A, A_phi, A_s_UNUSED);
	}
	
	//Free charges (right-hand side of Poisson-Boltzmann charge):
	ScalarFieldTilde rhoFreeTilde = rhoExplicitTilde + rhoLiquidTilde0;
	double Acoulomb = dot(rhoFreeTilde, O(phiTot));
	
	if(grad)
	{	ScalarFieldTilde A_phiTilde = -divergence(J(A_Dphi)); //dielectric part
		if(w1) A_phiTilde = w1 * A_phiTilde; //CANON nonlocality
		if(A_phi) A_phiTilde += J(A_phi); //ionic part
		A_phiTilde += KinvPhi; //vacuum part
		*grad = O(A_phiTilde - rhoFreeTilde);
		if(Kgrad)
		{	*Kgrad = (*preconditioner) * (*grad);
		}
	}
	return A0 + Avac + integral(A) - Acoulomb;
}

bool NonlinearPCM::report(int iter)
{	iterLast = iter;
	return false;
}


void NonlinearPCM::set_internal(const ScalarFieldTilde& rhoExplicitTilde, const ScalarFieldTilde& nCavityTilde)
{	//Store the explicit system charge:
	this->rhoExplicitTilde = rhoExplicitTilde; zeroNyquist(this->rhoExplicitTilde);
	A0 = 0.5 * dot(rhoExplicitTilde, O(coulomb(rhoExplicitTilde)));
	
	//Update cavity:
	nCavity = I(nCavityTilde + getFullCore());
	updateCavity();
	
	//Built-in charge for CANON:
	if(isNonlocal)
		rhoLiquidTilde0 = Nw0 * J(shape[0]);
	
	//Initialize the state if it hasn't been loaded:
	if(!phiTot) nullToZero(phiTot, gInfo);
}

double NonlinearPCM::get_Adiel_and_grad_internal(ScalarFieldTilde& Adiel_rhoExplicitTilde, ScalarFieldTilde& Adiel_nCavityTilde, IonicGradient* extraForces, matrix3<>* Adiel_RRT) const
{
	EnergyComponents& Adiel = ((NonlinearPCM*)this)->Adiel;
	ScalarFieldArray F_shape; nullToZero(F_shape, gInfo, shape.size());
	ScalarField F; nullToZero(F, gInfo); //Hessian energy contributions

	//Compute vacuum contributions:
	double minusFvac = dot(phiTot, L(phiTot)) * (1.0/(8*M_PI));
	if(Adiel_RRT) *Adiel_RRT += Lstress(phiTot, phiTot) * (1.0/(8*M_PI));
	
	//Compute dielectric contributions:
	{	VectorField Dphi = I(gradient(w1 ? w1*phiTot : phiTot)), F_Dphi;
		if(Adiel_RRT) nullToZero(F_Dphi, gInfo); //only needed for stress
		(*dielectricEval)(dielEnergyLookup, shape[0], Dphi, F, F_Dphi, F_shape[0]);
		if(Adiel_RRT) *Adiel_RRT += gInfo.dV * dotOuter(F_Dphi, Dphi);
	}
	
	//Compute screening contributions
	if(screeningEval)
	{	ScalarField phi = I(phiTot), F_phi_UNUSED;
		(*screeningEval)(ionEnergyLookup, shape.back(), phi, F, F_phi_UNUSED, F_shape.back());
	}
	
	//Compute the energy:
	ScalarFieldTilde phiExplicitTilde = coulomb(rhoExplicitTilde), phiLiquid0;
	ScalarFieldTilde rhoFreeTilde = rhoExplicitTilde + rhoLiquidTilde0;
	Adiel["Electrostatic"] = minusFvac - integral(F)
		+ dot(phiTot, O(rhoFreeTilde))
		- 0.5*dot(phiExplicitTilde, O(rhoExplicitTilde));
	if(isNonlocal)
	{	phiLiquid0 = coulomb(rhoLiquidTilde0);
		Adiel["Electrostatic"] -= 0.5 * dot(phiLiquid0, O(rhoLiquidTilde0));
	}
	if(Adiel_RRT)
	{	*Adiel_RRT += Adiel["Electrostatic"] * matrix3<>(1,1,1) //volume contribution
			- 0.5*coulombStress(rhoExplicitTilde, rhoExplicitTilde); //through coulomb in phiExt
	}
	
	//Collect cavity shape derivatives:
	ScalarFieldArray Adiel_shape = -1.0 * F_shape; //since energy contribution is -F
	if(isNonlocal)
	{	ScalarFieldTilde Adiel_rhoLiquidTilde0 = phiTot - phiLiquid0;
		Adiel_shape[0] += I(Nw0 * Adiel_rhoLiquidTilde0); 
	}
	
	//Propagate to derivatives w.r.t electronic charge and density:
	Adiel_rhoExplicitTilde = phiTot - phiExplicitTilde;
	ScalarField Adiel_nCavity;
	propagateCavityGradients(Adiel_shape, Adiel_nCavity, Adiel_rhoExplicitTilde, extraForces, Adiel_RRT);
	Adiel_nCavityTilde = J(Adiel_nCavity);
	
	accumExtraForces(extraForces, Adiel_nCavityTilde);
	return Adiel;
}

void NonlinearPCM::dumpDensities(const char* filenamePattern) const
{	PCM::dumpDensities(filenamePattern);

	//Output dielectric bound charge:
	string filename;
	{	ScalarField A, A_s_UNUSED; nullToZero(A, gInfo);
		VectorField Dphi = I(gradient(phiTot));
		VectorField& A_Dphi = Dphi; //Retrieve dA/dDphi = -P in place
		(*dielectricEval)(dielEnergyLookup, shape[0], Dphi, A, A_Dphi, A_s_UNUSED);
		ScalarField rhoDiel = divergence(A_Dphi); //bound charge due to dielectric (since A_Dphi = -P)
		FLUID_DUMP(rhoDiel, "RhoDiel");
	}
	
	//Output ionic bound charge (if any):
	if(screeningEval)
	{	ScalarField A, A_s_UNUSED, A_phi;
		ScalarField phi = I(phiTot);
		nullToZero(A, gInfo);
		nullToZero(A_phi, gInfo); //Retrieve dA/dphi = -rho_ion
		(*screeningEval)(ionEnergyLookup, shape.back(), phi, A, A_phi, A_s_UNUSED);
		FLUID_DUMP(-A_phi, "RhoIon");
	}
}
