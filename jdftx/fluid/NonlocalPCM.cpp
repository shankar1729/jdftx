/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

#include <core/DataIO.h>
#include <core/DataMultiplet.h>
#include <fluid/NonlocalPCM.h>
#include <fluid/PCM_internal.h>
#include <electronic/Everything.h>
#include <electronic/SphericalHarmonics.h>
#include <electronic/VanDerWaals.h>
#include <electronic/operators.h>
#include <cstring>

//Electrostatic kernel (gaussian convolved with l=1 Bessel)
inline double wDielTilde(double G, double sigma, double eta)
{	double sigmaG = sigma*G;
	double sigmaGsq = sigmaG*sigmaG;
	double etaG = eta*G;
	return exp(-0.5*sigmaGsq) * (etaG ? 3.*bessel_jl(1,etaG)/etaG : 1.);
}
//derivative if wDielTilde w.r.t eta:
inline double wDielTilde_eta(double G, double sigma, double eta)
{	double sigmaG = sigma*G;
	double sigmaGsq = sigmaG*sigmaG;
	double etaG = eta*G;
	return exp(-0.5*sigmaGsq) * (etaG ? 3.*(sin(etaG)-3*bessel_jl(1,etaG))/(etaG*eta) : 0.);
}

//Initialize Kkernel to inverse square-root of the Hessian in the uniform fluid limit
inline double setPreconditionerKernel(double G, double sigma, double eta, double epsBulk, double k2factor)
{	const double wDiel = wDielTilde(G, sigma, eta);
	double diagH = G*G + wDiel*wDiel*(G*G*(epsBulk-1.) + k2factor);
	return diagH>1e-12 ? 1./sqrt(diagH) : 0.;
}

NonlocalPCM::NonlocalPCM(const Everything& e, const FluidSolverParams& fsp) : PCM(e, fsp), wCavity(Sf[0])
{	//Compute the gaussian width parameter from Rvdw:
	sigmaVdw = 1.;
	for(int iter=0; iter<50; iter++) //-- solve (Ztot wCavity * Ztot wCavity)(2 Rvdw) = nc by fixed-point (Picard) iteration
	{	double sigmaVdwNew = (2*fsp.solvents[0]->Rvdw) / sqrt(-4. * log(fsp.nc * pow(2*sqrt(M_PI)*sigmaVdw, 3) / pow(fsp.Ztot,2)));
		if(fabs(sigmaVdwNew/sigmaVdw - 1.) < 1e-12) break;
		sigmaVdw = sigmaVdwNew;
	}
	//Initialize kernels:
	wCavity.init(0, e.gInfo.dGradial, e.gInfo.GmaxGrid, RadialFunctionG::gaussTilde, 1., sigmaVdw);
	wDiel.init(0, e.gInfo.dGradial, e.gInfo.GmaxGrid, wDielTilde, sigmaVdw, fsp.eta_wDiel);
	Kkernel.init(0, e.gInfo.dGradial, e.gInfo.GmaxGrid, setPreconditionerKernel, sigmaVdw, fsp.eta_wDiel, epsBulk, k2factor);
	logPrintf("   NonlocalPCM weight functions with sigma = %lg bohr and eta = %lg\n", sigmaVdw, fsp.eta_wDiel);
}

NonlocalPCM::~NonlocalPCM()
{	wDiel.free();
	Kkernel.free();
}


DataGptr NonlocalPCM::hessian(const DataGptr& phiTilde) const
{	DataGptr wphiTilde = wDiel * phiTilde;
	DataGptr wrhoTilde = ((epsBulk-1.)/(4*M_PI)) * divergence(J(shape * I(gradient(wphiTilde))));  //dielectric bound charge (except for convolution by wDiel)
	if(k2factor)
		wrhoTilde -= (k2factor/(4*M_PI)) * J(shape*I(wphiTilde)); //ionic bound charge (except for convolution by wDiel)
	return (-1./(4*M_PI*e.gInfo.detR)) * L(phiTilde) - wDiel * wrhoTilde; //return free charge = total charge - bound charge
}

DataGptr NonlocalPCM::precondition(const DataGptr& rTilde) const
{	return Kkernel*(J(epsInv*I(Kkernel*rTilde)));
}


void NonlocalPCM::set(const DataGptr& rhoExplicitTilde, const DataGptr& nCavityTilde)
{
	this->rhoExplicitTilde = clone(rhoExplicitTilde); zeroNyquist(this->rhoExplicitTilde);
	
	//Compute cavity shape function (0 to 1)
	nCavity = fsp.Ztot * I(wCavity * nCavityTilde);
	updateCavity();

	logPrintf("\tNonlocalPCM fluid occupying %lf of unit cell:", integral(shape)/e.gInfo.detR);
	logFlush();

	//Update the inhomogeneity factor of the preconditioner
	epsInv = inv(1. + (epsBulk-1.)*shape);
	
	//Initialize the state if it hasn't been loaded:
	if(!state) nullToZero(state, e.gInfo);
}


void NonlocalPCM::minimizeFluid()
{
	fprintf(e.fluidMinParams.fpLog, "\n\tWill stop at %d iterations, or sqrt(|r.z|)<%le\n", e.fluidMinParams.nIterations, e.fluidMinParams.knormThreshold);
	int nIter = solve(rhoExplicitTilde, e.fluidMinParams);
	logPrintf("\tCompleted after %d iterations.\n", nIter);
}

double NonlocalPCM::get_Adiel_and_grad(DataGptr& Adiel_rhoExplicitTilde, DataGptr& Adiel_nCavityTilde, IonicGradient& extraForces) const
{
	EnergyComponents& Adiel = ((NonlocalPCM*)this)->Adiel;
	const DataGptr& phi = state; // that's what we solved for in minimize

	//First-order correct estimate of electrostatic energy:
	DataGptr phiExt = coulomb(rhoExplicitTilde);
	Adiel["Electrostatic"] = -0.5*dot(phi, O(hessian(phi))) + dot(phi - 0.5*phiExt, O(rhoExplicitTilde));
	
	//Gradient w.r.t rhoExplicitTilde:
	Adiel_rhoExplicitTilde = phi - phiExt;

	//Compute gradient w.r.t shape function:
	DataGptr wphi = wDiel * phi;
	//--- Dielectric contributions:
	DataRptr Adiel_shape = (-(epsBulk-1)/(8*M_PI)) * lengthSquared(I(gradient(wphi)));
	//--- Screening contributions:
	if(k2factor)
	{	DataRptr Iwphi = I(wphi);
		Adiel_shape -= (k2factor/(8*M_PI)) * (Iwphi*Iwphi);
	}
	
	//Propagate shape gradients to A_nCavity:
	DataRptr Adiel_nCavity;
	propagateCavityGradients(Adiel_shape, Adiel_nCavity);
	Adiel_nCavityTilde = wCavity * J(fsp.Ztot*Adiel_nCavity);
	
	if(vdwForces) extraForces = *vdwForces;
	return Adiel;
}

void NonlocalPCM::loadState(const char* filename)
{	DataRptr Istate(DataR::alloc(e.gInfo));
	loadRawBinary(Istate, filename); //saved data is in real space
	state = J(Istate);
}

void NonlocalPCM::saveState(const char* filename) const
{	if(mpiUtil->isHead()) saveRawBinary(I(state), filename); //saved data is in real space
}

void NonlocalPCM::printDebug(FILE* fp) const
{	//Print derivative w.r.t electrostatic weight function parameter:
	//--- compute derivative w.r.t w*phi
	const DataGptr& phi = state;
	DataGptr wphi = wDiel*phi;
	DataGptr A_wphi = ((epsBulk-1)/(4*M_PI)) * divergence(J(shape * I(gradient(wphi))));
	if(k2factor) A_wphi -= (k2factor/(4*M_PI)) * J(shape * I(wphi));
	//--- propagate to derivative w.r.t eta
	RadialFunctionG wDiel_eta;
	wDiel_eta.init(0, e.gInfo.dGradial, e.gInfo.GmaxGrid, wDielTilde_eta, sigmaVdw, fsp.eta_wDiel);
	double A_eta = integral(I(A_wphi) * I(wDiel_eta*phi));
	wDiel_eta.free();
	fprintf(fp, "   E_wDiel_eta = %.15lg\n", A_eta);
}
