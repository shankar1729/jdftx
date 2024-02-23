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
#include <core/Util.h>


//Initialize Kkernel to square-root of the inverse kinetic operator

inline void setPreconditioner(int i, double Gsq, double* preconditioner, double epsBulk, double kappaSq)
{	preconditioner[i] = (Gsq || kappaSq) ? 1./(epsBulk*Gsq + kappaSq) : 0.;
}


NonlinearPCM::NonlinearPCM(const Everything& e, const FluidSolverParams& fsp)
: PCM(e, fsp), NonlinearCommon(fsp, epsBulk)
{
	//Initialize preconditioner:
	preconditioner = std::make_shared<RealKernel>(gInfo);
	applyFuncGsq(gInfo, setPreconditioner, preconditioner->data(), epsBulk, k2factor);
}

NonlinearPCM::~NonlinearPCM()
{
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
	ScalarField A, A_s_UNUSED; nullToZero(A, gInfo);
	VectorField Dphi = I(gradient(phiTot)), A_Dphi_null;
	VectorField& A_Dphi = grad ? Dphi : A_Dphi_null; //Retrieve gradient in place (since Dphi no longer needed)
	(*dielectricEval)(dielEnergyLookup, shape[0], Dphi, A, A_Dphi, A_s_UNUSED);
	
	ScalarField A_phi;
	if(screeningEval)
	{	ScalarField phi = I(phiTot);
		if(grad) nullToZero(A_phi, gInfo);
		(*screeningEval)(ionEnergyLookup, shape.back(), phi, A, A_phi, A_s_UNUSED);
	}
	
	if(grad)
	{	ScalarFieldTilde A_phiTilde = -divergence(J(A_Dphi));
		if(A_phi) A_phiTilde += J(A_phi);
		*grad = O(A_phiTilde - rhoExplicitTilde);
		if(Kgrad)
		{	*Kgrad = (*preconditioner) * (*grad);
		}
	}
	return A0 + integral(A) - dot(rhoExplicitTilde, O(phiTot));
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
	this->nCavity = I(nCavityTilde + getFullCore());
	updateCavity();
	
	//Initialize the state if it hasn't been loaded:
	if(!phiTot) nullToZero(phiTot, gInfo);
}

double NonlinearPCM::get_Adiel_and_grad_internal(ScalarFieldTilde& Adiel_rhoExplicitTilde, ScalarFieldTilde& Adiel_nCavityTilde, IonicGradient* extraForces, matrix3<>* Adiel_RRT) const
{
	EnergyComponents& Adiel = ((NonlinearPCM*)this)->Adiel;
	ScalarFieldArray Adiel_shape; nullToZero(Adiel_shape, gInfo, shape.size());
	ScalarField F; nullToZero(F, gInfo); //Hessian energy contributions

	//Compute dielectric contributions:
	{	VectorField Dphi = I(gradient(phiTot)), F_Dphi;
		if(Adiel_RRT) nullToZero(F_Dphi, gInfo); //only needed for stress
		(*dielectricEval)(dielEnergyLookup, shape[0], Dphi, F, F_Dphi, Adiel_shape[0]);
		if(Adiel_RRT) *Adiel_RRT += gInfo.dV * dotOuter(F_Dphi, Dphi);
	}

	//Compute screening contributions
	if(screeningEval)
	{	ScalarField phi = I(phiTot), F_phi_UNUSED;
		(*screeningEval)(ionEnergyLookup, shape.back(), phi, F, F_phi_UNUSED, Adiel_shape.back());
	}

	//Compute the energy:
	ScalarFieldTilde phiExplicitTilde = coulomb(rhoExplicitTilde);
	Adiel["Electrostatic"] = -integral(F) + dot(phiTot - 0.5*phiExplicitTilde, O(rhoExplicitTilde));
	if(Adiel_RRT)
	{	*Adiel_RRT += Adiel["Electrostatic"] * matrix3<>(1,1,1) //volume contribution
			- 0.5*coulombStress(rhoExplicitTilde, rhoExplicitTilde); //through coulomb in phiExt
	}
	
	//Derivatives w.r.t electronic charge and density:
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
