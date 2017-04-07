/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver

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
#include <fluid/LinearPCM.h>
#include <fluid/PCM_internal.h>
#include <electronic/operators.h>
#include <core/VectorField.h>
#include <core/ScalarFieldIO.h>
#include <core/Thread.h>

LinearPCM::LinearPCM(const Everything& e, const FluidSolverParams& fsp)
: PCM(e, fsp)
{
}

LinearPCM::~LinearPCM()
{	Kkernel.free();
}

ScalarFieldTilde LinearPCM::hessian(const ScalarFieldTilde& phiTilde) const
{	//Dielectric term:
	ScalarField epsilon = epsilonOverride ? epsilonOverride : 1. + (epsBulk-1.) * shape;
	ScalarFieldTilde rhoTilde = divergence(J(epsilon * I(gradient(phiTilde))));
	//Screening term:
	if(k2factor)
	{	ScalarField kappaSq = kappaSqOverride ? kappaSqOverride : k2factor * shape;
		rhoTilde -= J(kappaSq * I(phiTilde));
	}
	return (-1./(4*M_PI)) * rhoTilde;
}

ScalarFieldTilde LinearPCM::chi(const ScalarFieldTilde& phiTilde) const
{	return (-1./(4*M_PI*gInfo.detR)) * L(phiTilde) - hessian(phiTilde);
}

ScalarFieldTilde LinearPCM::precondition(const ScalarFieldTilde& rTilde) const
{	return Kkernel*(J(epsInv*I(Kkernel*rTilde)));
}

//Initialize Kkernel to square-root of the inverse kinetic operator
inline double setPreconditionerKernel(double G, double epsMean, double kRMS)
{	return (G || kRMS) ? 1./(epsMean*hypot(G, kRMS)) : 0.;
}

void LinearPCM::set_internal(const ScalarFieldTilde& rhoExplicitTilde, const ScalarFieldTilde& nCavityTilde)
{	//Store the explicit system charge:
	this->rhoExplicitTilde = rhoExplicitTilde; zeroNyquist(this->rhoExplicitTilde);
	
	//Update cavity:
	this->nCavity = I(nCavityTilde + getFullCore());
	updateCavity();

	//Update the preconditioner
	ScalarField epsilon = 1 + (epsBulk-1)*shape;
	ScalarField kappaSq = k2factor ? k2factor*shape : 0; //set kappaSq to null pointer if no screening
	updatePreconditioner(epsilon, kappaSq);
	
	//Initialize the state if it hasn't been loaded:
	if(!state) nullToZero(state, gInfo);
}

void LinearPCM::updatePreconditioner(const ScalarField& epsilon, const ScalarField& kappaSq)
{	epsInv = inv(epsilon);
	double epsMean = sum(epsilon) / gInfo.nr;
	double kappaSqMean = (kappaSq ? sum(kappaSq) : 0.) / gInfo.nr;
	Kkernel.init(0, 0.02, gInfo.GmaxGrid, setPreconditionerKernel, epsMean, sqrt(kappaSqMean/epsMean));
}

void LinearPCM::override(const ScalarField& epsilon, const ScalarField& kappaSq)
{	epsilonOverride = epsilon;
	kappaSqOverride = kappaSq;
	updatePreconditioner(epsilon, kappaSq);
}

void LinearPCM::minimizeFluid()
{	//Info:
	logPrintf("\tLinear fluid (dielectric constant: %g", epsBulk);
	if(k2factor) logPrintf(", screening length: %g Bohr", sqrt(epsBulk/k2factor));
	logPrintf(") occupying %lf of unit cell:", integral(shape)/gInfo.detR); logFlush();
	//Minimize:
	fprintf(e.fluidMinParams.fpLog, "\n\tWill stop at %d iterations, or sqrt(|r.z|)<%le\n", e.fluidMinParams.nIterations, e.fluidMinParams.knormThreshold);
	int nIter = solve(rhoExplicitTilde, e.fluidMinParams);
	logPrintf("\tCompleted after %d iterations at t[s]: %9.2lf\n", nIter, clock_sec());
}

double LinearPCM::get_Adiel_and_grad_internal(ScalarFieldTilde& Adiel_rhoExplicitTilde, ScalarFieldTilde& Adiel_nCavityTilde, IonicGradient* extraForces, bool electricOnly) const
{
	EnergyComponents& Adiel = ((LinearPCM*)this)->Adiel;
	ScalarFieldTilde phi = clone(state); // that's what we solved for in minimize
	
	//Neutrality constraint (and derivatives):
	double phi0_Qexp = 0.; ScalarField phi0_shape;
	if(k2factor)
	{	double chiPrefac = k2factor/(4*M_PI);
		ScalarField Iphi = I(phi);
		double Qexp = integral(rhoExplicitTilde);
		double sPhiInt = integral(shape * Iphi);
		double sInt = integral(shape);
		phi0_Qexp = 1./(sInt*chiPrefac);
		double phi0_sPhiInt = -1./sInt;
		double phi0 = phi0_Qexp * Qexp + phi0_sPhiInt * sPhiInt;
		double phi0_sInt = -phi0/sInt;
		phi0_shape = phi0_sInt + phi0_sPhiInt * Iphi;
		//Fix the G=0 of phi to include the constraint correction:
		phi->setGzero(phi->getGzero() + phi0);
	}
	
	//Calculate Coulomb contribution and gradients:
	ScalarFieldTilde Adiel_rhoFluid;
	{	ScalarFieldTilde rhoFluid = chi(phi);
		ScalarFieldTilde phiFluid = coulomb(rhoFluid);
		Adiel["Coulomb"] = dot(phiFluid, O(0.5*rhoFluid+rhoExplicitTilde));
		Adiel_rhoExplicitTilde = phiFluid;
		Adiel_rhoFluid = coulomb(rhoFluid+rhoExplicitTilde);
	}
	
	//Calculate internal energy and shape function gradients:
	ScalarField Adiel_shape;
	//--- dielectric
	{	double chiPrefac = (epsBulk-1)/(4*M_PI);
		VectorField IgradPhi = I(gradient(phi));
		ScalarField Aeps_shape = (0.5*chiPrefac) * lengthSquared(IgradPhi);
		Adiel["Aeps"] = integral(shape * Aeps_shape);
		Adiel_shape += Aeps_shape //internal energy contribution
			- chiPrefac * dotElemwise(IgradPhi, I(gradient(Adiel_rhoFluid))); //propagate Coulomb contribution
	}
	//--- screening
	if(k2factor)
	{	double chiPrefac = k2factor/(4*M_PI);
		ScalarField Iphi = I(phi), IAdiel_rhoFluid = I(Adiel_rhoFluid);
		ScalarField Akappa_shape = (0.5*chiPrefac) * (Iphi*Iphi);
		Adiel["Akappa"] = integral(shape * Akappa_shape);
		Adiel_shape += Akappa_shape //internal energy contribution
			- chiPrefac * (Iphi * IAdiel_rhoFluid); //propagate Coulomb contribution
		//Propagate neutrality constraint gradient contributions
		double Adiel_phi0 = chiPrefac * integral(shape * (Iphi - IAdiel_rhoFluid));
		Adiel_shape += Adiel_phi0 * phi0_shape;
		Adiel_rhoExplicitTilde->setGzero(Adiel_phi0 * phi0_Qexp);
	}
	
	//Propagate shape gradients to A_nCavity:
	ScalarField Adiel_nCavity;
	propagateCavityGradients(Adiel_shape, Adiel_nCavity, Adiel_rhoExplicitTilde, electricOnly);
	Adiel_nCavityTilde = J(Adiel_nCavity);
	
	setExtraForces(extraForces, Adiel_nCavityTilde);
	return Adiel;
}

void LinearPCM::loadState(const char* filename)
{	ScalarField Istate(ScalarFieldData::alloc(gInfo));
	loadRawBinary(Istate, filename); //saved data is in real space
	state = J(Istate);
}

void LinearPCM::saveState(const char* filename) const
{	if(mpiUtil->isHead()) saveRawBinary(I(state), filename); //saved data is in real space
}
