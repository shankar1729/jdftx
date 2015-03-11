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

	//Info:
	logPrintf("\tLinear fluid (dielectric constant: %g", epsBulk);
	if(k2factor) logPrintf(", screening length: %g Bohr", sqrt(epsBulk/k2factor));
	logPrintf(") occupying %lf of unit cell:", integral(shape)/gInfo.detR); logFlush();

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
{
	fprintf(e.fluidMinParams.fpLog, "\n\tWill stop at %d iterations, or sqrt(|r.z|)<%le\n", e.fluidMinParams.nIterations, e.fluidMinParams.knormThreshold);
	int nIter = solve(rhoExplicitTilde, e.fluidMinParams);
	logPrintf("\tCompleted after %d iterations.\n", nIter);
}

double LinearPCM::get_Adiel_and_grad_internal(ScalarFieldTilde& Adiel_rhoExplicitTilde, ScalarFieldTilde& Adiel_nCavityTilde, IonicGradient* extraForces, bool electricOnly) const
{
	EnergyComponents& Adiel = ((LinearPCM*)this)->Adiel;
	const ScalarFieldTilde& phi = state; // that's what we solved for in minimize

	//First-order correct estimate of electrostatic energy:
	ScalarFieldTilde phiExt = coulomb(rhoExplicitTilde);
	Adiel["Electrostatic"] = -0.5*dot(phi, O(hessian(phi))) + dot(phi - 0.5*phiExt, O(rhoExplicitTilde));
	
	//Gradient w.r.t rhoExplicitTilde:
	Adiel_rhoExplicitTilde = phi - phiExt;
	
	//Compute gradient w.r.t shape function:
	ScalarField Adiel_shape = (-(epsBulk-1)/(8*M_PI)) * lengthSquared(I(gradient(phi))); //dielectric contributions
	if(k2factor) Adiel_shape -= (k2factor/(8*M_PI)) * pow(I(phi),2); //ionic contributions
	
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
