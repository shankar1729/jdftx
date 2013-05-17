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
#include <electronic/LinearPCM.h>
#include <electronic/PCM_internal.h>
#include <electronic/operators.h>
#include <core/DataMultiplet.h>
#include <core/DataIO.h>
#include <core/Thread.h>

LinearPCM::LinearPCM(const Everything& e, const FluidSolverParams& fsp)
: PCM(e, fsp)
{
}

LinearPCM::~LinearPCM()
{	Kkernel.free();
}


DataGptr LinearPCM::hessian(const DataGptr& phiTilde) const
{	DataRptr epsilon = 1. + (epsBulk-1.) * shape;
	DataGptr rhoTilde = divergence(J(epsilon * I(gradient(phiTilde))));  //dielectric term
	if(k2factor) rhoTilde -= k2factor * J(shape*I(phiTilde)); // screening term
	return (-1./(4*M_PI)) * rhoTilde;
}

DataGptr LinearPCM::precondition(const DataGptr& rTilde) const
{	return Kkernel*(J(epsInv*I(Kkernel*rTilde)));
}

//Initialize Kkernel to square-root of the inverse kinetic operator
inline double setPreconditionerKernel(double G, double kRMS)
{	return (G || kRMS) ? 1./hypot(G, kRMS) : 0.;
}

void LinearPCM::set(const DataGptr& rhoExplicitTilde, const DataGptr& nCavityTilde)
{
	this->rhoExplicitTilde = clone(rhoExplicitTilde); zeroNyquist(this->rhoExplicitTilde);
	this->nCavity = I(nCavityTilde);

	updateCavity();
	
	//Info:
	logPrintf("\tLinear fluid (dielectric constant: %g", epsBulk);
	if(k2factor) logPrintf(", screening length: %g Bohr", sqrt(epsBulk/k2factor));
	logPrintf(") occupying %lf of unit cell:", integral(shape)/e.gInfo.detR); logFlush();

	//Update the preconditioner
	DataRptr epsilon = 1 + (epsBulk-1)*shape;
	DataRptr kappaSq = k2factor ? k2factor*shape : 0; //set kappaSq to null pointer if no screening
	epsInv = inv(epsilon);
	double kRMS = (kappaSq ? sqrt(sum(kappaSq)/sum(epsilon)) : 0.0);
	Kkernel.init(0, 0.02, e.gInfo.GmaxGrid, setPreconditionerKernel, kRMS);
	
	//Initialize the state if it hasn't been loaded:
	if(!state) nullToZero(state, e.gInfo);
}


void LinearPCM::minimizeFluid()
{
	fprintf(e.fluidMinParams.fpLog, "\n\tWill stop at %d iterations, or sqrt(|r.z|)<%le\n", e.fluidMinParams.nIterations, e.fluidMinParams.knormThreshold);
	int nIter = solve(rhoExplicitTilde, e.fluidMinParams);
	logPrintf("\tCompleted after %d iterations.\n", nIter);
}

double LinearPCM::get_Adiel_and_grad(DataGptr& Adiel_rhoExplicitTilde, DataGptr& Adiel_nCavityTilde, IonicGradient& extraForces) const
{
	EnergyComponents& Adiel = ((LinearPCM*)this)->Adiel;
	const DataGptr& phi = state; // that's what we solved for in minimize

	//The "electrostatic" gradient is the potential due to the bound charge alone:
	Adiel_rhoExplicitTilde = phi - (-4*M_PI)*Linv(O(rhoExplicitTilde));
	Adiel["Electrostatic"] = 0.5*dot(Adiel_rhoExplicitTilde, O(rhoExplicitTilde)) //True energy if phi was an exact solution
		+ 0.5*dot(O(phi), rhoExplicitTilde - hessian(phi)); //First order residual correction (remaining error is second order)
	
	//Compute gradient w.r.t shape function:
	//--- Dielectric contributions:
	DataRptr Adiel_shape = (-(epsBulk-1)/(8*M_PI)) * lengthSquared(I(gradient(phi)));
	//--- Screening contributions:
	if(k2factor)
	{	DataRptr Iphi = I(phi);
		Adiel_shape -= (k2factor/(8*M_PI)) * (Iphi*Iphi);
	}
	
	//Propagate shape gradients to A_nCavity:
	DataRptr Adiel_nCavity;
	propagateCavityGradients(Adiel_shape, Adiel_nCavity);
	Adiel_nCavityTilde = J(Adiel_nCavity);
	
	if(vdwForces) extraForces = *vdwForces;
	return Adiel;
}

void LinearPCM::loadState(const char* filename)
{	DataRptr Istate(DataR::alloc(e.gInfo));
	loadRawBinary(Istate, filename); //saved data is in real space
	state = J(Istate);
}

void LinearPCM::saveState(const char* filename) const
{	saveRawBinary(I(state), filename); //saved data is in real space
}
