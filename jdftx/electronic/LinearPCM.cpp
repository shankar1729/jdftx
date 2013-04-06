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
#include <core/DataMultiplet.h>
#include <core/DataIO.h>
#include <core/Thread.h>

LinearPCM::LinearPCM(const Everything& e, const FluidSolverParams& fsp)
: PCM(e, fsp), Kkernel(e.gInfo)
{
	logPrintf("   Cavity determined by nc: %lg and sigma: %lg\n", fsp.nc, fsp.sigma);
	Cavitation::print(fsp);
}

DataGptr LinearPCM::hessian(const DataGptr& phiTilde)
{	DataRptr epsilon = 1. + (fsp.epsBulk-1.) * shape;
	DataGptr rhoTilde = divergence(J(epsilon * I(gradient(phiTilde))));  //dielectric term
	if(k2factor) rhoTilde -= k2factor * J(shape*I(phiTilde)); // screening term
	return rhoTilde;
}

DataGptr LinearPCM::precondition(const DataGptr& rTilde)
{	return Kkernel*(J(epsInv*I(Kkernel*rTilde)));
}

//Initialize Kkernel to square-root of the inverse kinetic operator
inline void setPreconditionerKernel(int i, double G2, double* Kkernel, double kRMS)
{	if(i==0) Kkernel[i] = kRMS ? 1.0/kRMS : 0.0;
	else Kkernel[i] = 1.0/sqrt(G2 + pow(kRMS,2));
}

void LinearPCM::set(const DataGptr& rhoExplicitTilde, const DataGptr& nCavityTilde)
{
	this->rhoExplicitTilde = clone(rhoExplicitTilde); zeroNyquist(this->rhoExplicitTilde);
	this->nCavity = I(nCavityTilde);

	updateCavity();
	
	//Info:
	logPrintf("\tLinear fluid (dielectric constant: %g", fsp.epsBulk);
	if(fsp.ionicConcentration) logPrintf(", screening length: %g Bohr", sqrt(fsp.epsBulk/k2factor));
	logPrintf(") occupying %lf of unit cell:", integral(shape)/e.gInfo.detR); logFlush();

	//Update the preconditioner
	DataRptr epsilon = 1 + (fsp.epsBulk-1)*shape;
	DataRptr kappaSq = fsp.ionicConcentration ? k2factor*shape : 0; //set kappaSq to null pointer if no screening
	epsInv = inv(epsilon);
	double kRMS = (kappaSq ? sqrt(sum(kappaSq)/sum(epsilon)) : 0.0);
	applyFuncGsq(Kkernel.gInfo, setPreconditionerKernel, Kkernel.data, kRMS);
	Kkernel.set();
	
	//Initialize the state if it hasn't been loaded:
	if(!state) nullToZero(state, e.gInfo);
}


void LinearPCM::minimizeFluid()
{
	fprintf(e.fluidMinParams.fpLog, "\n\tWill stop at %d iterations, or sqrt(|r.z|)<%le\n", e.fluidMinParams.nIterations, e.fluidMinParams.knormThreshold);
	int nIter = solve((-4*M_PI)*rhoExplicitTilde, e.fluidMinParams);
	logPrintf("\tCompleted after %d iterations.\n", nIter);
}

double LinearPCM::get_Adiel_and_grad(DataGptr& Adiel_rhoExplicitTilde, DataGptr& Adiel_nCavityTilde, IonicGradient& extraForces) const
{
	EnergyComponents& Adiel = ((LinearPCM*)this)->Adiel;
	const DataGptr& phi = state; // that's what we solved for in minimize

	//The "electrostatic" gradient is the potential due to the bound charge alone:
	Adiel_rhoExplicitTilde = phi - (-4*M_PI)*Linv(O(rhoExplicitTilde));
	Adiel["Electrostatic"] = 0.5*dot(Adiel_rhoExplicitTilde, O(rhoExplicitTilde));
	
	//Compute gradient w.r.t shape function:
	//--- Dielectric contributions:
	DataRptrVec gradPhi = I(gradient(phi));
	DataRptr gradPhiSq = gradPhi[0]*gradPhi[0] + gradPhi[1]*gradPhi[1] + gradPhi[2]*gradPhi[2];
	DataRptr Adiel_shape = (-(fsp.epsBulk-1)/(8*M_PI)) * gradPhiSq; //dielectric part
	//--- Screening contributions:
	if(fsp.ionicConcentration)
	{	DataRptr Iphi = I(phi); //potential in real space
		Adiel_shape += (k2factor/(8*M_PI)) * (Iphi*Iphi);
	}
	
	//Propagate shape gradients to A_nCavity:
	DataRptr Adiel_nCavity;
	propagateCavityGradients(Adiel_shape, Adiel_nCavity);
	Adiel_nCavityTilde = J(Adiel_nCavity);
	
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

void LinearPCM::dumpDensities(const char* filenamePattern) const
{	string filename(filenamePattern);
	filename.replace(filename.find("%s"), 2, "Shape");
	logPrintf("Dumping '%s'... ", filename.c_str());  logFlush();
	saveRawBinary(shape, filename.c_str());
	logPrintf("done.\n"); logFlush();
}

void LinearPCM::dumpDebug(const char* filenamePattern) const
{
	// Prepares to dump
	string filename(filenamePattern);
	filename.replace(filename.find("%s"), 2, "Debug");
	logPrintf("Dumping '%s'... \t", filename.c_str());  logFlush();

	FILE* fp = fopen(filename.c_str(), "w");
	if(!fp) die("Error opening %s for writing.\n", filename.c_str());	
	
    PCM::dumpDebug(fp);
	
	fprintf(fp, "\n\nGradients wrt fit parameters:");
	DataRptrVec shape_x = gradient(shape);
	DataRptr surfaceDensity = sqrt(shape_x[0]*shape_x[0] + shape_x[1]*shape_x[1] + shape_x[2]*shape_x[2]);
	DataGptr Adiel_nCavityTilde;
	DataGptr Adiel_rhoExplicitTilde;
	IonicGradient temp; get_Adiel_and_grad(Adiel_rhoExplicitTilde, Adiel_nCavityTilde, temp);
	fprintf(fp, "\nE_nc = %f", integral(I(Adiel_nCavityTilde)*(-(1./fsp.nc)*nCavity)));
	fprintf(fp, "\nE_t = %f", integral(surfaceDensity));
	
	fclose(fp);
	logPrintf("done\n"); logFlush();
}
