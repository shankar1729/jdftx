/*-------------------------------------------------------------------
Copyright 2013 Deniz Gunceler, Ravishankar Sundararaman

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

#include <electronic/SCF.h>
#include <electronic/ElecMinimizer.h>

SCF::SCF(Everything& e): e(e)
{	//Set up the caching size:
	SCFparams& sp = e.scfParams;
	if(sp.vectorExtrapolation == SCFparams::VE_Plain) //No history needed for plain mixing
		sp.history = 1;
	overlap.init(sp.history, sp.history);
	
	eigenShiftInit();
}

void SCF::minimize()
{	
	ElecInfo& eInfo = e.eInfo;
	ElecVars& eVars = e.eVars;
	SCFparams sp = e.scfParams;
	
	logPrintf("\tWill mix electronic %s%s at each iteration.\n",
		(e.exCorr.needsKEdensity() ? "and kinetic " : ""),
		(sp.mixedVariable==SCFparams::MV_Density ? "density" : "potential"));

	bool subspaceRotation=false;
	std::swap(eInfo.subspaceRotation, subspaceRotation); //Switch off subspace rotation for SCF
	
	//Compute energy for the initial guess
	if(e.eInfo.fillingsUpdate != ElecInfo::ConstantFillings) //Compute Hsub and update fillings
	{	eVars.elecEnergyAndGrad(e.ener, 0, 0, true);
		updateFillings();
	}
	double E = eVars.elecEnergyAndGrad(e.ener); //Compute energy
	e.iInfo.augmentDensityGridGrad(e.eVars.Vscloc); //update Vscloc atom projections for ultrasoft psp's 
	
	double Eprev = 0.;
	for(int scfCounter=0; scfCounter<e.scfParams.nIterations; scfCounter++)
	{
		//If history is full, remove oldest member
		if(((int)pastResiduals.size() >= sp.history) or ((int)pastVariables.size() >= sp.history))
		{	size_t ndim = pastResiduals.size();
			if(ndim>1) overlap.set(0,ndim-1, 0,ndim-1, overlap(1,ndim, 1,ndim));
			pastVariables.erase(pastVariables.begin());
			pastResiduals.erase(pastResiduals.begin());
		}
		
		//Cache the old energy and variables
		Eprev = E;
		pastVariables.push_back(clone(getVariable()));
		
		//Band-structure minimize:
		if(not sp.verbose) { logSuspend(); e.elecMinParams.fpLog = nullLog; } // Silence eigensolver output
		bandMinimize(e);
		if(not sp.verbose) { logResume(); e.elecMinParams.fpLog = globalLog; }  // Resume output
		
		//Compute new density and energy
		if(e.eInfo.fillingsUpdate != ElecInfo::ConstantFillings) // Update fillings
			updateFillings();
		E = eVars.elecEnergyAndGrad(e.ener);
		const char* Elabel = relevantFreeEnergyName(e);
		if(e.exCorr.orbitalDep) //total E meaningless, use Eband instead
		{	E = e.ener.Eband;
			Elabel = "Eband";
		}
		
		// Debug print
		if(e.cntrl.shouldPrintEigsFillings)
			print_Hsub_eigs(e);
		if(e.cntrl.shouldPrintEcomponents)
		{	logPrintf("\nEcomponents:\n");
			e.ener.print();
			logPrintf("\n");
		}
		logFlush();

		//Calculate and cache residual:
		double residualNorm = 0.;
		{	DataRptrCollection residual = getVariable() - pastVariables.back();
			pastResiduals.push_back(residual);
			residualNorm = e.gInfo.dV * sqrt(dot(residual,residual));
		}
		logPrintf("SCF: Cycle: %2i   %s: %.15e   |Residual|: %.3e   dE: %.3e\n", scfCounter, Elabel, E, residualNorm, E-Eprev);
		
		//Check for convergence and update variable:
		if(fabs(E - Eprev) < sp.energyDiffThreshold) { logPrintf("SCF: Converged (|Delta E|<%le).\n\n", sp.energyDiffThreshold); break; }
		else if(residualNorm < sp.residualThreshold) { logPrintf("SCF: Converged (|Residual|<%le).\n\n", sp.residualThreshold); break; }
		else
		{	switch(sp.vectorExtrapolation)
			{	case SCFparams::VE_Plain: mixPlain(); break;
				case SCFparams::VE_DIIS: mixDIIS(); break;
			}
		}
		logFlush();
	}
	
	std::swap(eInfo.subspaceRotation, subspaceRotation); //Restore subspaceRotation to its original state
}

void SCF::mixPlain()
{	setVariable(pastVariables.back() + (1.-e.scfParams.damping)*pastResiduals.back());
}

inline double preconditionerKernel(double G, double f, double slope)
{	return (G ? f + slope*G : 1.);
}

DataRptrCollection precondition(DataRptrCollection& n, double Gmax, double f = 20.)
{
	// Set up preconditioner kernel for density/potential overlaps
	RadialFunctionG preconditioner;
	double slope = (1.-f)/Gmax;
	preconditioner.init(0, 0.02, Gmax, preconditionerKernel, f, slope);
	
	DataGptrCollection preconditioned(n.size());
	for(size_t i=0; i<n.size(); i++)
		preconditioned[i] = preconditioner * J(n[i]);
	
	return I(preconditioned);
}

void SCF::mixDIIS()
{	size_t ndim = pastResiduals.size();
	
	//Update the overlap matrix
	for(size_t j=0; j<ndim; j++)
	{	double thisOverlap = dot(pastResiduals[j], pastResiduals.back());
		overlap.set(j, ndim-1, thisOverlap);
		if(j!=ndim-1) overlap.set(ndim-1, j, thisOverlap);
	}
	
	//Plain mixing on the first few steps:
	if(ndim<3) { mixPlain(); return; }
	
	//Invert the residual overlap matrix to get the minimum of residual
	matrix cOverlap(ndim+1, ndim+1); // Add row and column to enforce normalization constraint
	cOverlap.set(0, ndim, 0, ndim, overlap(0, ndim, 0, ndim));
	for(size_t j=0; j<ndim; j++)
	{	cOverlap.set(j, ndim, 1);
		cOverlap.set(ndim, j, 1);
	}
	cOverlap.set(ndim, ndim, 0);
	matrix cOverlap_inv = inv(cOverlap);
	
	//Preconditioner:
	//TODO
	
	//Update variable:
	complex* coefs = cOverlap_inv.data();
	DataRptrCollection variable(pastVariables.back().size()); //zero
	for(size_t j=0; j<ndim; j++)
		variable += coefs[cOverlap_inv.index(j, ndim)].real() * (pastVariables[j] + (1.-e.scfParams.damping)*pastResiduals[j]);
	setVariable(variable);
	logPrintf("\tDIIS acceleration, lagrange multiplier is %.3e\n", coefs[cOverlap_inv.index(ndim, ndim)].real());
}

DataRptrCollection SCF::getVariable() const
{	bool mixDensity = (e.scfParams.mixedVariable==SCFparams::MV_Density);
	DataRptrCollection variable = mixDensity ? e.eVars.n : e.eVars.Vscloc;
	if(e.exCorr.needsKEdensity()) //append the relevant KE variable:
		for(DataRptr v: (mixDensity ? e.eVars.tau : e.eVars.Vtau)) variable.push_back(v);
	return variable;
}

void SCF::setVariable(DataRptrCollection variable)
{	bool mixDensity = (e.scfParams.mixedVariable==SCFparams::MV_Density);
	size_t iVar=0;
	for(DataRptr& v: (mixDensity ? e.eVars.n : e.eVars.Vscloc)) v = variable[iVar++];
	if(e.exCorr.needsKEdensity())
		for(DataRptr& v: (mixDensity ? e.eVars.tau : e.eVars.Vtau)) v = variable[iVar++];
	assert(iVar == variable.size());
	if(mixDensity) e.eVars.EdensityAndVscloc(e.ener); //Recompute Vscloc if mixing density
	e.iInfo.augmentDensityGridGrad(e.eVars.Vscloc); //update Vscloc atom projections for ultrasoft psp's 
}


void SCF::updateFillings()
{
	ElecInfo& eInfo = e.eInfo; 
	ElecVars& eVars = e.eVars;
	
	//Update nElectrons from mu, or mu from nElectrons as appropriate:
	eigenShiftApply(false);
	double mu; // Electron chemical potential
	if(std::isnan(eInfo.mu)) mu = eInfo.findMu(eVars.Hsub_eigs, eInfo.nElectrons);
	else
	{	mu = eInfo.mu; 
		((ElecInfo&)eInfo).nElectrons = eInfo.nElectronsFermi(mu, eVars.Hsub_eigs); 
	}
	//Compute fillings from aux hamiltonian eigenvalues:
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		eVars.F[q] = eInfo.fermi(mu, eVars.Hsub_eigs[q]);
	//Update TS and muN:
	eInfo.updateFillingsEnergies(e.eVars.F, e.ener);
	eigenShiftApply(true);
	
	// Print filling update information
	if(e.eInfo.fillingsUpdate)
	{	logPrintf("\tFillingsUpdate:  mu: %.15le  nElectrons: %.15le", mu, e.eInfo.nElectrons);
		if(e.eInfo.spinType == SpinZ)
		{	double spinPol = integral(e.eVars.n[0] - e.eVars.n[1]);
			logPrintf("  magneticMoment: %.5f", spinPol);
		}
		logPrintf("\n"); logFlush();
	}
}

void SCF::eigenShiftInit()
{	for(SCFparams::EigenShift& es: e.scfParams.eigenShifts)
	{	// Correct for the 0
		if(es.fromHOMO)
		{	es.n += e.eInfo.findHOMO(es.q);
			es.fromHOMO = false; // Needed if SCF is called multiple times, e.g. from an ionic minimize
		}
		// Check for a meaningful q and n
		if(es.q < 0) die("Eigenshift quantum number (q) must be greater than 0!\n");
		if(es.q > e.eInfo.nStates) die("Eigenshift quantum number (q) must be less than nStates=%i!\n", e.eInfo.nStates);
		if(es.n < 0) die("Eigenshift band index (n) must be greater than 0!\n");
		if(es.n > e.eInfo.nBands) die("Eigenshift band index (n) must be less than nBands=%i!\n", e.eInfo.nBands);		
	}
}

void SCF::eigenShiftApply(bool reverse)
{	int sign = reverse ? -1 : +1;
	for(const SCFparams::EigenShift& es: e.scfParams.eigenShifts)
		if(e.eInfo.isMine(es.q))
			e.eVars.Hsub_eigs[es.q][es.n] += sign * es.shift;
}
