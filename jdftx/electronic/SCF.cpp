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

SCF::SCF(Everything& e): e(e)
{	//Set up the caching size:
	SCFparams& sp = e.scfParams;
	if(sp.vectorExtrapolation == SCFparams::VE_Plain) //No history needed for plain mixing
		sp.history = 1;
	overlap.init(sp.history, sp.history);
	
	eigenShiftInit();
}

#define ifTau(command) if(e.exCorr.needsKEdensity()) command;

#define overlapResiduals(r1, r2) \
	(integral(r1[0]*r2[0]) + (r1.size() == 2 ? integral(r1[1]*r2[1]) : 0.))
	
#define	cacheResidual(var) \
	pastResiduals_ ## var.push_back(variable_ ## var - pastVariables_ ## var.back());

void SCF::minimize()
{	
	ElecInfo& eInfo = e.eInfo;
	ElecVars& eVars = e.eVars;
	SCFparams sp = e.scfParams;
	
	bool subspaceRotation=false;
	std::swap(eInfo.subspaceRotation, subspaceRotation); //Switch off subspace rotation for SCF
	
	//Compute energy for the initial guess
	if(e.eInfo.fillingsUpdate != ElecInfo::ConstantFillings) //Compute Hsub and update fillings
	{	eVars.elecEnergyAndGrad(e.ener, 0, 0, true);
		updateFillings();
	}
	double E = eVars.elecEnergyAndGrad(e.ener); //Compute energy
	e.iInfo.augmentDensityGridGrad(e.eVars.Vscloc); //update Vscloc atom projections for ultrasoft psp's 
	
	//Initialize the variable that defines the single particle (Kohn-Sham) Hamiltonian
	//It can be either the densities (electronic and KE) or the potentials (Vscloc and Vtau)
	DataRptrCollection& variable_n = (sp.mixedVariable == SCFparams::MV_Density ? eVars.n : eVars.Vscloc);
	DataRptrCollection& variable_tau = (sp.mixedVariable == SCFparams::MV_Density ? eVars.tau : eVars.Vtau);
	logPrintf("\nWill mix electronic and kinetic %s at each iteration.\n", (sp.mixedVariable==SCFparams::MV_Density ? "density" : "potential"));
	
	//Set up variable history for vector extrapolation
	std::vector<DataRptrCollection> pastVariables_n, pastVariables_tau, pastResiduals_n, pastResiduals_tau;
	
	double Eprev = 0.;
	
	logPrintf("\n------------------- SCF Cycle ---------------------\n");
	for(int scfCounter=0; scfCounter<e.scfParams.nIterations; scfCounter++)
	{
		//Clear history if full
		if(((int)pastResiduals_n.size() >= sp.history) or ((int)pastVariables_n.size() >= sp.history))
		{	double ndim = pastResiduals_n.size();
			if((sp.vectorExtrapolation != SCFparams::VE_Plain) and (sp.history!=1))
				overlap.set(0, ndim-1, 0, ndim-1, overlap(1, ndim, 1, ndim));
			pastVariables_n.erase(pastVariables_n.begin());
			ifTau(pastVariables_tau.erase(pastVariables_tau.begin()))
			if(sp.vectorExtrapolation != SCFparams::VE_Plain)
			{	pastResiduals_n.erase(pastResiduals_n.begin());
				ifTau(pastResiduals_tau.erase(pastResiduals_tau.begin()););
			}
		}
		
		//Cache the old energy and variables
		Eprev = E;
		pastVariables_n.push_back(clone(variable_n));
		ifTau(pastVariables_tau.push_back(clone(variable_tau)))
		
		//Solve at fixed hamiltonian
		e.cntrl.fixed_n = true;
		if(not sp.verbose) { logSuspend(); e.elecMinParams.fpLog = nullLog; } // Silence eigensolver output
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{	BandMinimizer bmin(e, q, true);
			bmin.minimize(e.elecMinParams);
		}
		e.eVars.setEigenvectors();
		if(not sp.verbose) { logResume(); e.elecMinParams.fpLog = globalLog; }  // Resume output
		e.cntrl.fixed_n = false;
		
		//Compute new density and energy
		if(e.eInfo.fillingsUpdate != ElecInfo::ConstantFillings) // Update fillings
			updateFillings();
		E = eVars.elecEnergyAndGrad(e.ener);
		
		// Debug print
		if(e.cntrl.shouldPrintEigsFillings)
			print_Hsub_eigs(e);
		if(e.cntrl.shouldPrintEcomponents)
		{	logPrintf("\nEcomponents:\n");
			e.ener.print();
			logPrintf("\n");
		}
		logFlush();

		//Calculate residual
		DataRptrCollection variableResidual = variable_n - pastVariables_n.back();
		double residual = e.gInfo.dV*sqrt(dot(variableResidual, variableResidual));
		
		/// PRINT ///
		logPrintf("%sResidualMinimize: Iter:\t%i\tEtot: %.15e\tResidual:%.3e\tdE: %.3e\n", (sp.verbose ? "\t" : ""), scfCounter, E, residual, E-Eprev);
		
		// Check for convergence, mix density or potential if otherwise
		if(fabs(E-Eprev) < sp.energyDiffThreshold)
		{	logPrintf("Residual Minimization Converged (|Delta E|<%le).\n\n", sp.energyDiffThreshold);
			break;
		}
		else
		{	
			if((sp.vectorExtrapolation == SCFparams::VE_Plain))
				mixPlain(variable_n, variable_tau, pastVariables_n.back(), pastVariables_tau.back(), 1.-sp.damping);
			else if(sp.vectorExtrapolation == SCFparams::VE_DIIS)
			{	cacheResidual(n)
				ifTau(cacheResidual(tau))
				mixDIIS(variable_n, variable_tau, pastVariables_n, pastVariables_tau, pastResiduals_n, pastResiduals_tau);
			}
		
			// Recompute Vscloc if mixing density
			if(e.scfParams.mixedVariable == SCFparams::MV_Density)
				e.eVars.EdensityAndVscloc(e.ener);
			e.iInfo.augmentDensityGridGrad(e.eVars.Vscloc); //update Vscloc atom projections for ultrasoft psp's 
		}
		
		logFlush();
	}
	
	std::swap(eInfo.subspaceRotation, subspaceRotation); //Restore subspaceRotation to its original state
}

void SCF::mixPlain(DataRptrCollection& variable_n, DataRptrCollection& variable_tau, 
					DataRptrCollection& mixingVariable_n, DataRptrCollection& mixingVariable_tau, double mixFraction)
{	//Mix old and new variable
	variable_n = mixFraction*variable_n + (1.-mixFraction)*mixingVariable_n;
	ifTau(variable_tau = mixFraction*variable_tau + (1.-mixFraction)*mixingVariable_tau)
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

void SCF::mixDIIS(DataRptrCollection& variable_n, DataRptrCollection& variable_tau, 
				  std::vector< DataRptrCollection >& pastVariables_n, std::vector< DataRptrCollection >& pastVariables_tau, 
				  std::vector< DataRptrCollection >& pastResiduals_n, std::vector< DataRptrCollection >& pastResiduals_tau)
{
	size_t ndim = pastResiduals_n.size();
	
	// Update the overlap matrix
	for(size_t j=0; j<ndim; j++)
	{	double thisOverlap = dot(pastResiduals_n[j], pastResiduals_n.back());
		overlap.set(j, ndim-1, thisOverlap);
		if(j != ndim-1) overlap.set(ndim-1, j, thisOverlap);
	}
	
	// If the number of dimensions is less than 0, do a plain mixing step
	if(ndim < 3)
	{	mixPlain(variable_n, variable_tau, pastVariables_n.back(), pastVariables_tau.back());
		return;
	}
	
	// diagonalize the residual overlap matrix to get the minimum of residual
	matrix cOverlap(ndim+1, ndim+1); // Add row and column to enforce normalization constraint
	cOverlap.set(0, ndim, 0, ndim, overlap(0, ndim, 0, ndim));
	for(size_t j=0; j<ndim; j++)
	{	cOverlap.set(j, ndim, 1);
		cOverlap.set(ndim, j, 1);
	}
	cOverlap.set(ndim, ndim, 0);
	
	matrix cOverlap_inv = inv(cOverlap);
	
	// Zero variables
	initZero(variable_n);
	ifTau(initZero(variable_tau))
	
	double damping = 1. - e.scfParams.damping;
	
	// Variables for the preconditioner
	//double Gmax = e.gInfo.GmaxGrid;
	
	// Accumulate
	complex* coefs = cOverlap_inv.data();
	for(size_t j=0; j<ndim; j++)
	{	variable_n += coefs[cOverlap_inv.index(j, ndim)].real() * ((1.-damping)*pastResiduals_n[j] + pastVariables_n[j]);
		ifTau(variable_tau +=coefs[cOverlap_inv.index(j, ndim)].real() * ((1.-damping)*pastResiduals_tau[j] + pastVariables_tau[j]))		
	}

	logPrintf("\tDIIS acceleration, lagrange multiplier is %.3e\n", coefs[cOverlap_inv.index(ndim, ndim)].real());

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
	{	logPrintf("\tFillingsMix:  mu: %.15le  nElectrons: %.15le", mu, e.eInfo.nElectrons);
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
