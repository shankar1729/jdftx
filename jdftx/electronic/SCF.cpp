/*-------------------------------------------------------------------
Copyright 2013 Deniz Gunceler

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
{
	// set up the cacheing size
	ResidualMinimizerParams& rp = e.residualMinimizerParams;
	if(rp.vectorExtrapolation == plain) // No history is needed for plain mixing
		rp.history = 1;	
	
	overlap.init(e.residualMinimizerParams.history, e.residualMinimizerParams.history);
	
	// Subspace rotation make no sense for residual minimize
	if(e.eInfo.fillingsUpdate != ElecInfo::ConstantFillings)
		e.eInfo.subspaceRotation = false;

	// Check eigenshifts and adjust index for 0
	for(size_t j=0; j<e.residualMinimizerParams.eigenShifts.size(); j++)
	{	
		EigenShift& eigenShift = e.residualMinimizerParams.eigenShifts[j];
		// Correct for the 0
		if(eigenShift.fromHOMO)
		{	eigenShift.n += e.eInfo.findHOMO(eigenShift.q);
			eigenShift.fromHOMO = false; // Needed if SCF is called multiple times, e.g. from an ionic minimize
		}
		// Check for a meaningful q and n
		if(eigenShift.q < 0)
			die("Eigenshift quantum number (q) must be greater than 0!\n");
		if(eigenShift.q > e.eInfo.nStates)
			die("Eigenshift quantum number (q) must be less than nStates=%i!\n", e.eInfo.nStates);
		if(eigenShift.n < 0)
			die("Eigenshift band index (n) must be greater than 0!\n");
		if(eigenShift.n > e.eInfo.nBands)
			die("Eigenshift band index (n) must be less than nBands=%i!\n", e.eInfo.nBands);		
		
	}
	
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
	ResidualMinimizerParams rp = e.residualMinimizerParams;
	
	// Compute energy for the initial guess
	e.cntrl.fixed_n = false;
	e.ener = Energies();
	e.iInfo.update(e.ener);
	double E;
	if(e.eInfo.fillingsUpdate != ElecInfo::ConstantFillings) // Compute Hsub and update fillings
	{	eVars.elecEnergyAndGrad(e.ener, 0, 0, true);
		updateFillings();
	}
	E = eVars.elecEnergyAndGrad(e.ener, 0, 0, 0);	// Compute energy

	e.ener = Energies(); 
	e.cntrl.fixed_n = true;
	
	// Initialize the variable that defines the single particle (Kohn-Sham) Hamiltonian
	// It can be either the densities (electronic and KE) or the potentials (Vscloc and Vtau)
	DataRptrCollection& variable_n = (rp.mixedVariable == density ? eVars.n : eVars.Vscloc);
	DataRptrCollection& variable_tau = (rp.mixedVariable == density ? eVars.tau : eVars.Vtau);
	logPrintf("\nWill mix electronic and kinetic %s at each iteration.\n", (rp.mixedVariable==density ? "density" : "potential"));
	
	// Set up variable history for vector extrapolation
	std::vector<DataRptrCollection> pastVariables_n, pastVariables_tau, pastResiduals_n, pastResiduals_tau;
	
	double Eprev = 0.;
	
	logPrintf("\n------------------- SCF Cycle ---------------------\n");
	for(int scfCounter=0; scfCounter<e.residualMinimizerParams.nIterations; scfCounter++)
	{	
		// Clear history if full
		if(((int)pastResiduals_n.size() >= rp.history) or ((int)pastVariables_n.size() >= rp.history))
		{	double ndim = pastResiduals_n.size();
			if((rp.vectorExtrapolation != plain) and (rp.history!=1))
				overlap.set(0, ndim-1, 0, ndim-1, overlap(1, ndim, 1, ndim));
			pastVariables_n.erase(pastVariables_n.begin());
			ifTau(pastVariables_tau.erase(pastVariables_tau.begin()))
			if(rp.vectorExtrapolation != plain)
			{	pastResiduals_n.erase(pastResiduals_n.begin());
				ifTau(pastResiduals_tau.erase(pastResiduals_tau.begin()););
			}
		}
		
		// Cache the old energy and variables
		Eprev = E;
		pastVariables_n.push_back(clone(variable_n));
		ifTau(pastVariables_tau.push_back(clone(variable_tau)))
		
		/// Solve at fixed hamiltonian ///
		e.cntrl.fixed_n = true; e.ener = Energies();
		if(not rp.verbose) { logSuspend(); e.elecMinParams.fpLog = nullLog; } // Silence eigensolver output
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{	BandMinimizer bmin(e, q, true);
			bmin.minimize(e.elecMinParams);
		}
		e.eVars.setEigenvectors();
		if(not rp.verbose) { logResume(); e.elecMinParams.fpLog = globalLog; }  // Resume output
		e.cntrl.fixed_n = false; e.ener = Energies();
		/// ///////////////////////// ///
		
		// Compute new density and energy
		if(e.eInfo.fillingsUpdate != ElecInfo::ConstantFillings) // Update fillings
			updateFillings();
		e.iInfo.update(e.ener);
		E = eVars.elecEnergyAndGrad(e.ener, 0, 0, 0);
		
		// Debug print
		if(e.cntrl.shouldPrintEigsFillings)
			print_Hsub_eigs(e);
		if(e.cntrl.shouldPrintEcomponents)
		{	logPrintf("\nEcomponents:\n");
			e.ener.print();
			logPrintf("\n");
		}
		logFlush();

		// Calculate residual
		DataRptrCollection variableResidual = variable_n - pastVariables_n.back();
		double residual = e.gInfo.dV*sqrt(dot(variableResidual, variableResidual));
		
		/// PRINT ///
		logPrintf("%sResidualMinimize: Iter:\t%i\tEtot: %.15e\tResidual:%.3e\tdE: %.3e\n", (rp.verbose ? "\t" : ""), scfCounter, E, residual, E-Eprev);
		
		// Check for convergence, mix density or potential if otherwise
		if(fabs(E-Eprev) < rp.energyDiffThreshold)
		{	logPrintf("Residual Minimization Converged (|Delta E|<%le).\n\n", rp.energyDiffThreshold);
			break;
		}
		else
		{	
			if((rp.vectorExtrapolation == plain))
				mixPlain(variable_n, variable_tau, pastVariables_n.back(), pastVariables_tau.back(), 1.-rp.damping);
			else if(rp.vectorExtrapolation == DIIS)
			{	cacheResidual(n)
				ifTau(cacheResidual(tau))
				mixDIIS(variable_n, variable_tau, pastVariables_n, pastVariables_tau, pastResiduals_n, pastResiduals_tau);
			}
		
			// Recompute Vscloc if mixing density
			if(e.residualMinimizerParams.mixedVariable == density)
				e.eVars.EdensityAndVscloc(e.ener);
			e.iInfo.augmentDensityGridGrad(e.eVars.Vscloc); //update Vscloc atom projections for ultrasoft psp's 
		}
		
		logFlush();
	}
}

void SCF::mixPlain(DataRptrCollection& variable_n, DataRptrCollection& variable_tau, 
					DataRptrCollection& mixingVariable_n, DataRptrCollection& mixingVariable_tau, double mixFraction)
{		// Mix old and new variable
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
	
	double damping = 1. - e.residualMinimizerParams.damping;
	
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
	ResidualMinimizerParams& rp = e.residualMinimizerParams;
	
	// Apply eigenshifts, if any
	for(size_t j=0; j<rp.eigenShifts.size(); j++)
		if(eInfo.isMine(rp.eigenShifts[j].q))
			e.eVars.Hsub_eigs[rp.eigenShifts[j].q][rp.eigenShifts[j].n] += rp.eigenShifts[j].shift;
	
	double mu; // Electron chemical potential
		
	//Update nElectrons from mu, or mu from nElectrons as appropriate:
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
	
	// Undo eigenshifts, if any
	for(size_t j=0; j<rp.eigenShifts.size(); j++)
		if(eInfo.isMine(rp.eigenShifts[j].q))
			e.eVars.Hsub_eigs[rp.eigenShifts[j].q][rp.eigenShifts[j].n] -= rp.eigenShifts[j].shift;
	
	// Print filling update information
	if(e.eInfo.fillingsUpdate)
	{	double muMix = e.eInfo.findMu(eVars.Hsub_eigs, e.eInfo.nElectrons);
		logPrintf("\tFillingsMix:  mu: %.15le  nElectrons: %.15le", muMix, e.eInfo.nElectrons);
		if(e.eInfo.spinType == SpinZ)
		{	double spinPol = integral(e.eVars.n[0] - e.eVars.n[1]);
			logPrintf("  magneticMoment: %.5f", spinPol);
		}
		logPrintf("\n");
		logFlush();
	}
}


