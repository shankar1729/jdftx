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
}

void SCF::minimize()
{
	ElecInfo& eInfo = e.eInfo;
	ElecVars& eVars = e.eVars;
	ResidualMinimizerParams& rp = e.residualMinimizerParams;
	
	// Compute energy for the initial guess
	e.cntrl.fixed_n = false;
	e.ener = Energies();
	e.iInfo.update(e.ener);
	eVars.elecEnergyAndGrad(e.ener, 0, 0, 0);	
	e.ener = Energies(); 
	e.cntrl.fixed_n = true;
	
	// Initialize the variable that defines the single particle (Kohn-Sham) Hamiltonian
	// It can be either the densities (electronic and KE) or the potentials (Vscloc and Vtau)
	DataRptrCollection& variable_n = (rp.mixedVariable == density ? eVars.n : eVars.Vscloc);
	DataRptrCollection& variable_tau = (rp.mixedVariable == density ? eVars.tau : eVars.Vtau);
	logPrintf("\nWill mix electronic and kinetic potential %s at each iteration.\n", (rp.mixedVariable==density ? "density" : "potential"));
	
	// Set up variable history for vector extrapolation
	std::vector<DataRptrCollection> pastVariables_n, pastVariables_tau, pastResiduals;
	
	#define ifTau(command) if(e.exCorr.needsKEdensity()) command;
	
	double Eprev = 0, E;
	
	logPrintf("\n------------------- SCF Cycle ---------------------\n");
	for(int scfCounter=0; scfCounter<e.residualMinimizerParams.nIterations; scfCounter++)
	{	
		// Clear history if full
		if((pastResiduals.size() >= rp.history) or (pastResiduals.size() >= rp.history))
		{	pastVariables_n.clear();
			pastVariables_tau.clear();
			ifTau(pastResiduals.clear())
		}
		
		// Cache the old energy and variables
		Eprev = E;
		pastVariables_n.push_back(clone(variable_n));
		ifTau(pastVariables_tau.push_back(clone(variable_tau)))
		
		// Solve at fixed hamiltonian
		e.cntrl.fixed_n = true; e.ener = Energies();
		logSuspend();
		for(int q = 0; q < eInfo.nStates; q++)
		{	
			BandMinimizer bmin(e, q, true);
			bmin.minimize(e.elecMinParams);
		}
		logResume();
		e.cntrl.fixed_n = false; e.ener = Energies();
		
		// Compute new density and energy
		e.iInfo.update(e.ener);
		E = eVars.elecEnergyAndGrad(e.ener, 0, 0, 0);
				
		logPrintf("SCF Iter: %i\tEprev: %f\tdE: %.2e\tEtot: %f\n\n", scfCounter, Eprev, fabs(E-Eprev), E);
		
		// Check for convergence, mix density or potential if otherwise
		if(fabs(E-Eprev) < rp.energyDiffThreshold)
		{	logPrintf("Residual Minimization Converged (|Delta E|<%le).\n", rp.energyDiffThreshold);
			break;
		}
		else
		{	mixHamiltonian(variable_n, variable_tau, pastVariables_n.back(), pastVariables_tau.back());
		}
	}
	
}

void SCF::mixHamiltonian(DataRptrCollection& variable_n, DataRptrCollection& variable_tau, 
						 DataRptrCollection& prevVariable_n, DataRptrCollection& prevVariable_tau, double mixFraction)
{
		// Mix old and new variable
		for(size_t s=0; s<e.eVars.n.size(); s++)
		{	variable_n[s] = mixFraction*variable_n[s] + (1.-mixFraction)*prevVariable_n[s];
			if(e.exCorr.needsKEdensity())
				e.eVars.tau[s] = mixFraction*variable_tau[s] + (1.-mixFraction)*prevVariable_tau[s];
		}

		// Recompute Vscloc if mixing density
		if(e.residualMinimizerParams.mixedVariable == density)
			e.eVars.EdensityAndVscloc(e.ener);

}