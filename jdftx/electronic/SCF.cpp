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
	int maxIter = 10;
	double energyDiffThreshold = 1e-5;
	
	ElecInfo& eInfo = e.eInfo;
	ElecVars& eVars = e.eVars;
	
	e.cntrl.fixed_n = false;
	e.ener = Energies();
	e.iInfo.update(e.ener);
	eVars.elecEnergyAndGrad(e.ener, 0, 0, 0);
	e.ener.print();			
	e.ener = Energies(); 
	e.cntrl.fixed_n = true;
	
	DataRptrCollection n_or_Vscloc_prev(eVars.n.size());
	DataRptrCollection tau_or_Vtau_prev(eVars.n.size());
	double Eprev = 0, E;
	
	logPrintf("\n------------------- SCF Cycle ---------------------\n");
	for(int scfCounter=0; scfCounter<maxIter; scfCounter++)
	{	// Cache the old energy and density
		Eprev = E;
		for(size_t s=0; s<eVars.n.size(); s++)
		{	n_or_Vscloc_prev[s] = eVars.n[s];
			if(e.exCorr.needsKEdensity())
				tau_or_Vtau_prev[s] = e.eVars.tau[s];
		}
		
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
		eVars.elecEnergyAndGrad(e.ener, 0, 0, 0);
		E = relevantFreeEnergy(e);
		
		// Mix old and new density
		double mixFraction = 0.5;
		for(size_t s=0; s<eVars.n.size(); s++)
		{	eVars.n[s] = mixFraction*eVars.n[s] + (1.-mixFraction)*n_or_Vscloc_prev[s];
			if(e.exCorr.needsKEdensity())
				eVars.tau[s] = mixFraction*eVars.tau[s] + (1.-mixFraction)*tau_or_Vtau_prev[s];
		}

		// Recompute Vscloc using mixed density
		eVars.EdensityAndVscloc(e.ener);
		e.ener.print();
			
		logPrintf("SCF Iter: %i\tmixFraction: %f\tEprev: %f\tdE: %.2e\tEtot: %f\n\n", scfCounter, mixFraction, Eprev, fabs(E-Eprev), E);
		
		// Check for convergence
		if(fabs(E-Eprev) < energyDiffThreshold)
			break;
	}
	
}