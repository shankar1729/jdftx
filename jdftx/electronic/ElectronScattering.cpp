/*-------------------------------------------------------------------
Copyright 2015 Ravishankar Sundararaman

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

#include <electronic/ElectronScattering.h>
#include <electronic/Everything.h>

ElectronScattering::ElectronScattering()
: eta(0.), Ecut(0.), fCut(1e-6), omegaMax(0.)
{
}

void ElectronScattering::dump(const Everything& e)
{
	logPrintf("\n----- Electron-electron scattering Im(Sigma) -----\n"); logFlush();

	//Update default parameters:
	if(!eta)
	{	eta = e.eInfo.kT;
		if(!eta) die("eta must be specified explicitly since electronic temperature is zero.\n");
	}
	if(!Ecut) Ecut = e.cntrl.Ecut;
	if(!omegaMax)
	{	//Determine maximum possible energy transfer:
		double oMin = DBL_MAX, oMax = -DBL_MAX; //occupied energy range
		double uMin = DBL_MAX, uMax = -DBL_MAX; //unoccupied energy range
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
			for(int b=0; b<e.eInfo.nBands; b++)
			{	double E = e.eVars.Hsub_eigs[q][b];
				double f = e.eVars.F[q][b];
				if(f > fCut) //sufficiently occupied
				{	oMin = std::min(oMin, E);
					oMax = std::max(oMax, E);
				}
				if(f < 1.-fCut) //sufficiently unoccupied
				{	uMin = std::min(uMin, E);
					uMax = std::max(uMax, E);
				}
			}
		mpiUtil->allReduce(oMin, MPIUtil::ReduceMin);
		mpiUtil->allReduce(oMax, MPIUtil::ReduceMax);
		mpiUtil->allReduce(uMin, MPIUtil::ReduceMin);
		mpiUtil->allReduce(uMax, MPIUtil::ReduceMax);
		omegaMax = std::max(uMax-uMin, oMax-oMin);
	}
	//--- print selected values after fixing defaults:
	logPrintf("Frequency resolution:    %lg\n", eta);
	logPrintf("Dielectric matrix Ecut:  %lg\n", Ecut);
	logPrintf("Maximum energy transfer: %lg\n", omegaMax);
	
	
	die("Not yet implemented.\n");
	logPrintf("\n"); logFlush();
}
