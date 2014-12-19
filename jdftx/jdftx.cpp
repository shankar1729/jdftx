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
#include <electronic/ColumnBundle.h>
#include <electronic/matrix.h>
#include <electronic/Dump.h>
#include <electronic/ElecMinimizer.h>
#include <electronic/LatticeMinimizer.h>
#include <electronic/InverseKohnSham.h>
#include <electronic/Vibrations.h>
#include <fluid/FluidSolver.h>
#include <core/Util.h>
#include <commands/parser.h>

//Program entry point
int main(int argc, char** argv)
{	//Parse command line, initialize system and logs:
	string inputFilename; bool dryRun, printDefaults;
	initSystemCmdline(argc, argv, "Performs Joint Density Functional Theory calculations.", inputFilename, dryRun, printDefaults);
	
	//Parse input file and setup
	Everything e; //the parent data structure for, well, everything
	ElecVars& eVars = e.eVars;
	parse(readInputFile(inputFilename), e, printDefaults);
	e.setup();
	Citations::print();
	if(dryRun)
	{	logPrintf("Dry run successful: commands are valid and initialization succeeded.\n");
		finalizeSystem();
		return 0;
	}
	
	if(e.cntrl.dumpOnly)
	{	//Single energy calculation so that all dependent quantities have been initialized:
		logPrintf("\n----------- Energy evaluation at fixed state -------------\n"); logFlush();
		eVars.elecEnergyAndGrad(e.ener, 0, 0, true); //calculate Hsub so that eigenvalues are available (used by many dumps)
		logPrintf("# Energy components:\n"); e.ener.print(); logPrintf("\n");
	}
	else if(e.cntrl.fixed_H)
	{	//Band structure calculation - ion and fluid minimization need to be handled differently
		if(eVars.nFilenamePattern.length())
		{	//If starting from density, compute potential:
			eVars.EdensityAndVscloc(e.ener);
			if(eVars.fluidSolver && eVars.fluidSolver->needsGummel())
			{	//Relies on the gummel loop, so EdensityAndVscloc would not have invoked minimize
				eVars.fluidSolver->minimizeFluid();
				eVars.EdensityAndVscloc(e.ener); //update Vscloc
			}
		}
		e.iInfo.augmentDensityGridGrad(eVars.Vscloc); //update Vscloc atom projections for ultrasoft psp's 
		if(e.cntrl.invertKS) //Inverse Kohn-Sham problem (sequence of band structure calculations)
		{	InverseKohnSham inverseKS(e);
			inverseKS.minimize(e.inverseKSminParams);
		}
		else //Single band structure calculation
		{	logPrintf("\n----------- Band structure minimization -------------\n"); logFlush();
			bandMinimize(e); // Do the band-structure minimization
		}
	}
	else if(e.vibrations) //Bypasses ionic/lattice minimization, calls electron/fluid minimization loops at various ionic configurations
	{	e.vibrations->calculate();
	}
	else if(e.latticeMinParams.nIterations)
	{	//Lattice minimization loop (which invokes the ionic minimization loop)
		LatticeMinimizer lmin(e);
		lmin.minimize(e.latticeMinParams);
	}
	else
	{	//Ionic minimization loop (which calls electron/fluid minimization loops)
		IonicMinimizer imin(e);
		imin.minimize(e.ionicMinParams);
	}

	//Final dump:
	e.dump(DumpFreq_End, 0);
	
	finalizeSystem();
	return 0;
}
