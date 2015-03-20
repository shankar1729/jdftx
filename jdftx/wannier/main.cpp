/*-------------------------------------------------------------------
Copyright 2014 Ravishankar Sundararaman

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
#include <core/Util.h>
#include <commands/parser.h>
#include <wannier/Wannier.h>

//Program entry point
int main(int argc, char** argv)
{	//Parse command line, initialize system and logs:
	string inputFilename; bool dryRun, printDefaults;
	initSystemCmdline(argc, argv, "Compute maximally-localized Wannier functions.", inputFilename, dryRun, printDefaults);

	//Parse input file:
	WannierEverything e;
	parse(readInputFile(inputFilename), e, printDefaults);
	
	//Set initial filenames and prevent unnecessary setup below:
	e.eVars.wfnsFilename = e.wannier.getFilename(Wannier::FilenameInit, "wfns");
	e.eVars.eigsFilename = e.wannier.getFilename(Wannier::FilenameInit, "eigenvals");
	e.eVars.nFilenamePattern.clear();
	e.eVars.VFilenamePattern.clear();
	e.eVars.HauxFilename.clear();
	e.eVars.fluidParams.fluidType = FluidNone;
	
	//Setup:
	e.setup();
	Citations::print();
	if(dryRun)
	{	logPrintf("Dry run successful: commands are valid and initialization succeeded.\n");
		finalizeSystem();
		return 0;
	}
	
	if(e.wannier.saveMomenta && e.eInfo.hasU)
	{	//Calculate U_rho needed for the DFT+U correction to the [r,H] momentum matrix elements:
		e.iInfo.rhoAtom_initZero(e.eVars.rhoAtom);
		e.iInfo.rhoAtom_initZero(e.eVars.U_rhoAtom);
		e.iInfo.rhoAtom_calc(e.eVars.F, e.eVars.C, e.eVars.rhoAtom);
		e.iInfo.rhoAtom_computeU(e.eVars.rhoAtom, e.eVars.U_rhoAtom);
	}
	
	e.wannier.saveMLWF();
	
	finalizeSystem();
	return 0;
}
