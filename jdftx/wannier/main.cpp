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
#include <core/matrix.h>
#include <core/Util.h>
#include <commands/parser.h>
#include <wannier/Wannier.h>

//Program entry point
int main(int argc, char** argv)
{	//Parse command line, initialize system and logs:
	InitParams ip("Compute maximally-localized Wannier functions.");
	initSystemCmdline(argc, argv, ip);

	//Parse input file:
	WannierEverything e;
	parse(readInputFile(ip.inputFilename), e, ip.printDefaults);
	
	//Set initial filenames and prevent unnecessary setup below:
	e.eVars.wfnsFilename = e.wannier.getFilename(Wannier::FilenameInit, "wfns");
	e.eVars.eigsFilename = e.wannier.eigsFilename.length()
		?	e.wannier.eigsFilename
		:	e.wannier.getFilename(Wannier::FilenameInit, "eigenvals");
	e.eVars.nFilenamePattern.clear();
	e.eVars.VFilenamePattern.clear();
	e.eVars.fluidParams.fluidType = FluidNone;
	if(ip.dryRun) e.eVars.skipWfnsInit = true;
	//--- fillings (if available)
	e.eInfo.initialFillingsFilename = e.wannier.getFilename(Wannier::FilenameInit, "fillings");
	if(fileSize(e.eInfo.initialFillingsFilename.c_str()) <= 0)
		e.eInfo.initialFillingsFilename.clear();
	
	//Setup:
	e.setup();
	Citations::print();
	if(ip.dryRun)
	{	logPrintf("Dry run successful: commands are valid and initialization succeeded.\n");
		finalizeSystem();
		return 0;
	}
	
	e.wannier.saveMLWF();
	
	finalizeSystem();
	return 0;
}
