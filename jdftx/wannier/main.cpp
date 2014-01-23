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

extern Wannier wannier; //defined in wannier/commands.cpp so that it can be initialized from inputs

//Program entry point
int main(int argc, char** argv)
{	//Parse command line, initialize system and logs:
	string inputFilename; bool dryRun, printDefaults;
	initSystemCmdline(argc, argv, "Compute maximally-localized Wannier functions.", inputFilename, dryRun, printDefaults);

	//Parse input file and setup
	Everything e;
	parse(inputFilename.c_str(), e, printDefaults);
	e.setup();
	wannier.setup(e);
	Citations::print();
	if(dryRun)
	{	logPrintf("Dry run successful: commands are valid and initialization succeeded.\n");
		finalizeSystem();
		return 0;
	}
	
	wannier.saveMLWF();
	
	finalizeSystem();
	return 0;
}
