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
#include <phonon/Phonon.h>

//Program entry point
int main(int argc, char** argv)
{	//Parse command line, initialize system and logs:
	string inputFilename; bool dryRun, printDefaults;
	initSystemCmdline(argc, argv, "Compute maximally-localized Phonon functions.", inputFilename, dryRun, printDefaults);

	Phonon phonon;
	phonon.input = readInputFile(inputFilename);
	
	//Read input file and setup unit cell:
	phonon.setup(printDefaults);
	Citations::print();
	if(dryRun)
	{	logPrintf("Dry run successful: commands are valid and initialization succeeded.\n");
		finalizeSystem();
		return 0;
	}
	
	//Calculate:
	phonon.dump();
	
	finalizeSystem();
	return 0;
}
