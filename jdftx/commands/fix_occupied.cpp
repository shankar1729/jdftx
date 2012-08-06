/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

#include <commands/command.h>
#include <electronic/Everything.h>

struct CommandFixOccupied : public Command
{
	CommandFixOccupied() : Command("fix-occupied")
	{
		format = "[<fThreshold>=0]";
		comments = "Fix orbitals with fillings larger than <fThreshold> in band-structure calculations\n"
			"The occupied orbitals must therefore be read in using the wavefunction command.\n";
		
		require("wavefunction");
		require("fix-electron-density");
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.cntrl.occupiedThreshold, 0., "fThreshold");
		if(e.cntrl.occupiedThreshold<0) throw string("fThreshold must be >= 0");
		if(!e.eVars.wfnsFilename.length()) throw string("wavefunctions must be read in from file for fixed occupied states.");
		e.cntrl.fixOccupied = true;
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg", e.cntrl.occupiedThreshold);
	}
}
commandFixOccupied;
