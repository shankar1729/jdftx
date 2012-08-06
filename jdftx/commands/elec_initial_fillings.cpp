/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

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

struct CommandElecInitialFillings : public Command
{
	CommandElecInitialFillings() : Command("elec-initial-fillings")
	{
		format = "read <filename> [<nBandsOld>]";
		comments =
			"Initial electronic fillings are read from <filename> instead of computed automatically (default)\n"
			"\t<nBandsOld> specifies the number of bands in the file being read, if different from current nBands\n"
			"\textra new bands are filled with 0s, fillings of extra old bands are accumulated into other bands\n"
			"\tdefault: 0 => fillings file has same number of bands as this run";
		
		require("elec-n-bands");
		forbid("initial-state");
	}

	void process(ParamList& pl, Everything& e)
	{	string key; pl.get(key, string(), "read", true);
		if(key!=string("read")) throw "First parameter must be 'read', encountered " + key;
		pl.get(e.eInfo.initialFillingsFilename, string(), "filename", true);
		pl.get(e.eInfo.nBandsOld, 0, "nBandsOld");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("read %s %d", e.eInfo.initialFillingsFilename.c_str(), e.eInfo.nBandsOld);
	}
}
commandElecInitialFillings;

