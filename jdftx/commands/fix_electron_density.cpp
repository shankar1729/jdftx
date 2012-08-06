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

struct CommandFixElectronDensity : public Command
{
	CommandFixElectronDensity() : Command("fix-electron-density")
	{
		format = "<nFilename> | <nupFilename> <ndnFilename>";
		comments = "Perform band structure calculations at fixed electron density\n"
			"(or spin densities) read from the specified file(s)";
		
		require("spintype");
		forbid("elec-fermi-fillings");
	}

	void process(ParamList& pl, Everything& e)
	{	std::vector<string>& nFilename = e.eVars.nFilename;
		nFilename.resize(e.eInfo.spinType==SpinNone ? 1 : 2);
		for(unsigned s=0; s<nFilename.size(); s++)
		{	const char* componentName = nFilename.size()==1 ? "nFilename" : (s==0 ? "nupFilename" : "ndnFilename");
			pl.get(nFilename[s], string(), componentName, true);
		}
		e.cntrl.fixed_n = true;
	}

	void printStatus(Everything& e, int iRep)
	{	std::vector<string>& nFilename = e.eVars.nFilename;
		for(unsigned s=0; s<nFilename.size(); s++)
			logPrintf("%s ", e.eVars.nFilename[s].c_str());
	}
}
commandFixElectronDensity;
