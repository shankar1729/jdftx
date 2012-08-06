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

struct CommandTargetMu : public Command
{
	CommandTargetMu() : Command("target-mu")
	{
		format = "<mu> [<Cinitial>=1.0] [<dnMix>=0.7]";
		comments =
			"Fixed chemical potential <mu> (instead of fixed charge)\n"
			"\tCinitial: Initial capacitance, affects only first few mixing steps\n"
			"\tdnMix: Scale the ideal step in n by this factor";
		hasDefault = false;

		require("fluid");
		require("ionic-screening");
		require("elec-fermi-fillings");
		forbid("elec-initial-charge");
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.eInfo.mu, 0.0, "mu", true);
		pl.get(e.eInfo.Cmeasured, 1.0, "Cinitial");
		pl.get(e.eInfo.dnMixFraction, 0.7, "dnMix");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg %lg %lg\n", e.eInfo.mu, e.eInfo.Cmeasured, e.eInfo.dnMixFraction);
	}
}
commandTargetMu;
