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
#include <electronic/SpeciesInfo.h>

struct CommandChargeball : public Command
{
	CommandChargeball() : Command("chargeball")
	{
		format = "<species-id> <unusedInteger> <width> [<norm>=2.0]";
		comments =
			"Gaussian chargeball of norm <norm> and width <width> for species <id>\n"
			"(<unusedInteger> used to be the unused <Z_core>, and has been retained for compatibility)\n"
			"This feature is deprecated, use a pseudopotential with partial core correction instead.";
		allowMultiple = true;

		require("ion-species");
		require("lattice");
		require("coords-type");
	}

	void process(ParamList& pl, Everything& e)
	{	string id;
		pl.get(id, string(), "species-id", true);
		for(auto sp: e.iInfo.species)
			if(sp->name == id)
			{	int unusedInteger;
				pl.get(unusedInteger, 0, "unusedInteger", true);
				pl.get(sp->width_chargeball, 0.0, "width", true);
				pl.get(sp->Z_chargeball, 2.0, "norm");
				return;
			}
		throw string("Species "+id+" has not been defined");
	}

	void printStatus(Everything& e, int iRep)
	{	if(unsigned(iRep) < e.iInfo.species.size())
		{	const SpeciesInfo& sp = *(e.iInfo.species[iRep]);
			logPrintf("%s -1 %lg %lg", sp.name.c_str(), sp.width_chargeball, sp.Z_chargeball);
		}
	}
}
commandChargeball;
