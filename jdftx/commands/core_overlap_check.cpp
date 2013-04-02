/*-------------------------------------------------------------------
Copyright 2011 Deniz Gunceler

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
#include <core/Units.h>

EnumStringMap<coreOverlapCheck> overlapCheckMap
(	additive, "additive",
	vector, "vector",
	none, "none"
);

struct CommandCoreOverlapCheck : public Command
{
	CommandCoreOverlapCheck() : Command("core-overlap-check")
	{
		format = "<condition> = vectorial";
		comments = "Checks for core overlaps between ionic pseudopotentials.\n"
				   "\tadditive -> checks for interatomic distance < (R1 + R2)\n"
				   "\tvector   -> checks for interatomic distance < sqrt(R1^2 + R2^2)   (default)\n"
				   "\tnone\n";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.iInfo.coreOverlapCondition, vector, overlapCheckMap, "overlap check");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s", overlapCheckMap.getString(e.iInfo.coreOverlapCondition));
	}
}
commandCoreOverlapCheck;
