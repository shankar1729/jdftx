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

struct CommandForcesOutputCoords : public Command
{
	CommandForcesOutputCoords() : Command("forces-output-coords")
	{
		format = "<coords>=" + forcesOutputCoordsMap.optionList();
		comments =
			"Coordinate system to use for force output in logPrintf-file as well as dump:\n"
			"\tPositions: Use the same coordinate system as ionic position input (selected by coords-type) [default].\n"
			"\tLattice:   Use (covariant) lattice coordinates\n"
			"\tCartesian: Use cartesian coordinates\n"
			"\tContravariant: Use contravariant lattice coordinates (covariant multiplied by inv(RT.R))";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.iInfo.forcesOutputCoords, ForcesCoordsPositions, forcesOutputCoordsMap, "coords");
	}

	void printStatus(Everything& e, int iRep)
	{	fputs(forcesOutputCoordsMap.getString(e.iInfo.forcesOutputCoords), globalLog);
	}
}
commandForcesOutputCoords;
