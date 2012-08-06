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

EnumStringMap<CoordsType> coordsMap(
	CoordsLattice, "lattice",
	CoordsCartesian, "cartesian" );

struct CommandCoordsType : public Command
{
	CommandCoordsType() : Command("coords-type")
	{
		format = "<coords>=" + coordsMap.optionList();
		comments = "Coordinate system used in specifying ion positions (default: lattice)";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.iInfo.coordsType, CoordsLattice, coordsMap, "coords");
	}

	void printStatus(Everything& e, int iRep)
	{	fputs(coordsMap.getString(e.iInfo.coordsType), globalLog);
	}
}
commandCoordsType;
