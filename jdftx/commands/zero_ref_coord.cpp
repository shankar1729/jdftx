/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver

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
/*
struct CommandZeroRefCoord : public Command
{
    CommandZeroRefCoord() : Command("zero-ref-coord")
	{
		format = "<x0> <x1> <x2>";
		comments = "Specify zero reference of potential (in coordinates specified by coords-type)\n"
				   "Default zero reference at center of cell.";
		
		//Dependencies due to coordinate system option:
		require("latt-scale");
		require("coords-type");
	}

	void process(ParamList& pl, Everything& e)
	{
		vector3<> pos;
		for(int k=0; k<3; k++)
				{	ostringstream oss; oss << "x" << k;
					pl.get(pos[k], 0.0, oss.str(), true);
				}
		if(e.iInfo.coordsType == CoordsCartesian) pos = inv(e.gInfo.R)*pos;
		e.eVars.fluidParams.zeroRef = pos;
	}

	void printStatus(Everything& e, int iRep)
	{
			for (int k=0; k < 3; k++) logPrintf("%20.15le ",e.eVars.fluidParams.zeroRef[k]);
	}
}
CommandZeroRefCoord;
*/