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

struct CommandIon : public Command
{
	CommandIon() : Command("ion")
	{
		format = "<species-id> <x0> <x1> <x2> <moveScale>";
		comments =
			"Add an atom of species <species-id> at coordinates (<x0>,<x1>,<x2>).\n"
			"<moveScale> preconditions the motion of this ion (set 0 to hold fixed)";
		allowMultiple = true;

		require("ion-species");
		require("lattice");
		require("latt-scale");
		require("coords-type");
	}

	void process(ParamList& pl, Everything& e)
	{	string id;
		pl.get(id, string(), "species-id", true);
		for(auto sp: e.iInfo.species)
			if(sp->name == id)
			{	vector3<> pos;
				for(int k=0; k<3; k++)
				{	ostringstream oss; oss << "x" << k;
					pl.get(pos[k], 0.0, oss.str(), true);
				}
				double moveScale;
				pl.get(moveScale, 0.0, "moveScale", true);
				//Transform coordinates if necessary
				if(e.iInfo.coordsType == CoordsCartesian)
					pos = inv(e.gInfo.R) * pos;
				//Add position and move scale factor to list:
				sp->atpos.push_back(pos);
				sp->moveScale.push_back(moveScale);
				return;
			}
		throw string("Species "+id+" has not been defined");
	}

	void printStatus(Everything& e, int iRep)
	{	int iIon=0;
		for(auto sp: e.iInfo.species)
			for(unsigned at=0; at<sp->atpos.size(); at++)
			{	if(iIon==iRep)
				{	vector3<> pos = sp->atpos[at];
					if(e.iInfo.coordsType == CoordsCartesian)
						pos = e.gInfo.R * pos; //report cartesian positions
					logPrintf("%s %18.14lf %18.14lf %18.14lf %lg", sp->name.c_str(),
						pos[0], pos[1], pos[2], sp->moveScale[at]);
				}
				iIon++;
			}
	}
}
commandIon;
