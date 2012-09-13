/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Deniz Gunceler

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
		format = "<species-id> <x0> <x1> <x2> <moveScale> [<constraint type> <x0> <x1> <x2>]";
		comments =
			"Add an atom of species <species-id> at coordinates (<x0>,<x1>,<x2>).\n"
			"<moveScale> preconditions the motion of this ion (set 0 to hold fixed)\n"
			"\nYou can specify a constraint to the direction which the ions can relax in.\n"
			"There are 2 types of constraints, planar and linear.\nThese options constrain the motion "
			"of the ions to the input direction, or the input plane specified by its normal.  "
			"\nIn both cases, the coordinates are the same as the ionic coordinates.\n"
			"If the symmetries are turned on, consistency of these constraints wrt to symmetries will be checked.";
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
				//Transform coordinates if necessary
				if(e.iInfo.coordsType == CoordsCartesian)
					pos = inv(e.gInfo.R) * pos;
				//Add position and move scale factor to list:
				sp->atpos.push_back(pos);
				
				//Look for constraints
				SpeciesInfo::Constraint constraint;
				pl.get(constraint.moveScale, 0.0, "moveScale", true);
				if(constraint.moveScale < 0.) throw string("You can't have a negative moveScale for an ion!")
				pl.get(constraint.type, SpeciesInfo::Constraint::None, constraintTypeMap, "Type");
				if(constraint.type != SpeciesInfo::Constraint::None)
				{	if(!constraint.moveScale)
						throw string("Constraint specified after moveScale = 0");
					pl.get(constraint.x[0], 0.0, "x0", true);				  
					pl.get(constraint.x[1], 0.0, "x1", true);
					pl.get(constraint.x[2], 0.0, "x2", true);
					if(not constraint.x.length_squared())
						throw string("Constraint vector must be non-null");
					if(e.iInfo.coordsType == CoordsLattice)
						constraint.x = ~inv(e.gInfo.R) * constraint.x;
				}
				sp->constraints.push_back(constraint);
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
					logPrintf("%s %19.15lf %19.15lf %19.15lf %lg", sp->name.c_str(),
						pos[0], pos[1], pos[2], sp->constraints[at].moveScale);
					if(sp->constraints[at].type != SpeciesInfo::Constraint::None)
						sp->constraints[at].print(globalLog);
				}
				iIon++;
			}
	}
}
commandIon;
