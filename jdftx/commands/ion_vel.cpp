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
#include <fluid/Euler.h>

struct CommandIonVelocity : public Command
{
	CommandIonVelocity() : Command("ion-vel", "jdftx/Ionic/Dynamics")
	{
		format = "<species-id> <x0> <x1> <x2> <v0> <v1> <v2>";
		comments =
			"Add an atom of species <species-id> at coordinates (<x0>,<x1>,<x2>)\n"
			"with an initial velocity of (<v0>,<v1>,<v2)>\n"
			"\n";
		allowMultiple = true;

		require("ion-species");
		//Dependencies due to coordinate system option:
		require("latt-scale");
		require("coords-type");
		require("ionic-dynamics");
		forbid("ion");
	}

	void process(ParamList& pl, Everything& e)
	{	//Find species:
		string id; pl.get(id, string(), "species-id", true);
		auto sp = findSpecies(id, e);
		if(!sp) throw string("Species "+id+" has not been defined");
		//Read coordinates:
		vector3<> pos,vel;
		for(int k=0; k<3; k++)
		{	ostringstream oss; oss << "x" << k;
			pl.get(pos[k], 0.0, oss.str(), true);
		}
		bool velocitiesGiven = true;
		for(int k=0; k<3; k++)
		{	ostringstream oss; oss << "v" << k;
			pl.get(vel[k], (double)NAN, oss.str(), false);
			if (std::isnan(vel[k]))
			      velocitiesGiven = (velocitiesGiven && false); // if any of the velocity components is missing, then don't push_back
		}
		//Transform coordinates if necessary
		if(e.iInfo.coordsType == CoordsCartesian)
		{	pos = inv(e.gInfo.R) * pos;
			vel = inv(e.gInfo.R) * vel;
		}
		//Add position to list:
		sp->atpos.push_back(pos);
		if (velocitiesGiven)
			sp->velocities.push_back(vel);
		
		//No need for constraints
		SpeciesInfo::Constraint constraint;
		constraint.moveScale = 1.0;
		constraint.type = SpeciesInfo::Constraint::None;
		sp->constraints.push_back(constraint);
	}

	void printStatus(Everything& e, int iRep)
	{	int iIon=0;
		for(auto sp: e.iInfo.species)
			for(unsigned at=0; at<sp->atpos.size(); at++)
			{	if(iIon==iRep)
				{	vector3<> pos = sp->atpos[at];
					vector3<> vel = sp->velocities[at];
					if(e.iInfo.coordsType == CoordsCartesian)
					{	pos = e.gInfo.R * pos; //report cartesian positions
						vel = e.gInfo.R * vel;
					}
					logPrintf("%s %19.15lf %19.15lf %19.15lf %19.15lf %19.15lf %19.15lf", sp->name.c_str(),
						pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]);
				}
				iIon++;
			}
	}
}
commandIonVelocity;
