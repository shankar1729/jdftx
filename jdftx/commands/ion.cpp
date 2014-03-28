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
		format = "<species-id> <x0> <x1> <x2> <moveScale> [<constraint type>="
			+ constraintTypeMap.optionList() + " <d0> <d1> <d2>]";
		comments =
			"Add an atom of species <species-id> at coordinates (<x0>,<x1>,<x2>).\n"
			"<moveScale> preconditions the motion of this ion (set 0 to hold fixed)\n"
			"In addition, the ion may be constrained to a line or a plane with line\n"
			"direction or plane normal equal to (<d0>,<d1>,<d2>) in the coordinate\n"
			"system selected by command coords-type. Note that the constraints must\n"
			"be consistent with respect to symmetries (if enabled).";
		allowMultiple = true;

		require("ion-species");
		//Dependencies due to coordinate system option:
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
				if(constraint.moveScale < 0.)   // Check for negative moveScales
					throw string("moveScale cannot be negative");
				pl.get(constraint.type, SpeciesInfo::Constraint::None, constraintTypeMap, "Type");
				if(constraint.type != SpeciesInfo::Constraint::None)
				{	if(!constraint.moveScale)
						throw string("Constraint specified after moveScale = 0");
					pl.get(constraint.d[0], 0.0, "d0", true);				  
					pl.get(constraint.d[1], 0.0, "d1", true);
					pl.get(constraint.d[2], 0.0, "d2", true);
					if(not constraint.d.length_squared())
						throw string("Constraint vector must be non-null");
					if(e.iInfo.coordsType == CoordsLattice) //Constraints transform like forces:
						constraint.d = ~inv(e.gInfo.R) * constraint.d; 
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
						sp->constraints[at].print(globalLog, e);
				}
				iIon++;
			}
	}
}
commandIon;

//-------------------------------------------------------------------------------------------------

struct CommandInitialMagneticMoments : public Command
{
	CommandInitialMagneticMoments() : Command("initial-magnetic-moments")
	{
		format = "<species> <M1> <M2> ... <Mn> [<species2> ...]";
		comments =
			"Specify initial magnetic moments, defined as the difference between\n"
			"up and down electron counts, on each atom of one or more species.\n"
			"For each species, the initial magnetic moments are applied to the\n"
			"atoms in the order of ion commands for that species.\n"
			"This may be used to construct a spin-polarized reference density\n"
			"for LCAO initialization of the Kohn-Sham orbitals.";
		
		require("ion");
		require("spintype");
	}

	void process(ParamList& pl, Everything& e)
	{	if(e.eInfo.spinType == SpinNone)
			throw string("Cannot specify magnetic moments in an unpolarized calculation");
		string id;
		pl.get(id, string(), "species", true);
		while(id.length())
		{	bool spFound = false;
			for(auto sp: e.iInfo.species)
				if(sp->name == id)
				{	sp->initialMagneticMoments.resize(sp->atpos.size());
					for(unsigned a=0; a<sp->atpos.size(); a++)
					{	ostringstream oss; oss << "M" << (a+1);
						pl.get(sp->initialMagneticMoments[a], 0., oss.str(), true);
					}
					spFound = true;
					break;
				}
			if(!spFound) throw string("Species "+id+" has not been defined");
			//Check for additional species:
			pl.get(id, string(), "species");
		}
	}

	void printStatus(Everything& e, int iRep)
	{	for(auto sp: e.iInfo.species)
			if(sp->initialMagneticMoments.size())
			{	logPrintf(" \\\n\t%s", sp->name.c_str());
				for(double M: sp->initialMagneticMoments) logPrintf(" %lg", M);
			}
	}
}
commandInitialMagneticMoments;

//-------------------------------------------------------------------------------------------------

struct CommandInitialOxidationState : public Command
{
	CommandInitialOxidationState() : Command("initial-oxidation-state")
	{
		format = "<species> <oxState> [<species2> ...]";
		comments =
			"Specify initial oxidation state assumed for each species in LCAO.\n"
			"This may significantly improve LCAO convergence in some cases.";
		
		require("ion-species");
	}

	void process(ParamList& pl, Everything& e)
	{	string id;
		pl.get(id, string(), "species", true);
		while(id.length())
		{	bool spFound = false;
			for(auto sp: e.iInfo.species)
				if(sp->name == id)
				{	pl.get(sp->initialOxidationState, 0., "oxState", true);
					spFound = true;
					break;
				}
			if(!spFound) throw string("Species "+id+" has not been defined");
			//Check for additional species:
			pl.get(id, string(), "species");
		}
	}

	void printStatus(Everything& e, int iRep)
	{	for(auto sp: e.iInfo.species)
			if(sp->initialOxidationState)
				logPrintf(" \\\n\t%s %lg", sp->name.c_str(), sp->initialOxidationState);
	}
}
commandInitialOxidationState;

