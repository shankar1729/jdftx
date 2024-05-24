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

struct CommandIon : public Command
{
	CommandIon() : Command("ion", "jdftx/Ionic/Geometry")
	{
		format = "<species-id> <x0> <x1> <x2> [v <vx0> <vx1> <vx2>] <moveScale> [<constraint type>="
			+ constraintTypeMap.optionList() + " <d0> <d1> <d2> [<group> [HyperPlane <d0> ...]]]";
		comments =
			"Add an atom of species <species-id> at coordinates (<x0>,<x1>,<x2>).\n"
			"\n"
			"Optionally, for dynamics, specify ion velocity <v0>,<v1>,<v2> after keyword 'v'.\n"
			"\n"
			"<moveScale> preconditions the motion of this ion (set 0 to hold fixed)\n"
			"\n"
			"In addition, the ion may be constrained to a line or a plane with line\n"
			"direction or plane normal equal to (<d0>,<d1>,<d2>) in the coordinate\n"
			"system selected by command coords-type. Note that the constraints must\n"
			"be consistent with respect to symmetries (if enabled).\n"
			"\n"
			"The HyperPlane constraint allows constraining collective motion of many\n"
			"ions by restricting their motion to a hyperplane with normal specified\n"
			"by (<d0>,<d1>,<d2>) for all ions specifying a hyperplane constraint.\n"
			"By default, all hyperplane-constrained ions are included in a single\n"
			"group; use optional <group> label to specify multiple hyper-planes.\n"
			"Multiple hyperplane constraints may also be added to each atom,\n"
			"but this requires an explicit group label for each hyperplane.\n"
			"\n"
			"Note that when coords-type is lattice, the constraints are in covariant\n"
			"lattice coordinates (like direction indices) for line constraints, but\n"
			"contravariant coordinates (like plane indices) for plane constraints.";
		
		allowMultiple = true;

		require("ion-species");
		//Dependencies due to coordinate system option:
		require("latt-scale");
		require("coords-type");
	}

	void process(ParamList& pl, Everything& e)
	{	//Find species:
		string id; pl.get(id, string(), "species-id", true);
		auto sp = findSpecies(id, e, true);
		if(!sp) throw string("Species "+id+" has not been defined");
		//Read coordinates:
		vector3<> pos;
		for(int k=0; k<3; k++)
		{	ostringstream oss; oss << "x" << k;
			pl.get(pos[k], 0., oss.str(), true);
		}
		//Transform coordinates if necessary
		if(e.iInfo.coordsType == CoordsCartesian)
			pos = inv(e.gInfo.R) * pos;
		//Add position to list:
		sp->atpos.push_back(pos);
		sp->velocities.push_back(vector3<>(NAN,NAN,NAN));
		
		//Look for optional velocity and get moveScale:
		string key;
		pl.get(key, string(), "moveScale|v", true);
		SpeciesInfo::Constraint constraint;
		if(key == "v")
		{	//Read velocity (update the NAN added above):
			vector3<>& vel = sp->velocities.back();
			for(int k=0; k<3; k++)
			{	ostringstream oss; oss << "v" << k;
				pl.get(vel[k], 0., oss.str(), true);
			}
			//Transform coordinates if necessary
			if(e.iInfo.coordsType == CoordsCartesian)
				vel = inv(e.gInfo.R) * vel;
			//Get moveScale from command line beyond the velocity:
			pl.get(constraint.moveScale, 0., "moveScale", true);
		}
		else
		{	//No velocity: get moveScale from key
			ParamList(key+" ").get(constraint.moveScale, 0., "moveScale", true);
		}
		
		//Parse constraints
		if(constraint.moveScale < 0.)   // Check for negative moveScales
			throw string("moveScale cannot be negative");
		pl.get(constraint.type, SpeciesInfo::Constraint::None, constraintTypeMap, "Type");
		if(constraint.type != SpeciesInfo::Constraint::None)
		{	if(!constraint.moveScale)
				throw string("Constraint specified after moveScale = 0");
			pl.get(constraint.d[0], 0.0, "d0", true);				  
			pl.get(constraint.d[1], 0.0, "d1", true);
			pl.get(constraint.d[2], 0.0, "d2", true);
			if(e.iInfo.coordsType == CoordsLattice) //Transform to Cartesian (taking care of covariant/contravariant for line/plane directions)
				switch(constraint.type)
				{	case SpeciesInfo::Constraint::Linear:       constraint.d = e.gInfo.R * constraint.d; break;
					case SpeciesInfo::Constraint::Planar:       constraint.d = ~inv(e.gInfo.R) * constraint.d; break;
					case SpeciesInfo::Constraint::HyperPlane:   constraint.d = ~inv(e.gInfo.R) * constraint.d; break;
					default: break;
				}
			if(constraint.type == SpeciesInfo::Constraint::HyperPlane)
			{	string groupLabel;
				pl.get(groupLabel, string(), "group"); //optional group label for hyperplane constraint
				constraint.hyperplane.assign(1, std::make_pair(constraint.d, groupLabel));
				constraint.d = vector3<>(); //zero out the constraint d to prevent incorrect symmetry checks
				//Look for additional hyperplane constraints:
				while(true)
				{	string key; pl.get(key, string(), "Type");
					if(not key.length()) break; //End of constraint list
					if(key != "HyperPlane") throw string("Additional constraints must be of type HyperPlane");
					vector3<> d;
					pl.get(d[0], 0.0, "d0", true);				  
					pl.get(d[1], 0.0, "d1", true);
					pl.get(d[2], 0.0, "d2", true);
					pl.get(groupLabel, string(), "group", true); //group label required for additional hyperplanes
					if(e.iInfo.coordsType == CoordsLattice) d = ~inv(e.gInfo.R) * d; //transform to Cartesian
					constraint.hyperplane.push_back(std::make_pair(d, groupLabel));
				}
			}
			else //Note d can be null for a guiven atom with HyperPlane
			{	if(not constraint.d.length_squared())
					throw string("Constraint vector must be non-null");
			}
		}
		sp->constraints.push_back(constraint);
	}

	void printStatus(Everything& e, int iRep)
	{	int iIon=0;
		for(auto sp: e.iInfo.species)
			for(unsigned at=0; at<sp->atpos.size(); at++)
			{	if(iIon==iRep)
				{	//Species and position:
					vector3<> pos = sp->atpos[at];
					if(e.iInfo.coordsType == CoordsCartesian)
						pos = e.gInfo.R * pos; //report cartesian positions
					logPrintf("%s %19.15lf %19.15lf %19.15lf", sp->name.c_str(), pos[0], pos[1], pos[2]);
					//Optional velocity:
					vector3<> vel = sp->velocities[at];
					if(not std::isnan(vel.length_squared()))
					{	if(e.iInfo.coordsType == CoordsCartesian)
							vel = e.gInfo.R * vel; //report cartesian velocities
						logPrintf(" v %19.15lf %19.15lf %19.15lf", vel[0], vel[1], vel[2]);
					}
					logPrintf(" %lg", sp->constraints[at].moveScale);
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
	CommandInitialMagneticMoments() : Command("initial-magnetic-moments", "jdftx/Initialization")
	{
		format = "<species> <M1> <M2> ... <Mn> [<species2> ...]\n"
			"      | <species> <M1> <theta1> <phi1> ... <Mn> <thetan> <phin> [<species2> ...]";
		comments =
			"Specify initial magnetic moments, defined as the difference between\n"
			"up and down electron counts, on each atom of one or more species.\n"
			"\n"
			"The second syntax with polar angles (in degrees) must be used\n"
			"for noncollinear  magnetism calculations.\n"
			"\n"
			"For each species, the initial magnetic moments are applied\n"
			"to the atoms in the order of ion commands for that species.\n"
			"This may be used to construct a spin-polarized reference density\n"
			"for LCAO initialization of the Kohn-Sham orbitals.";
		
		require("ion");
		require("spintype");
	}

	void process(ParamList& pl, Everything& e)
	{	if(e.eInfo.spinType==SpinNone || e.eInfo.spinType==SpinOrbit)
			throw string("Cannot specify magnetic moments in an unpolarized calculation");
		string id;
		pl.get(id, string(), "species", true);
		while(id.length())
		{	auto sp = findSpecies(id, e);
			if(!sp) throw string("Species "+id+" has not been defined");
			sp->initialMagneticMoments.resize(sp->atpos.size());
			for(unsigned a=0; a<sp->atpos.size(); a++)
			{	double M=0., theta=0., phi=0.;
				pl.get(M, 0., "M", true);
				if(e.eInfo.spinType == SpinVector)
				{	pl.get(theta, 0., "theta", true);
					pl.get(phi, 0., "phi", true);
				}
				//Store as Cartesian:
				sp->initialMagneticMoments[a] = M * polarUnitVector(phi*M_PI/180, theta*M_PI/180);
			}
			//Check for additional species:
			pl.get(id, string(), "species");
		}
	}

	void printStatus(Everything& e, int iRep)
	{	for(auto sp: e.iInfo.species)
			if(sp->initialMagneticMoments.size())
			{	logPrintf(" \\\n\t%s", sp->name.c_str());
				for(const vector3<>& M: sp->initialMagneticMoments)
				{	if(e.eInfo.spinType == SpinVector)
					{	vector3<> euler;
						if(M.length()) getEulerAxis(M, euler);
						euler *= 180./M_PI; //convert to degrees
						logPrintf(" %lg %lg %lg ", M.length(), euler[1], euler[0]);
					}
					else logPrintf(" %lg", M[2]); //SpinZ
				}
			}
	}
}
commandInitialMagneticMoments;

//-------------------------------------------------------------------------------------------------

struct CommandInitialOxidationState : public Command
{
	CommandInitialOxidationState() : Command("initial-oxidation-state", "jdftx/Initialization")
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
		{	auto sp = findSpecies(id, e);
			if(!sp) throw string("Species "+id+" has not been defined");
			pl.get(sp->initialOxidationState, 0., "oxState", true);
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

EnumStringMap<coreOverlapCheck> overlapCheckMap
(	additive, "additive",
	vector, "vector",
	none, "none"
);


struct CommandCoreOverlapCheck : public Command
{
	CommandCoreOverlapCheck() : Command("core-overlap-check", "jdftx/Ionic/Optimization")
	{
		format = "<condition>";
		comments = "Checks for core overlaps between ionic pseudopotentials based on <condition>:\n"
				   "+ additive: checks for interatomic distance < (R1 + R2)\n"
				   "+ vector: checks for interatomic distance < sqrt(R1^2 + R2^2) (default)\n"
				   "+ none";
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
