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
#include <electronic/DOS.h>

EnumStringMap<DOS::Weight::Type> weightTypeMap
(	DOS::Weight::Total, "Total",
	DOS::Weight::Slice, "Slice",
	DOS::Weight::Sphere, "Sphere",
	DOS::Weight::AtomSlice, "AtomSlice",
	DOS::Weight::AtomSphere, "AtomSphere",
	DOS::Weight::File, "File"
);

struct CommandDensityOfStates : public Command
{
	CommandDensityOfStates() : Command("density-of-states")
	{
		format = "[<occupied>=" + boolMap.optionList() + "] [<Etol>=1e-6] [<weightType1> ...] [<weightType2> ...] ...";
		comments =
			"Compute density of states. The results are printed to a text file\n"
			"with name corresponding to variable name 'dos' (see dump-name).\n"
			"(Spin polarized calculations output variables 'dosUp' and 'dosDn'.)\n"
			"Density of states with different weight functions may be computed\n"
			"simultaneously, and they are all output as columns in the same file\n"
			"in the same order that they appear in this command, with the energy\n"
			"in the first column. The energy is in Hartrees, and the density of\n"
			"states is in electrons/UnitCell/Hartree.\n"
			"   Compute occupied density of states if <occupied> = yes, and total\n"
			"density of states otherwise. <Etol> specifies the resolution in energy\n"
			"within which eigenvalues are identified, and is used as the band width\n"
			"for Gamma-point only calculations. The remaining arguments list the\n"
			"weight functions, and each function is specified by a mode followed\n"
			"by arguments specific to that mode. The available modes and their\n"
			"arguments are:\n"
			"   Total\n"
			"      Compute the total density of states (no arguments)\n"
			"   Slice  <c0> <c1> <c2>   <r>   <i0> <i1> <i2>\n"
			"      Density of states in a planar slab centered at (<c0>,<c1>,<c2>)\n"
			"      in the coordinate system selected by coords-type, parallel to\n"
			"      the lattice plane specified by Miller indices (<i0>,<i1>,<i2>),\n"
			"      with half-width <r> bohrs normal to the lattice plane.\n"
			"   Sphere  <c0> <c1> <c2>   <r>\n"
			"      Density of states in a sphere of radius <r> bohrs centered at\n"
			"      (<c0>,<c1>,<c2>) in the coordinate system selected by coords-type.\n"
			"   AtomSlice  <species> <atomIndex>   <r>   <i0> <i1> <i2>\n"
			"      Like Slice mode, with center located at atom number <atomIndex>\n"
			"      (1-based index, in input file order) of species name <species>.\n"
			"   AtomSphere  <species> <atomIndex>   <r>\n"
			"      Like Sphere mode, with center located at atom number <atomIndex>\n"
			"      (1-based index, in input file order) of species name <species>.\n"
			"   File <filename>\n"
			"      Arbitrary real-space weight function read from file <filename>.\n"
			"      (double-precision binary, same format as electron density output)\n"
			"      A file with all 1.0's would yield the same result as mode Total.\n"
			"Any number of weight functions may be specified; only the total density\n"
			"of states is output if no weight functions are specified. By default,\n"
			"dumping occurs only at end of calculation, but this may be altered by\n"
			"explicitly including DOS in a dump command of appropriate frequency.\n"
			"(see command dump).";
		hasDefault = false;
		
		require("ion"); //This ensures that this command is processed after all ion commands 
		// (which in turn are processed after lattice and all ion-species commands)
	}

	void process(ParamList& pl, Everything& e)
	{	e.dump.dos = std::make_shared<DOS>();
		DOS& dos = *(e.dump.dos);
		pl.get(dos.occupied, false, boolMap, "occupied");
		pl.get(dos.Etol, 1e-6, "Etol");
		//Check for weight functions:
		while(true)
		{	DOS::Weight weight;
			pl.get(weight.type, DOS::Weight::Delim, weightTypeMap, "weightType");
			if(weight.type==DOS::Weight::Delim) break; //end of input
			//Get center coordinates for Slice or Sphere mode:
			if(weight.type==DOS::Weight::Slice || weight.type==DOS::Weight::Sphere)
			{	vector3<> center;
				pl.get(center[0], 0., "c0", true);
				pl.get(center[1], 0., "c1", true);
				pl.get(center[2], 0., "c2", true);
				weight.center = (e.iInfo.coordsType==CoordsCartesian) ? inv(e.gInfo.R) * center : center; //internally store in lattice coordinates
			}
			//Get species and atom index for AtomSlice or AtomSphere mode
			if(weight.type==DOS::Weight::AtomSlice || weight.type==DOS::Weight::AtomSphere)
			{	//Find specie:
				string spName; pl.get(spName, string(), "species", true);
				bool spFound = false;
				for(size_t sp=0; sp<e.iInfo.species.size(); sp++)
					if(e.iInfo.species[sp]->name == spName)
					{	weight.specieIndex = sp;
						spFound = true;
						break;
					}
				if(!spFound)
					throw "Could not find species with name '" + spName + "'";
				//Get atom index:
				pl.get(weight.atomIndex, size_t(0), "atomIndex", true);
				if(!weight.atomIndex)
					throw string("Atom index should be a positive integer");
				if(weight.atomIndex > e.iInfo.species[weight.specieIndex]->atpos.size())
					throw "Atom index exceeds number of atoms for species '" + spName + "'";
				weight.atomIndex--; //store 0-based index internally
			}
			//Get radius / half-width for all sphere and slice modes
			if(weight.type==DOS::Weight::Slice || weight.type==DOS::Weight::Sphere
				|| weight.type==DOS::Weight::AtomSlice || weight.type==DOS::Weight::AtomSphere)
			{	pl.get(weight.radius, 0., "r", true);
				if(weight.radius <= 0)
					throw string("Radius / half-width of weight function must be > 0");
			}
			//Get lattice plane direction for slice modes
			if(weight.type==DOS::Weight::Slice || weight.type==DOS::Weight::AtomSlice)
			{	pl.get(weight.direction[0], 0, "i0", true);
				pl.get(weight.direction[1], 0, "i1", true);
				pl.get(weight.direction[2], 0, "i2", true);
				if(!weight.direction.length_squared())
					throw string("Lattice plane direction (0,0,0) is invalid");
			}
			//Get filename for File mode:
			if(weight.type==DOS::Weight::File)
			{	pl.get(weight.filename, string(), "filename", true);
				//Check if file exists and is readable:
				FILE* fp = fopen(weight.filename.c_str(), "r");
				if(!fp) throw "File '"+weight.filename+"' cannot be opened for reading.\n";
				fclose(fp);
			}
			dos.weights.push_back(weight);
		}
		e.dump.insert(std::make_pair(DumpFreq_End, DumpDOS)); //dump at end
	}

	void printStatus(Everything& e, int iRep)
	{	assert(e.dump.dos);
		DOS& dos = *(e.dump.dos);
		logPrintf("%s", boolMap.getString(dos.occupied));
		for(const DOS::Weight& weight: dos.weights)
		{	logPrintf(" \\\n\t%s", weightTypeMap.getString(weight.type));
			switch(weight.type)
			{	case DOS::Weight::Total:
				case DOS::Weight::Delim: //never encountered, just to suppress compiler warning
					break; //no arguments
				case DOS::Weight::Slice:
				case DOS::Weight::Sphere:
				{	vector3<> center = (e.iInfo.coordsType==CoordsCartesian) ? e.gInfo.R * weight.center : weight.center;
					logPrintf(" %lg %lg %lg   %lg", center[0], center[1], center[2], weight.radius);
					if(weight.type == DOS::Weight::Slice)
						logPrintf("   %d %d %d", weight.direction[0], weight.direction[1], weight.direction[2]);
					break;
				}
				case DOS::Weight::AtomSlice:
				case DOS::Weight::AtomSphere:
				{	logPrintf(" %s %lu   %lg", e.iInfo.species[weight.specieIndex]->name.c_str(), weight.atomIndex+1, weight.radius);
					if(weight.type == DOS::Weight::AtomSlice)
						logPrintf("   %d %d %d", weight.direction[0], weight.direction[1], weight.direction[2]);
					break;
				}
				case DOS::Weight::File:
					logPrintf(" %s", weight.filename.c_str());
					break;
			}
		}
	}
}
commandDensityOfStates;
