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

EnumStringMap<int> lCodeMap(
	0, "s",
	1, "p",
	2, "d",
	3, "f",
	-1, "c" //custom channel / last channel in file
);

struct CommandIonSpecies : public Command
{
	CommandIonSpecies() : Command("ion-species")
	{
		format = "<id> <Z> <mass> <pot-file> <pulayfile> <projector>\n"
			"\t | <id> <Z> <mass> <cpi-file> <pulayfile> fhi [<unused> <local-l>]\n"
			"\t | [<path>/]<id>.fhi [<pulayfile>]\n"
			"\t | [<path>/]<id>.uspp [<startBytes>=4] [<endBytes>=4] [<pulayfile>]";
		comments =
			"Define ion species identified by <id> with valence charge <Z> and mass <mass>.\n"
			"  <pulayfile>: filename for pulay info (\"none\" to disable)\n"
			"  <projector>: If \"projector\", pot file has a projectors\n"
			"  <local-l> = "+lCodeMap.optionList()+", select local channel (for fhi psp's)\n"
			"In the .fhi/.uspp options (recommended), all the info is obtained from the file,\n"
			"including the <id>, which is obtained from the filename (minus path/extension).\n"
			"If <pulayfile> is not specified and [<path>/]<id>.pulay exists, pulay data will\n"
			"be read from that file; set <pulayfile> to \"none\" to ignore such files.\n"
			"The USPP format is a fortran sequential binary file with start and end markers\n"
			"with fortran-compiler dependent lengths. The defaults are valid for gfortran";
		allowMultiple = true;
	}

	void process(ParamList& pl, Everything& e)
	{	std::shared_ptr<SpeciesInfo> specie(new SpeciesInfo);
		pl.get(specie->name, string(), "id", true);
		
		//Check for .fhi and .uspp "smart" pseudopotentials
		bool smart = false;
		if(specie->name.length()>4 && specie->name.substr(specie->name.length()-4)==".fhi")
		{	specie->pspFormat = SpeciesInfo::Fhi;
			smart = true;
		}
		if(specie->name.length()>5 && specie->name.substr(specie->name.length()-5)==".uspp")
		{	specie->pspFormat = SpeciesInfo::Uspp;
			smart = true;
			//Get the optional fortran record lengths:
			pl.get(specie->recStartLen, 4, "startBytes");
			pl.get(specie->recStopLen, 4, "endBytes");
		}
		
		if(smart)
		{	specie->potfilename = specie->name;
			//Drop the extension from the name:
			specie->name = specie->name.substr(0, specie->name.find_last_of("."));
			//Drop the path from the name:
			size_t lastSlash = specie->name.find_last_of("\\/");
			if(lastSlash != string::npos)
				specie->name = specie->name.substr(lastSlash+1);
			//Check for a pulay file:
			pl.get(specie->pulayfilename, string(), "pulayfile"); //manual pulay file override
			if(!specie->pulayfilename.length()) //check for default pulay file:
			{	specie->pulayfilename =
					specie->potfilename.substr(0, specie->potfilename.find_last_of("."))
					+ ".pulay";
				FILE* fpPulay = fopen(specie->pulayfilename.c_str(), "r");
				if(!fpPulay) specie->pulayfilename = "none"; //disable if such a file does not exist
			}
			//Everything else is obtained from the smart pseudopotential during SpeciesInfo::setup()
		}
		else
		{	//Get the required quantities from the command parameters:
			pl.get(specie->Z, 0.0, "Z", true);
			pl.get(specie->mass, 0.0, "mass", true);
			pl.get(specie->potfilename, string(), "psp-file", true);
			pl.get(specie->pulayfilename, string(), "pulayfile", true);
			//Determine the format using the projector/fhi specification:
			string projfhi;
			pl.get(projfhi, string(), "noprojector");
			if(projfhi=="fhi")
			{	specie->pspFormat = SpeciesInfo::Cpi;
				specie->readProjectors = false;
				string lMax; pl.get(lMax, string("unused"), "unused"); //unused, lmax is determined from file
				pl.get(specie->lLocCpi, -1, lCodeMap, "local-l");
			}
			else
			{	specie->pspFormat = SpeciesInfo::Pot;
				if(projfhi=="projector") specie->readProjectors = true;
			}
		}
		//Check for duplicates, add to the list:
		for(auto sp: e.iInfo.species)
			if(specie->name==sp->name)
				throw string("Ion species "+specie->name+" had been defined more than once");
		e.iInfo.species.push_back(specie);
	}

	void printStatus(Everything& e, int iRep)
	{	const SpeciesInfo& specie = *(e.iInfo.species[iRep]);
		if(specie.pspFormat==SpeciesInfo::Fhi)
			logPrintf("%s", specie.potfilename.c_str());
		else if(specie.pspFormat==SpeciesInfo::Uspp)
			logPrintf("%s %d %d", specie.potfilename.c_str(), specie.recStartLen, specie.recStopLen);
		else
		{	logPrintf("%s %lg %lg %s %s %s", specie.name.c_str(), specie.Z, specie.mass,
				specie.potfilename.c_str(), specie.pulayfilename.c_str(), specie.readProjectors ? "projector" : "");
			if(specie.pspFormat==SpeciesInfo::Cpi) logPrintf("fhi <unused> %s", lCodeMap.getString(specie.lLocCpi));
		}
	}
}
commandIonSpecies;
