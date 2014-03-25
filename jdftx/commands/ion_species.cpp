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
		format = "<id> <Z> <mass> <pot-file> <pulayfile>\n"
			"\t | <id> <Z> <mass> <cpi-file> <pulayfile> fhi [<unused> <local-l>]\n"
			"\t | [<path>/]<id>.fhi [<pulayfile>]\n"
			"\t | [<path>/]<id>.uspp [<startBytes>=4] [<endBytes>=4] [<pulayfile>]";
		comments =
			"Define ion species identified by <id> with valence charge <Z> and mass <mass>.\n"
			"  <pulayfile>: filename for pulay info (\"none\" to disable)\n"
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
			string fhi;
			pl.get(fhi, string(), "");
			if(fhi=="fhi")
			{	specie->pspFormat = SpeciesInfo::Cpi;
				string lMax; pl.get(lMax, string("unused"), "unused"); //unused, lmax is determined from file
				pl.get(specie->lLocCpi, -1, lCodeMap, "local-l");
			}
			else specie->pspFormat = SpeciesInfo::Pot;
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
		{	logPrintf("%s %lg %lg %s %s", specie.name.c_str(), specie.Z, specie.mass,
				specie.potfilename.c_str(), specie.pulayfilename.c_str());
			if(specie.pspFormat==SpeciesInfo::Cpi) logPrintf(" fhi <unused> %s", lCodeMap.getString(specie.lLocCpi));
		}
	}
}
commandIonSpecies;


struct CommandChargeball : public Command
{
	CommandChargeball() : Command("chargeball")
	{
		format = "<species-id> <unusedInteger> <width> [<norm>=2.0]";
		comments =
			"Gaussian chargeball of norm <norm> and width <width> for species <id>\n"
			"(<unusedInteger> used to be the unused <Z_core>, and has been retained for compatibility)\n"
			"This feature is deprecated, use a pseudopotential with partial core correction instead.";
		allowMultiple = true;

		require("ion-species");
	}

	void process(ParamList& pl, Everything& e)
	{	string id;
		pl.get(id, string(), "species-id", true);
		for(auto sp: e.iInfo.species)
			if(sp->name == id)
			{	int unusedInteger;
				pl.get(unusedInteger, 0, "unusedInteger", true);
				pl.get(sp->width_chargeball, 0.0, "width", true);
				pl.get(sp->Z_chargeball, 2.0, "norm");
				return;
			}
		throw string("Species "+id+" has not been defined");
	}

	void printStatus(Everything& e, int iRep)
	{	if(unsigned(iRep) < e.iInfo.species.size())
		{	const SpeciesInfo& sp = *(e.iInfo.species[iRep]);
			logPrintf("%s -1 %lg %lg", sp.name.c_str(), sp.width_chargeball, sp.Z_chargeball);
		}
	}
}
commandChargeball;


struct CommandTauCore : public Command
{
	CommandTauCore() : Command("tau-core")
	{
		format = "<species-id> [<rCut>=0] [<plot>=yes|no]";
		comments =
			"Control generation of kinetic energy core correction for species <id>.\n"
			"The core KE density is set to the Thomas-Fermi + von-Weisacker functional\n"
			"of the core electron density (if any), and is pseudized inside within <rCut>\n"
			"If <rCut>=0, it is chosen to be 1.5 times the location of the first radial\n"
			"maximum in the TF+VW KE density. Optionally, if <plot>=yes, the resulting\n"
			"core KE density (and electron density) are output to a gnuplot-friendly file.";
		allowMultiple = true;

		require("ion-species");
	}

	void process(ParamList& pl, Everything& e)
	{	string id;
		pl.get(id, string(), "species-id", true);
		for(auto sp: e.iInfo.species)
			if(sp->name == id)
			{	pl.get(sp->tauCore_rCut, 0., "rCut", true);
				pl.get(sp->tauCorePlot, false, boolMap, "plot");
				return;
			}
		throw string("Species "+id+" has not been defined");
	}

	void printStatus(Everything& e, int iRep)
	{	if(unsigned(iRep) < e.iInfo.species.size())
		{	const SpeciesInfo& sp = *(e.iInfo.species[iRep]);
			logPrintf("%s %lg %s", sp.name.c_str(), sp.tauCore_rCut, boolMap.getString(sp.tauCorePlot));
		}
	}
}
commandTauCore;


struct CommandAddU : public Command
{
	CommandAddU() : Command("add-U")
	{	format = "<species> <orbDesc> <UminusJ> [ <species2> ... ]";
		comments =
			"Add U correction (DFT+U) to specified species and orbitals, in the simplified\n"
			"rotationally-invariant scheme of [Dudarev et al, Phys. Rev. B 57, 1505], where\n"
			"the correction depends only on U - J.\n"
			"  <species> is one of the pseudopotential identifiers.\n"
			"  <orbDesc> is one of s,p,d or f.\n"
			"  <UminusJ> = U-J is the on-site correction energy in hartrees.\n"
			"Repeat the sequence for corrections to multiple species.\n"
			"If pseudoatom has multiple shells of same angular momentum, prefix <orbDesc>\n"
			"with a number e.g. 1p or 2p to select the first or second p shell respectively.";
		
		require("ion-species");
	}
	
	void process(ParamList& pl, Everything& e)
	{	e.eInfo.hasU = false;
		string id;
		pl.get(id, string(), "species", true);
		while(id.length())
		{	bool spFound = false;
			for(auto sp: e.iInfo.species)
				if(sp->name == id)
				{	SpeciesInfo::PlusU plusU;
					//Get the orbital description:
					string orbCode; pl.get(orbCode, string(), "orbDesc", true);
					size_t lPos = string("spdf").find_first_of(orbCode.back());
					if(lPos==string::npos) throw  "'" + orbCode + "' is not a valid orbital code.";
					plusU.l = int(lPos);
					plusU.n = 0;
					if(orbCode.length() > 1)
					{	plusU.n = atoi(orbCode.substr(0, orbCode.length()-1).c_str()) - 1;
						if(plusU.n<0) throw string("Principal quantum number in orbital description must be a positive integer");
					}
					//Get the value:
					pl.get(plusU.UminusJ, 0., "UminusJ", true);
					//Add U decsriptor to species:
					sp->plusU.push_back(plusU);
					spFound = true;
					e.eInfo.hasU = true;
					break;
				}
			if(!spFound) throw string("Species "+id+" has not been defined");
			//Check for additional species:
			pl.get(id, string(), "species");
		}
		Citations::add("Simplified rotationally-invariant DFT+U", "S. L. Dudarev et al., Phys. Rev. B 57, 1505 (1998)");
	}
	
	void printStatus(Everything& e, int iRep)
	{	bool first = true;
		for(auto sp: e.iInfo.species)
			for(auto plusU: sp->plusU)
			{	if(!first) logPrintf(" \\\n"); first = false;
				ostringstream oss;
				if(plusU.n) oss << (plusU.n + 1);
				oss << string("spdf")[plusU.l];
				logPrintf("\t%s %s %lg", sp->name.c_str(), oss.str().c_str(), plusU.UminusJ);
			}
	}
}
commandAddU;
