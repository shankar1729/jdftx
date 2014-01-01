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

struct CommandVexternal : public Command
{
	CommandVexternal() : Command("Vexternal")
	{
		format = "<filename> | <filenameUp> <filenameDn>";
		comments =
			"Include an external potential (in hartrees) for the electrons\n"
			"(real space binary). Specify two files if V is spin-polarized.";
	}

	void process(ParamList& pl, Everything& e)
	{	e.eVars.VexternalFilename.resize(1);
		pl.get(e.eVars.VexternalFilename[0], string(), "filename", true);
		//Check if a second file has been specified:
		string filenameDn;
		pl.get(filenameDn, string(), "filenameDn");
		if(filenameDn.length())
			e.eVars.VexternalFilename.push_back(filenameDn);
	}

	void printStatus(Everything& e, int iRep)
	{	for(const string& filename: e.eVars.VexternalFilename)
			logPrintf("%s ", filename.c_str());
	}
}
commandVexternal;


struct CommandBoxPot : public Command
{
	CommandBoxPot() : Command("box-potential")
	{
		format = "xmin xmax ymin ymax zmin zmax Vin Vout [<convolve_radius>=0.1]";
		comments =
			"Include an step-function shaped external potential (in hartrees) for the electrons";
	
		allowMultiple = true;
	}

	void process(ParamList& pl, Everything& e)
	{	ElecVars::BoxPotential bP;
		const char* dirNames[3] = { "x", "y", "z" };
		for(int k=0; k<3; k++)
		{	pl.get(bP.min[k], 0., dirNames[k]+string("min"), true);
			pl.get(bP.max[k], 0., dirNames[k]+string("max"), true);
			if(bP.max[k]<bP.min[k])
				throw(string("max must be smaller than min for each dimension"));
		}
		pl.get(bP.Vin, 0., "Vin", true);
		pl.get(bP.Vout, 0., "Vout", true);
		pl.get(bP.convolve_radius, 0.2, "convolve_radius", false);
		e.eVars.boxPot.push_back(bP);
	}

	void printStatus(Everything& e, int iRep)
	{	const ElecVars::BoxPotential& bP = e.eVars.boxPot[iRep];
		logPrintf("%.5g %.5g %.5g %.5g %.5g %.5g    %.5g %.5g  %.5g",
			 bP.min[0], bP.max[0], bP.min[1], bP.max[1], bP.min[2], bP.max[2],
			 bP.Vin, bP.Vout, bP.convolve_radius);
	}
}
commandBoxPot;

