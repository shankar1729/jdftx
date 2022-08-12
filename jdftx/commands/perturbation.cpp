/*
 * perturbation.cpp
 *
 *  Created on: Jul 22, 2022
 *      Author: brandon
 */

#include <commands/minimize.h>

struct CommandSolvePerturbation : public CommandMinimize
{
	CommandSolvePerturbation() : CommandMinimize("solve-perturbation", "jdftx/Electronic/Optimization")
	{	emptyParamError =
			"   Note: nIterations defaults to 0 for lattice minimization,\n"
			"      and must be set manually to enable this feature.";
	}

    MinimizeParams& target(Everything& e) { return e.vptParams; }
    void process(ParamList& pl, Everything& e)
    {	e.vptParams.nIterations = 0; //override default value (100) in MinimizeParams.h
    	e.vptParams.energyDiffThreshold = 1e-6;
    	e.vptParams.dirUpdateScheme = MinimizeParams::FletcherReeves;
    	CommandMinimize::process(pl, e);
	}
}
CommandSolvePerturbation;


struct CommandDVexternal : public Command
{
	CommandDVexternal() : Command("dVexternal", "jdftx/Electronic/Parameters")
	{
		format = "<filename> | <filenameUp> <filenameDn>";
		comments =
			"Include an external potential (in hartrees) for the electrons\n"
			"(real space binary). Specify two files if V is spin-polarized.";
	}

	void process(ParamList& pl, Everything& e)
	{	e.vptInfo.dVexternalFilename.resize(1);
		pl.get(e.vptInfo.dVexternalFilename[0], string(), "filename", true);
		//Check if a second file has been specified:
		//string filenameDn;
		//pl.get(filenameDn, string(), "filenameDn");
		//if(filenameDn.length())
		//	e.eVars.VexternalFilename.push_back(filenameDn);
	}

	void printStatus(Everything& e, int iRep)
	{	for(const string& filename: e.vptInfo.dVexternalFilename)
			logPrintf("%s ", filename.c_str());
	}
}
commandDVexternal;
