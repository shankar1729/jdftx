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
		//TODO: Check if a second file has been specified:
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
CommandDVexternal;

struct CommandSetQpoint : public Command
{
	CommandSetQpoint() : Command("set-qpoint", "jdftx/Electronic/Parameters")
	{
		format = "<q0> <q1> <q2>";
		comments =
			"Set the bloch phase of the perturbation potential to exp(i q*r) where q = (<q0>, <q1>, <q2>)\n";
	}

    void process(ParamList& pl, Everything& e)
    {
		pl.get(e.vptInfo.qvec[0], 0.0, "q0");
		pl.get(e.vptInfo.qvec[1], 0.0, "q1");
		pl.get(e.vptInfo.qvec[2], 0.0, "q2");

		e.vptInfo.incommensurate = e.vptInfo.qvec.isNonzero();
		//TODO: Band structure in case q=0?
	}

	void printStatus(Everything& e, int iRep)
	{
		logPrintf("%f %f %f", e.vptInfo.qvec[0], e.vptInfo.qvec[1], e.vptInfo.qvec[2]);
	}
}
CommandSetQpoint;

struct CommandReadOffsetWfns : public Command
{
	CommandReadOffsetWfns() : Command("read-offset-wfns", "jdftx/Electronic/Parameters")
	{
		format = "<filename>";
		comments =
			"Reads in results of band structure calculations for incommensurate perturbations\n";
	}

    void process(ParamList& pl, Everything& e)
    {
		pl.get(e.vptInfo.wfnsFilename, string(), "filename", true);
	}

	void printStatus(Everything& e, int iRep)
	{
		logPrintf("%s", e.vptInfo.wfnsFilename.c_str());
	}
}
CommandReadOffsetWfns;

struct CommandTestPerturbationOperators : public Command
{
	CommandTestPerturbationOperators() : Command("test-perturbation-ops", "miscellaneous")
	{
		comments =
			"Perform finite difference tests of perturbation operators.\n";
	}

    void process(ParamList& pl, Everything& e)
    {
		e.vptInfo.testing = true;
	}

	void printStatus(Everything& e, int iRep)
	{
		return;
	}
}
CommandTestPerturbationOperators;

