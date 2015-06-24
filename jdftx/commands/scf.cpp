/*-------------------------------------------------------------------
Copyright 2013 Deniz Gunceler, Ravishankar Sundararaman

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

//--------- PulayParams parser -----

enum PulayParamsMember
{	PPM_nIterations,
	PPM_energyDiffThreshold,
	PPM_residualThreshold,
	PPM_mixFraction,
	PPM_qKerker,
	PPM_qMetric,
	PPM_history
};

EnumStringMap<PulayParamsMember> pulayParamsMap
(	PPM_nIterations, "nIterations",
	PPM_energyDiffThreshold, "energyDiffThreshold",
	PPM_residualThreshold, "residualThreshold",
	PPM_mixFraction, "mixFraction",
	PPM_qKerker, "qKerker",
	PPM_qMetric, "qMetric",
	PPM_history, "history"
);

EnumStringMap<PulayParamsMember> pulayParamsDescMap
(	PPM_nIterations, "maximum iterations (single point calculation if 0)",
	PPM_energyDiffThreshold, "convergence threshold for energy difference between successive iterations",
	PPM_residualThreshold, "convergence threshold for the residual in the mixed variable",
	PPM_mixFraction, "mix fraction (default 0.5)",
	PPM_qKerker, "wavevector controlling Kerker preconditioning (default: 0.8 bohr^-1)",
	PPM_qMetric, "wavevector controlling the metric for overlaps (default: 0.8 bohr^-1)",
	PPM_history, "number of past residuals that are cached and used for mixing"
);

//Base class for pulay-mixing commands
struct CommandPulay : public Command
{
	CommandPulay(string name) : Command(name)
	{
		format = "<key1> <value1> <key2> <value2> ...";
		//Derived classes must set all other necessary parameters
	}
	
protected:
	void process(ParamList& pl, Everything& e, PulayParams& pp)
	{	while(true)
		{	string keyStr;
			pl.get(keyStr, string(), "key", false);
			if(!keyStr.length()) return; //End of input
			PulayParamsMember key;
			if(pulayParamsMap.getEnum(keyStr.c_str(), key))
			{	//Process base-class parameters:
				switch(key)
				{	case PPM_energyDiffThreshold: pl.get(pp.energyDiffThreshold, 1e-8, "energyDiffThreshold", true); break;
					case PPM_residualThreshold: pl.get(pp.residualThreshold, 1e-7, "residualThreshold", true); break;
					case PPM_nIterations: pl.get(pp.nIterations, 50, "nIterations", true); break;
					case PPM_mixFraction: pl.get(pp.mixFraction, 0.5, "mixFraction", true); break;
					case PPM_qKerker: pl.get(pp.qKerker, 0.8, "qKerker", true); break;
					case PPM_qMetric: pl.get(pp.qMetric, 0.8, "qMetric", true); break;
					case PPM_history: pl.get(pp.history, 10, "history", true); if(pp.history<1) throw string("<history> must be >= 1"); break;
				}
			}
			else process_sub(keyStr, pl, e);
		}
	}
	
	void printStatus(const PulayParams& pp)
	{
		#define PRINT(param,format) logPrintf(" \\\n\t" #param "\t" #format, pp.param);
		PRINT(nIterations, %i)
		PRINT(energyDiffThreshold, %lg)
		PRINT(residualThreshold, %lg)
		PRINT(mixFraction, %lg)
		PRINT(qKerker, %lg)
		PRINT(qMetric, %lg)
		PRINT(history, %d)
		#undef PRINT
	}
	
	//Derived class should handle keys other than those in PulayParams, and throw error if key is not recognized
	virtual void process_sub(string keyStr, ParamList& pl, Everything& e) = 0;
};

//-------- CommandElectronicSCF ------

enum SCFparamsMember
{
	SCFpm_nEigSteps,
	SCFpm_eigDiffThreshold,
	SCFpm_mixedVariable,
	SCFpm_verbose,
	SCFpm_mixFractionMag,
	SCFpm_maxOverlap,
};

EnumStringMap<SCFparamsMember> scfParamsMap
(	SCFpm_nEigSteps, "nEigSteps",
	SCFpm_eigDiffThreshold, "eigDiffThreshold",
	SCFpm_mixedVariable, "mixedVariable",
	SCFpm_verbose, "verbose",
	SCFpm_mixFractionMag, "mixFractionMag",
	SCFpm_maxOverlap, "maximum-overlap-method"
);
EnumStringMap<SCFparamsMember> scfParamsDescMap
(	SCFpm_nEigSteps, "number of eigenvalue steps per iteration (if 0, limited by electronic-minimize nIterations)",
	SCFpm_eigDiffThreshold, "convergence threshold for the RMS difference in KS eigenvalues between successive iterations",
	SCFpm_mixedVariable, "whether density or potential will be mixed at each step",
	SCFpm_verbose, "whether the inner eigenvalue solver will print or not",
	SCFpm_mixFractionMag, "mix fraction for magnetization density / potential (default 1.5)",
	SCFpm_maxOverlap, "uses the maximum-overlap method to determine fillings, discards elec-fermi-fillings"
);

EnumStringMap<SCFparams::MixedVariable> scfMixing
(	SCFparams::MV_Density, "Density",
	SCFparams::MV_Potential, "Potential"
);

struct CommandElectronicScf: public CommandPulay
{
	CommandElectronicScf() : CommandPulay("electronic-scf")
	{	
		comments =
			"Enables self-consistent field optimization of electronic state.\n"
			"Possible keys and value types to control SCF optimization:"
			+ addDescriptions(pulayParamsMap.optionList(), linkDescription(pulayParamsMap, pulayParamsDescMap))
			+ addDescriptions(scfParamsMap.optionList(), linkDescription(scfParamsMap, scfParamsDescMap))
			+ "\n\nAny number of these key-value pairs may be specified in any order.";
		hasDefault = false;
		forbid("fix-electron-density");
		forbid("fix-electron-potential");
		require("elec-eigen-algo");
	}
	
	void process(ParamList& pl, Everything& e)
	{	e.cntrl.scf = true;
		SCFparams& sp = e.scfParams;
		sp.nEigSteps = (e.cntrl.elecEigenAlgo==ElecEigenCG) ? 40 : 2; //default eigenvalue steps based on algo
		CommandPulay::process(pl, e, sp);
	}
	
	void process_sub(string keyStr, ParamList& pl, Everything& e)
	{	SCFparams& sp = e.scfParams;
		SCFparamsMember key;
		if(scfParamsMap.getEnum(keyStr.c_str(), key))
		{	switch(key)
			{	case SCFpm_nEigSteps: pl.get(sp.nEigSteps, 0, "nEigSteps", true); break;
				case SCFpm_eigDiffThreshold: pl.get(sp.eigDiffThreshold, 1e-8, "eigDiffThreshold", true); break;
				case SCFpm_mixedVariable: pl.get(sp.mixedVariable, SCFparams::MV_Potential, scfMixing, "mixedVariable", true); break;
				case SCFpm_verbose: pl.get(sp.verbose, false, boolMap, "verbose", true); break;
				case SCFpm_mixFractionMag: pl.get(sp.mixFractionMag, 1.5, "mixFractionMag", true); break;
				case SCFpm_maxOverlap: e.eInfo.fillingsUpdate = ElecInfo::MaximumOverlapMethod; break;
			}
		}
		else throw string("Parameter <key> must be one of " + pulayParamsMap.optionList() + "|" + scfParamsMap.optionList());
	}
	
	void printStatus(Everything& e, int iRep)
	{	const SCFparams& sp = e.scfParams;
		CommandPulay::printStatus(sp); //base class parameters
		#define PRINT(param,format) logPrintf(" \\\n\t" #param "\t" #format, sp.param);
		PRINT(nEigSteps, %i)
		PRINT(eigDiffThreshold, %lg)
		logPrintf(" \\\n\tmixedVariable\t%s", scfMixing.getString(sp.mixedVariable));
		logPrintf(" \\\n\tverbose\t%s", boolMap.getString(sp.verbose));
		PRINT(mixFractionMag, %lg)
		if(e.eInfo.fillingsUpdate == ElecInfo::MaximumOverlapMethod) logPrintf(" \\\n\tmaximum-overlap-method");
		#undef PRINT
	}
}
commandElectronicScf;


struct CommandEigenShift : public Command
{
	CommandEigenShift() : Command("eigen-shift")
	{
		format = "<qnum> <band> <shift> [<fromHOMO>=yes]";
		comments =
			"When fillings are computed during electronic-scf, shifts Kohn-Sham\n"
			"eigenvalue of (qnum, band) by shift. Band index is calculated\n"
			"relative to HOMO if <fromHOMO> = yes and is absolute otherwise.\n";
		allowMultiple = true;
		require("electronic-scf");
		require("elec-fermi-fillings");
		forbid("custom-filling");
	}

	void process(ParamList& pl, Everything& e)
	{	SCFparams::EigenShift es;
		pl.get(es.q, 0, "qnum", true);
		pl.get(es.n, 0, "band", true);
		pl.get(es.shift, 0., "filling", true);
		pl.get(es.fromHOMO, true, boolMap, "fromHOMO", false);
		e.scfParams.eigenShifts.push_back(es);
	}

	void printStatus(Everything& e, int iRep)
	{	const SCFparams::EigenShift& es =  e.scfParams.eigenShifts[iRep];
		logPrintf("%i %i %.5e %s", es.q, es.n, es.shift, boolMap.getString(es.fromHOMO));
	}
}
commandEigenShift;


struct CommandPcmNonlinearScf: public CommandPulay
{
	CommandPcmNonlinearScf() : CommandPulay("pcm-nonlinear-scf")
	{	
		comments =
			"Enables self-consistent field optimization for nonlinear PCM fluids.\n"
			"Possible keys and value types to control SCF optimization:"
			+ addDescriptions(pulayParamsMap.optionList(), linkDescription(pulayParamsMap, pulayParamsDescMap))
			+ "\n\nAny number of these key-value pairs may be specified in any order.";
		hasDefault = false;
	}
	
	void process(ParamList& pl, Everything& e)
	{	e.eVars.fluidParams.nonlinearSCF = true;
		PulayParams& pp = e.eVars.fluidParams.scfParams;
		pp.energyLabel = "Adiel";
		pp.linePrefix = "NonlinearFluidSCF: ";
		pp.energyFormat = "%+.15lf";
		pp.energyDiffThreshold = 1e-7;
		pp.nIterations = 20;
		pp.fpLog = globalLog;
		CommandPulay::process(pl, e, pp); //only base class parameters
	}
	
	void process_sub(string keyStr, ParamList& pl, Everything& e)
	{	throw string("Parameter <key> must be one of " + pulayParamsMap.optionList());
	}
	
	void printStatus(Everything& e, int iRep)
	{	CommandPulay::printStatus(e.eVars.fluidParams.scfParams); //only base class parameters
	}
}
commandPcmNonlinearScf;
