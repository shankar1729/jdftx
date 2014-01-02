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

enum SCFparamsMember
{
	nIterations, 
	energyDiffThreshold,
	residualThreshold,
	mixedVariable,
	vectorExtrapolation,
	verbose,
	damping,
	history,
	//Delimiter used in parsing:
	scfDelim
};

EnumStringMap<SCFparamsMember> scfParamsMap
(	nIterations, "nIterations",
	energyDiffThreshold, "energyDiffThreshold",
	residualThreshold, "residualThreshold",
	mixedVariable, "mixedVariable",
	vectorExtrapolation, "vectorExtrapolation",
	verbose, "verbose",
	damping, "damping",
	history, "history"
);
EnumStringMap<SCFparamsMember> scfParamsDescMap
(	nIterations, "maximum iterations (single point calculation if 0)",
	energyDiffThreshold, "convergence threshold for energy difference between successive iterations",
	residualThreshold, "convergence threshold for the residual in the mixed variable (density or potential)",
    mixedVariable, "whether density or potential will be mixed at each step",
	vectorExtrapolation, "algorithm to use in vector extrapolation: plainMixing, DIIS",
	verbose, "whether the inner eigenvalue solver will print or not",
	damping, "damping parameter for vector extrapolation (default 0.5).",
	history, "Number of past residuals and vectors are kept cached and used in DIIS"
);

EnumStringMap<SCFparams::MixedVariable> scfMixing
(	SCFparams::MV_Density, "Density",
	SCFparams::MV_Potential, "Potential"
);

EnumStringMap<SCFparams::VectorExtrapolation> scfExtrapolation
(	SCFparams::VE_Plain, "Plain",
	SCFparams::VE_DIIS, "DIIS"
);

struct CommandsScfParams: public Command
{
	CommandsScfParams() : Command("electronic-scf")
	{	
		format = "<key1> <value1> <key2> <value2> ...";
		comments = "Enables self-consistent residual minimization.  If provided, keys adjust SCF parameters. Possible keys and value types are:"
			+ addDescriptions(scfParamsMap.optionList(), linkDescription(scfParamsMap, scfParamsDescMap))
			+ "\nAny number of these key-value pairs may be specified in any order.";
		hasDefault = false;
		forbid("fix-electron-density");
		forbid("spin-restricted");
	}
	
	void process(ParamList& pl, Everything& e)
	{	e.cntrl.scf = true;
	
		while(true)
		{	SCFparamsMember key;
			pl.get(key, scfDelim, scfParamsMap, "key", false);
			
			switch(key)
			{	case energyDiffThreshold: pl.get(e.scfParams.energyDiffThreshold, 1e-8, "energyDiffThreshold", true); break;
				case residualThreshold: pl.get(e.scfParams.residualThreshold, 1e-7, "residualThreshold", true); break;
				case mixedVariable: pl.get(e.scfParams.mixedVariable, SCFparams::MV_Potential, scfMixing, "mixedVariable", true); break;
				case nIterations: pl.get(e.scfParams.nIterations, 20, "nIterations", true); break;
				case vectorExtrapolation: pl.get(e.scfParams.vectorExtrapolation, SCFparams::VE_DIIS, scfExtrapolation, "vectorExtrapolation", true); break;
				case verbose: pl.get(e.scfParams.verbose, false, boolMap, "verbose", true); break;
				case damping: pl.get(e.scfParams.damping, 0.5, "damping", true); break;
				case history: pl.get(e.scfParams.history, 15, "history", true); break;
				case scfDelim: return; //end of input
			}
		}
	}
	
	void printStatus(Everything& e, int iRep)
	{	
		#define PRINT(param,format) logPrintf(" \\\n\t" #param "\t" #format, e.scfParams.param);
		PRINT(nIterations, %i)
		PRINT(energyDiffThreshold, %lg)
		PRINT(residualThreshold, %lg)
		PRINT(nIterations, %d)
		logPrintf(" \\\n\tmixedVariable\t%s", scfMixing.getString(e.scfParams.mixedVariable));
		logPrintf(" \\\n\tvectorExtrapolation\t%s", scfExtrapolation.getString(e.scfParams.vectorExtrapolation));
		if(e.scfParams.verbose) logPrintf(" \\\n\tverbose");
		PRINT(damping, %lg)
		PRINT(history, %d)
	} 
	
} commandsScfParams;


struct CommandEigenShift : public Command
{
	CommandEigenShift() : Command("eigen-shift")
	{
		format = "<qnum> <band> <shift> [<fromHOMO>]";
		comments = "When fillings are computed during electronic-scf, shifts KS eigenvalue of (qnum, band) by shift\n"
					"Band index is calculated from HOMO, unless set to false\n";
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
