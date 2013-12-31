/*-------------------------------------------------------------------
Copyright 2013 Deniz Gunceler

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

enum  ResMinParameter
{
	nIterations, 
	energyDiffThreshold,
	mixedVariable,
	vectorExtrapolation,
	verbose,
	damping,
	history,
	//Delimiter used in parsing:
	scfDelim
};

EnumStringMap<ResMinParameter> ResMinParameterMap
(	nIterations, "nIterations",
	energyDiffThreshold, "energyDiffThreshold",
    mixedVariable, "mixedVariable",
    vectorExtrapolation, "vectorExtrapolation",
	verbose, "verbose",
	damping, "damping",
	history, "history"
);
EnumStringMap<ResMinParameter> resMinParameterDescrMap
(	nIterations, "maximum iterations (single point calculation if 0)",
	energyDiffThreshold, "convergence threshold for energy difference between successive iterations",
    mixedVariable, "whether density or potential will be mixed at each step",
	vectorExtrapolation, "algorithm to use in vector extrapolation: plainMixing, DIIS",
	verbose, "whether the inner eigenvalue solver will print or not",
	damping, "damping parameter for vector extrapolation (default 0.5).",
	history, "Number of past residuals and vectors are kept cached and used in DIIS"
);

EnumStringMap<MixedVariable> resMinMixing
(	density, "density",
	potential, "potential"
);

EnumStringMap<VectorExtrapolation> resMinExtrapolation
(	plain, "plain",
	DIIS, "DIIS"
);

struct CommandsScfParams: public Command
{
	CommandsScfParams() : Command("residual-minimize")
	{	
		format = "<key1> <value1> <key2> <value2> ...";
		comments = "Enables self-consistent residual minimization.  If provided, keys adjust SCF parameters. Possible keys and value types are:"
			+ addDescriptions(ResMinParameterMap.optionList(), linkDescription(ResMinParameterMap, resMinParameterDescrMap))
			+ "\nAny number of these key-value pairs may be specified in any order.";
		hasDefault = false;
		forbid("fix-electron-density");
		forbid("spin-restricted");
	}
	
	void process(ParamList& pl, Everything& e)
	{	e.cntrl.minimisingResidual = true;
	
		while(true)
		{	ResMinParameter key;
			pl.get(key, scfDelim, ResMinParameterMap, "key", false);
			
			switch(key)
			{	case energyDiffThreshold:
					pl.get(e.residualMinimizerParams.energyDiffThreshold, 1e-6, "energyDiffThreshold", true);
					break;
				case mixedVariable:
					pl.get(e.residualMinimizerParams.mixedVariable, density, resMinMixing, "mixedVariable", true);
					break;
				case nIterations:
					pl.get(e.residualMinimizerParams.nIterations, 10, "nIterations", true);
					break;
				case vectorExtrapolation:
					pl.get(e.residualMinimizerParams.vectorExtrapolation, plain, resMinExtrapolation, "vectorExtrapolation", true);					
					break;
				case verbose:
					e.residualMinimizerParams.verbose = true;				
					break;
				case damping:
					pl.get(e.residualMinimizerParams.damping, 0.5, "damping", true);
					break;
				case history:
					pl.get(e.residualMinimizerParams.history, 15, "history", true);
					break;
				case scfDelim: return; //end of input
			}
			
		}	
	}
	
	void printStatus(Everything& e, int iRep)
	{	
		#define PRINT(param,format) logPrintf(" \\\n\t" #param "\t" #format, e.residualMinimizerParams.param);
		
		PRINT(nIterations, %i)
		PRINT(energyDiffThreshold, %lg)
		logPrintf(" \\\n\tmixedVariable\t%s", resMinMixing.getString(e.residualMinimizerParams.mixedVariable));
		logPrintf(" \\\n\tvectorExtrapolation\t%s", resMinExtrapolation.getString(e.residualMinimizerParams.vectorExtrapolation));
		logPrintf(" \\\n\thistory\t%i", e.residualMinimizerParams.history);
		if(e.residualMinimizerParams.verbose) logPrintf(" \\\n\tverbose");
		PRINT(damping, %lg)
	} 
	
} commandsScfParams;


struct CommandEigenShift : public Command
{
	CommandEigenShift() : Command("eigen-shift")
	{
		format = "<qnum> <band> <shift> [<fromHOMO>]";
		comments = "When fillings are computed during residual-minimize, shifts KS eigenvalue of (qnum, band) by shift\n"
					"Band index is calculated from HOMO, unless set to false\n";
		allowMultiple = true;
		require("residual-minimize");
		require("elec-fermi-fillings");
		forbid("custom-filling");
	}

	void process(ParamList& pl, Everything& e)
	{	
		int q, n;
		double shift;
		bool fromHOMO;
		
		pl.get(q, 0, "qnum", true);
		pl.get(n, 0, "band", true);
		pl.get(shift, 0., "filling", true);
		pl.get(fromHOMO, true, boolMap, "fromHOMO", false);

		e.residualMinimizerParams.eigenShifts.push_back(EigenShift(q,n,shift, fromHOMO));
		
	}

	void printStatus(Everything& e, int iRep)
	{	std::vector<EigenShift>& eigenShifts =  e.residualMinimizerParams.eigenShifts;
		logPrintf("%i %i %.5e", eigenShifts[iRep].q, eigenShifts[iRep].n, eigenShifts[iRep].shift);
	}
}
commandEigenShift;