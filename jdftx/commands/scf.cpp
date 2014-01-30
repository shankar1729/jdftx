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
	SCFpm_nIterations, 
	SCFpm_energyDiffThreshold,
	SCFpm_eigDiffThreshold,
	SCFpm_residualThreshold,
	SCFpm_mixedVariable,
	SCFpm_vectorExtrapolation,
	SCFpm_verbose,
	SCFpm_mixFraction,
	SCFpm_qKerker,
	SCFpm_qMetric,
	SCFpm_history,
	SCFpm_single_particle_constraint,
	//Delimiter used in parsing:
	SCFpm_scfDelim
};

EnumStringMap<SCFparamsMember> scfParamsMap
(	SCFpm_nIterations, "nIterations",
	SCFpm_energyDiffThreshold, "energyDiffThreshold",
	SCFpm_eigDiffThreshold, "eigDiffThreshold",
	SCFpm_residualThreshold, "residualThreshold",
	SCFpm_mixedVariable, "mixedVariable",
	SCFpm_vectorExtrapolation, "vectorExtrapolation",
	SCFpm_verbose, "verbose",
	SCFpm_mixFraction, "mixFraction",
	SCFpm_qKerker, "qKerker",
	SCFpm_qMetric, "qMetric",
	SCFpm_history, "history",
	SCFpm_single_particle_constraint, "single-particle-constraint"
);
EnumStringMap<SCFparamsMember> scfParamsDescMap
(	SCFpm_nIterations, "maximum iterations (single point calculation if 0)",
	SCFpm_energyDiffThreshold, "convergence threshold for energy difference between successive iterations",
	SCFpm_eigDiffThreshold, "convergence threshold for the RMS difference in KS eigenvalues between successive iterations",
	SCFpm_residualThreshold, "convergence threshold for the residual in the mixed variable (density or potential)",
	SCFpm_mixedVariable, "whether density or potential will be mixed at each step",
	SCFpm_vectorExtrapolation, "algorithm to use in vector extrapolation: plainMixing, DIIS",
	SCFpm_verbose, "whether the inner eigenvalue solver will print or not",
	SCFpm_mixFraction, "maximum fraction of new variable mixed in at each step (default 0.5)",
	SCFpm_qKerker, "wavevector controlling Kerker preconditioning (default: auto-set to Gmin)",
	SCFpm_qMetric, "wavevector controlling the DIIS metric (default: auto-set to Gmin)",
	SCFpm_history, "number of past residuals and vectors are kept cached and used in DIIS",
	SCFpm_single_particle_constraint, "enforces the single particle constraint on the exchange"
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
		forbid("fix-electron-potential");
		forbid("spin-restricted");
	}
	
	void process(ParamList& pl, Everything& e)
	{	e.cntrl.scf = true;
	
		while(true)
		{	SCFparamsMember key;
			pl.get(key, SCFpm_scfDelim, scfParamsMap, "key", false);
			
			switch(key)
			{	case SCFpm_energyDiffThreshold: pl.get(e.scfParams.energyDiffThreshold, 1e-8, "energyDiffThreshold", true); break;
				case SCFpm_eigDiffThreshold: pl.get(e.scfParams.eigDiffThreshold, 1e-8, "eigDiffThreshold", true); break;
				case SCFpm_residualThreshold: pl.get(e.scfParams.residualThreshold, 1e-7, "residualThreshold", true); break;
				case SCFpm_mixedVariable: pl.get(e.scfParams.mixedVariable, SCFparams::MV_Potential, scfMixing, "mixedVariable", true); break;
				case SCFpm_nIterations: pl.get(e.scfParams.nIterations, 20, "nIterations", true); break;
				case SCFpm_vectorExtrapolation: pl.get(e.scfParams.vectorExtrapolation, SCFparams::VE_DIIS, scfExtrapolation, "vectorExtrapolation", true); break;
				case SCFpm_verbose: pl.get(e.scfParams.verbose, false, boolMap, "verbose", true); break;
				case SCFpm_mixFraction: pl.get(e.scfParams.mixFraction, 0.5, "mixFraction", true); break;
				case SCFpm_qKerker: pl.get(e.scfParams.qKerker, -1., "qKerker", true); break;
				case SCFpm_qMetric: pl.get(e.scfParams.qMetric, -1., "qMetric", true); break;
				case SCFpm_history: pl.get(e.scfParams.history, 10, "history", true); break;
				case SCFpm_single_particle_constraint: pl.get(e.scfParams.sp_constraint, 0., "single-particle-constraint", true); break;
				case SCFpm_scfDelim: return; //end of input
			}
		}
	}
	
	void printStatus(Everything& e, int iRep)
	{	
		#define PRINT(param,format) logPrintf(" \\\n\t" #param "\t" #format, e.scfParams.param);
		PRINT(nIterations, %i)
		PRINT(energyDiffThreshold, %lg)
		PRINT(eigDiffThreshold, %lg)
		PRINT(residualThreshold, %lg)
		logPrintf(" \\\n\tmixedVariable\t%s", scfMixing.getString(e.scfParams.mixedVariable));
		logPrintf(" \\\n\tvectorExtrapolation\t%s", scfExtrapolation.getString(e.scfParams.vectorExtrapolation));
		if(e.scfParams.verbose) logPrintf(" \\\n\tverbose");
		PRINT(mixFraction, %lg)
		PRINT(qKerker, %lg)
		PRINT(qMetric, %lg)
		PRINT(history, %d)
		#undef PRINT
	}
	
} commandsScfParams;


struct CommandEigenShift : public Command
{
	CommandEigenShift() : Command("eigen-shift")
	{
		format = "<qnum> <band> <shift> [<fromHOMO>=yes]";
		comments = "When fillings are computed during electronic-scf, shifts KS eigenvalue of (qnum, band) by shift\n"
					"Band index is calculated from HOMO, unless fromHOMO is set to no\n";
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
