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
	SCFpm_nEigSteps,
	SCFpm_energyDiffThreshold,
	SCFpm_eigDiffThreshold,
	SCFpm_residualThreshold,
	SCFpm_mixedVariable,
	SCFpm_verbose,
	SCFpm_mixFraction,
	SCFpm_mixFractionMag,
	SCFpm_qKerker,
	SCFpm_qMetric,
	SCFpm_history,
	SCFpm_single_particle_constraint,
	SCFpm_maxOverlap,
	//Delimiter used in parsing:
	SCFpm_scfDelim
};

EnumStringMap<SCFparamsMember> scfParamsMap
(	SCFpm_nIterations, "nIterations",
	SCFpm_nEigSteps, "nEigSteps",
	SCFpm_energyDiffThreshold, "energyDiffThreshold",
	SCFpm_eigDiffThreshold, "eigDiffThreshold",
	SCFpm_residualThreshold, "residualThreshold",
	SCFpm_mixedVariable, "mixedVariable",
	SCFpm_verbose, "verbose",
	SCFpm_mixFraction, "mixFraction",
	SCFpm_mixFractionMag, "mixFractionMag",
	SCFpm_qKerker, "qKerker",
	SCFpm_qMetric, "qMetric",
	SCFpm_history, "history",
	SCFpm_single_particle_constraint, "single-particle-constraint",
	SCFpm_maxOverlap, "maximum-overlap-method"
);
EnumStringMap<SCFparamsMember> scfParamsDescMap
(	SCFpm_nIterations, "maximum iterations (single point calculation if 0)",
	SCFpm_nEigSteps, "number of eigenvalue steps per iteration (if 0, limited by electronic-minimize nIterations)",
	SCFpm_energyDiffThreshold, "convergence threshold for energy difference between successive iterations",
	SCFpm_eigDiffThreshold, "convergence threshold for the RMS difference in KS eigenvalues between successive iterations",
	SCFpm_residualThreshold, "convergence threshold for the residual in the mixed variable (density or potential)",
	SCFpm_mixedVariable, "whether density or potential will be mixed at each step",
	SCFpm_verbose, "whether the inner eigenvalue solver will print or not",
	SCFpm_mixFraction, "mix fraction for total density / potential (default 0.5)",
	SCFpm_mixFractionMag, "mix fraction for magnetization density / potential (default 1.5)",
	SCFpm_qKerker, "wavevector controlling Kerker preconditioning (default: 0.8 bohr^-1)",
	SCFpm_qMetric, "wavevector controlling the DIIS metric (default: 0.8 bohr^-1)",
	SCFpm_history, "number of past residuals and vectors are kept cached and used in DIIS",
	SCFpm_single_particle_constraint, "enforces the single particle constraint on the exchange",
	SCFpm_maxOverlap, "uses the maximum-overlap method to determine fillings, discards elec-fermi-fillings"
);

EnumStringMap<SCFparams::MixedVariable> scfMixing
(	SCFparams::MV_Density, "Density",
	SCFparams::MV_Potential, "Potential"
);

struct CommandElectronicScf: public Command
{
	CommandElectronicScf() : Command("electronic-scf")
	{	
		format = "<key1> <value1> <key2> <value2> ...";
		comments = "Enables self-consistent residual minimization.  If provided, keys adjust SCF parameters. Possible keys and value types are:"
			+ addDescriptions(scfParamsMap.optionList(), linkDescription(scfParamsMap, scfParamsDescMap))
			+ "\n\nAny number of these key-value pairs may be specified in any order.";
		hasDefault = false;
		forbid("fix-electron-density");
		forbid("fix-electron-potential");
		forbid("spin-restricted");
		require("elec-eigen-algo");
	}
	
	void process(ParamList& pl, Everything& e)
	{	e.cntrl.scf = true;
		SCFparams& sp = e.scfParams;
		sp.nEigSteps = (e.cntrl.elecEigenAlgo==ElecEigenCG) ? 40 : 2; //default eigenvalue steps based on algo
		while(true)
		{	SCFparamsMember key;
			pl.get(key, SCFpm_scfDelim, scfParamsMap, "key", false);
			
			switch(key)
			{	case SCFpm_energyDiffThreshold: pl.get(sp.energyDiffThreshold, 1e-8, "energyDiffThreshold", true); break;
				case SCFpm_eigDiffThreshold: pl.get(sp.eigDiffThreshold, 1e-8, "eigDiffThreshold", true); break;
				case SCFpm_residualThreshold: pl.get(sp.residualThreshold, 1e-7, "residualThreshold", true); break;
				case SCFpm_mixedVariable: pl.get(sp.mixedVariable, SCFparams::MV_Potential, scfMixing, "mixedVariable", true); break;
				case SCFpm_nIterations: pl.get(sp.nIterations, 50, "nIterations", true); break;
				case SCFpm_nEigSteps: pl.get(sp.nEigSteps, 0, "nEigSteps", true); break;
				case SCFpm_verbose: pl.get(sp.verbose, false, boolMap, "verbose", true); break;
				case SCFpm_mixFraction: pl.get(sp.mixFraction, 0.5, "mixFraction", true); break;
				case SCFpm_mixFractionMag: pl.get(sp.mixFractionMag, 1.5, "mixFractionMag", true); break;
				case SCFpm_qKerker: pl.get(sp.qKerker, 0.8, "qKerker", true); break;
				case SCFpm_qMetric: pl.get(sp.qMetric, 0.8, "qMetric", true); break;
				case SCFpm_history: pl.get(sp.history, 10, "history", true); if(sp.history<1) throw string("<history> must be >= 1"); break;
				case SCFpm_single_particle_constraint: pl.get(sp.sp_constraint, 0., "single-particle-constraint", true); break;
				case SCFpm_maxOverlap: sp.MOMenabled=true; break;
				case SCFpm_scfDelim: return; //end of input
			}
		}
	}
	
	void printStatus(Everything& e, int iRep)
	{	const SCFparams& sp = e.scfParams;
		#define PRINT(param,format) logPrintf(" \\\n\t" #param "\t" #format, sp.param);
		PRINT(nIterations, %i)
		PRINT(nEigSteps, %i)
		PRINT(energyDiffThreshold, %lg)
		PRINT(eigDiffThreshold, %lg)
		PRINT(residualThreshold, %lg)
		logPrintf(" \\\n\tmixedVariable\t%s", scfMixing.getString(sp.mixedVariable));
		if(sp.verbose) logPrintf(" \\\n\tverbose");
		PRINT(mixFraction, %lg)
		PRINT(mixFractionMag, %lg)
		PRINT(qKerker, %lg)
		PRINT(qMetric, %lg)
		PRINT(history, %d)
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
