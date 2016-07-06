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

//! @file elec_misc.cpp Commands controlling electronic fillings algorithms

struct CommandElecFermiFillings : public Command
{
	CommandElecFermiFillings() : Command("elec-fermi-fillings", "jdftx/Electronic/Parameters")
	{
		format = "0 <kT>";
		comments =
			"Fermi-Dirac fillings at a temperature <kT> (in Hartrees).\n"
			"(The first unused argument, which used to be <mixInterval>,\n"
			" is present only for backward compatibility; set it to zero.)";
		
		require("lcao-params");
		forbid("fix-electron-density");
		forbid("fix-electron-potential");
	}

	void process(ParamList& pl, Everything& e)
	{	ElecInfo& eInfo = e.eInfo;
		string unusedParam;
		pl.get(unusedParam, string(), "unusedParam", true);
		eInfo.fillingsUpdate = ElecInfo::FillingsHsub;
		pl.get(eInfo.kT, 0.0, "kT", true);
	}

	void printStatus(Everything& e, int iRep)
	{	const ElecInfo& eInfo = e.eInfo;
		logPrintf("0 %lg", eInfo.kT);
	}
}
commandElecFermiFillings;

//-------------------------------------------------------------------------------------------------

struct CommandTargetMu : public Command
{
	CommandTargetMu() : Command("target-mu", "jdftx/Electronic/Parameters")
	{
		format = "<mu>";
		comments =
			"Fixed chemical potential <mu> (instead of fixed charge).\n"
			"Note that <mu> is absolute (relative to vacuum level) and in Hartrees.\n"
			"For example, potential V (in Volts) relative to SHE corresponds to\n"
			"mu = -(Vref + V)/27.2114, where Vref is the absolute SHE potential\n"
			"in Volts below vacuum; you could set Vref = 4.44 based on experiment\n"
			"or use the value calibrated using potentials of zero charge with\n"
			"the solvation model in use.";

		require("fluid-cation");
		require("fluid-anion");
		require("elec-fermi-fillings");
		forbid("elec-initial-charge");
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.eInfo.mu, 0.0, "mu", true);
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg\n", e.eInfo.mu);
	}
}
commandTargetMu;

//-------------------------------------------------------------------------------------------------

struct CommandElecInitialCharge : public Command
{
    CommandElecInitialCharge() : Command("elec-initial-charge", "jdftx/Initialization")
	{
		format = "<QNet>";
		comments =
			"Initialize a system with <QNet> excess electrons compared to a neutral system.";
		
		forbid("target-mu");
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.eInfo.Qinitial, 0., "QNet", true);
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lf", e.eInfo.Qinitial);
	}
}
CommandElecInitialCharge;

//-------------------------------------------------------------------------------------------------

struct CommandElecInitialMagnetization : public Command
{
    CommandElecInitialMagnetization() : Command("elec-initial-magnetization", "jdftx/Initialization")
	{
		format = "<M> <constrain>=yes|no";
		comments =
			"Initialize system with total magnetization <M> (= Nup - Ndn).\n"
			"With Fermi fillings, the magnetization will be constrained if <constrain>=yes,\n"
			"and will equilibriate otherwise. Without Fermi fillings, the magnetization will\n"
			"remain constrained regardless. Only valid for spintype z-spin.";
		
		require("spintype");
	}

	void process(ParamList& pl, Everything& e)
	{	if(e.eInfo.spinType != SpinZ)
			throw("Total magnetization can only be specified for spintype z-spin");
		pl.get(e.eInfo.Minitial, 0., "M", true);
		pl.get(e.eInfo.Mconstrain, true, boolMap, "constrain", true);
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lf %s", e.eInfo.Minitial, boolMap.getString(e.eInfo.Mconstrain));
	}
}
CommandElecInitialMagnetization;

//-------------------------------------------------------------------------------------------------

struct CommandElecInitialFillings : public Command
{
	CommandElecInitialFillings() : Command("elec-initial-fillings", "jdftx/Initialization")
	{
		format = "read <filename> [<nBandsOld>]";
		comments =
			"Initial electronic fillings are read from <filename> instead of computed automatically (default).\n"
			"+ <nBandsOld> specifies the number of bands in the file being read, if different from current nBands.\n"
			"  Extra new bands are filled with 0s, fillings of extra old bands are accumulated into other bands.\n"
			"  Default: 0 => fillings file has same number of bands as this run.";
		
		forbid("initial-state");
	}

	void process(ParamList& pl, Everything& e)
	{	string key; pl.get(key, string(), "read", true);
		if(key!=string("read")) throw "First parameter must be 'read', encountered " + key;
		pl.get(e.eInfo.initialFillingsFilename, string(), "filename", true);
		pl.get(e.eInfo.nBandsOld, 0, "nBandsOld");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("read %s %d", e.eInfo.initialFillingsFilename.c_str(), e.eInfo.nBandsOld);
	}
}
commandElecInitialFillings;

//-------------------------------------------------------------------------------------------------

struct CommandElecInitialEigenvals : public Command
{
	CommandElecInitialEigenvals() : Command("elec-initial-eigenvals", "jdftx/Initialization")
	{
		format = "<filename>";
		comments = "Read the initial eigenvalues for variable fillings (default: derive from subspace hamiltonian)\n";
		
		forbid("initial-state");
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.eVars.eigsFilename, string(), "filename", true);
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s", e.eVars.eigsFilename.c_str());
	}
}
commandElecInitialEigenvals;

//-------------------------------------------------------------------------------------------------

struct CommandSubspaceRotationFactor : public Command
{
	CommandSubspaceRotationFactor() : Command("subspace-rotation-factor", "jdftx/Electronic/Optimization")
	{	format = "<factor> <adjust>=yes|no";
		comments = 
			"Preconditioning factor for subspace rotations generated by the\n"
			"auxiliary Hamiltonian used for minimization with variable fillings.\n"
			"With <adjust> = yes (default), the factor is heuristically adjusted\n"
			"during the minimization to optimize convergence.";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.cntrl.subspaceRotationFactor, 30.0, "factor");
		pl.get(e.cntrl.subspaceRotationAdjust, true, boolMap, "adjust");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg %s", e.cntrl.subspaceRotationFactor, boolMap.getString(e.cntrl.subspaceRotationAdjust));
	}
}
commandSubspaceRotationFactor;
