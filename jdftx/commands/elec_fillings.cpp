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

EnumStringMap<ElecInfo::SmearingType> smearingTypeMap(
	ElecInfo::SmearingFermi, "Fermi",
	ElecInfo::SmearingGauss, "Gauss",
	ElecInfo::SmearingCold, "Cold" );

EnumStringMap<ElecInfo::SmearingType> smearingTypeDescMap(
	ElecInfo::SmearingFermi, "Use a Fermi-Dirac function for fillings",
	ElecInfo::SmearingGauss, "Use a gaussian-based (erfc) function for fillings",
	ElecInfo::SmearingCold, "Use the cold smearing function \\cite ColdSmearing to approximate zero temperature" );

struct CommandElecSmearing : public Command
{
	CommandElecSmearing() : Command("elec-smearing", "jdftx/Electronic/Parameters")
	{
		format = "<smearingType>="+smearingTypeMap.optionList()+" <smearingWidth>";
		comments =
			"Use variable electronic fillings using a smearing function selected by <smearingType>:"
			+ addDescriptions(smearingTypeMap.optionList(), linkDescription(smearingTypeMap, smearingTypeDescMap))
			+ "\n\nwith width set by <smearingWidth> in Hartrees.\n"
			"The width corresponds to kT (electronic temperature) for Fermi smearing,\n"
			"and sigma/2 for the Gauss and Cold smearing options: this convention\n"
			"results in roughly the same rate of k-point convergence for all three\n"
			"methods using the same width. However, the entropy contribution at the\n"
			"same width will follow the order Fermi > Gauss >> Cold.";
		
		require("lcao-params");
	}

	void process(ParamList& pl, Everything& e)
	{	ElecInfo& eInfo = e.eInfo;
		pl.get(eInfo.smearingType, ElecInfo::SmearingFermi, smearingTypeMap, "smearingType", true);
		pl.get(eInfo.smearingWidth, 0.0, "smearingWidth", true);
		if(e.eInfo.smearingWidth<=0) throw string("<smearingWidth> must be positive.\n");
		eInfo.fillingsUpdate = ElecInfo::FillingsHsub;
	}

	void printStatus(Everything& e, int iRep)
	{	const ElecInfo& eInfo = e.eInfo;
		logPrintf("%s %lg", smearingTypeMap.getString(eInfo.smearingType), eInfo.smearingWidth);
	}
}
commandElecFermiFillings;

struct DeprecatedCommandElecFermiFillings : public DeprecatedCommand
{	DeprecatedCommandElecFermiFillings() : DeprecatedCommand("elec-fermi-fillings") { }
	
	std::pair<string,string> replace(ParamList& pl) const
	{	string unusedParam; double kT;
		pl.get(unusedParam, string(), "unusedParam", true);
		pl.get(kT, 0.0, "kT", true);
		ostringstream oss; oss << " Fermi " << kT;
		return std::make_pair(string("elec-smearing"), oss.str());
	}
}
deprecatedCommandElecFermiFillings;

//-------------------------------------------------------------------------------------------------

struct CommandTargetMu : public Command
{
	CommandTargetMu() : Command("target-mu", "jdftx/Electronic/Parameters")
	{
		format = "<mu> [<outerLoop>=no]";
		comments =
			"Fixed chemical potential <mu> (instead of fixed charge).\n"
			"Note that <mu> is absolute (relative to vacuum level) and in Hartrees.\n"
			"For example, potential V (in Volts) relative to SHE corresponds to\n"
			"mu = -(Vref + V)/27.2114, where Vref is the absolute SHE potential\n"
			"in Volts below vacuum; you could set Vref = 4.44 based on experiment\n"
			"or use the value calibrated using potentials of zero charge with\n"
			"the solvation model in use.\n"
			"\n"
			"The default setting <outerLoop>=no directly performs variational minimization\n"
			"or SCF in the grand canonical ensemble: keeping mu fixed throughout, letting\n"
			"the number of electrons adjust continuously \\cite GC-DFT.\n"
			"\n"
			"Setting <outerLoop>=yes instead performs a sequence of conventional fixed-charge\n"
			"optimizations, adjusting mu in an outer loop using the secant method.\n"
			"This is usually much slower, and is only recommended if the default\n"
			"direct grand canonical method fails.";

		require("fluid-cation");
		require("fluid-anion");
		require("elec-smearing");
		forbid("elec-initial-charge");
		forbid("fix-electron-density");
		forbid("fix-electron-potential");
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.eInfo.mu, 0.0, "mu", true);
		pl.get(e.eInfo.muLoop, false, boolMap, "outerLoop");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg %s\n", e.eInfo.mu, boolMap.getString(e.eInfo.muLoop));
	}
}
commandTargetMu;

//-------------------------------------------------------------------------------------------------

struct CommandTargetBz : public Command
{
	CommandTargetBz() : Command("target-Bz", "jdftx/Electronic/Parameters")
	{
		format = "<Bz>";
		comments =
			"Fixed magnetic field <Bz> (instead of magnetization).\n"
			"Note that <Bz> is in atomic units (1 T is approximately 4.26E-6 a.u.)\n"
			"and in an electron-is-positive convention (Bz > 0 favors up spins).\n"
			"Requires smearing and is only valid for spintype z-spin.";

		require("elec-smearing");
		require("elec-initial-magnetization");
		forbid("fix-electron-density");
		forbid("fix-electron-potential");
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.eInfo.Bz, 0., "Bz", true);
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg\n", e.eInfo.Bz);
	}
}
commandTargetBz;

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
			"remain constrained regardless. Note: target-Bz will override <constrain> to no.\n"
			"Only valid for spintype z-spin.";
		
		require("spintype");
	}

	void process(ParamList& pl, Everything& e)
	{	if(e.eInfo.spinType != SpinZ)
			throw(string("Total magnetization can only be specified for spintype z-spin"));
		pl.get(e.eInfo.Minitial, 0., "M", true);
		//Constraint parameter:
		bool Mconstrain = true;
		pl.get(Mconstrain, true, boolMap, "constrain", true);
		e.eInfo.Bz = Mconstrain ? NAN : 0.; //nan B signifies M constrained, M free otherwise
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lf %s", e.eInfo.Minitial, boolMap.getString(std::isnan(e.eInfo.Bz)));
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
	{	pl.get(e.cntrl.subspaceRotationFactor, 1., "factor");
		pl.get(e.cntrl.subspaceRotationAdjust, true, boolMap, "adjust");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg %s", e.cntrl.subspaceRotationFactor, boolMap.getString(e.cntrl.subspaceRotationAdjust));
	}
}
commandSubspaceRotationFactor;
