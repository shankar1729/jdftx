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
	CommandElecFermiFillings() : Command("elec-fermi-fillings")
	{
		format = "<mixInterval> <kT> [<alpha>=0.5]";
		comments =
			"Fermi-Dirac fillings at a temperature <kT> (in Hartrees).\n"
			"  If <mixInterval> is zero, use the auxilliary hamiltonian method (recommended).\n"
			"  Else, mix fermi functions every <mixInterval> iterations with mixing fraction <alpha>.";
		
		require("lcao-params");
		forbid("fix-electron-density");
		forbid("fix-electron-potential");
	}

	void process(ParamList& pl, Everything& e)
	{	ElecInfo& eInfo = e.eInfo;
		pl.get(eInfo.mixInterval, 0, "mixInterval", true);
		//Determine algorithm based on mixInterval:
		if(eInfo.mixInterval<0) throw string("<mixInterval> must be positive");
		eInfo.fillingsUpdate = eInfo.mixInterval ? ElecInfo::FermiFillingsMix : ElecInfo::FermiFillingsAux;
		pl.get(eInfo.kT, 0.0, "kT", true);
		if(eInfo.mixInterval) pl.get(eInfo.fillingMixFraction, 0.5, "alpha");
	}

	void printStatus(Everything& e, int iRep)
	{	const ElecInfo& eInfo = e.eInfo;
		logPrintf("%d %lg", eInfo.mixInterval, eInfo.kT);
		if(eInfo.mixInterval) logPrintf(" %lg", eInfo.fillingMixFraction);
	}
}
commandElecFermiFillings;

//-------------------------------------------------------------------------------------------------

struct CommandTargetMu : public Command
{
	CommandTargetMu() : Command("target-mu")
	{
		format = "<mu> [<Cinitial>=1.0] [<dnMix>=0.7]";
		comments =
			"Fixed chemical potential <mu> (instead of fixed charge)\n"
			"When using elec-fermi-fillings with non-zero mixInterval (deprecated),\n"
			"the following parameters control the convergence:\n"
			"+ Cinitial: Initial capacitance, affects only first few mixing steps\n"
			"+ dnMix: Scale the ideal step in n by this factor";
		hasDefault = false;

		require("fluid-cation");
		require("fluid-anion");
		require("elec-fermi-fillings");
		forbid("elec-initial-charge");
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.eInfo.mu, 0.0, "mu", true);
		pl.get(e.eInfo.Cmeasured, 1.0, "Cinitial");
		pl.get(e.eInfo.dnMixFraction, 0.7, "dnMix");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg %lg %lg\n", e.eInfo.mu, e.eInfo.Cmeasured, e.eInfo.dnMixFraction);
	}
}
commandTargetMu;

//-------------------------------------------------------------------------------------------------

struct CommandElecInitialCharge : public Command
{
    CommandElecInitialCharge() : Command("elec-initial-charge")
	{
		format = "<QNet> | <QupNet> <QdnNet>";
		comments =
			"Initialize a system with <QNet> excess unpolarized electrons\n"
			"or with <QupNet> and <QdnNet> excess electrons in each spin channel.\n"
			"Both versions are valid irrespective of whether system is spin-polarized.";

		forbid("target-mu");
	}

	void process(ParamList& pl, Everything& e)
	{	double qTemp;
		e.eInfo.qNet.clear();
		//First q (required):
		pl.get(qTemp, 0., "QNet", true);
		e.eInfo.qNet.push_back(qTemp);
		//Second q (optional):
		pl.get(qTemp, std::numeric_limits<double>::quiet_NaN(), "QNet");
		if(!std::isnan(qTemp)) e.eInfo.qNet.push_back(qTemp);
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%le", e.eInfo.qNet[0]);
		if(e.eInfo.qNet.size()==2)
			logPrintf(" %le", e.eInfo.qNet[1]);
	}
}
CommandElecInitialCharge;

//-------------------------------------------------------------------------------------------------

struct CommandElecInitialFillings : public Command
{
	CommandElecInitialFillings() : Command("elec-initial-fillings")
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

struct CommandElecInitialHaux : public Command
{
	CommandElecInitialHaux() : Command("elec-initial-Haux")
	{
		format = "<filename>";
		comments = "Read the auxilliary hamiltonian for direct fillings (default: set to subspace hamiltonian)\n";
		
		forbid("initial-state");
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.eVars.HauxFilename, string(), "filename", true);
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s", e.eVars.HauxFilename.c_str());
	}
}
commandElecInitialHaux;

//-------------------------------------------------------------------------------------------------

struct CommandSubspaceRotationFactor : public Command
{
	CommandSubspaceRotationFactor() : Command("subspace-rotation-factor")
	{	format = "<factor>";
		comments = "preconditioning factor for subspace rotations";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.eVars.subspaceRotationFactor, 30.0, "factor");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg", e.eVars.subspaceRotationFactor);
	}
}
commandSubspaceRotationFactor;
