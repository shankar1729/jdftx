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

//! @file elec_misc.cpp Miscellaneous properties of the electronic system

struct CommandElecCutoff : public Command
{
	CommandElecCutoff() : Command("elec-cutoff")
	{
		format = "<Ecut>";
		comments = "Electronic planewave cutoff in Hartree";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.cntrl.Ecut, 20.0, "Ecut");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg", e.cntrl.Ecut);
	}
}
commandElecCutoff;

//-------------------------------------------------------------------------------------------------

struct CommandElecNbands : public Command
{
	CommandElecNbands() : Command("elec-n-bands")
	{
		format = "<n>";
		comments = "Manually specify the number of bands (Default: set nBands assuming insulator)";
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.eInfo.nBands, 0, "n", true);
		if(e.eInfo.nBands<=0) throw string("<n> must be positive.\n");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%d", e.eInfo.nBands);
	}
}
commandElecNbands;

//-------------------------------------------------------------------------------------------------

EnumStringMap<SpinType> spinMap
(	SpinNone, "no-spin",
	SpinZ, "z-spin"
);

struct CommandSpinType : public Command
{
	CommandSpinType() : Command("spintype")
	{
		format = "<type>=" + spinMap.optionList();
		comments = "Select spin-polarization type";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.eInfo.spinType, SpinNone, spinMap, "type");
	}

	void printStatus(Everything& e, int iRep)
	{	fputs(spinMap.getString(e.eInfo.spinType), globalLog);
	}
}
commandSpinType;

//-------------------------------------------------------------------------------------------------

struct CommandFixElectronDensity : public Command
{
	CommandFixElectronDensity() : Command("fix-electron-density")
	{
		format = "<nFilename> | <nupFilename> <ndnFilename>";
		comments = "Perform band structure calculations at fixed electron density\n"
			"(or spin densities) read from the specified file(s)";
		
		require("spintype");
		forbid("elec-fermi-fillings");
		forbid("elec-ex-corr-compare");
		forbid("residual-minimize");
	}

	void process(ParamList& pl, Everything& e)
	{	std::vector<string>& nFilename = e.eVars.nFilename;
		nFilename.resize(e.eInfo.spinType==SpinNone ? 1 : 2);
		for(unsigned s=0; s<nFilename.size(); s++)
		{	const char* componentName = nFilename.size()==1 ? "nFilename" : (s==0 ? "nupFilename" : "ndnFilename");
			pl.get(nFilename[s], string(), componentName, true);
		}
		e.cntrl.fixed_n = true;
	}

	void printStatus(Everything& e, int iRep)
	{	std::vector<string>& nFilename = e.eVars.nFilename;
		for(unsigned s=0; s<nFilename.size(); s++)
			logPrintf("%s ", e.eVars.nFilename[s].c_str());
	}
}
commandFixElectronDensity;

//-------------------------------------------------------------------------------------------------

struct CommandFixOccupied : public Command
{
	CommandFixOccupied() : Command("fix-occupied")
	{
		format = "[<fThreshold>=0]";
		comments = "Fix orbitals with fillings larger than <fThreshold> in band-structure calculations\n"
			"The occupied orbitals must be read in using the wavefunction / initial-state commands.\n";
		
		require("fix-electron-density");
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.cntrl.occupiedThreshold, 0., "fThreshold");
		if(e.cntrl.occupiedThreshold<0) throw string("fThreshold must be >= 0");
		e.cntrl.fixOccupied = true;
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg", e.cntrl.occupiedThreshold);
	}
}
commandFixOccupied;

//-------------------------------------------------------------------------------------------------

struct CommandReorthogonalizeOrbitals : public Command
{
	CommandReorthogonalizeOrbitals() : Command("reorthogonalize-orbitals")
	{
		format = "[<interval=20> [<threshold>=1.5]";
		comments =
			"Every <interval> electronic steps, re-orthogonalize analytically-continued\n"
			"orbitals if the condition number of their overlap matrix crosses <threshold>.\n"
			"Set <interval> = 0 to disable this check.";
		
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.cntrl.overlapCheckInterval, 20, "interval");
		pl.get(e.cntrl.overlapConditionThreshold, 1.5, "threshold");
		if(e.cntrl.overlapCheckInterval<0) throw string("<interval> must be non-negative");
		if(e.cntrl.overlapConditionThreshold<=1.) throw string("<threshold> must be > 1");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%d %lg", e.cntrl.overlapCheckInterval, e.cntrl.overlapConditionThreshold);
	}
}
commandReorthogonalizeOrbitals;
