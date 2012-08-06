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

enum WfnsInit { WfnsLCAO, WfnsRandom, WfnsRead, WfnsReadRS };

EnumStringMap<WfnsInit> wfnsInitMap(
	WfnsLCAO, "lcao",
	WfnsRandom, "random",
	WfnsRead, "read",
	WfnsReadRS, "read-rs" );

struct CommandWavefunction : public Command
{
	CommandWavefunction() : Command("wavefunction")
	{
		format =
			"lcao\n"
			"           | random\n"
			"           | read <filename> [<nBandsOld>] [<EcutOld>] [<kdepOld>]\n"
			"           | read-rs <filename-pattern> [<nBandsOld>] [<NxOld>] [<NyOld>] [<NzOld>]";
		comments =
			"Wavefunction initialization: use atomic orbitals (default), randomize or read from files:\n"
			"  read expects <filename> to point to a single file with fourier-space G-sphere wavefunctions.\n"
			"  read-rs expects <filename> to be a printf format with 2 %%d's, the first for state index and\n"
			"     the second for band. Each 'column' will be loaded from a separate file accordingly.\n"
			"  <nBandsOld> can be used to specify a wavefunction which has different bands\n"
			"     extra bands will be discarded, unspecified bands will be randomized and orthogonalized\n"
			"     Reminder: nBandsOlds for fillings file is specified separately in elec-initial-fillings\n"
			"     default: 0 => old and current nBands must match exactly.\n"
			"  <EcutOld> can be used to specify a wavefunction with different planewave cutoff\n"
			"     the wavefunction will be appropriately up/down-sampled in Foruier space\n"
			"     default: 0.0 => old and current Ecut must match exactly\n"
			"  <kdepOld>="+kdepMap.optionList()+" specifies whether the input wavefunction has a k-point dependent basis\n"
			"     default: kpoint-dependent\n"
			"  <N*old> specifies fftbox size of the input data when reading real-space wavefunctions\n"
			"     the wavefunction will be appropriately up/down-sampled in Foruier space\n"
			"     default: 0 => old and current fftbox must match exactly";
		hasDefault = false;
		
		forbid("initial-state");
	}

	void process(ParamList& pl, Everything& e)
	{	WfnsInit wfnsInit; pl.get(wfnsInit, WfnsRandom, wfnsInitMap, "option", true);
		switch(wfnsInit)
		{	case WfnsLCAO:
				//ElecVars constructor defaults to initLCAO=true with no wfnsFilename
				break;
			case WfnsRandom:
				e.eVars.initLCAO = false;
				break;
			case WfnsRead:
				e.eVars.readWfnsRealspace = false;
				pl.get(e.eVars.wfnsFilename, string(), "filename", true);
				pl.get(e.eVars.nBandsOld, 0, "nBandsOld");
				pl.get(e.eVars.EcutOld, 0.0, "EcutOld");
				pl.get(e.eVars.kdepOld, BasisKpointDep, kdepMap, "kdepOld");
				break;
			case WfnsReadRS:
				e.eVars.readWfnsRealspace = true;
				pl.get(e.eVars.wfnsFilename, string(), "filename", true);
				pl.get(e.eVars.nBandsOld, 0, "nBandsOld");
				pl.get(e.eVars.NxOld, 0, "NxOld");
				pl.get(e.eVars.NyOld, 0, "NyOld");
				pl.get(e.eVars.NzOld, 0, "NzOld");
				break;
		}
	}

	void printStatus(Everything& e, int iRep)
	{	if(!e.eVars.wfnsFilename.length())
			logPrintf(e.eVars.initLCAO ? "lcao" : "random");
		else if(!e.eVars.readWfnsRealspace)
			logPrintf("read %s %d %lf %s", e.eVars.wfnsFilename.c_str(), e.eVars.nBandsOld,
				e.eVars.EcutOld, kdepMap.getString(e.eVars.kdepOld));
		else
			logPrintf("read-rs %s %d %d %d %d", e.eVars.wfnsFilename.c_str(), e.eVars.nBandsOld,
				e.eVars.NxOld, e.eVars.NyOld, e.eVars.NzOld);
	}
}
commandWavefunction;
