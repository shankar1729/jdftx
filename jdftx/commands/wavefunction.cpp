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
#include <electronic/ColumnBundle.h>

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
			"           | read <filename> [<nBandsOld>] [<EcutOld>]\n"
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
			"  <N*old> specifies fftbox size of the input data when reading real-space wavefunctions\n"
			"     the wavefunction will be appropriately up/down-sampled in Fourier space\n"
			"     default: 0 => old and current fftbox must match exactly";
		hasDefault = false;
		
		forbid("initial-state");
	}

	void process(ParamList& pl, Everything& e)
	{	WfnsInit wfnsInit; pl.get(wfnsInit, WfnsLCAO, wfnsInitMap, "option", true);
		switch(wfnsInit)
		{	case WfnsLCAO:
				e.eVars.initLCAO = true;
				break;
			case WfnsRandom:
				//ElecVars constructor defaults to initLCAO=false with no wfnsFilename
				e.eVars.initLCAO = false;
				break;
			case WfnsRead:
			{	pl.get(e.eVars.wfnsFilename, string(), "filename", true);
				auto conversion = std::make_shared<ColumnBundleReadConversion>();
				conversion-> realSpace = false;
				pl.get(conversion->nBandsOld, -1, "nBandsOld");
				pl.get(conversion->EcutOld, 0.0, "EcutOld");
				if(conversion->nBandsOld>=0) e.eVars.readConversion = conversion;
				break;
			}
			case WfnsReadRS:
			{	pl.get(e.eVars.wfnsFilename, string(), "filename", true);
				auto conversion = std::make_shared<ColumnBundleReadConversion>();
				conversion-> realSpace = true;
				pl.get(conversion->nBandsOld, 0, "nBandsOld");
				pl.get(conversion->S_old[0], 0, "NxOld");
				pl.get(conversion->S_old[1], 0, "NyOld");
				pl.get(conversion->S_old[2], 0, "NzOld");
				e.eVars.readConversion = conversion;
				break;
			}
		}
	}

	void printStatus(Everything& e, int iRep)
	{	if(!e.eVars.wfnsFilename.length())
			logPrintf(e.eVars.initLCAO ? "lcao" : "random");
		else if(!e.eVars.readConversion)
			logPrintf("read %s", e.eVars.wfnsFilename.c_str());
		else if(!e.eVars.readConversion->realSpace)
			logPrintf("read %s %d %lf", e.eVars.wfnsFilename.c_str(),
				e.eVars.readConversion->nBandsOld, e.eVars.readConversion->EcutOld);
		else
			logPrintf("read-rs %s %d %d %d %d", e.eVars.wfnsFilename.c_str(), e.eVars.readConversion->nBandsOld,
				e.eVars.readConversion->S_old[0], e.eVars.readConversion->S_old[1], e.eVars.readConversion->S_old[2]);
	}
}
commandWavefunction;
