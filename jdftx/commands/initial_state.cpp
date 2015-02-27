/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver

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

struct CommandInitialState : public Command
{
	CommandInitialState() : Command("initial-state")
	{
		format = "<filename-pattern>";
		comments = "Initialize state from a filename pattern which contains a $VAR,\n"
			"equivalent to invoking the following commands:\n"
			"+ wavefunction          read  <filename-pattern>/$VAR/wfns\n"
			"+ elec-initial-fillings read  <filename-pattern>/$VAR/fillings\n"
			"+ elec-initial-Haux           <filename-pattern>/$VAR/Haux\n"
			"+ fluid-initial-state         <filename-pattern>/$VAR/fluidState\n"
			"\n"
			"(where A/x/y is sed for 'find x in A and replace it with y'.)\n"
			"This command will invoke the read only for those state variables for which\n"
			"the corresponding files exist, leaving the rest with default initialization.\n"
			"When using SCF, this will also read scfHistory and eigenvalues if available.";
		
		forbid("wavefunction");
		forbid("elec-initial-fillings");
		forbid("elec-initial-Haux");
		forbid("fluid-initial-state");
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(filenamePattern, string(), "filename-pattern", true);
		if(filenamePattern.find("$VAR")==string::npos)
			throw "<filename-pattern> = " + filenamePattern + " doesn't contain '$VAR'";
		processPattern("wfns", e.eVars.wfnsFilename);
		processPattern("fillings", e.eInfo.initialFillingsFilename);
		if(!e.eInfo.initialFillingsFilename.length())
			processPattern("fill", e.eInfo.initialFillingsFilename); //alternate naming convention
		processPattern("Haux", e.eVars.HauxFilename);
		processPattern("fluidState", e.eVars.fluidInitialStateFilename);
		if(!e.eVars.fluidInitialStateFilename.length())
			processPattern("fS", e.eVars.fluidInitialStateFilename); //alternate naming convention
		processPattern("scfHistory", e.scfParams.historyFilename);
		processPattern("eigenvals", e.eVars.eigsFilename);
	}

	void printStatus(Everything& e, int iRep)
	{	fputs(filenamePattern.c_str(), globalLog);
	}

private:
	string filenamePattern;
	
	void processPattern(string varName, string& target) const
	{	string filename = filenamePattern;
		filename.replace(filename.find("$VAR"),4, varName);
		if(isReadable(filename))
			target = filename; //file exists and is readable
	}
}
commandInitialState;

//-----------------------------------------------------------------------

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
			"+ read expects <filename> to point to a single file with fourier-space G-sphere wavefunctions.\n"
			"+ read-rs expects <filename> to be a printf format with 2 %%d's, the first for state index and\n"
			"   the second for band. Each 'column' will be loaded from a separate file accordingly.\n"
			"   For spinor wavefunctions, each spinor component has a separate second index, so that\n"
			"   the first band is read from 0 and 1, the second one from 2 and 3 and so on.\n"
			"+ <nBandsOld> can be used to specify a wavefunction which has different bands\n"
			"   extra bands will be discarded, unspecified bands will be randomized and orthogonalized.\n"
			"   Reminder: nBandsOlds for fillings file is specified separately in elec-initial-fillings.\n"
			"   Default: 0 => old and current nBands must match exactly.\n"
			"+ <EcutOld> can be used to specify a wavefunction with different planewave cutoff.\n"
			"   The wavefunction will be appropriately up/down-sampled in Fourier space.\n"
			"   Default: 0.0 => old and current Ecut must match exactly.\n"
			"+ <N*old> specify fftbox dimensions of the input data when reading real-space wavefunctions.\n"
			"   The wavefunction will be appropriately up/down-sampled in Fourier space.\n"
			"   Default: 0 => old and current fftbox must match exactly.";
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
