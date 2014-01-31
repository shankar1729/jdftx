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

struct CommandInitialState : public Command
{
	CommandInitialState() : Command("initial-state")
	{
		format = "<filename-pattern>";
		comments = "Initialize state from a filename pattern which contains a $VAR,\n"
			"equivalent to invoking the following commands:\n"
			"   wavefunction          read  <filename-pattern>/$VAR/wfns\n"
			"   elec-initial-fillings read  <filename-pattern>/$VAR/fillings\n"
			"   elec-initial-Haux           <filename-pattern>/$VAR/Haux\n"
			"   fluid-initial-state         <filename-pattern>/$VAR/fluidState\n"
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
		FILE* fp = fopen(filename.c_str(), "r");
		if(fp) target = filename; //file exists and is readable
	}
}
commandInitialState;
