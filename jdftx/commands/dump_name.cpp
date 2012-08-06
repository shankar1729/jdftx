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


struct CommandDumpName : public Command
{
	CommandDumpName() : Command("dump-name")
	{
		format = "<format>";
		comments = 
			"  Control the filename pattern for dump output:\n"
			"    <format> is an arbitrary format string that will be substituted according to:\n"
			"       $VAR -> name of the variable being dumped (this must be present somewhere in the string)\n"
			"       $STAMP -> time-stamp at the start of dump";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.dump.format, string("$STAMP.$VAR"), "format");
		if(e.dump.format.find("$VAR")==string::npos)
			throw "<format> = " + e.dump.format + " doesn't contain the pattern $VAR";
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s", e.dump.format.c_str());
	}
}
commandDumpName;
