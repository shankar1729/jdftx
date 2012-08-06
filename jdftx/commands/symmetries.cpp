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

EnumStringMap<SymmetryMode> symmMap(
	SymmetriesNone, "none",
	SymmetriesAutomatic, "automatic",
	SymmetriesManual, "manual" );

struct CommandSymmetries : public Command
{
	CommandSymmetries() : Command("symmetries")
	{
		format = "<symm>=" + symmMap.optionList() + " <moveAtoms>=" + boolMap.optionList();
		comments = "\tnone: symmetries are off\n"
			"\tautomatic: automatic calculation of symmetries (default)\n"
			"\tmanual: symmetries specified using symmetry-matrix command\n"
			"Check if symmetries could be increased by translating all the atoms\n"
			"and quit after suggesting the translated positions if <moveAtoms>=yes.";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.symm.mode, SymmetriesAutomatic, symmMap, "symm");
		pl.get(e.symm.shouldMoveAtoms, false, boolMap, "moveAtoms");
	}

	void printStatus(Everything& e, int iRep)
	{	fputs(symmMap.getString(e.symm.mode), globalLog);
	}
}
commandSymmetries;
