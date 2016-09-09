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
#include <core/LatticeUtils.h>

EnumStringMap<SymmetryMode> symmMap(
	SymmetriesNone, "none",
	SymmetriesAutomatic, "automatic",
	SymmetriesManual, "manual" );

struct CommandSymmetries : public Command
{
	CommandSymmetries() : Command("symmetries", "jdftx/Miscellaneous")
	{
		format = "<symm>=" + symmMap.optionList();
		comments = "+ none: symmetries are off\n"
			"+ automatic: automatic calculation of symmetries (default)\n"
			"+ manual: symmetries specified using symmetry-matrix command.";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.symm.mode, SymmetriesAutomatic, symmMap, "symm");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s", symmMap.getString(e.symm.mode));
	}
}
commandSymmetries;


struct CommandSymmetryMatrix : public Command
{
	CommandSymmetryMatrix() : Command("symmetry-matrix", "jdftx/Miscellaneous")
	{
		format = " \\\n\t<s00> <s01> <s02> \\\n\t<s10> <s11> <s12> \\\n\t<s20> <s21> <s22> \\\n\t<a0> <a1> <a2>";
		comments = "Specify symmetry operator matrices explicitly. The top 3 x 3 block\n"
			"contains the integer rotation matrix in lattice coordinates, while the\n"
			"final row contains the subsequent translation in lattice coordinates.\n"
			"Requires symmetries command to be called with manual argument.";
		allowMultiple = true;

		require("symmetries");
	}

	void process(ParamList& pl, Everything& e)
	{	if(e.symm.mode != SymmetriesManual)
			throw string("symmetry-matrix needs symmetries to be called with \"manual\"");

		SpaceGroupOp op;
		for(int j=0; j<3; j++) for(int k=0; k<3; k++)
		{	ostringstream oss; oss << "s" << j << k;
			pl.get(op.rot(j,k), 0, oss.str(), true);
		}
		for(int k=0; k<3; k++)
		{	ostringstream oss; oss << "a" << k;
			pl.get(op.a[k], 0., oss.str(), true);
		}
		e.symm.sym.push_back(op);
	}

	void printStatus(Everything& e, int iRep)
	{	for (int j=0; j < 3; j++)
		{	logPrintf(" \\\n\t");
			for(int k=0; k<3; k++) logPrintf("%d ",e.symm.sym[iRep].rot(j,k));
		}
		logPrintf(" \\\n\t");
		for(int k=0; k<3; k++) logPrintf("%lg ",e.symm.sym[iRep].a[k]);
	}
}
commandSymmetryMatrix;


struct CommandSymmetryThreshold : public Command
{
	CommandSymmetryThreshold() : Command("symmetry-threshold", "jdftx/Miscellaneous")
	{
		format = "[<threshold>=1e-4]";
		comments = "Relative threshold parameter used for symmetry detection.";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(symmThreshold, 1e-4, "threshold");
		symmThresholdSq = symmThreshold*symmThreshold;
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg", symmThreshold);
	}
}
commandSymmetryThreshold;

