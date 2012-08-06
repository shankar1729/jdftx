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

struct CommandSymmetryMatrix : public Command
{
	CommandSymmetryMatrix() : Command("symmetry-matrix")
	{
		format = "<s00> <s01> <s02> <s10> <s11> <s12> <s20> <s21> <s22>";
		comments = "Specify symmetry operator matrices explicitly.\n"
			"Requires symmetries command to be called with manual argument";
		allowMultiple = true;

		require("symmetries");
	}

	void process(ParamList& pl, Everything& e)
	{	if(e.symm.mode != SymmetriesManual)
			throw string("symmetry-matrix needs symmetries to be called with \"manual\"");

		matrix3<int> m;
		for(int j=0; j<3; j++) for(int k=0; k<3; k++)
		{	ostringstream oss; oss << "s" << j << k;
			pl.get(m(j,k), 0, oss.str(), true);
		}
		e.symm.sym.push_back(m);
	}

	void printStatus(Everything& e, int iRep)
	{	for (int j=0; j < 3; j++)
		{	logPrintf(" \\\n\t");
			for (int k=0; k < 3; k++) logPrintf("%d ",e.symm.sym[iRep](j,k));
		}
	}
}
commandSymmetryMatrix;
