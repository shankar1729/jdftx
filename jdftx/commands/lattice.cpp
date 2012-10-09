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

struct CommandLattice : public Command
{
	CommandLattice() : Command("lattice")
	{
		format = " \\\n\t<R00> <R01> <R02> \\\n\t<R10> <R11> <R12> \\\n\t<R20> <R21> <R22>";
		comments = "Specify lattice vectors (in columns)";
	}

	void process(ParamList& pl, Everything& e)
	{	for(int j=0; j<3; j++) for(int k=0; k<3; k++)
		{	ostringstream oss; oss << "R" << j << k;
			pl.get(e.gInfo.R(j,k), 0.0, oss.str(), true);
		}
	}

	void printStatus(Everything& e, int iRep)
	{	for (int j=0; j < 3; j++)
		{	logPrintf(" \\\n\t");
			for (int k=0; k < 3; k++) logPrintf("%20.15le ",e.gInfo.R(j,k));
		}
	}
}
commandLattice;


struct CommandLattScale : public Command
{
	CommandLattScale() : Command("latt-scale")
	{
		format = "<s0> <s1> <s2>";
		comments = "Scale lattice vector i by factor <si>";
		hasDefault = true;

		require("lattice");
	}

	void process(ParamList& pl, Everything& e)
	{	for(int k=0; k<3; k++)
		{	ostringstream oss; oss << "s" << k;
			double scale;
			pl.get(scale, 1.0, oss.str());
			e.gInfo.R.set_col(k, scale * e.gInfo.R.column(k));
		}
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("1 1 1  # Scale has been absorbed into lattice");
	}
}
commandLattScale;


struct CommandLattMoveScale : public Command
{
	CommandLattMoveScale() : Command("latt-move-scale")
	{
		format = "<s0> <s1> <s2>";
		comments = "Preconditioning factor for each lattice vector (must be commensurate with symmetries)";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	vector3<>& s = e.cntrl.lattMoveScale;
		pl.get(s[0], 1., "s0");
		pl.get(s[1], 1., "s1");
		pl.get(s[2], 1., "s2");
	}

	void printStatus(Everything& e, int iRep)
	{	const vector3<>& s = e.cntrl.lattMoveScale;
		logPrintf("%lg %lg %lg", s[0], s[1], s[2]);
	}
}
commandLattMoveScale;
