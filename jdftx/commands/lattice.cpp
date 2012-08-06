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

