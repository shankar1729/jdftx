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
