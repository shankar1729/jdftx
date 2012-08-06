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

struct CommandKpointFolding : public Command
{
	CommandKpointFolding() : Command("kpoint-folding")
	{
		format = "[<n0>=1] [<n1>=1] [<n2>=1]";
		comments = "Fold k-points in direction i by factor <ni> (for i=0,1,2)";
		hasDefault = true;

		require("kpoint");
	}

	void process(ParamList& pl, Everything& e)
	{	for(int k=0; k<3; k++)
		{	ostringstream oss; oss << "n" << k;
			string paramName = oss.str();

			pl.get(e.eInfo.kfold[k], 1, paramName);
			if(e.eInfo.kfold[k]<1)
				logPrintf("<%s> must be a positive integer.", paramName.c_str());
		}
	}

	void printStatus(Everything& e, int iRep)
	{	for(int k=0; k<3; k++) logPrintf("%d ", e.eInfo.kfold[k]);
	}
}
commandKpointFolding;
