/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

//An enum entry for each configurable option of NonlocalPCMparams
enum NonlocalPCMparamsMember
{	NPM_lMax,
	NPM_Delim //!< delimiter to detect end of input
};

EnumStringMap<NonlocalPCMparamsMember> npmMap
(	NPM_lMax, "lMax"
);
EnumStringMap<NonlocalPCMparamsMember> npmDescMap
(	NPM_lMax, "maximum angular momentum to include in response"
);


//! Minimization parameters for ions, electrons or fluid
struct CommandNonlocalPCM : public Command
{
	CommandNonlocalPCM() : Command("nonlocal-pcm")
	{
		format = "<key1> <value1> <key2> <value2> ...";
		comments = "where possible keys and value types are:"
			+ addDescriptions(npmMap.optionList(), linkDescription(npmMap, npmDescMap))
			+ "\nAny number of these key-value pairs may be specified in any order.";
	}

	void process(ParamList& pl, Everything& e)
	{	NonlocalPCMparams& np = e.eVars.fluidParams.npcmParams;
		while(true)
		{	NonlocalPCMparamsMember key;
			pl.get(key, NPM_Delim, npmMap, "key");
			switch(key)
			{	case NPM_lMax: pl.get(np.lMax, 0, "lMax", true); break;
				case NPM_Delim: return; //end of input
			}
		}
	}

	void printStatus(Everything& e, int iRep)
	{	NonlocalPCMparams& np = e.eVars.fluidParams.npcmParams;
		logPrintf(" \\\n\tlMax      %d", np.lMax);
	}
} commandNonlocalPCM;
