/*-------------------------------------------------------------------
Copyright 2013 Deniz Gunceler

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

enum  SCFparameter
{
	nIterations, //! maximum iterations (single point calculation if 0)
	energyDiffThreshold, //! convergence threshold for energy difference between successive iterations
	mixingMode, //! whether density or potentials
	//Delimiter used in parsing:
	scfDelim
};

EnumStringMap<SCFparameter> scfParamMap
(	nIterations, "nIterations",
	energyDiffThreshold, "energyDiffThreshold",
    mixingMode, "mixingMode"
);
EnumStringMap<SCFparameter> scfParamDescMap
(	nIterations, "maximum iterations (single point calculation if 0)",
	energyDiffThreshold, "convergence threshold for energy difference between successive iterations",
    mixingMode, "whether density or potential will be mixed at each SCF step"
);

enum SCFmixing
{
	density, //! Mix electron density (n) and kinetic energy density (tau)
	potential //! Mix the local electronic potential (Vscloc) and the kinetic energy potential (Vtau)
};

EnumStringMap<SCFmixing> scfMixing
(	density, "density",
	potential, "potential"
);

struct CommandsScfParams: public Command
{
	CommandsScfParams() : Command("residual-minimization")
	{	
		format = "<key1> <value1> <key2> <value2> ...";
		comments = "Enables self-consistent residual minimization.  If provided, keys adjust SCF parameters. Possible keys and value types are:"
			+ addDescriptions(scfParamMap.optionList(), linkDescription(scfParamMap, scfParamDescMap))
			+ "\nAny number of these key-value pairs may be specified in any order.";
		hasDefault = false;
		forbid("fix-electron-density");
	}
	
	void process(ParamList& pl, Everything& e)
	{	e.cntrl.minimisingResidual = true;
	
	/*	while(true)
		{	SCFparameter key;
			pl.get(key, scfDelim, scfParamMap, "key");
			#define READ_AND_CHECK(param, op, val) \
				case ##param: \
					pl.get(fsp.param, val, #param, true); \
					if(!(fsp.param op val)) throw string(#param " must be " #op " " #val); \
					break;
			switch(key)
			{	READ_AND_CHECK(nIter, >, 0.)
				READ_AND_CHECK(energyDiffThreshold, >, 0.)
				case mixingMode:
					pl.get(e., val, #param, true);
				case scfDelim: return; //end of input
			}
			#undef READ_AND_CHECK
		}*/
		
	}
	
	void printStatus(Everything& e, int iRep)
	{	fputs("scf", globalLog);
	} 
	
} commandsScfParams;
