/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman, Kendra Letchworth Weaver

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
#include <core/Units.h>
#include <electronic/FluidSolver.h>

/*
struct CommandFluidIonMixing : public Command
{
	CommandFluidIonMixing() : Command("fluid-ion-mixing")
	{
		format = "<fluid-ion-id> <rSolv> <ESolv=0.004>";
		comments =
			"Add a mixing functional between <fluid-ion-id> and water with interaction strength <ESolv>\n"
			"and gaussian kernel radius <rSolv> within a joint-density-functional theory fluid.";
		allowMultiple = true;

		require("fluid-ion");
	}
	
	void process(ParamList& pl, Everything& e)
	{
		string id;
		pl.get(id, string(), "fluid-ion-id", true);
		if((e.eVars.fluidParams.fluidType==FluidFittedCorrelations)||(e.eVars.fluidParams.fluidType==FluidScalarEOS)
			||(e.eVars.fluidParams.fluidType==FluidBondedVoids))			
		{
			for(int i=0; i<int(e.eVars.fluidParams.hSIons.size()); i++)
			{
				if( e.eVars.fluidParams.hSIons[i].name == id)
				{
					e.eVars.fluidParams.hSIons[i].MixFunc = true;
					pl.get(e.eVars.fluidParams.hSIons[i].Rsolv, 0.0, "rSolv", true);
					pl.get(e.eVars.fluidParams.hSIons[i].Esolv, 0.004, "ESolv", false);
					return;
				}
			}
			throw string("Fluid ion "+id+" has not been defined.");
		}
		throw string("fluid-ion-mixing must be used with a JDFT3.0 water functional.\n");
	}

	void printStatus(Everything& e, int iRep)
	{	
		int iIon = 0;
		int iMix = 0;
		while (iIon<int(e.eVars.fluidParams.hSIons.size()))
		{
			if((iMix==iRep)&&(e.eVars.fluidParams.hSIons[iIon].MixFunc))
			{
				logPrintf("%s %16.10lf %16.10lf", e.eVars.fluidParams.hSIons[iIon].name.c_str(),
					e.eVars.fluidParams.hSIons[iIon].Rsolv, e.eVars.fluidParams.hSIons[iIon].Esolv);
				return;
			}
			else if(e.eVars.fluidParams.hSIons[iIon].MixFunc)
			{
				iMix++;
				iIon++;
			}
			else
			{
				iIon++;
			}
		}		
	}
}
CommandFluidIonMixing;
*/