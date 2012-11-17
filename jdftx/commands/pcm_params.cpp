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


struct CommandPcmParams : public Command
{
	CommandPcmParams() : Command("pcm-params")
	{
		format = "[<nc>] [<sigma>]";
		comments =  "Parameters for the PCM Cavity\n"
					"Control the critical density <nc> and smoothing parameter <sigma> for PCM cavities\n"
					"Defaults are set by the fluid chosen, but can be manually overwritten";
		hasDefault = true;
		require("fluid");
	}

	void process(ParamList& pl, Everything& e)
	{	FluidSolverParams& fsp = e.eVars.fluidParams;
		double nc=7e-4, sigma=0.6; //defaults for FluidLinear
		switch(e.eVars.fluidType)
		{	case FluidLinear: //defaults set above
				break;
			case FluidLinearPCM:
				nc = 0.000413567616685;
				sigma = 0.6;
				break;
			case FluidNonlinearPCM:
				nc = 0.0013516934742;
				sigma = 0.6;
				break;
			default: //Other fluids do not use these parameters
				break;
		}
		pl.get(fsp.nc, nc, "nc");
		pl.get(fsp.sigma, sigma, "sigma");
		
		// Set the cavitation terms to 0. These can be overwritten by the cavitation command
		fsp.cavityTension = 0;
		fsp.cavityPressure = 0;
	}

	void printStatus(Everything& e, int iRep)
	{	const FluidSolverParams& fsp = e.eVars.fluidParams;
		logPrintf("%lg %lg", fsp.nc, fsp.sigma);
	}
}
commandPcmParams;

struct CommandCavitation : public Command
{
	CommandCavitation() : Command("cavitation")
	{	format = "<cavity tension> <cavity pressure>";
		comments =  "Overwrites the non-electrostatic (cavitation + van der waals) terms for solvation.\n"
					"If not explicitly called, defaults are set by the fluid type chosen.";
		hasDefault = true;
		require("pcm-params");
	}
	
	void process(ParamList& pl, Everything& e)
	{	FluidSolverParams& fsp = e.eVars.fluidParams;
		double cavityTension=0., cavityPressure=0.; // Defaults for fluid linear
		
		switch(e.eVars.fluidType)
		{	case FluidLinear: //defaults set above
				break;
			case FluidLinearPCM:
				cavityTension = 1.96639422455e-05;
				cavityPressure = -4.33957737175e-06;
				break;
			case FluidNonlinearPCM:
				cavityTension = 1.77987130105e-05;
				cavityPressure = -3.23501563062e-06;
				break;
			default: //Other fluids do not use these parameters
				break;
		}
		
		pl.get(fsp.cavityTension, cavityTension, "cavityTension");
		pl.get(fsp.cavityPressure, cavityPressure, "cavityPressure");
	}
	
	void printStatus(Everything& e, int iRep)
	{	const FluidSolverParams& fsp = e.eVars.fluidParams;
		logPrintf("%lg %lg", fsp.cavityTension, fsp.cavityPressure);
	}
	
}
commandCavitation;