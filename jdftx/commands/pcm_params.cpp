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
	CommandPcmParams() : Command("jdft1-shape") //TODO: Change name to pcm-params and update documentation
	{
		format = "[<nc>=7e-4] [<sigma>=0.6]";
		comments = "Control the critical density <nc> and smoothing parameter <sigma> for PCM cavities";
		hasDefault = true;
		require("fluid");
	}

	void process(ParamList& pl, Everything& e)
	{	FluidSolverParams& fsp = e.eVars.fluidParams;
		double nc=7e-4, sigma=0.6, cavityPressure=0., cavityTension=0.; //defaults for FluidLinear
		switch(e.eVars.fluidType)
		{	case FluidLinear: //defaults set above
				break;
			case FluidLinearPCM:
				//TODO
				break;
			case FluidNonlinearPCM:
				//TODO
				break;
			default: //Other fluids do not use these parameters
				break;
		}
		pl.get(fsp.nc, nc, "nc");
		pl.get(fsp.sigma, sigma, "sigma");
		pl.get(fsp.cavityTension, cavityTension, "cavityTension");
		pl.get(fsp.cavityPressure, cavityPressure, "cavityPressure");
	}

	void printStatus(Everything& e, int iRep)
	{	const FluidSolverParams& fsp = e.eVars.fluidParams;
		logPrintf("%lg %lg %lg %lg", fsp.nc, fsp.sigma, fsp.cavityTension, fsp.cavityPressure);
	}
}
commandPcmParams;
