/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver

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

EnumStringMap<FluidType> fluidTypeMap
(	FluidNone, "None",
	FluidLinearPCM, "LinearPCM",
	FluidNonlinearPCM, "NonlinearPCM",
	FluidNonlocalPCM, "NonlocalPCM",
	FluidFittedCorrelations, "FittedCorrelations",
	FluidScalarEOS, "ScalarEOS",
	FluidScalarEOSCustom, "ScalarEOSCustom",
	FluidBondedVoids, "BondedVoids",
	FluidHSIonic, "HSIonic"
);

struct CommandFluid : public Command
{
	CommandFluid() : Command("fluid")
	{
		format = "[<type>=None] [<Temperature>=298K]";
		comments = "Enable joint density functional theory\n"
			"\t<type> = " + fluidTypeMap.optionList() + "\n"
			"\t<Temperature> in Kelvin";
		hasDefault = true;
		require("coulomb-interaction");
	}

	void process(ParamList& pl, Everything& e)
	{	FluidSolverParams& fsp = e.eVars.fluidParams;
		pl.get(fsp.fluidType, FluidNone, fluidTypeMap, "type");
		if((e.coulombParams.geometry != CoulombParams::Periodic)
			&& (fsp.fluidType != FluidNone))
			throw string("Fluids cannot be used with a truncated coulomb interaction");
		pl.get(fsp.T, 298.0, "Temperature"); fsp.T *= Kelvin; //convert to atomic units
	}

	void printStatus(Everything& e, int iRep)
	{	const FluidSolverParams& fsp = e.eVars.fluidParams;
		logPrintf("%s", fluidTypeMap.getString(fsp.fluidType));
		if(fsp.fluidType != FluidNone)
			logPrintf(" %lf", fsp.T/Kelvin);
	}
}
commandFluid;
