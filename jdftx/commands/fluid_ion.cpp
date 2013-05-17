

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
/*
struct CommandFluidIon : public Command
{
	CommandFluidIon() : Command("fluid-ion")
	{
		format = "<fluid-ion-id> <rHS> <Z=0.0> <Concentration=1.0 M>";
		comments =
			"Add a monoatomic hard sphere of type <fluid-ion-id> with charge <Z> in electrons\n"
			"and hard sphere radius <rHS> in bohr to the joint-density-functional theory fluid\n"
		    "with <Concentration> in mol/L.";
		allowMultiple = true;

		require("fluid");
	}

	void process(ParamList& pl, Everything& e)
	{	const FluidType& fluidType = e.eVars.fluidParams.fluidType;
		if((fluidType==FluidHSIonic)||(fluidType==FluidFittedCorrelations)
			||(fluidType==FluidScalarEOS)||(fluidType==FluidBondedVoids))			
		{
			HardSphereIon ion;
			pl.get(ion.name, string(), "fluid-ion-id", true);
			pl.get(ion.rHS, 0.0, "rHS", true);
			pl.get(ion.Z, 0.0, "Z", false);
			double conc;
			pl.get(conc, 1.0, "Concentration", false);
			ion.Concentration = conc * mol/liter;
			ion.MixFunc = false;
			ion.convCouplingModel = ConvCouplingNone; //Sets default to no Convolution Coupling
			e.eVars.fluidParams.hSIons.push_back(ion);
			return;
		}
		throw string("fluid-ion must be used with a JDFT3.0 Fluid\n");
	}

	void printStatus(Everything& e, int iRep)
	{	
		for(int iIon=0; iIon<int(e.eVars.fluidParams.hSIons.size()); iIon++)
			{	if(iIon==iRep)
				{	logPrintf("%s %16.10lf %16.10lf %10.4lf", e.eVars.fluidParams.hSIons[iIon].name.c_str(),
							  e.eVars.fluidParams.hSIons[iIon].rHS, -1.0*e.eVars.fluidParams.hSIons[iIon].Z,
							  e.eVars.fluidParams.hSIons[iIon].Concentration * liter/mol);
				}
			}
	}
}
CommandFluidIon;
*/