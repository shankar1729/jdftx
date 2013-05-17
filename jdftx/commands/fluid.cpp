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
	FluidClassicalDFT, "ClassicalDFT"
);

struct CommandFluid : public Command
{
	CommandFluid() : Command("fluid")
	{
		format = "[<type>=None] [<Temperature>=298K] [<Pressure>=1.01325bar]";
		comments = "Enable joint density functional theory\n"
			"\t<type> = " + fluidTypeMap.optionList() + "\n"
			"\t<Temperature> in Kelvin and <Pressure> in bars.";
		hasDefault = true;
		require("coulomb-interaction");
	}

	void process(ParamList& pl, Everything& e)
	{	FluidSolverParams& fsp = e.eVars.fluidParams;
		pl.get(fsp.fluidType, FluidNone, fluidTypeMap, "type");
		if((e.coulombParams.geometry != CoulombParams::Periodic)
			&& (fsp.fluidType != FluidNone))
			throw string("Fluids cannot be used with a truncated coulomb interaction");
		pl.get(fsp.T, 298., "Temperature"); fsp.T *= Kelvin; //convert to atomic units
		pl.get(fsp.P, 1.01325, "Pressure"); fsp.P *= Bar; //convert to atomic units
	}

	void printStatus(Everything& e, int iRep)
	{	const FluidSolverParams& fsp = e.eVars.fluidParams;
		logPrintf("%s", fluidTypeMap.getString(fsp.fluidType));
		if(fsp.fluidType != FluidNone)
			logPrintf(" %lf %lf", fsp.T/Kelvin, fsp.P/Bar);
	}
}
commandFluid;


EnumStringMap<FluidComponent::Name> solventMap
(	FluidComponent::H2O, "H2O",
	FluidComponent::CHCl3, "CHCl3",
	FluidComponent::CCl4, "CCl4",
	FluidComponent::DMC, "DMC",
	FluidComponent::EC, "EC",
	FluidComponent::PC, "PC",
	FluidComponent::DMF, "DMF",
	FluidComponent::THF, "THF"
);
EnumStringMap<FluidComponent::Name> cationMap
(	FluidComponent::Sodium, "Sodium"
);
EnumStringMap<FluidComponent::Name> anionMap
(	FluidComponent::Chloride, "Chloride"
);
EnumStringMap<FluidComponent::Functional> functionalMap
(	FluidComponent::ScalarEOS, "ScalarEOS",
	FluidComponent::FittedCorrelations, "FittedCorrelations",
	FluidComponent::BondedVoids, "BondedVoids",
	FluidComponent::MeanFieldLJ, "MeanFieldLJ"
);

struct CommandFluidSolvent : public Command
{
    CommandFluidSolvent() : Command("fluid-solvent")
	{
		format = "[<name>=H2O] [<functional>=ScalarEOS] [<concentration>=bulk]";
		comments = "Add solvent component " + solventMap.optionList() + " to fluid.\n"
			"For classical DFT fluids, the excess functional may be one of " + functionalMap.optionList() + "\n"
			"Optionally override <concentration> in mol/liter\n";
		
		require("fluid");
		require("pcm-variant");
		hasDefault = true;
		allowMultiple = true;
	}
	
	void process(ParamList& pl, Everything& e)
	{	FluidComponent::Name name; pl.get(name, FluidComponent::H2O, solventMap, "name");
		FluidComponent::Functional functional; pl.get(functional, FluidComponent::ScalarEOS, functionalMap, "functional");
		auto c = std::make_shared<FluidComponent>(name, e.eVars.fluidParams.T, functional);
		//Optionally override concentration:
		pl.get(c->Nbulk, c->Nbulk/(mol/liter), "concentration");
		c->Nbulk *= (mol/liter);
		if(c->Nbulk <= 0.) throw string("<concentration> must be positive");
		//Add to component list:
		e.eVars.fluidParams.addComponent(c);
		//Set PCM parameters if necessary:
		switch(e.eVars.fluidParams.fluidType)
		{	case FluidLinearPCM:
			case FluidNonlinearPCM:
			case FluidNonlocalPCM:
				if(e.eVars.fluidParams.solvents.size()>1)
					throw string("PCMs require exactly one solvent component - more than one specified.");
				e.eVars.fluidParams.setPCMparams();
				break;
			default:;
		}
	}
	
	void printStatus(Everything& e, int iRep)
	{	const FluidComponent& c = *(e.eVars.fluidParams.solvents[iRep]);
		logPrintf("%s %s %lg", solventMap.getString(c.name), functionalMap.getString(c.functional), c.Nbulk/(mol/liter));
	}
}
commandFluidSolvent;

