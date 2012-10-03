/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

This file is part of JDFTx.cd 

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
#include <electronic/FluidJDFTx.h>


EnumStringMap<ConvolutionCouplingSiteModel> couplingMap(
	ConvCouplingExponential, "Exponential",
	ConvCouplingExpCuspless, "ExpCuspless",
	ConvCouplingRadialFunction, "ReadRadialFunction",
	ConvCouplingBinaryKernel, "ReadBinaryKernel" );

struct CommandFluidCoupling : public Command
{
	CommandFluidCoupling() : Command("fluid-coupling")
	{
		format = "Exponential/ExpCuspless [<oxygenWidth>=0.35003] [<hydrogenWidth>=0.26342] [oxygenFilename] [hydrogenFilename] \n"
				 "| ReadRadialFunction/ReadBinaryKernel <oxygenFilename> <hydrogenFilename> ";
		comments = "Specify electron density model for coupling (default: ExpCuspless)\n"
					"	If Exponential or ExpCuspless, the oxygen <oxygenWidth> and hydrogen <hydrogenWidth>\n"
					"	exponential widths may also be specified. Radial functions contained in\n"
					"   [oxygenFilename] and [hydrogenFilename] may be added to the exponential models.\n"
					"	If ReadRadialFunction/ReadBinaryKernel, the filenames containing the \n" 
					"	oxygen and hydrogen electron density models must be specified.";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.eVars.fluidParams.convCouplingH2OModel, ConvCouplingExpCuspless, couplingMap, "couplingType");
		if (e.eVars.fluidParams.convCouplingH2OModel == ConvCouplingExponential || e.eVars.fluidParams.convCouplingH2OModel == ConvCouplingExpCuspless )
		{	
			//pl.get(e.eVars.fluidParams.oxygenWidth, 0.39285, "oxygenWidth");
			//pl.get(e.eVars.fluidParams.hydrogenWidth, 0.46108, "hydrogenWidth");
			pl.get(e.eVars.fluidParams.oxygenWidth, 0.35003, "oxygenWidth"); // nonlinear fits to Wannier functions in cubic Ice constrained to SPCE site charges
			pl.get(e.eVars.fluidParams.hydrogenWidth, 0.26342, "hydrogenWidth");
			pl.get(e.eVars.fluidParams.oxygenFilename, string(""), "oxygenFilename");
			pl.get(e.eVars.fluidParams.hydrogenFilename, string(""), "hydrogenFilename");
			return;
		}
		if ((e.eVars.fluidParams.convCouplingH2OModel == ConvCouplingRadialFunction)
			||(e.eVars.fluidParams.convCouplingH2OModel == ConvCouplingBinaryKernel))
		{
			pl.get(e.eVars.fluidParams.oxygenFilename, string(""), "oxygenFilename", true);
			pl.get(e.eVars.fluidParams.hydrogenFilename, string(""), "hydrogenFilename", true);
			return;
		}
	}

	void printStatus(Everything& e, int iRep)
	{	fputs(couplingMap.getString(e.eVars.fluidParams.convCouplingH2OModel), globalLog);
		if(e.eVars.fluidParams.convCouplingH2OModel == ConvCouplingExponential || e.eVars.fluidParams.convCouplingH2OModel == ConvCouplingExpCuspless )
			logPrintf(" %lg %lg %s %s", e.eVars.fluidParams.oxygenWidth, e.eVars.fluidParams.hydrogenWidth,
					  e.eVars.fluidParams.oxygenFilename.c_str(), e.eVars.fluidParams.hydrogenFilename.c_str());
		if ((e.eVars.fluidParams.convCouplingH2OModel == ConvCouplingRadialFunction)
			||(e.eVars.fluidParams.convCouplingH2OModel == ConvCouplingBinaryKernel))
			logPrintf(" %s %s", e.eVars.fluidParams.oxygenFilename.c_str(), e.eVars.fluidParams.hydrogenFilename.c_str());
		
	}
	
	
}
commandFluidCoupling;

struct CommandFluidCharge : public Command
{
	CommandFluidCharge() : Command("fluid-charge")
	{
		format = "[<ZO>=0.8476] [<ZH>=-0.4238] [<Znuc_O>=8.0] [<Znuc_H>=1.0]";
		comments = "For water-based fluids, sets oxygen and hydrogen site charges (<ZO> and <ZH>)\n"
					"and nuclear charges (<Znuc_O> and <Znuc_H>) for convolution coupling purposes only.\n"
					"Eventually, classical DFT site charges should match coupling site charges.";
		hasDefault = true;
	}
	
	void process(ParamList& pl, Everything& e)
	{
		pl.get(e.eVars.fluidParams.oxygenSiteCharge, 0.8476, "ZO");
		pl.get(e.eVars.fluidParams.hydrogenSiteCharge, -0.4238, "ZH");
		pl.get(e.eVars.fluidParams.oxygenZnuc, 8.0, "Znuc_O");
		pl.get(e.eVars.fluidParams.hydrogenZnuc, 1.0, "Znuc_H");
	}
	
	void printStatus(Everything& e, int iRep)
	{	
			logPrintf(" %lg %lg %lg %lg", e.eVars.fluidParams.oxygenSiteCharge, e.eVars.fluidParams.hydrogenSiteCharge,
					  e.eVars.fluidParams.oxygenZnuc, e.eVars.fluidParams.hydrogenZnuc);	
	}
}
commandFluidCharge;


struct CommandFluidIonCoupling : public Command
{
	CommandFluidIonCoupling() : Command("fluid-ion-coupling")
	{
						
		format = "<fluid-ion-id>  Exponential/ExpCuspless <width> <Znuc> [filename]\n"
				 "| ReadRadialFunction/ReadBinaryKernel <filename>";
		comments =
			"Add a convolution coupling functional between <fluid-ion-id> and electrons by specifying an \n"
			"ion radial electron density model.\n"
			"	If Exponential or ExpCuspless, the  exponential width <width> and nuclear charge <Znuc> must be specified.\n"
			"	Optionally, a radial function contained in [filename] may be added to the exponential model.\n"
			"	If ReadRadialFunction/ReadBinaryKernel, the filename containing the \n" 
			"	electron density model must be specified.\n"
			"	The nuclear charge of the ion is determined by <Znuc> \n" 
			"	or the integral of the electron density model specified by <filename> \n" 
			"	and the charge <Z> specified by fluid-ion (default is neutral).";
		allowMultiple = true;

		require("fluid-ion");
	}
	
	void process(ParamList& pl, Everything& e)
	{
		string id;
		pl.get(id, string(), "fluid-ion-id", true);
		if((e.eVars.fluidType==FluidLischner10)||(e.eVars.fluidType==FluidScalarEOS)
			||(e.eVars.fluidType==FluidBondedVoids)||(e.eVars.fluidType==FluidHSIonic))			
		{
			for(int i=0; i<int(e.eVars.fluidParams.hSIons.size()); i++)
			{
				HardSphereIon& hSIon = e.eVars.fluidParams.hSIons[i];
				if(hSIon.name == id)
				{
					pl.get(hSIon.convCouplingModel, ConvCouplingExpCuspless, couplingMap, "couplingType", true);
					if (hSIon.convCouplingModel == ConvCouplingExponential || e.eVars.fluidParams.convCouplingH2OModel == ConvCouplingExpCuspless )
					{
						pl.get(hSIon.CouplingWidth, 0.0, "width", true);
						pl.get(hSIon.Znuc, 0.0, "Znuc", true);
						pl.get(hSIon.CouplingFilename, string(""), "filename", false);
						return;
					}
					if ((hSIon.convCouplingModel == ConvCouplingRadialFunction)
						||(hSIon.convCouplingModel == ConvCouplingBinaryKernel))
					{
						pl.get(hSIon.CouplingFilename, string(""), "filename", true);
						return;
					}
				}
			}
			throw string("Fluid ion "+id+" has not been defined.");
		}
		throw string("fluid-ion-coupling must be used with a JDFT3.0 liquid functional.\n");
	}

	void printStatus(Everything& e, int iRep)
	{	
		int iIon = 0;
		int iCpl = 0;
		while (iIon<int(e.eVars.fluidParams.hSIons.size()))
		{
			HardSphereIon& hSIon = e.eVars.fluidParams.hSIons[iIon];
			if((iCpl==iRep)&&(hSIon.convCouplingModel==ConvCouplingExponential || e.eVars.fluidParams.convCouplingH2OModel == ConvCouplingExpCuspless))
			{
				logPrintf("%s ", hSIon.name.c_str()),
				fputs(couplingMap.getString(hSIon.convCouplingModel), globalLog);
				logPrintf(" %16.10lf %16.10lf %s",	hSIon.CouplingWidth, hSIon.Znuc, hSIon.CouplingFilename.c_str());
				return;
			}
			else if((iCpl==iRep)&&((hSIon.convCouplingModel==ConvCouplingRadialFunction)||(hSIon.convCouplingModel==ConvCouplingBinaryKernel)))
			{
				logPrintf("%s ", hSIon.name.c_str()),
				fputs(couplingMap.getString(hSIon.convCouplingModel), globalLog);
				logPrintf(" %s", hSIon.CouplingFilename.c_str());
				return;
			}
			else if(hSIon.convCouplingModel!=ConvCouplingNone)
			{
				iCpl++;
				iIon++;
			}
			else
			{
				iIon++;
			}
		}		
	}
}
CommandFluidIonCoupling;
