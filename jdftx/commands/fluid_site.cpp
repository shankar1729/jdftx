/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver

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
#include <fluid/FluidSolver.h>
#include <fluid/FluidComponent.h>
#include <fluid/Fex_ScalarEOS.h>

enum FluidSiteParameter
{	
	FSp_Znuc, //!<magnitude of the nuclear charge (positive) 
	FSp_sigmaNuc, //!<gaussian width of the nuclear charge (positive)
	FSp_Zelec, //!<magnitude of electron charge (positive)
	FSp_aElec, //!<cuspless-exponential width of electron charge distribution
	FSp_elecFilename, //!<filename to read in additional radial realspace electron charge distribution,
	FSp_elecFilenameG, //!<filename to read in additional radial Gspace electron charge distribution
	FSp_alpha, //!<isotropic polarizability
	FSp_aPol, //!<cuspless-exponential width of polarizability
	FSp_Delim //!< Delimiter used in parsing:
};

EnumStringMap<FluidSiteParameter> FSParamMap(
	FSp_Znuc, "Znuc", 
	FSp_sigmaNuc, "sigmaNuc",
	FSp_Zelec, "Zelec",
	FSp_aElec, "aElec",
	FSp_alpha, "alpha",
	FSp_aPol, "aPol",
	FSp_elecFilename, "elecFilename",
	FSp_elecFilenameG, "elecFilenameG"	);

EnumStringMap<FluidSiteParameter> FSParamDescMap(
	FSp_Znuc, "magnitude of the nuclear charge (positive)", 
	FSp_sigmaNuc, "gaussian width of the nuclear charge (positive)",
	FSp_Zelec, "magnitude of electron charge (positive)",
	FSp_aElec, "cuspless-exponential width of electron charge distribution",
	FSp_elecFilename, "filename to read in additional radial realspace electron charge distribution",
	FSp_elecFilenameG, "filename to read in additional radial Gspace electron charge distribution",
	FSp_alpha, "isotropic polarizability",
	FSp_aPol, "cuspless-exponential width of polarizability"	);

//Kendra: list of fluid components supported in JDFT3 at present
EnumStringMap<FluidComponent::Name> fluidComponentMap(
	//solvents
	FluidComponent::H2O, "H2O",
	FluidComponent::CHCl3, "CHCl3",
	FluidComponent::CCl4, "CCl4"/*,
	FluidComponent::DMC, "DMC",
	FluidComponent::EC, "EC",
	FluidComponent::PC, "PC",
	FluidComponent::DMF, "DMF",
	FluidComponent::THF, "THF",
	FluidComponent::EthylEther, "EthylEther",
  	FluidComponent::Chlorobenzene, "Chlorobenzene",
	FluidComponent::Isobutanol, "Isobutanol",
	FluidComponent::CarbonDisulfide, "CarbonDisulfide",
	FluidComponent::CustomSolvent, "CustomSolvent",
	//cations
	FluidComponent::Sodium, "Na+",
	FluidComponent::CustomCation, "CustomCation",
	//anions
	FluidComponent::Chloride, "Cl-", 
	FluidComponent::CustomAnion, "CustomAnion"*/ );


struct CommandFluidSiteParams : public Command
{
	

	CommandFluidSiteParams() : Command("fluid-site-params")
	{	
		format = " <solvent> <siteName> <key1> <value1> <key2> <value2> ...";
		comments = "Set parameters of <siteName> site for solvent component <solvent>=" + fluidComponentMap.optionList() 
			+ ".\nChanges default parameters of existing site. \nPossible keys and value types are:"
			+ addDescriptions(FSParamMap.optionList(), linkDescription(FSParamMap, FSParamDescMap))
			+ "\nAny number of these key-value pairs may be specified in any order.";
		
		require("fluid-solvent");
		hasDefault = true;
		allowMultiple = true;
	}
	
void process(ParamList& pl, Everything& e)
	{	
		if(e.eVars.fluidParams.fluidType == FluidNonlinearPCM || e.eVars.fluidParams.fluidType == FluidLinearPCM)
			return;
		FluidSolverParams& fsp = e.eVars.fluidParams;

		//Read in and check name of the solvent, get index of the solvent in FluidComponent
		FluidComponent::Name solventName; 
		pl.get(solventName, fsp.components[0]->name, fluidComponentMap, "solvent", false);
		std::shared_ptr<FluidComponent> FC;
		for (const auto& c : fsp.components)
		{
			if (solventName == c->name) 
				FC=c;
		}
		if (!FC)
			throw string("Choice of <solvent> is not valid.\n Hint: Issue fluid-solvent first");
		
		//Read in name of the site
		string siteName;
		pl.get(siteName, FC->molecule.sites[0]->name, "siteName", false);
		std::shared_ptr<Molecule::Site> site;
		for (const auto& s : FC->molecule.sites)
		{
			if(siteName == s->name) 
				site=s;
		}
		if (!site)	
			throw string("Choice of <siteName> is not valid.");

		//Read parameters:
		while(true)
		{	FluidSiteParameter key;
			pl.get(key, FSp_Delim, FSParamMap, "key");
			#define READ_AND_CHECK(param, op, val) \
				case FSp_##param: \
					pl.get(site->param, val, #param, true); \
					if(!(site->param op val)) throw string(#param " must be " #op " " #val); \
					break;
			switch(key)
			{	
				READ_AND_CHECK(Znuc,>,0.) 
				READ_AND_CHECK(sigmaNuc,>,0.)
				READ_AND_CHECK(Zelec,>=,0.) 
				READ_AND_CHECK(aElec,>,0.) 
				READ_AND_CHECK(elecFilename,!=,string("")) 
				READ_AND_CHECK(elecFilenameG,!=,string("")) 
				READ_AND_CHECK(alpha,>,0.)
				READ_AND_CHECK(aPol,>,0.)
				case FSp_Delim: return; //end of input
			}
			#undef READ_AND_CHECK
		}		
	}
	
void printStatus(Everything& e, int iRep)
	{	
		//prints all the sites and parameters, even if the default is unchanged
		#define PRINT(param) logPrintf(" \\\n\t" #param " %lg", s->param);
		
		if(e.eVars.fluidParams.fluidType == FluidNonlinearPCM || e.eVars.fluidParams.fluidType == FluidLinearPCM)
			return;
		if(iRep==0)
		{
			int counter=0;
			const FluidSolverParams& fsp = e.eVars.fluidParams;
			for (const auto& c : fsp.components)
			{
				string cName = fluidComponentMap.getString(c->name);
				for (const auto& s : c->molecule.sites)
				{	
					string sName = s->name;
					if(counter) 
						logPrintf("\nfluid-site-params ");
					logPrintf("%s %s",cName.c_str(),sName.c_str()),
					#define PRINT(param) logPrintf(" \\\n\t" #param " %lg", s->param);
					PRINT(Znuc)
					PRINT(sigmaNuc)
					PRINT(Zelec)
					PRINT(aElec)
					PRINT(alpha)
					PRINT(aPol)
					logPrintf(" \\\n\telecFilename ");
					if (s->elecFilename.length())
						logPrintf("%s", s->elecFilename.c_str());
					logPrintf(" \\\n\telecFilenameG ");
					if (s->elecFilenameG.length())
						logPrintf("%s", s->elecFilenameG.c_str());
					#undef PRINT
					
					counter++;	
				}
			}
		}			
	}
}
commandFSParams;


/*
EnumStringMap<ConvolutionCouplingSiteModel> couplingMap(
	ConvCouplingExponential, "Exponential",
	ConvCouplingExpCuspless, "ExpCuspless",
	ConvCouplingRadialFunction, "ReadRadialFunction",
	ConvCouplingBinaryKernel, "ReadBinaryKernel" );

struct CommandFluidCouplingScale : public Command
{
	CommandFluidCouplingScale() : Command("fluid-coupling-scale")
	{
		format = "<scale=1/9>";
		comments = "Scale Von Weisacker correction to kinetic energy by a constant factor <scale>.\n"
					"	Default is 1/9, corresponding to the gradient expansion for homogeneous electron gas.\n";
		hasDefault = true; 
		allowMultiple = false;
	}
	
	void process(ParamList& pl, Everything& e)
	{
		pl.get(e.eVars.fluidParams.convCouplingScale, 1.0/9.0, "scale");
		if(e.eVars.fluidParams.convCouplingScale == 0.0)
			die("Must specify nonzero scale for Von Weisacker correction.\nTo use bare Thomas-Fermi use command fluid-ex-corr instead.\n"); 
	}
	
	void printStatus(Everything& e, int iRep)
	{	
		logPrintf(" %lg", e.eVars.fluidParams.convCouplingScale);	
	}
	
}commandFluidCouplingScale;

struct CommandFluidSitePosition : public Command
{
	CommandFluidSitePosition() : Command("fluid-site-position")
	{
		format = "<fluid-site-id> <x0> <x1> <x2>";
		comments = "Specify a site of type <fluid-site-id> at location <x0> <x1> <x2> in Bohr.\n"; 
		hasDefault = false; 
		allowMultiple = true;
		require("fluid-site");
	}
	
	void process(ParamList& pl, Everything& e)
	{
		string id;
		pl.get(id, string(), "fluid-site-id", true);
		if((e.eVars.fluidParams.fluidType==FluidScalarEOS)||(e.eVars.fluidParams.fluidType==FluidScalarEOSCustom))			
		{
			for(int i=0; i<int(e.eVars.fluidParams.H2OSites.size()); i++)
			{
				H2OSite& site = e.eVars.fluidParams.H2OSites[i];
				if(site.name == id)
				{
						vector3<> pos;
						pl.get(pos[0], 0.0, "x0");
						pl.get(pos[1], 0.0, "x1");
						pl.get(pos[2], 0.0, "x2");
						site.Positions.push_back(pos);
						return;	
				}
			}
			throw string("fluid site "+id+" has not been defined.");
		}
		throw string("fluid-site-position must be used with a ScalarEOS Fluid\n");
	}
	
	void printStatus(Everything& e, int iRep)
	{	
		int iSite = 0;
		int counter = 0;
		while (iSite<int(e.eVars.fluidParams.H2OSites.size()))
		{
			int iPos = 0;
			H2OSite& site = e.eVars.fluidParams.H2OSites[iSite];
			while (iPos<int(site.Positions.size()))
			{
				vector3<> pos = site.Positions[iPos];
				if (iRep == counter)
				{
					logPrintf("%s %19.15lf %19.15lf %19.15lf",site.name.c_str(),pos[0],pos[1],pos[2]);
					return;
				}
				iPos++;
				counter++;
			}
			iSite++;
		}
	}
}commandFluidSitePosition;

struct CommandFluidSite : public Command
{
	CommandFluidSite() : Command("fluid-site")
	{
		format = "<fluid-site-id> <Z_nuc> <Z_site> <Z_at> <Exponential/ExpCuspless <width> [filename]\n"
				 "| ReadRadialFunction/ReadBinaryKernel <filename>";
		comments = "Specify site <fluid-site-id> with nuclear charge <Z_nuc>, site charge <Z_site>,\n"
					"   atomic number <Z_at>, and coupling electron density model specified by\n" 
					"   Exponential|ExpCuspless|ReadRadialFunction|ReadBinaryKernel\n"
					"	If Exponential or ExpCuspless, the exponential width must be specified and\n"
					"	a radial function contained in [filename] may be added to the exponential model.\n"
					"	If ReadRadialFunction/ReadBinaryKernel, the filename containing the \n" 
					"	electron density model must be specified.\n"
					"	<Z_at> specifies the atomic number of the element used for van Der Waals corrections between the fluid and explicit system.\n"
					"	If Z_at = 0, there are no van der Waals corrections between a site of this type and the explicit system.\n"
					"	Note that any dipole must be along the z-direction.\n";
		hasDefault = false;
		allowMultiple = true;
		require("fluid");
	}
	
	void process(ParamList& pl, Everything& e)
	{
		if((e.eVars.fluidParams.fluidType==FluidScalarEOS)||(e.eVars.fluidParams.fluidType==FluidScalarEOSCustom))			
		{
			H2OSite site;
			pl.get(site.name, string(), "fluid-site-id", true);
			//later make special cases here for O and H
			pl.get(site.Znuc, 0.0, "Z_nuc",true);
			pl.get(site.Z, 0.0, "Z_site",true);
			pl.get(site.atomicNumber, 0, "Z_at",true);
			pl.get(site.ccSiteModel, ConvCouplingExpCuspless, couplingMap, "couplingType");
			if (site.ccSiteModel == ConvCouplingExponential || site.ccSiteModel == ConvCouplingExpCuspless )
			{	
				pl.get(site.CouplingWidth, 0.0, "width", true);  
				pl.get(site.CouplingFilename, string(""), "filename");
			}
			else if ((site.ccSiteModel == ConvCouplingRadialFunction)
				||(site.ccSiteModel == ConvCouplingBinaryKernel))
			{
				pl.get(site.CouplingFilename, string(""), "filename", true);
			}
			e.eVars.fluidParams.H2OSites.push_back(site);
			e.eVars.fluidParams.fluidType = FluidScalarEOSCustom;
			return;
		}
		throw string("fluid-site must be used with a ScalarEOS Fluid\n");
	}
	
	void printStatus(Everything& e, int iRep)
	{	
		for(int iSite=0; iSite<int(e.eVars.fluidParams.H2OSites.size()); iSite++)
			{	if(iSite==iRep)
				{	
					H2OSite &site = e.eVars.fluidParams.H2OSites[iSite];
					logPrintf("%s %16.10lf %16.10lf %d ", site.name.c_str(), site.Znuc, site.Z, site.atomicNumber);
					fputs(couplingMap.getString(site.ccSiteModel), globalLog);
					if(site.ccSiteModel == ConvCouplingExponential || site.ccSiteModel == ConvCouplingExpCuspless )
						logPrintf(" %lg %s", site.CouplingWidth, site.CouplingFilename.c_str());
					if ((site.ccSiteModel == ConvCouplingRadialFunction) || (site.ccSiteModel == ConvCouplingBinaryKernel))
						logPrintf(" %s", site.CouplingFilename.c_str());
				}
			}
	}
	
}commandFluidSite;


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
		if((e.eVars.fluidParams.fluidType==FluidFittedCorrelations)||(e.eVars.fluidParams.cavityTension==FluidScalarEOS)
			||(e.eVars.fluidParams.fluidType==FluidBondedVoids)||(e.eVars.fluidParams.fluidType==FluidHSIonic))			
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
