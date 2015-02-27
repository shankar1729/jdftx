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
	FluidSaLSA, "SaLSA",
	FluidClassicalDFT, "ClassicalDFT"
);

struct CommandFluid : public Command
{
	CommandFluid() : Command("fluid")
	{
		format = "[<type>=None] [<Temperature>=298K] [<Pressure>=1.01325bar]";
		comments = "Perform joint density functional theory with fluid of <type>:\n"
			"\n+ None:\n\n"
			"   Standard vacuum DFT calculation with no solvation model.\n"
			"\n+ LinearPCM:\n\n"
			"   Use a solvation model that includes linear dielectric (and/or ionic)\n"
			"    response. Select a specific linear solvation model using pcm-variant.\n"
			"\n+ NonlinearPCM:\n\n"
			"   Use a solvation model that includes nonlinear dielectric (and/or ionic)\n"
			"   response, and accounts for dielectric saturation effects.\n"
			"   Select a specific nonlinear solvation model using pcm-variant.\n"
			"\n+ SaLSA:\n\n"
			"   Use the non-empirical nonlocal-response solvation model based on the\n"
			"   Spherically-averaged Liquid Susceptibility Ansatz.\n"
			"\n+ ClassicalDFT:\n\n"
			"   Full joint density-functional theory with a classical density-functional\n"
			"   description of the solvent. See fluid-solvent, fluid-cation, fluid-anion\n"
			"   and related commands for controlling the classical density-functional theory.\n"
			"\n"
			"Optionally adjust the fluid <Temperature> (in Kelvin) and <Pressure> (in bars).";
		hasDefault = true;
		require("coulomb-interaction");
	}

	void process(ParamList& pl, Everything& e)
	{	FluidSolverParams& fsp = e.eVars.fluidParams;
		pl.get(fsp.fluidType, FluidNone, fluidTypeMap, "type");
		if((e.coulombParams.geometry != CoulombParams::Periodic) && (fsp.fluidType != FluidNone))
			e.coulombParams.embedFluidMode = true; //require embedding in fluid mode (periodic Coulomb kernels in larger box)
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


struct CommandFluidGummelLoop : public Command
{
	CommandFluidGummelLoop() : Command("fluid-gummel-loop")
	{
		format = "[<maxIterations>=10] [<Atol>=1e-5]";
		comments =
			"Settings for the fluid <--> electron self-consistency loop:\n"
			"+ <maxIterations>: Max number of electron and fluid minimization pairs\n"
			"+ <Atol>: Free energy convergence criterion for this outer loop";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.cntrl.fluidGummel_nIterations, 10, "maxIterations");
		pl.get(e.cntrl.fluidGummel_Atol, 1e-5, "Atol");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%d %le", e.cntrl.fluidGummel_nIterations, e.cntrl.fluidGummel_Atol);
	}
}
commandFluidGummelLoop;


struct CommandFluidInitialState : public Command
{
	CommandFluidInitialState() : Command("fluid-initial-state")
	{
		format = "<filename>";
		comments = "Read initial state of a fluid (compatible with *.fluidState from dump End State)";
		
		forbid("initial-state");
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.eVars.fluidInitialStateFilename, string(), "filename", true);
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s", e.eVars.fluidInitialStateFilename.c_str());
	}
}
commandFluidInitialState;


struct CommandFluidVdwScale : public Command
{
	CommandFluidVdwScale() : Command("fluid-vdwScale")
	{
		format = "<scale=0.75>";
		comments = "Scale van der Waals interactions between fluid and explicit system by a constant factor <scale>.\n\n"
			"Default is fluid specific and ranges between 0.4 to 1.3.\n"
			"Set to 0 to use the prefactor corresponding to fluid exchange-correlation.";
		require("fluid-solvent");
	}
	
	void process(ParamList& pl, Everything& e)
	{
		pl.get(e.eVars.fluidParams.vdwScale, 0.75, "scale");
	}
	
	void printStatus(Everything& e, int iRep)
	{	
		logPrintf(" %lg", e.eVars.fluidParams.vdwScale);	
	}
}
commandFluidVDWscale;


EnumStringMap<FluidComponent::Name> solventMap
(	FluidComponent::H2O, "H2O",
	FluidComponent::CHCl3, "CHCl3",
	FluidComponent::CCl4, "CCl4",
	FluidComponent::CH3CN, "CH3CN",
	FluidComponent::DMC, "DMC",
	FluidComponent::EC, "EC",
	FluidComponent::PC, "PC",
	FluidComponent::DMF, "DMF",
	FluidComponent::THF, "THF",
	FluidComponent::DMSO, "DMSO",
	FluidComponent::CH2Cl2, "CH2Cl2",
	FluidComponent::Ethanol, "Ethanol",
	FluidComponent::Methanol, "Methanol",
	FluidComponent::Octanol, "Octanol",
	FluidComponent::EthylEther, "EthylEther",
  	FluidComponent::Chlorobenzene, "Chlorobenzene",
	FluidComponent::Isobutanol, "Isobutanol",
	FluidComponent::CarbonDisulfide, "CarbonDisulfide",
	FluidComponent::Glyme, "Glyme",
	FluidComponent::EthyleneGlycol, "EthyleneGlycol"
);
EnumStringMap<FluidComponent::Name> cationMap
(	FluidComponent::Sodium, "Na+",
	FluidComponent::HydratedSodium, "Na(6H2O)+",
	FluidComponent::Potassium, "K+",
        FluidComponent::HydratedPotassium, "K(6H2O)+",
	FluidComponent::Hydronium, "H3O+",
	FluidComponent::HydratedHydronium, "H3O(4H2O)+"
);
EnumStringMap<FluidComponent::Name> anionMap
(	FluidComponent::Chloride, "Cl-",
	FluidComponent::Fluoride, "F-",
	FluidComponent::Perchlorate, "ClO4-",
	FluidComponent::Hydroxide, "OH-",
	FluidComponent::HydratedHydroxide, "OH(4H2O)-"
);
EnumStringMap<FluidComponent::Functional> functionalMap
(	FluidComponent::ScalarEOS, "ScalarEOS",
	FluidComponent::FittedCorrelations, "FittedCorrelations",
	FluidComponent::BondedVoids, "BondedVoids",
	FluidComponent::MeanFieldLJ, "MeanFieldLJ"
);

EnumStringMap<FluidComponent::TranslationMode> translationModeMap
(	FluidComponent::ConstantSpline, "ConstantSpline",
	FluidComponent::LinearSpline, "LinearSpline",
	FluidComponent::Fourier, "Fourier"
);

EnumStringMap<FluidComponent::Representation> representationMap
(	FluidComponent::MuEps, "MuEps",
	FluidComponent::Pomega, "Pomega",
	FluidComponent::PsiAlpha, "PsiAlpha"
);

const EnumStringMap<S2quadType>& s2quadTypeMap = S2quadTypeMap;

enum FluidComponentMember
{	FCM_epsBulk, //!< bulk dielectric constant
	FCM_pMol, //!< dipole moment of each molecule in e-bohr
	FCM_epsInf, //!< optical-frequency dielectric constant
	FCM_Pvap, //!< vapor pressure in Eh/bohr^3
	FCM_sigmaBulk, //!< bulk surface tension in Eh/bohr^2
	FCM_Rvdw, //!< effective van der Waals radius of the fluid (derived from equation of state) in bohrs
	FCM_Res, //!< electrostatic radius of solvent (derived from nonlocal response) in bohrs
	//Extras for ClassicalDFT:
	FCM_epsLJ, //!< Lennard-Jones well depth for Mean-Field LJ excess functional
	FCM_representation, //!< ideal gas representation
	FCM_s2quadType, //!< orientation quadrature type
	FCM_quad_nBeta, //!< number of beta samples for Euler quadrature
	FCM_quad_nAlpha, //!< number of alpha samples for Euler quadrature
	FCM_quad_nGamma, //!< number of gamma samples for Euler quadrature
	FCM_translationMode, //!< translation operator type
	FCM_Nnorm, //!< unit cell molecule count constraint
	//Delimiter used in parsing:
	FCM_Delim
};
EnumStringMap<FluidComponentMember> fcmMap
(	FCM_epsBulk,       "epsBulk",
	FCM_pMol,          "pMol",
	FCM_epsInf,        "epsInf",
	FCM_Pvap,          "Pvap",
	FCM_sigmaBulk,     "sigmaBulk",
	FCM_Rvdw,          "Rvdw",
	FCM_Res,           "Res",
	FCM_epsLJ,          "epsLJ",
	FCM_representation, "representation",
	FCM_s2quadType,     "s2quadType",
	FCM_quad_nBeta,     "quad_nBeta",
	FCM_quad_nAlpha,    "quad_nAlpha",
	FCM_quad_nGamma,    "quad_nGamma",
	FCM_translationMode,"translation",
	FCM_Nnorm,          "Nnorm"
);
EnumStringMap<FluidComponentMember> fcmDescMap
(	FCM_epsBulk, "bulk dielectric constant",
	FCM_pMol, "dipole moment of each molecule in e-bohr",
	FCM_epsInf, "optical-frequency dielectric constant",
	FCM_Pvap, "vapor pressure in Eh/bohr^3",
	FCM_sigmaBulk, "bulk surface tension in Eh/bohr^2",
	FCM_Rvdw, "effective van der Waals radius of the fluid (derived from equation of state) in bohrs",
	FCM_Res, "electrostatic radius of solvent (derived from nonlocal response) in bohrs",
	FCM_epsLJ, "Lennard-Jones well depth for Mean-Field LJ excess functional",
	FCM_representation, "ideal gas representation: " + addDescriptions(representationMap.optionList(), nullDescription, "\n   - "),
	FCM_s2quadType, "orientation quadrature type:" + addDescriptions(s2quadTypeMap.optionList(), nullDescription, "\n   - "),
	FCM_quad_nBeta, "number of beta samples for Euler quadrature",
	FCM_quad_nAlpha, "number of alpha samples for Euler quadrature",
	FCM_quad_nGamma, "number of gamma samples for Euler quadrature",
	FCM_translationMode, "translation operator type: " + addDescriptions(translationModeMap.optionList(), nullDescription, "\n   - "),
	FCM_Nnorm, "unit cell molecule count constraint"
);

//Abstract base class for fluid-solvent fluid-cation and fluid-anion
struct CommandFluidComponent : public Command
{
private:
	const EnumStringMap<FluidComponent::Name>& nameMap;
	FluidComponent::Name defaultName;
	FluidComponent::Functional defaultFunctional;
	bool defaultEnabled;
	
protected:
	CommandFluidComponent(string suffix, const EnumStringMap<FluidComponent::Name>& nameMap, FluidComponent::Name defaultName, FluidComponent::Functional defaultFunctional, bool defaultEnabled)
	: Command("fluid-"+suffix), nameMap(nameMap), defaultName(defaultName), defaultFunctional(defaultFunctional), defaultEnabled(defaultEnabled)
	{
		format = (defaultEnabled ? ("[<name>=" + string(nameMap.getString(defaultName)) +"] [<concentration>=bulk]") : "<name> <concentration>")
			+ " [<functional>=" + string(functionalMap.getString(defaultFunctional)) + "]"
			" [<key1> <value1>] ...";
		comments = "Add " + suffix + " component to fluid, where <name> may be one of:"
			+ addDescriptions(nameMap.optionList(), nullDescription) + "\n\n"
			+ string(defaultEnabled
				? "The concentration may be specified explicitly in mol/liter or set to 'bulk' to use bulk fluid value."
				: "The concentration must be specified explicitly in mol/liter.") + "\n"
			"For classical DFT fluids, the excess functional may be one of " + functionalMap.optionList() + ".\n"
			"\n"
			"Optional component properties may be overridden using an arbitrary number of <key> <value> pairs.\n"
			+ (suffix=="solvent"
				? ("Possible keys and value types are:"
					+ addDescriptions(fcmMap.optionList(), linkDescription(fcmMap, fcmDescMap))
					+ "\n\nAny number of these key-value pairs may be specified in any order.")
				: string("See command fluid-solvent for a description of adjustable properties."));

		hasDefault = defaultEnabled;
		allowMultiple = true;
	}
	
public:
	void process(ParamList& pl, Everything& e)
	{	//Read parameters:
		FluidComponent::Name name; pl.get(name, defaultName, nameMap, "name", !defaultEnabled);
		string keyNbulk; pl.get(keyNbulk, string("bulk"), "concentration", !defaultEnabled);
		FluidComponent::Functional functional; pl.get(functional, defaultFunctional, functionalMap, "functional");
		//Make component and add to list:
		auto c = std::make_shared<FluidComponent>(name, e.eVars.fluidParams.T, functional);
		e.eVars.fluidParams.addComponent(c);
		//Check and optionally override concentration:
		if((!defaultEnabled) || (keyNbulk != "bulk"))
		{	istringstream iss(keyNbulk);
			double Nbulk = 0.; iss >> Nbulk;
			if(iss.fail() || !iss.eof()) throw string("Conversion of parameter <concentration> failed");
			if(Nbulk <= 0.) throw string("<concentration> must be positive");
			c->Nbulk = Nbulk * (mol/liter);
		}
		//Optional properties
		while(true)
		{	FluidComponentMember key;
			pl.get(key, FCM_Delim, fcmMap, "key");
			#define READ_AND_CHECK(param, op, val) \
				case FCM_##param: \
					pl.get(c->param, val, #param, true); \
					if(!(c->param op val)) throw string(#param " must be " #op " " #val); \
					break;
			#define READ_ENUM(param, paramDefault) \
				case FCM_##param: \
					pl.get(c->param, paramDefault, param##Map, #param, true); \
					break;
			switch(key)
			{	READ_AND_CHECK(epsBulk, >, 1.)
				READ_AND_CHECK(pMol, >=, 0.)
				READ_AND_CHECK(epsInf, >=, 1.)
				READ_AND_CHECK(Pvap, >, 0.)
				READ_AND_CHECK(sigmaBulk, >, 0.)
				READ_AND_CHECK(Rvdw, >, 0.)
				READ_AND_CHECK(Res, >, 0.)
				READ_AND_CHECK(epsLJ, >, 0.)
				READ_ENUM(representation, FluidComponent::MuEps)
				READ_ENUM(s2quadType, QuadOctahedron)
				READ_AND_CHECK(quad_nBeta, >, 0u)
				READ_AND_CHECK(quad_nAlpha, >=, 0u)
				READ_AND_CHECK(quad_nGamma, >=, 0u)
				READ_ENUM(translationMode, FluidComponent::LinearSpline)
				READ_AND_CHECK(Nnorm, >=, 0.)
				case FCM_Delim: return; //end of input
			}
			#undef READ_AND_CHECK
			#undef READ_ENUM
		}
	}
	
	void print(const Everything& e, const FluidComponent& c)
	{	logPrintf("%s %lg %s", nameMap.getString(c.name), c.Nbulk/(mol/liter), functionalMap.getString(c.functional));
		#define PRINT(param) logPrintf(" \\\n\t" #param " %lg", c.param);
		#define PRINT_UINT(param) logPrintf(" \\\n\t" #param " %u", c.param);
		#define PRINT_ENUM(param) logPrintf(" \\\n\t" #param " %s", param##Map.getString(c.param));
		PRINT(epsBulk)
		PRINT(pMol)
		PRINT(epsInf)
		PRINT(Pvap)
		PRINT(sigmaBulk)
		PRINT(Rvdw)
		PRINT(Res)
		if(e.eVars.fluidParams.fluidType == FluidClassicalDFT)
		{	PRINT(epsLJ)
			PRINT_ENUM(representation)
			PRINT_ENUM(s2quadType)
			PRINT_UINT(quad_nBeta)
			PRINT_UINT(quad_nAlpha)
			PRINT_UINT(quad_nGamma)
			PRINT_ENUM(translationMode)
			PRINT(Nnorm)
		}
		#undef PRINT
		#undef PRINT_INT
		#undef PRINT_ENUM
	}
};

struct CommandFluidSolvent : public CommandFluidComponent
{
	CommandFluidSolvent() : CommandFluidComponent("solvent", solventMap, FluidComponent::H2O, FluidComponent::ScalarEOS, true)
    {	require("pcm-variant"); //which in turn requires fluid
	}
	
	void process(ParamList& pl, Everything& e)
	{	CommandFluidComponent::process(pl, e);
		//Set PCM parameters if necessary:
		switch(e.eVars.fluidParams.fluidType)
		{	case FluidLinearPCM:
			case FluidNonlinearPCM:
			case FluidSaLSA:
				if(e.eVars.fluidParams.solvents.size()>1)
					throw string("PCMs require exactly one solvent component - more than one specified.");
				e.eVars.fluidParams.setPCMparams();
				break;
				
			case FluidClassicalDFT:
				e.eVars.fluidParams.setCDFTparams();
				break;
			default:;
		}
	}
	
	void printStatus(Everything& e, int iRep)
	{	print(e, *(e.eVars.fluidParams.solvents[iRep]));
	}
}
commandFluidSolvent;

struct CommandFluidCation : public CommandFluidComponent
{
	CommandFluidCation() : CommandFluidComponent("cation", cationMap, FluidComponent::Sodium, FluidComponent::MeanFieldLJ, false)
    {	require("fluid-solvent"); //which in turn requires fluid indirectly
	}
	
	void process(ParamList& pl, Everything& e)
	{	CommandFluidComponent::process(pl, e);
		//Ions not yet supported in this version of ClassicalDFT
		//if(e.eVars.fluidParams.fluidType == FluidClassicalDFT)
		//	throw string("Ions not yet supported in ClassicalDFT fluid.");
	}

	void printStatus(Everything& e, int iRep)
	{	print(e, *(e.eVars.fluidParams.cations[iRep]));
	}
}
commandFluidCation;

struct CommandFluidAnion : public CommandFluidComponent
{
	CommandFluidAnion() : CommandFluidComponent("anion", anionMap, FluidComponent::Chloride, FluidComponent::MeanFieldLJ, false)
    {	require("fluid-solvent"); //which in turn requires fluid indirectly
	}
	
	void process(ParamList& pl, Everything& e)
	{	CommandFluidComponent::process(pl, e);
		//Ions not yet supported in this version of ClassicalDFT
		//if(e.eVars.fluidParams.fluidType == FluidClassicalDFT)
		//	throw string("Ions not yet supported in ClassicalDFT fluid.");
	}

	void printStatus(Everything& e, int iRep)
	{	print(e, *(e.eVars.fluidParams.anions[iRep]));
	}
}
commandFluidAnion;

enum FluidSiteParameter
{	
	FSp_Znuc, //!<magnitude of the nuclear charge (positive) 
	FSp_sigmaNuc, //!<gaussian width of the nuclear charge (positive)
	FSp_Zelec, //!<magnitude of electron charge (positive)
	FSp_aElec, //!<exponential decay width of electron charge distribution
	FSp_sigmaElec, //!<width of peak in electron charge distribution
	FSp_rcElec, //!<Location of peak in electron charge distribution	
	FSp_elecFilename, //!<filename to read in additional radial realspace electron charge distribution,
	FSp_elecFilenameG, //!<filename to read in additional radial Gspace electron charge distribution
	FSp_alpha, //!<isotropic polarizability
	FSp_aPol, //!<cuspless-exponential width of polarizability
	FSp_Rhs, //!< hard sphere radius
	FSp_Delim //!< Delimiter used in parsing:
};

EnumStringMap<FluidSiteParameter> FSParamMap(
	FSp_Znuc, "Znuc", 
	FSp_sigmaNuc, "sigmaNuc",
	FSp_Zelec, "Zelec",
	FSp_aElec, "aElec",
	FSp_sigmaElec, "sigmaElec",
	FSp_rcElec,"rcElec",
	FSp_alpha, "alpha",
	FSp_aPol, "aPol",
	FSp_Rhs, "Rhs",
	FSp_elecFilename, "elecFilename",
	FSp_elecFilenameG, "elecFilenameG"	);

EnumStringMap<FluidSiteParameter> FSParamDescMap(
	FSp_Znuc, "magnitude of the nuclear charge (positive)", 
	FSp_sigmaNuc, "gaussian width of the nuclear charge (positive)",
	FSp_Zelec, "magnitude of electron charge (positive)",
	FSp_aElec, "exponential decay width of electron charge distribution",
	FSp_sigmaElec, "width of peak in electron charge distribution",
	FSp_rcElec, "location of peak in electron charge distribution",
	FSp_elecFilename, "filename to read in additional radial realspace electron charge distribution",
	FSp_elecFilenameG, "filename to read in additional radial Gspace electron charge distribution",
	FSp_alpha, "isotropic polarizability",
	FSp_aPol, "cuspless-exponential width of polarizability",
	FSp_Rhs, "hard sphere radius for use in FMT"						);

//Kendra: list of fluid components supported 
EnumStringMap<FluidComponent::Name> fluidComponentMap(
	//solvents
	FluidComponent::H2O, "H2O",
	FluidComponent::CHCl3, "CHCl3",
	FluidComponent::CCl4, "CCl4",
	FluidComponent::CH3CN, "CH3CN",/*
	FluidComponent::DMC, "DMC",
	FluidComponent::EC, "EC",
	FluidComponent::PC, "PC",
	FluidComponent::DMF, "DMF",
	FluidComponent::THF, "THF",
	FluidComponent::EthylEther, "EthylEther",
  	FluidComponent::Chlorobenzene, "Chlorobenzene",
	FluidComponent::Isobutanol, "Isobutanol",
	FluidComponent::CarbonDisulfide, "CarbonDisulfide",
	FluidComponent::CustomSolvent, "CustomSolvent",*/
	//cations
	FluidComponent::Sodium, "Na+",
	FluidComponent::HydratedSodium, "Na(H2O)4+",
	FluidComponent::CustomCation, "CustomCation",
	//anions
	FluidComponent::Chloride, "Cl-", 
	FluidComponent::Fluoride, "F-", 
	FluidComponent::Perchlorate, "ClO4-",
	FluidComponent::CustomAnion, "CustomAnion" );

struct CommandFluidSiteParams : public Command
{
	CommandFluidSiteParams() : Command("fluid-site-params")
	{	
		format = " <component> <siteName> <key1> <value1> <key2> <value2> ...";
		comments = "Set parameters of site <siteName> for fluid <component> which may be one of:"
			+ addDescriptions(fluidComponentMap.optionList(), nullDescription)
			+ "\n\nPossible keys and value types are:"
			+ addDescriptions(FSParamMap.optionList(), linkDescription(FSParamMap, FSParamDescMap))
			+ "\n\nAny number of these key-value pairs may be specified in any order.";
		
		require("fluid-solvent");
		hasDefault = false;
		allowMultiple = true;
	}
	
	void process(ParamList& pl, Everything& e)
	{
		if(e.eVars.fluidParams.fluidType == FluidNonlinearPCM || e.eVars.fluidParams.fluidType == FluidLinearPCM || e.eVars.fluidParams.fluidType == FluidNone)
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
				READ_AND_CHECK(Znuc,>=,0.) 
				READ_AND_CHECK(sigmaNuc,>=,0.)
				READ_AND_CHECK(Zelec,>=,0.) 
				READ_AND_CHECK(aElec,>,0.)
				READ_AND_CHECK(sigmaElec,>=,0.)
				READ_AND_CHECK(rcElec,>=,0.)
				READ_AND_CHECK(elecFilename,!=,string("")) 
				READ_AND_CHECK(elecFilenameG,!=,string("")) 
				READ_AND_CHECK(alpha,>,0.)
				READ_AND_CHECK(aPol,>,0.)
				READ_AND_CHECK(Rhs,>,0.)
				case FSp_Delim: return; //end of input
			}
			#undef READ_AND_CHECK
		}		
	}
	
void printStatus(Everything& e, int iRep)
	{	
		//prints all the sites and parameters, even if the default is unchanged
		#define PRINT(param) logPrintf(" \\\n\t" #param " %lg", s->param);
		
		if(e.eVars.fluidParams.fluidType == FluidNonlinearPCM || e.eVars.fluidParams.fluidType == FluidLinearPCM || e.eVars.fluidParams.fluidType == FluidNone)
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
					PRINT(sigmaElec)
				        PRINT(rcElec)
					PRINT(alpha)
					PRINT(aPol)
					PRINT(Rhs)
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


EnumStringMap<FMixFunctional> fMixMap
(	
	LJPotential, "LJPotential",
	GaussianKernel, "GaussianKernel"
);

struct CommandFluidMixingFunctional : public Command 
{
    CommandFluidMixingFunctional() : Command("fluid-mixing-functional")
	{
	  format = "<fluid1> <fluid2> <energyScale> [<lengthScale>] [<FMixType>=LJPotential]";
	  comments = 
		"Couple named fluids <fluid1> and <fluid2> which could each be one of:"
		+ addDescriptions(fluidComponentMap.optionList(), nullDescription)
		+ "\n\ntogether with a mixing functional of type:"
		+ addDescriptions(fMixMap.optionList(), nullDescription)
		+ "\n\nwith strength <energyScale> (in Eh) and range parameter <lengthScale> (in bohrs).\n"; 

	  require("fluid-solvent"); //which in turn requires fluid indirectly
	  allowMultiple = true;
	}
	
	void process(ParamList& pl, Everything& e)
	{	
	      FluidSolverParams& fsp = e.eVars.fluidParams;
	      FmixParams fmp;
	      
	      FluidComponent::Name name1; 
	      pl.get(name1, fsp.components[0]->name, fluidComponentMap, "fluid1", true);
	      
	      FluidComponent::Name name2; 
	      pl.get(name2, fsp.components[0]->name, fluidComponentMap, "fluid2", true);
	      
	
	      for(const std::shared_ptr<FluidComponent> c: fsp.components)
	      {
		if(c->name == name1)  fmp.fluid1 = c;
		if(c->name == name2)  fmp.fluid2 = c;
	      }
	      
	      if(!fmp.fluid1)
		throw string("Choice of <fluid1> = %s is not valid.\n Hint: Issue fluid-solvent first.",name1);
	      if(!fmp.fluid2)
		throw string("Choice of <fluid2> = %s is not valid.\n Hint: Issue fluid-solvent first.",name2);
	      if(fmp.fluid1->name == fmp.fluid2->name)
		throw string("<fluid1>=<fluid2> Cannot specify mixing functional for the same fluid.");
	     
	      double default_energyscale = sqrt(fmp.fluid1->epsLJ*fmp.fluid2->epsLJ);
	      pl.get(fmp.energyScale, default_energyscale,"energyScale", true);
	      
	      double default_lengthscale = fmp.fluid1->Rvdw + fmp.fluid2->Rvdw;
	      pl.get(fmp.lengthScale, default_lengthscale,"lengthScale");	     
	     
	      FMixFunctional defaultFunctional = LJPotential;
	      pl.get(fmp.FmixType, defaultFunctional, fMixMap, "FMixType");
	     
	      fsp.FmixList.push_back(fmp);
	}
	
	void printStatus(Everything& e, int iRep)
	{	
	    const FluidSolverParams& fsp = e.eVars.fluidParams;
	    const FmixParams& fmp = fsp.FmixList[iRep];
	    
	    string c1Name = fluidComponentMap.getString(fmp.fluid1->name);
	    string c2Name = fluidComponentMap.getString(fmp.fluid2->name);
	    string fmixName = fMixMap.getString(fmp.FmixType);
	    
	    logPrintf("%s %s %lg %lg %s",c1Name.c_str(),c2Name.c_str(),fmp.energyScale,fmp.lengthScale, fmixName.c_str());
	}
}
commandFluidMixingFunctional;

	
struct CommandFluidDielectricConstant : public Command
{
    CommandFluidDielectricConstant() : Command("fluid-dielectric-constant")
	{
		format = "[<epsBulkOverride>=0] [<epsInfOverride>=0]";
		comments = "Override bulk static or high frequency dieelctric constant of fluid (if non-zero values specified)";
	}
	
	void process(ParamList& pl, Everything& e)
	{	FluidSolverParams& fsp = e.eVars.fluidParams;
		pl.get(fsp.epsBulkOverride, 0., "epsBulkOverride");
		pl.get(fsp.epsInfOverride, 0., "epsInfOverride");
	}
	
	void printStatus(Everything& e, int iRep)
	{	const FluidSolverParams& fsp = e.eVars.fluidParams;
		logPrintf("%lg %lg", fsp.epsBulkOverride, fsp.epsInfOverride);
	}
}
commandFluidDielectricConstant;
