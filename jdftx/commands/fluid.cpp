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
		comments = "Enable joint density functional theory\n"
			"\t<type> = " + fluidTypeMap.optionList() + "\n"
			"\t<Temperature> in Kelvin and <Pressure> in bars.";
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
			"\t<maxIterations>: Max number of electron and fluid minimization pairs\n"
			"\t<Atol>: Free energy convergence criterion for this outer loop";
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


struct CommandFluidVDWscale : public Command
{
	CommandFluidVDWscale() : Command("fluid-vdwScale")
	{
		format = "<scale=0.75>";
		comments = "Scale van der Waals interactions between fluid and explicit system by a constant factor <scale>.\n"
					"   Default is scale=0.75 (scale factor for PBE).\n"
					"   Set to 0 to use prefactor corresponding to fluid exchange-correlation.";
					
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
	FluidComponent::EthylEther, "EthylEther",
  	FluidComponent::Chlorobenzene, "Chlorobenzene",
	FluidComponent::Isobutanol, "Isobutanol",
	FluidComponent::CarbonDisulfide, "CarbonDisulfide"
);
EnumStringMap<FluidComponent::Name> cationMap
(	FluidComponent::Sodium, "Na+"
);
EnumStringMap<FluidComponent::Name> anionMap
(	FluidComponent::Chloride, "Cl-"
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
	FCM_representation, "ideal gas representation = " + representationMap.optionList(),
	FCM_s2quadType, "orientation quadrature type = " + s2quadTypeMap.optionList(),
	FCM_quad_nBeta, "number of beta samples for Euler quadrature",
	FCM_quad_nAlpha, "number of alpha samples for Euler quadrature",
	FCM_quad_nGamma, "number of gamma samples for Euler quadrature",
	FCM_translationMode, "translation operator type = " + translationModeMap.optionList(),
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
		comments = "Add " + suffix + " component " + nameMap.optionList() + " to fluid.\n"
			+ string(defaultEnabled
				? "The concentration may be specified explicitly in mol/liter or set to 'bulk' to use bulk fluid value."
				: "The concentration must be specified explicitly in mol/liter.") + "\n"
			"For classical DFT fluids, the excess functional may be one of " + functionalMap.optionList() + "\n"
			"Optional component properties may be overridden using an arbitrary number of <key> <value> pairs.\n"
			+ (suffix=="solvent"
				? ("Possible keys and value types are:"
					+ addDescriptions(fcmMap.optionList(), linkDescription(fcmMap, fcmDescMap))
					+ "\nAny number of these key-value pairs may be specified in any order.")
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
		if(e.eVars.fluidParams.fluidType == FluidClassicalDFT)
			throw string("Ions not yet supported in ClassicalDFT fluid.");
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
		if(e.eVars.fluidParams.fluidType == FluidClassicalDFT)
			throw string("Ions not yet supported in ClassicalDFT fluid.");
	}

	void printStatus(Everything& e, int iRep)
	{	print(e, *(e.eVars.fluidParams.anions[iRep]));
	}
}
commandFluidAnion;
	
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
