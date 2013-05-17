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

EnumStringMap<PCMVariant> pcmVariantMap
(	PCM_SGA13,   "SGA13", 
	PCM_GLSSA13, "GLSSA13",
	PCM_LA12,    "LA12", 
	PCM_PRA05,   "PRA05"
);
EnumStringMap<PCMVariant> pcmVariantDescMap
(	PCM_SGA13,   "PCM with weighted-density cavitation and dispersion [R. Sundararaman, D. Gunceler and T.A. Arias, (under preparation)]", 
	PCM_GLSSA13, "PCM with empirical cavity tension [D. Gunceler, K. Letchworth-Weaver, R. Sundararaman, K.A. Schwarz and T.A. Arias, arXiv:1301.6189]",
	PCM_LA12,    "PCM with no cavitation/dispersion contributions [K. Letchworth-Weaver and T.A. Arias, Phys. Rev. B 86, 075140 (2012)]", 
	PCM_PRA05,   "PCM with no cavitation/dispersion contributions [S.A. Petrosyan SA, A.A. Rigos and T.A. Arias, J Phys Chem B. 109, 15436 (2005)]"
);

struct CommandPcmVariant : public Command
{
	CommandPcmVariant() : Command("pcm-variant")
	{
		format = "[<variant>=GLSSA13]";
		comments = "Select LinearPCM or NonlinearPCM <variant> amongst:"
			+ addDescriptions(pcmVariantMap.optionList(), linkDescription(pcmVariantMap, pcmVariantDescMap));
		hasDefault = true;
		require("fluid");
	}

	void process(ParamList& pl, Everything& e)
	{	FluidSolverParams& fsp = e.eVars.fluidParams;
		if(fsp.fluidType==FluidNonlocalPCM) fsp.pcmVariant = PCM_SLSA13; //only option for NonlocalPCM
		else pl.get(fsp.pcmVariant, PCM_GLSSA13, pcmVariantMap, "variant");
	}

	void printStatus(Everything& e, int iRep)
	{	const FluidSolverParams& fsp = e.eVars.fluidParams;
		if(fsp.fluidType==FluidNonlocalPCM) logPrintf("SLSA13"); //only option for NonlocalPCM
		else logPrintf("%s", pcmVariantMap.getString(fsp.pcmVariant));
	}
}
commandPcmVariant;

enum PCMparameter
{	PCMp_lMax, //!< angular momentum truncation in nonlocal PCM
	PCMp_nc, //!< critical density for the PCM cavity shape function
	PCMp_sigma, //!< smoothing factor for the PCM cavity shape function
	PCMp_cavityTension, //!< effective surface tension (including dispersion etc.) of the cavity (hartree per bohr^2)
	PCMp_Delim //!< Delimiter used in parsing:
};
EnumStringMap<PCMparameter> pcmParamMap
(	PCMp_lMax,          "lMax",
	PCMp_nc,            "nc",
	PCMp_sigma,         "sigma",
	PCMp_cavityTension, "cavityTension"
);
EnumStringMap<PCMparameter> pcmParamDescMap
(	PCMp_lMax, "angular momentum truncation in nonlocal PCM",
	PCMp_nc, "critical density for the PCM cavity shape function",
	PCMp_sigma, "smoothing factor for the PCM cavity shape function",
	PCMp_cavityTension, "effective surface tension (including dispersion etc.) of the cavity (hartree per bohr^2)"
);

struct CommandPcmParams : public Command
{
	CommandPcmParams() : Command("pcm-params")
	{	
		format = "<key1> <value1> <key2> <value2> ...";
		comments = "Adjust PCM solvent parameters. Possible keys and value types are:"
			+ addDescriptions(pcmParamMap.optionList(), linkDescription(pcmParamMap, pcmParamDescMap))
			+ "\nAny number of these key-value pairs may be specified in any order.";
		
		require("fluid-solvent");
	}

	void process(ParamList& pl, Everything& e)
	{	FluidSolverParams& fsp = e.eVars.fluidParams;
		while(true)
		{	PCMparameter key;
			pl.get(key, PCMp_Delim, pcmParamMap, "key");
			#define READ_AND_CHECK(param, op, val) \
				case PCMp_##param: \
					pl.get(fsp.param, val, #param, true); \
					if(!(fsp.param op val)) throw string(#param " must be " #op " " #val); \
					break;
			switch(key)
			{	READ_AND_CHECK(lMax, >=, 1)
				READ_AND_CHECK(nc, >, 0.)
				READ_AND_CHECK(sigma, >, 0.)
				READ_AND_CHECK(cavityTension, <, DBL_MAX)
				case PCMp_Delim: return; //end of input
			}
			#undef READ_AND_CHECK
		}
	}

	void printStatus(Everything& e, int iRep)
	{	const FluidSolverParams& fsp = e.eVars.fluidParams;
		#define PRINT(param) logPrintf(" \\\n\t" #param " %lg", fsp.param);
		logPrintf(" \\\n\t lMax %d", fsp.lMax);
		PRINT(nc)
		PRINT(sigma)
		PRINT(cavityTension)
		#undef PRINT
	}
}
commandPcmParams;


/*
EnumStringMap<FluidSolverParams::SolventName> pcmSolventMap
(
	FluidSolverParams::H2O,   "H2O",
	FluidSolverParams::CHCl3, "CHCl3",
	FluidSolverParams::CCl4,  "CCl4",
	FluidSolverParams::DMC,   "DMC",
	FluidSolverParams::EC,    "EC",
	FluidSolverParams::PC,    "PC",
	FluidSolverParams::DMF,   "DMF",
	FluidSolverParams::THF,   "THF"
);
EnumStringMap<FluidSolverParams::SolventName> pcmSolventDescMap
(	FluidSolverParams::H2O,   "Water",
	FluidSolverParams::CHCl3, "Chloroform",
	FluidSolverParams::CCl4,  "Carbon tetrachloride",
	FluidSolverParams::DMC,   "Dimethyl carbonate",
	FluidSolverParams::EC,    "Ethylene carbonate",
	FluidSolverParams::PC,    "Propylene carbonate",
	FluidSolverParams::DMF,   "Dimethylformamide",
	FluidSolverParams::THF,   "Tetrahydrofuran"
);

struct CommandPcmSolvent : public Command
{
	CommandPcmSolvent() : Command("pcm-solvent")
	{
		format = "[<solvent>=H2O]";
		comments = "Load default PCM parameters for <solvent>, which is one of:"
			+ addDescriptions(pcmSolventMap.optionList(), linkDescription(pcmSolventMap, pcmSolventDescMap));
		hasDefault = true;
		require("pcm-variant");
	}

	void process(ParamList& pl, Everything& e)
	{	FluidSolverParams& fsp = e.eVars.fluidParams;
		pl.get(fsp.solventName, FluidSolverParams::H2O, pcmSolventMap, "solvent");
		fsp.setPCMparams();
	}

	void printStatus(Everything& e, int iRep)
	{	const FluidSolverParams& fsp = e.eVars.fluidParams;
		logPrintf("%s", pcmSolventMap.getString(fsp.solventName));
	}
}
commandPcmSolvent;


enum PCMparameter
{	PCMp_lMax, //!< angular momentum truncation in nonlocal PCM
	//Fit parameters:
	PCMp_nc, //!< critical density for the PCM cavity shape function
	PCMp_sigma, //!< smoothing factor for the PCM cavity shape function
	PCMp_cavityTension, //!< effective surface tension (including dispersion etc.) of the cavity (hartree per bohr^2)
	//Physical parameters:
	PCMp_epsBulk, //!< bulk dielectric constant
	PCMp_Nbulk, //!< bulk number-density of molecules in bohr^-3
	PCMp_pMol, //!< dipole moment of each molecule in e-bohr
	PCMp_epsInf, //!< optical-frequency dielectric constant
	PCMp_Pvap, //!< vapor pressure in Eh/bohr^3
	PCMp_sigmaBulk, //!< bulk surface tension in Eh/bohr^2
	PCMp_Rvdw, //!< effective van der Waals radius of the fluid (derived from equation of state) in bohrs
	PCMp_Res, //!< electrostatic radius of solvent (derived from nonlocal response) in bohrs
	//Delimiter used in parsing:
	PCMp_Delim
};
EnumStringMap<PCMparameter> pcmParamMap
(	PCMp_lMax,          "lMax",
	PCMp_nc,            "nc",
	PCMp_sigma,         "sigma",
	PCMp_cavityTension, "cavityTension",
	PCMp_epsBulk,       "epsBulk",
	PCMp_Nbulk,         "Nbulk",
	PCMp_pMol,          "pMol",
	PCMp_epsInf,        "epsInf",
	PCMp_Pvap,          "Pvap",
	PCMp_sigmaBulk,     "sigmaBulk",
	PCMp_Rvdw,          "Rvdw",
	PCMp_Res,           "Res"
);
EnumStringMap<PCMparameter> pcmParamDescMap
(	PCMp_lMax, "angular momentum truncation in nonlocal PCM"
	PCMp_nc, "critical density for the PCM cavity shape function",
	PCMp_sigma, "smoothing factor for the PCM cavity shape function",
	PCMp_cavityTension, "effective surface tension (including dispersion etc.) of the cavity (hartree per bohr^2)",
	PCMp_epsBulk, "bulk dielectric constant",
	PCMp_Nbulk, "bulk number-density of molecules in bohr^-3",
	PCMp_pMol, "dipole moment of each molecule in e-bohr",
	PCMp_epsInf, "optical-frequency dielectric constant",
	PCMp_Pvap, "vapor pressure in Eh/bohr^3",
	PCMp_sigmaBulk, "bulk surface tension in Eh/bohr^2",
	PCMp_Rvdw, "effective van der Waals radius of the fluid (derived from equation of state) in bohrs",
	PCMp_Res, "electrostatic radius of solvent (derived from nonlocal response) in bohrs"
);

struct CommandPcmParams : public Command
{
	CommandPcmParams() : Command("pcm-params")
	{	
		format = "<key1> <value1> <key2> <value2> ...";
		comments = "Adjust PCM solvent parameters. Possible keys and value types are:"
			+ addDescriptions(pcmParamMap.optionList(), linkDescription(pcmParamMap, pcmParamDescMap))
			+ "\nAny number of these key-value pairs may be specified in any order.";
		require("pcm-solvent");
	}

	void process(ParamList& pl, Everything& e)
	{	FluidSolverParams& fsp = e.eVars.fluidParams;
		while(true)
		{	PCMparameter key;
			pl.get(key, PCMp_Delim, pcmParamMap, "key");
			#define READ_AND_CHECK(param, op, val) \
				case PCMp_##param: \
					pl.get(fsp.param, val, #param, true); \
					if(!(fsp.param op val)) throw string(#param " must be " #op " " #val); \
					break;
			switch(key)
			{	READ_AND_CHECK(lMax, >=, 1)
				READ_AND_CHECK(nc, >, 0.)
				READ_AND_CHECK(sigma, >, 0.)
				READ_AND_CHECK(cavityTension, <, DBL_MAX)
				READ_AND_CHECK(epsBulk, >, 1.)
				READ_AND_CHECK(Nbulk, >, 0.)
				READ_AND_CHECK(pMol, >=, 0.)
				READ_AND_CHECK(epsInf, >=, 1.)
				READ_AND_CHECK(Pvap, >, 0.)
				READ_AND_CHECK(sigmaBulk, >, 0.)
				READ_AND_CHECK(Rvdw, >, 0.)
				READ_AND_CHECK(Res, >, 0.)
				case PCMp_Delim: return; //end of input
			}
			#undef READ_AND_CHECK
		}
	}

	void printStatus(Everything& e, int iRep)
	{	const FluidSolverParams& fsp = e.eVars.fluidParams;
		#define PRINT(param) logPrintf(" \\\n\t" #param " %lg", fsp.param);
		PRINT(lMax)
		PRINT(nc)
		PRINT(sigma)
		PRINT(cavityTension)
		PRINT(epsBulk)
		PRINT(Nbulk)
		PRINT(pMol)
		PRINT(epsInf)
		PRINT(Pvap)
		PRINT(sigmaBulk)
		PRINT(Rvdw)
		PRINT(Res)
		#undef PRINT
	}
}
commandPcmParams;


struct CommandIonicScreening : public Command
{
	CommandIonicScreening() : Command("ionic-screening")
	{
		format = "<concentration> <Zelectrolyte> <linear> <Rcation> <Ranion>";
		comments =
			"\t<concentration>: molar concentration of ions (default: 0.0, which turns off ionic screening)\n"
			"\t<Zelectrolyte>: magnitude of charge of the cations and anions (assumed equal)\n"
			"\t<linear>: linearity of screening = " + boolMap.optionList() + " (default: no)\n"
			"\t<Rcation>: Cation ionic radius in Angstroms (default: 1.16 [Na+])\n"
			"\t<Ranion>: Anion ionic radius in Angstroms (default: 1.67 [Cl-])";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	FluidSolverParams& fsp = e.eVars.fluidParams;
		pl.get(fsp.ionicConcentration, 0.0, "concentration");
		pl.get(fsp.ionicZelectrolyte, 1, "Zelectrolyte");
		pl.get(fsp.ionicRadiusMinus, 1.16, "Rcation"); //Note 'minus' is wrt electron-positive convention
		pl.get(fsp.ionicRadiusPlus, 1.67, "Ranion"); //Note 'plus' is wrt electron-positive convention
		//Check parameters
		if(fsp.ionicConcentration < 0.) throw("Ionic concentration must be non-negative");
		if(fsp.ionicZelectrolyte <= 0.) throw("Ionic charge magnitude must be positive");
		if(fsp.ionicRadiusMinus <= 0.) throw("Cation radius must be positive");
		if(fsp.ionicRadiusPlus <= 0.) throw("Anion radius must be positive");
		//convert to atomic units
		fsp.ionicConcentration *= mol/liter;
		fsp.ionicRadiusPlus *= Angstrom;
		fsp.ionicRadiusMinus *= Angstrom;
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg %d %lg %lg",
			e.eVars.fluidParams.ionicConcentration/(mol/liter), //report back in mol/liter
			e.eVars.fluidParams.ionicZelectrolyte,
			e.eVars.fluidParams.ionicRadiusMinus/Angstrom,
			e.eVars.fluidParams.ionicRadiusPlus/Angstrom);
	}
}
commandIonicScreening;
*/

struct CommandPCMnonlinearDebug : public Command
{
    CommandPCMnonlinearDebug() : Command("pcm-nonlinear-debug")
	{
		format = "<linearDielectric>=" + boolMap.optionList() + " <linearScreening>=" + boolMap.optionList();
		comments = "Emulate linear response of the dielectric or screening within NonlinearPCM (for debugging purposes only)";
	}
	
	void process(ParamList& pl, Everything& e)
	{	FluidSolverParams& fsp = e.eVars.fluidParams;
		pl.get(fsp.linearDielectric, false, boolMap, "linearDielectric", true);
		pl.get(fsp.linearScreening, false, boolMap, "linearScreening", true);
	}
	
	void printStatus(Everything& e, int iRep)
	{	const FluidSolverParams& fsp = e.eVars.fluidParams;
		logPrintf("%s %s", boolMap.getString(fsp.linearDielectric), boolMap.getString(fsp.linearScreening));
	}
}
commandPCMnonlinearDebug;
