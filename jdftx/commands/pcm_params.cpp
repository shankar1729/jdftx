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
(	PCM_CANDLE,  "CANDLE",
	PCM_SGA13,   "SGA13", 
	PCM_GLSSA13, "GLSSA13",
	PCM_LA12,    "LA12",
	PCM_SoftSphere, "SoftSphere",
	PCM_SCCS_g09,      "SCCS_g09",
	PCM_SCCS_g03,      "SCCS_g03",
	PCM_SCCS_g03p,     "SCCS_g03p",
	PCM_SCCS_g09beta,  "SCCS_g09beta",
	PCM_SCCS_g03beta,  "SCCS_g03beta",
	PCM_SCCS_g03pbeta, "SCCS_g03pbeta",
	PCM_SCCS_cation,   "SCCS_cation",
	PCM_SCCS_anion,    "SCCS_anion"
);
EnumStringMap<PCMVariant> pcmVariantDescMap
(	PCM_CANDLE,  "Charge-asymmetry corrected, local-response, nonlocal-cavity solvation model \\cite CANDLE",
	PCM_SGA13,   "PCM with weighted-density cavitation and dispersion \\cite CavityWDA", 
	PCM_GLSSA13, "PCM with empirical cavity tension \\cite NonlinearPCM",
	PCM_LA12,    "PCM with no cavitation/dispersion contributions \\cite PCM-Kendra", 
	PCM_SoftSphere, "Soft-sphere continuum solvation model \\cite PCM-SoftSphere",
	PCM_SCCS_g09,      "g09 parametrization of SCCS local linear model for water \\cite PCM-SCCS",
	PCM_SCCS_g03,      "g03 parametrization of SCCS local linear model for water \\cite PCM-SCCS",
	PCM_SCCS_g03p,     "g03' parametrization of SCCS local linear model for water \\cite PCM-SCCS",
	PCM_SCCS_g09beta,  "g09+beta parametrization of SCCS local linear model for water \\cite PCM-SCCS",
	PCM_SCCS_g03beta,  "g03+beta parametrization of SCCS local linear model for water \\cite PCM-SCCS",
	PCM_SCCS_g03pbeta, "g03'+beta parametrization of SCCS local linear model for water \\cite PCM-SCCS",
	PCM_SCCS_cation,   "cations-only parametrization of SCCS local linear model for water \\cite PCM-SCCS-charged",
	PCM_SCCS_anion,    "anions-only parametrization of SCCS local linear model for water \\cite PCM-SCCS-charged"
);

struct CommandPcmVariant : public Command
{
	CommandPcmVariant() : Command("pcm-variant", "jdftx/Fluid/Parameters")
	{
		format = "[<variant>=GLSSA13]";
		comments = "Select <variant> of LinearPCM or NonlinearPCM that determines\n"
			"the cavity and related energies (cavitation, dispersion etc.).\n"
			"CANDLE and SCCS variants are only supported for LinearPCM.\n"
			"Here, <variant> must be one of:"
			+ addDescriptions(pcmVariantMap.optionList(), linkDescription(pcmVariantMap, pcmVariantDescMap));
		hasDefault = true;
		require("fluid");
	}

	void process(ParamList& pl, Everything& e)
	{	FluidSolverParams& fsp = e.eVars.fluidParams;
		pl.get(fsp.pcmVariant, PCM_GLSSA13, pcmVariantMap, "variant");
		if(fsp.fluidType==FluidSaLSA)
			 fsp.pcmVariant = PCM_SaLSA; //only option for SaLSA
		//Check variant compatibility with fluidType
		if(fsp.fluidType!=FluidNone)
		{	//check only when fluid is not None, so that you can switch any
			//fluid input file to vacuum simply by commenting out fluid line
			if(fsp.fluidType!=FluidLinearPCM && ( fsp.pcmVariant==PCM_CANDLE || isPCM_SCCS(fsp.pcmVariant) ) )
				throw string("CANDLE and SCCS variants can only be used with fluid LinearPCM");
		}
	}

	void printStatus(Everything& e, int iRep)
	{	const FluidSolverParams& fsp = e.eVars.fluidParams;
		if(fsp.fluidType==FluidSaLSA) logPrintf("SaLSA"); //only option for SaLSA
		else logPrintf("%s", pcmVariantMap.getString(fsp.pcmVariant));
	}
}
commandPcmVariant;

enum PCMparameter
{	PCMp_lMax, //!< angular momentum truncation in SaLSA
	PCMp_nc, //!< critical density for the PCM cavity shape function
	PCMp_sigma, //!< smoothing factor for the PCM cavity shape function
	PCMp_cavityTension, //!< effective surface tension (including dispersion etc.) of the cavity (hartree per bohr^2)
	PCMp_cavityPressure, //!< effective pressure on the cavity (hartree per bohr^3) for SCCS and SoftSphere models
	PCMp_cavityScale, //!< atomic radius scale factor for soft sphere solvation model
	PCMp_ionSpacing, //!< extra spacing from dielectric to ionic cavity in bohrs for soft sphere model
	PCMp_zMask0, //z center in lattice coordinates for cavity mask
	PCMp_zMaskH, //half-width in z lattice coordinates for cavity mask
	PCMp_zMaskIonH, //half-width in z lattice coordinates for ion cavity mask
	PCMp_zMaskSigma, //smoothness of z-mask in bohrs
	PCMp_rhoMin, //!< min electron density (bohr^-3) for SCCS cavity switching function
	PCMp_rhoMax, //!< max electron density (bohr^-3) for SCCS cavity switching function
	PCMp_rhoDelta, //!< electron density change (bohr^-3) for SCCS cavity area calculation
	PCMp_eta_wDiel, //!< fit parameter for dielectric cavity in CANDLE
	PCMp_sqrtC6eff, //!< sqrt(effective molecule C6 coefficient) for CANDLE
	PCMp_pCavity, //!< sensitivity of cavity to surface electric fields [e-a0/Eh] in CANDLE
	PCMp_Ztot, //! Total valence charge on the solvent, used by CANDLE
	PCMp_screenOverride, //! Overrides screening length
	PCMp_Delim //!< Delimiter used in parsing
};
EnumStringMap<PCMparameter> pcmParamMap
(	PCMp_lMax,          "lMax",
	PCMp_nc,            "nc",
	PCMp_sigma,         "sigma",
	PCMp_cavityTension, "cavityTension",
	PCMp_cavityPressure, "cavityPressure",
	PCMp_cavityScale, "cavityScale",
	PCMp_ionSpacing, "ionSpacing",
	PCMp_zMask0, "zMask0",
	PCMp_zMaskH, "zMaskH",
	PCMp_zMaskIonH, "zMaskIonH",
	PCMp_zMaskSigma, "zMaskSigma",
	PCMp_rhoMin, "rhoMin",
	PCMp_rhoMax, "rhoMin",
	PCMp_rhoDelta, "rhoDelta",
	PCMp_eta_wDiel, "eta_wDiel",
	PCMp_sqrtC6eff, "sqrtC6eff",
	PCMp_pCavity, "pCavity",
	PCMp_Ztot, "Ztot",
	PCMp_screenOverride, "screenOverride"
);
EnumStringMap<PCMparameter> pcmParamDescMap
(	PCMp_lMax, "angular momentum truncation in SaLSA",
	PCMp_nc, "critical density for the PCM cavity shape function",
	PCMp_sigma, "smoothing factor for the PCM cavity shape function",
	PCMp_cavityTension, "effective surface tension (including dispersion etc.) of the cavity (hartree per bohr^2)",
	PCMp_cavityPressure, "effective pressure on the cavity (hartree per bohr^3) for SCCS and soft sphere models",
	PCMp_cavityScale, "atomic radius scale factor for soft sphere model",
	PCMp_ionSpacing, "extra spacing from dielectric to ionic cavity in bohrs for soft sphere model",
	PCMp_zMask0, "center in z lattice coordinates for cavity mask (default: 0)",
	PCMp_zMaskIonH, "half-width in z lattice coordinates for ion cavity mask (default: 0 => disabled)",
	PCMp_zMaskH, "half-width in z lattice coordinates for cavity mask (default: 0 => disabled)",
	PCMp_zMaskSigma, "smoothness of z-mask in bohrs (default: 0.5)",
	PCMp_rhoMin, "min electron density (bohr^-3) for SCCS cavity switching function",
	PCMp_rhoMax, "max electron density (bohr^-3) for SCCS cavity switching function",
	PCMp_rhoDelta, "electron density change (bohr^-3) for SCCS cavity area calculation",
	PCMp_eta_wDiel, "fit parameter for dielectric cavity in CANDLE",
	PCMp_sqrtC6eff, "sqrt(effective molecule C6 coefficient) for CANDLE",
	PCMp_pCavity, "sensitivity of cavity to surface electric fields [a.u.] in CANDLE",
	PCMp_Ztot, "total valence charge on the solvent, used by CANDLE",
	PCMp_screenOverride, "overrides the screening length calculated from fluid-components"
);

struct CommandPcmParams : public Command
{
	CommandPcmParams() : Command("pcm-params", "jdftx/Fluid/Parameters")
	{	
		format = "<key1> <value1> <key2> <value2> ...";
		comments = "Adjust PCM solvent parameters. Possible keys and value types are:"
			+ addDescriptions(pcmParamMap.optionList(), linkDescription(pcmParamMap, pcmParamDescMap))
			+ "\n\nAny number of these key-value pairs may be specified in any order.";
		
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
			{	READ_AND_CHECK(lMax, >=, 0)
				READ_AND_CHECK(nc, >, 0.)
				READ_AND_CHECK(sigma, >, 0.)
				READ_AND_CHECK(cavityTension, <, DBL_MAX)
				READ_AND_CHECK(cavityPressure, <, DBL_MAX)
				READ_AND_CHECK(cavityScale, >, 0.)
				READ_AND_CHECK(ionSpacing, >=, 0.)
				READ_AND_CHECK(zMask0, <, DBL_MAX)
				READ_AND_CHECK(zMaskH, >=, 0.)
				READ_AND_CHECK(zMaskIonH, >=, 0.)
				READ_AND_CHECK(zMaskSigma, >, 0.)
				READ_AND_CHECK(rhoMin, >, 0.)
				READ_AND_CHECK(rhoMax, >, 0.)
				READ_AND_CHECK(rhoDelta, >, 0.)
				READ_AND_CHECK(eta_wDiel, >=, 0.)
				READ_AND_CHECK(sqrtC6eff, >=, 0.)
				READ_AND_CHECK(pCavity, <, DBL_MAX)
				READ_AND_CHECK(Ztot, >, 0.)
				READ_AND_CHECK(screenOverride, >, 0.)
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
		PRINT(cavityPressure)
		PRINT(cavityScale)
		PRINT(ionSpacing)
		PRINT(zMask0)
		PRINT(zMaskH)
		PRINT(zMaskIonH)
		PRINT(zMaskSigma)
		PRINT(rhoMin)
		PRINT(rhoMax)
		PRINT(rhoDelta)
		PRINT(eta_wDiel)
		PRINT(sqrtC6eff)
		PRINT(pCavity)
		PRINT(Ztot)
		PRINT(screenOverride)
		#undef PRINT
	}
}
commandPcmParams;


struct CommandPcmNonlinearDebug : public Command
{
    CommandPcmNonlinearDebug() : Command("pcm-nonlinear-debug", "jdftx/Fluid/Parameters")
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
commandPcmNonlinearDebug;



struct CommandIonWidth : public Command
{
	CommandIonWidth() : Command("ion-width", "jdftx/Fluid/Parameters")
	{
		format = "Ecut | fftbox | <width>";
		comments = "Manually specify width of gaussian representations of nuclear charge in bohr\n"
			"or set automatically based on either energy cut-off (Ecut) or grid spacing (fftbox).\n"
			"The widened charges are only used in the interaction with fluids, and does not\n"
			"affect the energy of the electronic system. Default is Ecut-based selection\n"
			"for LinearPCM/NonlinearPCM and 0 (no widening) for all other fluid types.";
		hasDefault = true;
		require("fluid");
	}

	void process(ParamList& pl, Everything& e)
	{	string key; pl.get(key, string(), "width");
		if(!key.length()) //default
		{	switch(e.eVars.fluidParams.fluidType)
			{	case FluidLinearPCM:
				case FluidNonlinearPCM:
					e.iInfo.ionWidthMethod = IonInfo::IonWidthEcut;
					break;
				default:
					e.iInfo.ionWidthMethod = IonInfo::IonWidthManual;
					e.iInfo.ionWidth = 0.;
			}
		}
		else if(key=="Ecut") e.iInfo.ionWidthMethod = IonInfo::IonWidthEcut;
		else if(key=="fftbox") e.iInfo.ionWidthMethod = IonInfo::IonWidthFFTbox;
		else
		{	istringstream iss(key);
			iss >> e.iInfo.ionWidth;
			if(iss.fail()) throw string("<width> must be Ecut, fftbox or a value in bohrs");
			e.iInfo.ionWidthMethod = IonInfo::IonWidthManual;
		}
	}

	void printStatus(Everything& e, int iRep)
	{	switch(e.iInfo.ionWidthMethod)
		{	case IonInfo::IonWidthFFTbox: logPrintf("fftbox"); break;
			case IonInfo::IonWidthEcut: logPrintf("Ecut"); break;
			case IonInfo::IonWidthManual: logPrintf("%lg", e.iInfo.ionWidth); break;
		}
	}
}
commandIonWidth;
