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
(	PCM_SG14NL,  "SG14NL",
	PCM_SGA13,   "SGA13", 
	PCM_GLSSA13, "GLSSA13",
	PCM_LA12,    "LA12", 
	PCM_PRA05,   "PRA05"
);
EnumStringMap<PCMVariant> pcmVariantDescMap
(	PCM_SG14NL,  "Charge-asymmetry corrected, local-response, nonlocal-cavity PCM with weighted-density cavitation and dispersion [EXPERIMENTAL]",
	PCM_SGA13,   "PCM with weighted-density cavitation and dispersion [R. Sundararaman, D. Gunceler and T.A. Arias, (under preparation)]", 
	PCM_GLSSA13, "PCM with empirical cavity tension [D. Gunceler, K. Letchworth-Weaver, R. Sundararaman, K.A. Schwarz and T.A. Arias, arXiv:1301.6189]",
	PCM_LA12,    "PCM with no cavitation/dispersion contributions [K. Letchworth-Weaver and T.A. Arias, Phys. Rev. B 86, 075140 (2012)]", 
	PCM_PRA05,   "PCM with no cavitation/dispersion contributions [S.A. Petrosyan SA, A.A. Rigos and T.A. Arias, J Phys Chem B. 109, 15436 (2005)]"
);

struct CommandPcmVariant : public Command
{
	CommandPcmVariant() : Command("pcm-variant", "Fluid parameters")
	{
		format = "[<variant>=GLSSA13]";
		comments = "Select LinearPCM or NonlinearPCM <variant> amongst:"
			+ addDescriptions(pcmVariantMap.optionList(), linkDescription(pcmVariantMap, pcmVariantDescMap));
		hasDefault = true;
		require("fluid");
	}

	void process(ParamList& pl, Everything& e)
	{	FluidSolverParams& fsp = e.eVars.fluidParams;
		if(fsp.fluidType==FluidSaLSA) fsp.pcmVariant = PCM_SaLSA; //only option for SaLSA
		else
		{	pl.get(fsp.pcmVariant, PCM_GLSSA13, pcmVariantMap, "variant");
			if(fsp.pcmVariant==PCM_SG14NL && fsp.fluidType!=FluidLinearPCM)
				throw string("SG14NL can only be used with fluid LinearPCM");
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
	PCMp_eta_wDiel, //!< fit parameter for dielectric cavity in SG14NL
	PCMp_sqrtC6eff, //!< sqrt(effective molecule C6 coefficient) for SG14NL
	PCMp_pCavity, //!< sensitivity of cavity to surface electric fields [e-a0/Eh] in SG14NL
	PCMp_Delim //!< Delimiter used in parsing
};
EnumStringMap<PCMparameter> pcmParamMap
(	PCMp_lMax,          "lMax",
	PCMp_nc,            "nc",
	PCMp_sigma,         "sigma",
	PCMp_cavityTension, "cavityTension",
	PCMp_eta_wDiel, "eta_wDiel",
	PCMp_sqrtC6eff, "sqrtC6eff",
	PCMp_pCavity, "pCavity"
);
EnumStringMap<PCMparameter> pcmParamDescMap
(	PCMp_lMax, "angular momentum truncation in SaLSA",
	PCMp_nc, "critical density for the PCM cavity shape function",
	PCMp_sigma, "smoothing factor for the PCM cavity shape function",
	PCMp_cavityTension, "effective surface tension (including dispersion etc.) of the cavity (hartree per bohr^2)",
	PCMp_eta_wDiel, "fit parameter for dielectric cavity in SG14NL",
	PCMp_sqrtC6eff, "sqrt(effective molecule C6 coefficient) for SG14NL",
	PCMp_pCavity, "sensitivity of cavity to surface electric fields [e-a0/Eh] in SG14NL"
);

struct CommandPcmParams : public Command
{
	CommandPcmParams() : Command("pcm-params", "Fluid parameters")
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
			{	READ_AND_CHECK(lMax, >=, 0)
				READ_AND_CHECK(nc, >, 0.)
				READ_AND_CHECK(sigma, >, 0.)
				READ_AND_CHECK(cavityTension, <, DBL_MAX)
				READ_AND_CHECK(eta_wDiel, >=, 0.)
				READ_AND_CHECK(sqrtC6eff, >=, 0.)
				READ_AND_CHECK(pCavity, <, DBL_MAX)
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
		PRINT(eta_wDiel)
		PRINT(sqrtC6eff)
		PRINT(pCavity)
		#undef PRINT
	}
}
commandPcmParams;


struct CommandPCMnonlinearDebug : public Command
{
    CommandPCMnonlinearDebug() : Command("pcm-nonlinear-debug", "Fluid parameters")
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



struct CommandIonWidth : public Command
{
	CommandIonWidth() : Command("ion-width", "Fluid parameters")
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
