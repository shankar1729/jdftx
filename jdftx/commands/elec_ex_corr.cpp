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


EnumStringMap<ExCorrType> exCorrTypeMap(
	ExCorrLDA_PZ, "lda", //Default LDA is PZ
	ExCorrLDA_PZ, "lda-PZ", 
	ExCorrLDA_PW, "lda-PW", 
	ExCorrLDA_PW_prec, "lda-PW-prec", 
	ExCorrLDA_VWN, "lda-VWN", 
	ExCorrLDA_Teter, "lda-Teter",
	ExCorrGGA_PBE, "gga", //Default GGA is PBE
	ExCorrGGA_PBE, "gga-PBE",
	ExCorrGGA_PBEsol, "gga-PBEsol",
	ExCorrGGA_PW91, "gga-PW91",
	ExCorrMGGA_TPSS, "mgga-TPSS",
	ExCorrMGGA_revTPSS, "mgga-revTPSS",
	ExCorrHYB_PBE0, "hyb-PBE0",
	ExCorrHYB_HSE06, "hyb-HSE06",
	ExCorrHF, "Hartree-Fock"
);
EnumStringMap<ExCorrType> exCorrDescriptionMap(
	ExCorrLDA_PZ, "Perdew-Zunger LDA", 
	ExCorrLDA_PW, "Perdew-Wang LDA", 
	ExCorrLDA_PW_prec, "Perdew-Wang LDA with extended precision (used by PBE)",
	ExCorrLDA_VWN, "Vosko-Wilk-Nusair LDA", 
	ExCorrLDA_Teter, "Teter93 LSDA",
	ExCorrGGA_PBE, "Perdew-Burke-Ernzerhof GGA",
	ExCorrGGA_PBEsol, "Perdew-Burke-Ernzerhof GGA reparametrized for solids",
	ExCorrGGA_PW91, "Perdew-Wang GGA",
	ExCorrMGGA_TPSS, "Tao-Perdew-Staroverov-Scuseria meta GGA",
	ExCorrMGGA_revTPSS, "revised Tao-Perdew-Staroverov-Scuseria meta GGA",
	ExCorrHYB_PBE0, "Hybrid PBE with 1/4 exact exchange",
	ExCorrHYB_HSE06, "HSE06 'wPBEh' hybrid with 1/4 screened exact exchange",
	ExCorrHF, "Full exact exchange with no correlation"
);

EnumStringMap<KineticType> kineticTypeMap(
	KineticTF, "lda-TF",
	KineticVW, "gga-vW", 
	KineticPW91, "gga-PW91k"
);
EnumStringMap<KineticType> kineticDescriptionMap(
	KineticTF, "Thomas-Fermi LDA kinetic energy", 
	KineticVW, "von Weisacker correction to LDA kinetic energy", 
	KineticPW91, "Perdew-Wang GGA kinetic energy parameterized by Lembarki and Chermette"
);


#ifdef LIBXC_ENABLED
#include <xc.h>

//! Exchange functionals
EnumStringMap<int> xcMap_X(
	XC_LDA_X_1D               , "lda-x-1d",
	XC_LDA_X_2D               , "lda-x-2d",
	XC_LDA_X                  , "lda-x",

	XC_GGA_X_DK87_R1          , "gga-x-dk87-r1",
	XC_GGA_X_DK87_R2          , "gga-x-dk87-r2",
	XC_GGA_X_FT97_A           , "gga-x-ft97-a",
	XC_GGA_X_FT97_B           , "gga-x-ft97-b",
	XC_GGA_X_G96              , "gga-x-g96",
	XC_GGA_X_LG93             , "gga-x-lg93",
	XC_GGA_X_MPBE             , "gga-x-mpbe",
	XC_GGA_X_mPW91            , "gga-x-mpw91",
	XC_GGA_X_OPTB88_VDW       , "gga-x-optb88-vdw",
	XC_GGA_X_OPTPBE_VDW       , "gga-x-optpbe-vdw",
	XC_GGA_X_OPTX             , "gga-x-optx",
	XC_GGA_X_PBEA             , "gga-x-pbea",
	XC_GGA_X_PBE              , "gga-x-pbe",
	XC_GGA_X_PBE_JSJR         , "gga-x-pbe-jsjr",
	XC_GGA_X_PBEK1_VDW        , "gga-x-pbek1-vdw",
	XC_GGA_X_PBE_R            , "gga-x-pbe-r",
	XC_GGA_X_PBE_SOL          , "gga-x-pbe-sol",
	XC_GGA_X_PW86             , "gga-x-pw86",
	XC_GGA_X_PW91             , "gga-x-pw91",
	XC_GGA_X_RGE2             , "gga-x-rge2",
	XC_GGA_X_RPBE             , "gga-x-rpbe",
	XC_GGA_X_WC               , "gga-x-wc",
	XC_GGA_X_XPBE             , "gga-x-xpbe",
	XC_GGA_X_2D_B86           , "gga-x-2d-b86",
	XC_GGA_X_2D_B86_MGC       , "gga-x-2d-b86-mgc",
	XC_GGA_X_2D_B88           , "gga-x-2d-b88",
	XC_GGA_X_2D_PBE           , "gga-x-2d-pbe",
	XC_GGA_X_AM05             , "gga-x-am05",
	XC_GGA_X_B86              , "gga-x-b86",
	XC_GGA_X_B86_MGC          , "gga-x-b86-mgc",
	XC_GGA_X_B88              , "gga-x-b88",
	XC_GGA_X_BAYESIAN         , "gga-x-bayesian",

	XC_MGGA_X_LTA             , "mgga-x-lta",
	XC_MGGA_X_TPSS            , "mgga-x-tpss",
	XC_MGGA_X_M06L            , "mgga-x-m06l",
	XC_MGGA_X_GVT4            , "mgga-x-gvt4",
	XC_MGGA_X_TAU_HCTH        , "mgga-x-tau-hcth",
	XC_MGGA_X_BR89            , "mgga-x-br89",
	XC_MGGA_X_BJ06            , "mgga-x-bj06",
	XC_MGGA_X_TB09            , "mgga-x-tb09",
	XC_MGGA_X_RPP09           , "mgga-x-rpp09"
);

//! Correlation functionals
EnumStringMap<int> xcMap_C(
	XC_LDA_C_1D_CSC           , "lda-c-1d-csc",
	XC_LDA_C_2D_AMGB          , "lda-c-2d-amgb",
	XC_LDA_C_2D_PRM           , "lda-c-2d-prm",
	XC_LDA_C_GL               , "lda-c-gl",
	XC_LDA_C_HL               , "lda-c-hl",
	XC_LDA_C_ML1              , "lda-c-ml1",
	XC_LDA_C_ML2              , "lda-c-ml2",
	XC_LDA_C_OB_PW            , "lda-c-ob-pw",
	XC_LDA_C_OB_PZ            , "lda-c-ob-pz",
	XC_LDA_C_PW               , "lda-c-pw",
	XC_LDA_C_PW_MOD           , "lda-c-pw-mod",
	XC_LDA_C_PZ               , "lda-c-pz",
	XC_LDA_C_PZ_MOD           , "lda-c-pz-mod",
	XC_LDA_C_RPA              , "lda-c-rpa",
	XC_LDA_C_vBH              , "lda-c-vbh",
	XC_LDA_C_VWN              , "lda-c-vwn",
	XC_LDA_C_VWN_RPA          , "lda-c-vwn-rpa",
	XC_LDA_C_WIGNER           , "lda-c-wigner",
	XC_LDA_C_XALPHA           , "lda-c-xalpha",

	XC_GGA_C_AM05             , "gga-c-am05",
	XC_GGA_C_LM               , "gga-c-lm",
	XC_GGA_C_LYP              , "gga-c-lyp",
	XC_GGA_C_P86              , "gga-c-p86",
	XC_GGA_C_PBE              , "gga-c-pbe",
	XC_GGA_C_PBE_JRGX         , "gga-c-pbe-jrgx",
	XC_GGA_C_PBE_SOL          , "gga-c-pbe-sol",
	XC_GGA_C_PW91             , "gga-c-pw91",
	XC_GGA_C_RGE2             , "gga-c-rge2",
	XC_GGA_C_XPBE             , "gga-c-xpbe",

	XC_MGGA_C_TPSS            , "mgga-c-tpss",
	XC_MGGA_C_VSXC            , "mgga-c-vsxc"
);

//! Combined exchange-correlation functionals
EnumStringMap<int> xcMap_XC(
	XC_LDA_XC_TETER93         , "lda-xc-teter93",

	XC_GGA_XC_B97_1           , "gga-xc-b97-1",
	XC_GGA_XC_B97_2           , "gga-xc-b97-2",
	XC_GGA_XC_B97_3           , "gga-xc-b97-3",
	XC_GGA_XC_B97_D           , "gga-xc-b97-d",
	XC_GGA_XC_B97             , "gga-xc-b97",
	XC_GGA_XC_B97_K           , "gga-xc-b97-k",
	XC_GGA_XC_EDF1            , "gga-xc-edf1",
	XC_GGA_XC_HCTH_120        , "gga-xc-hcth-120",
	XC_GGA_XC_HCTH_147        , "gga-xc-hcth-147",
	XC_GGA_XC_HCTH_407        , "gga-xc-hcth-407",
	XC_GGA_XC_HCTH_93         , "gga-xc-hcth-93",
	XC_GGA_XC_MPWLYP1W        , "gga-xc-mpwlyp1w",
	XC_GGA_XC_PBE1W           , "gga-xc-pbe1w",
	XC_GGA_XC_PBELYP1W        , "gga-xc-pbelyp1w",
	XC_GGA_XC_SB98_1a         , "gga-xc-sb98-1a",
	XC_GGA_XC_SB98_1b         , "gga-xc-sb98-1b",
	XC_GGA_XC_SB98_1c         , "gga-xc-sb98-1c",
	XC_GGA_XC_SB98_2a         , "gga-xc-sb98-2a",
	XC_GGA_XC_SB98_2b         , "gga-xc-sb98-2b",
	XC_GGA_XC_SB98_2c         , "gga-xc-sb98-2c",
	XC_GGA_XC_XLYP            , "gga-xc-xlyp",

	XC_HYB_GGA_XC_B3LYP       , "hyb-gga-xc-b3lyp",
	XC_HYB_GGA_XC_B3P86       , "hyb-gga-xc-b3p86",
	XC_HYB_GGA_XC_B3PW91      , "hyb-gga-xc-b3pw91",
	XC_HYB_GGA_XC_B97_1       , "hyb-gga-xc-b97-1",
	XC_HYB_GGA_XC_B97_2       , "hyb-gga-xc-b97-2",
	XC_HYB_GGA_XC_B97_3       , "hyb-gga-xc-b97-3",
	XC_HYB_GGA_XC_B97         , "hyb-gga-xc-b97",
	XC_HYB_GGA_XC_B97_K       , "hyb-gga-xc-b97-k",
	XC_HYB_GGA_XC_mPW3LYP     , "hyb-gga-xc-mpw3lyp",
	XC_HYB_GGA_XC_mPW3PW      , "hyb-gga-xc-mpw3pw",
	XC_HYB_GGA_XC_O3LYP       , "hyb-gga-xc-o3lyp",
	XC_HYB_GGA_XC_PBEH        , "hyb-gga-xc-pbeh",
	XC_HYB_GGA_XC_SB98_1a     , "hyb-gga-xc-sb98-1a",
	XC_HYB_GGA_XC_SB98_1b     , "hyb-gga-xc-sb98-1b",
	XC_HYB_GGA_XC_SB98_1c     , "hyb-gga-xc-sb98-1c",
	XC_HYB_GGA_XC_SB98_2a     , "hyb-gga-xc-sb98-2a",
	XC_HYB_GGA_XC_SB98_2b     , "hyb-gga-xc-sb98-2b",
	XC_HYB_GGA_XC_SB98_2c     , "hyb-gga-xc-sb98-2c",
	XC_HYB_GGA_XC_X3LYP       , "hyb-gga-xc-x3lyp"
);

//Get description by temporarily initializing functional:
string getLibXCdescription(const string& name, const EnumStringMap<int>& map)
{	int xcCode = 0;
	bool xcFound = map.getEnum(name.c_str(), xcCode);
	assert(xcFound && xcCode);
	xc_func_type func;
	if(xc_func_init(&func, xcCode, XC_UNPOLARIZED) != 0)
		die("Error obtaining description for LibXC functional %s.\n", name.c_str());
	string desc(func.info->name);
	xc_func_end(&func);
	return desc;
}
string getLibXCdescription_X(const string& name) { return getLibXCdescription(name, xcMap_X); }
string getLibXCdescription_C(const string& name) { return getLibXCdescription(name, xcMap_C); }
string getLibXCdescription_XC(const string& name) { return getLibXCdescription(name, xcMap_XC); }
#endif //LIBXC_ENABLED


struct CommandElecExCorr : public Command
{
	CommandElecExCorr(const char* cmdName = "elec-ex-corr") : Command(cmdName)
	{
		format = "<functional>";
		comments = "Specify the exchange-correlation functional, where <functional> is one of:"
			+ addDescriptions(exCorrTypeMap.optionList(), linkDescription(exCorrTypeMap, exCorrDescriptionMap))
			+ ".\nNote that lda is an alias for lda-pz, and gga for gga-pbe.\n";
		hasDefault = true;
		emptyParamError = "   eXchange/Correlation functional(s) must be specified.";
		
		#ifdef LIBXC_ENABLED
		format += "\n\t| <funcX> <funcC>\n\t| <funcXC>";
		comments +=
			"The second and third lines use eXchange/Correlation functionals from libXC\n"
			"(the internal threaded/gpu implementations above are usually much faster).\n"
			"Here, <funcX> is one of:"
			+ addDescriptions(xcMap_X.optionList(), getLibXCdescription_X)
			+ ",\n<funcC> is one of:"
			+ addDescriptions(xcMap_C.optionList(), getLibXCdescription_C)
			+ ",\nand <funcXC> is one of:"
			+ addDescriptions(xcMap_XC.optionList(), getLibXCdescription_XC)
			+ ".";
		#else
		comments += "Additional functionals can be enabled by compiling with LibXC support.";
		#endif
	}

	void process(ParamList& pl, Everything& e)
	{	process(pl, e.exCorr);
	}
	
	void printStatus(Everything& e, int iRep)
	{	printStatus(e.exCorr);
	}
	
protected:
	void process(ParamList& pl, ExCorr& exCorr)
	{	string key;
		pl.get(key, string(), "functional");
		if(!key.length()) return;
		else
		{	if(exCorrTypeMap.getEnum(key.c_str(), exCorr.exCorrType)) //Found internal ExCorr functional
			{	//Set the functional name:
				exCorr.xcName = string(exCorrTypeMap.getString(exCorr.exCorrType));
			}
			#ifdef LIBXC_ENABLED
			else if(xcMap_X.getEnum(key.c_str(), exCorr.xcExchange)) {} //Found LibXC Exchange functional
			else if(xcMap_XC.getEnum(key.c_str(), exCorr.xcExcorr)) {} //Found LibXC ExCorr functional
			#endif
			else throw key + " is not a recognized exchange or exchange-correlation functional";
		}
		#ifdef LIBXC_ENABLED
		if(exCorr.xcExchange || exCorr.xcExcorr)
		{	//LibXC will be used:
			exCorr.exCorrType = ExCorrLibXC;
			if(exCorr.xcExchange)
				pl.get(exCorr.xcCorr, 0, xcMap_C, "funcC", true); //required correlation functional
			//Set the short name:
			if(exCorr.xcExcorr)
				exCorr.xcName = string(xcMap_XC.getString(exCorr.xcExcorr));
			else
				exCorr.xcName = string(xcMap_X.getString(exCorr.xcExchange))
					+ ':' + string(xcMap_C.getString(exCorr.xcCorr));
		}
		#endif
	}

	void printStatus(ExCorr& exCorr)
	{	switch(exCorr.exCorrType)
		{
			#ifdef LIBXC_ENABLED
			case ExCorrLibXC:
			{	if(exCorr.xcExcorr)
					logPrintf("%s", xcMap_XC.getString(exCorr.xcExcorr));
				else
					logPrintf("%s %s", xcMap_X.getString(exCorr.xcExchange), xcMap_C.getString(exCorr.xcCorr));
				break;
			}
			#endif
			default:
				logPrintf("%s", exCorrTypeMap.getString(exCorr.exCorrType));
		}
	}
}
commandElecExCorr;


struct CommandElecExCorrCompare : public CommandElecExCorr
{
	CommandElecExCorrCompare() : CommandElecExCorr("elec-ex-corr-compare")
	{
		format = "<functional>";
		comments =
			"Compute total energies for other functionals at the final state for comparison.\n"
			"This command may be specified multiple times. It invokes 'dump End ExcCompare'\n"
			"automatically, but the compute frequency can be controlled using dump explicitly.\n"
			"The available options for each parameter are identical to elec-ex-corr";
		hasDefault = false;
		allowMultiple = true;
		emptyParamError = "   eXchange/Correlation functional(s) must be specified.";
		
		#ifdef LIBXC_ENABLED
		format += "\n\t| <funcX> <funcC>\n\t| <funcXC>";
		#endif
		forbid("fix-electron-density");
	}
	
	void process(ParamList& pl, Everything& e)
	{	e.exCorrDiff.push_back(std::shared_ptr<ExCorr>(new ExCorr));
		CommandElecExCorr::process(pl, *e.exCorrDiff.back());
		e.dump.insert(std::make_pair(DumpFreq_End, DumpExcCompare));
	}
	
	void printStatus(Everything& e, int iRep)
	{	CommandElecExCorr::printStatus(*e.exCorrDiff[iRep]);
	}
}
commandElecExCorrCompare;


struct CommandFluidExCorr : public CommandElecExCorr
{
	CommandFluidExCorr() : CommandElecExCorr("fluid-ex-corr")
	{
		format = "<kinetic> [<exchange-correlation>]";
		comments =
			"Kinetic energy functional for fluid convolution coupling where <kinetic> is one of:"
			+ addDescriptions(kineticTypeMap.optionList(), linkDescription(kineticTypeMap, kineticDescriptionMap)) +
			".\nThe available options for <exchange-correlation> are identical to elec-ex-corr\n"
			"and defaults to the functional set by elec-ex-corr.";
		hasDefault = true;
		emptyParamError = "   A kinetic energy functional must be specified.";
		require("elec-ex-corr");
	}
	
	void process(ParamList& pl, Everything& e)
	{	ExCorr& fluidExCorr = e.eVars.fluidParams.exCorr;
		pl.get(fluidExCorr.kineticType, KineticTF, kineticTypeMap, "kinetic");
		//Set default exchange-correlation to be that of elec-ex-corr
		fluidExCorr.exCorrType = e.exCorr.exCorrType;
		fluidExCorr.xcName = e.exCorr.xcName;
		#ifdef LIBXC_ENABLED
		fluidExCorr.xcExchange = e.exCorr.xcExchange;
		fluidExCorr.xcCorr = e.exCorr.xcCorr;
		fluidExCorr.xcExcorr = e.exCorr.xcExcorr;
		#endif
		CommandElecExCorr::process(pl, fluidExCorr);
	}
	
	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s ", kineticTypeMap.getString(e.eVars.fluidParams.exCorr.kineticType));
		CommandElecExCorr::printStatus(e.eVars.fluidParams.exCorr);
	}
}
commandFluidExCorr;
