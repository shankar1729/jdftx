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
	ExCorrORB_GLLBsc, "orb-GLLBsc",
	ExCorrPOT_LB94, "pot-LB94",
	ExCorrHYB_PBE0, "hyb-PBE0",
	ExCorrHYB_HSE06, "hyb-HSE06",
	ExCorrHYB_HSE12, "hyb-HSE12",
	ExCorrHYB_HSE12s, "hyb-HSE12s",
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
	ExCorrORB_GLLBsc, "Orbital-dependent GLLB-sc potential (no total energy)",
	ExCorrPOT_LB94, "van Leeuwen-Baerends model potential (no total energy)",
	ExCorrHYB_PBE0, "Hybrid PBE with 1/4 exact exchange",
	ExCorrHYB_HSE06, "HSE06 'wPBEh' hybrid with 1/4 screened exact exchange",
	ExCorrHYB_HSE12, "Reparametrized screened exchange functional for accuracy (w=0.185 A^-1 and a=0.313)",
	ExCorrHYB_HSE12s, "Reparametrized screened exchange functional for k-point convergence (w=0.408 A^-1 and a=0.425)",
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
EnumStringMap<int> xcMap_X
(	XC_LDA_X                      , "lda-x",
	XC_LDA_X_2D                   , "lda-x-2d",
	XC_LDA_X_1D                   , "lda-x-1d",

	XC_GGA_X_SSB_SW               , "gga-x-ssb-sw",
	XC_GGA_X_SSB                  , "gga-x-ssb",
	XC_GGA_X_SSB_D                , "gga-x-ssb-d",
	XC_GGA_X_BPCCAC               , "gga-x-bpccac",
	XC_GGA_X_PBE                  , "gga-x-pbe",
	XC_GGA_X_PBE_R                , "gga-x-pbe-r",
	XC_GGA_X_B86                  , "gga-x-b86",
	XC_GGA_X_HERMAN               , "gga-x-herman",
	XC_GGA_X_B86_MGC              , "gga-x-b86-mgc",
	XC_GGA_X_B88                  , "gga-x-b88",
	XC_GGA_X_G96                  , "gga-x-g96",
	XC_GGA_X_PW86                 , "gga-x-pw86",
	XC_GGA_X_PW91                 , "gga-x-pw91",
	XC_GGA_X_OPTX                 , "gga-x-optx",
	XC_GGA_X_DK87_R1              , "gga-x-dk87-r1",
	XC_GGA_X_DK87_R2              , "gga-x-dk87-r2",
	XC_GGA_X_LG93                 , "gga-x-lg93",
	XC_GGA_X_FT97_A               , "gga-x-ft97-a",
	XC_GGA_X_FT97_B               , "gga-x-ft97-b",
	XC_GGA_X_PBE_SOL              , "gga-x-pbe-sol",
	XC_GGA_X_RPBE                 , "gga-x-rpbe",
	XC_GGA_X_WC                   , "gga-x-wc",
	XC_GGA_X_MPW91                , "gga-x-mpw91",
	XC_GGA_X_AM05                 , "gga-x-am05",
	XC_GGA_X_PBEA                 , "gga-x-pbea",
	XC_GGA_X_MPBE                 , "gga-x-mpbe",
	XC_GGA_X_XPBE                 , "gga-x-xpbe",
	XC_GGA_X_2D_B86_MGC           , "gga-x-2d-b86-mgc",
	XC_GGA_X_BAYESIAN             , "gga-x-bayesian",
	XC_GGA_X_PBE_JSJR             , "gga-x-pbe-jsjr",
	XC_GGA_X_2D_B88               , "gga-x-2d-b88",
	XC_GGA_X_2D_B86               , "gga-x-2d-b86",
	XC_GGA_X_2D_PBE               , "gga-x-2d-pbe",
	XC_GGA_X_OPTB88_VDW           , "gga-x-optb88-vdw",
	XC_GGA_X_PBEK1_VDW            , "gga-x-pbek1-vdw",
	XC_GGA_X_OPTPBE_VDW           , "gga-x-optpbe-vdw",
	XC_GGA_X_RGE2                 , "gga-x-rge2",
	XC_GGA_X_RPW86                , "gga-x-rpw86",
	XC_GGA_X_KT1                  , "gga-x-kt1",
	XC_GGA_X_MB88                 , "gga-x-mb88",
	XC_GGA_X_SOGGA                , "gga-x-sogga",
	XC_GGA_X_SOGGA11              , "gga-x-sogga11",
	XC_GGA_X_C09X                 , "gga-x-c09x",
	XC_GGA_X_LB                   , "gga-x-lb",
	XC_GGA_X_LBM                  , "gga-x-lbm",
	XC_GGA_X_OL2                  , "gga-x-ol2",
	XC_GGA_X_APBE                 , "gga-x-apbe",
	XC_GGA_X_HTBS                 , "gga-x-htbs",
	XC_GGA_X_AIRY                 , "gga-x-airy",
	XC_GGA_X_LAG                  , "gga-x-lag",
	XC_GGA_X_WPBEH                , "gga-x-wpbeh",
	XC_GGA_X_HJS_PBE              , "gga-x-hjs-pbe",
	XC_GGA_X_HJS_PBE_SOL          , "gga-x-hjs-pbe-sol",
	XC_GGA_X_HJS_B88              , "gga-x-hjs-b88",
	XC_GGA_X_HJS_B97X             , "gga-x-hjs-b97x",
	XC_GGA_X_ITYH                 , "gga-x-ityh",

	XC_MGGA_X_LTA                 , "mgga-x-lta",
	XC_MGGA_X_TPSS                , "mgga-x-tpss",
	XC_MGGA_X_M06_L               , "mgga-x-m06-l",
	XC_MGGA_X_GVT4                , "mgga-x-gvt4",
	XC_MGGA_X_TAU_HCTH            , "mgga-x-tau-hcth",
	XC_MGGA_X_BR89                , "mgga-x-br89",
	XC_MGGA_X_BJ06                , "mgga-x-bj06",
	XC_MGGA_X_TB09                , "mgga-x-tb09",
	XC_MGGA_X_RPP09               , "mgga-x-rpp09",
	XC_MGGA_X_2D_PRHG07           , "mgga-x-2d-prhg07",
	XC_MGGA_X_2D_PRHG07_PRP10     , "mgga-x-2d-prhg07-prp10",
	XC_MGGA_X_REVTPSS             , "mgga-x-revtpss",
	XC_MGGA_X_PKZB                , "mgga-x-pkzb",
	XC_MGGA_X_M05                 , "mgga-x-m05",
	XC_MGGA_X_M05_2X              , "mgga-x-m05-2x",
	XC_MGGA_X_M06_HF              , "mgga-x-m06-hf",
	XC_MGGA_X_M06                 , "mgga-x-m06",
	XC_MGGA_X_M06_2X              , "mgga-x-m06-2x",
	XC_MGGA_X_M08_HX              , "mgga-x-m08-hx",
	XC_MGGA_X_M08_SO              , "mgga-x-m08-so"
);

//! Correlation functionals
EnumStringMap<int> xcMap_C
(	XC_LDA_C_WIGNER               , "lda-c-wigner",
	XC_LDA_C_RPA                  , "lda-c-rpa",
	XC_LDA_C_HL                   , "lda-c-hl",
	XC_LDA_C_GL                   , "lda-c-gl",
	XC_LDA_C_XALPHA               , "lda-c-xalpha",
	XC_LDA_C_VWN                  , "lda-c-vwn",
	XC_LDA_C_VWN_RPA              , "lda-c-vwn-rpa",
	XC_LDA_C_PZ                   , "lda-c-pz",
	XC_LDA_C_PZ_MOD               , "lda-c-pz-mod",
	XC_LDA_C_OB_PZ                , "lda-c-ob-pz",
	XC_LDA_C_PW                   , "lda-c-pw",
	XC_LDA_C_PW_MOD               , "lda-c-pw-mod",
	XC_LDA_C_OB_PW                , "lda-c-ob-pw",
	XC_LDA_C_2D_AMGB              , "lda-c-2d-amgb",
	XC_LDA_C_2D_PRM               , "lda-c-2d-prm",
	XC_LDA_C_vBH                  , "lda-c-vbh",
	XC_LDA_C_1D_CSC               , "lda-c-1d-csc",
	XC_LDA_C_ML1                  , "lda-c-ml1",
	XC_LDA_C_ML2                  , "lda-c-ml2",
	XC_LDA_C_GOMBAS               , "lda-c-gombas",
	XC_LDA_C_PW_RPA               , "lda-c-pw-rpa",
	XC_LDA_C_1D_LOOS              , "lda-c-1d-loos",
	XC_LDA_C_RC04                 , "lda-c-rc04",
	XC_LDA_C_VWN_1                , "lda-c-vwn-1",
	XC_LDA_C_VWN_2                , "lda-c-vwn-2",
	XC_LDA_C_VWN_3                , "lda-c-vwn-3",
	XC_LDA_C_VWN_4                , "lda-c-vwn-4",

	XC_GGA_C_OP_XALPHA            , "gga-c-op-xalpha",
	XC_GGA_C_OP_G96               , "gga-c-op-g96",
	XC_GGA_C_OP_PBE               , "gga-c-op-pbe",
	XC_GGA_C_OP_B88               , "gga-c-op-b88",
	XC_GGA_C_FT97                 , "gga-c-ft97",
	XC_GGA_C_SPBE                 , "gga-c-spbe",
	XC_GGA_C_REVTCA               , "gga-c-revtca",
	XC_GGA_C_TCA                  , "gga-c-tca",
	XC_GGA_C_PBE                  , "gga-c-pbe",
	XC_GGA_C_LYP                  , "gga-c-lyp",
	XC_GGA_C_P86                  , "gga-c-p86",
	XC_GGA_C_PBE_SOL              , "gga-c-pbe-sol",
	XC_GGA_C_PW91                 , "gga-c-pw91",
	XC_GGA_C_AM05                 , "gga-c-am05",
	XC_GGA_C_XPBE                 , "gga-c-xpbe",
	XC_GGA_C_LM                   , "gga-c-lm",
	XC_GGA_C_PBE_JRGX             , "gga-c-pbe-jrgx",
	XC_GGA_C_RGE2                 , "gga-c-rge2",
	XC_GGA_C_WL                   , "gga-c-wl",
	XC_GGA_C_WI                   , "gga-c-wi",
	XC_GGA_C_SOGGA11              , "gga-c-sogga11",
	XC_GGA_C_WI0                  , "gga-c-wi0",
	XC_GGA_C_SOGGA11_X            , "gga-c-sogga11-x",
	XC_GGA_C_APBE                 , "gga-c-apbe",
	XC_GGA_C_OPTC                 , "gga-c-optc",

	XC_MGGA_C_TPSS                , "mgga-c-tpss",
	XC_MGGA_C_VSXC                , "mgga-c-vsxc",
	XC_MGGA_C_M06_L               , "mgga-c-m06-l",
	XC_MGGA_C_M06_HF              , "mgga-c-m06-hf",
	XC_MGGA_C_M06                 , "mgga-c-m06",
	XC_MGGA_C_M06_2X              , "mgga-c-m06-2x",
	XC_MGGA_C_M05                 , "mgga-c-m05",
	XC_MGGA_C_M05_2X              , "mgga-c-m05-2x",
	XC_MGGA_C_PKZB                , "mgga-c-pkzb",
	XC_MGGA_C_BC95                , "mgga-c-bc95"
);

//! Combined exchange-correlation functionals
EnumStringMap<int> xcMap_XC
(	XC_LDA_XC_TETER93             , "lda-xc-teter93",

	XC_GGA_XC_HCTH_407P           , "gga-xc-hcth-407p",
	XC_GGA_XC_HCTH_P76            , "gga-xc-hcth-p76",
	XC_GGA_XC_HCTH_P14            , "gga-xc-hcth-p14",
	XC_GGA_XC_B97_GGA1            , "gga-xc-b97-gga1",
	XC_GGA_XC_HCTH_A              , "gga-xc-hcth-a",
	XC_GGA_XC_KT2                 , "gga-xc-kt2",
	XC_GGA_XC_TH1                 , "gga-xc-th1",
	XC_GGA_XC_TH2                 , "gga-xc-th2",
	XC_GGA_XC_TH3                 , "gga-xc-th3",
	XC_GGA_XC_TH4                 , "gga-xc-th4",
	XC_GGA_XC_HCTH_93             , "gga-xc-hcth-93",
	XC_GGA_XC_HCTH_120            , "gga-xc-hcth-120",
	XC_GGA_XC_HCTH_147            , "gga-xc-hcth-147",
	XC_GGA_XC_HCTH_407            , "gga-xc-hcth-407",
	XC_GGA_XC_EDF1                , "gga-xc-edf1",
	XC_GGA_XC_XLYP                , "gga-xc-xlyp",
	XC_GGA_XC_B97                 , "gga-xc-b97",
	XC_GGA_XC_B97_1               , "gga-xc-b97-1",
	XC_GGA_XC_B97_2               , "gga-xc-b97-2",
	XC_GGA_XC_B97_D               , "gga-xc-b97-d",
	XC_GGA_XC_B97_K               , "gga-xc-b97-k",
	XC_GGA_XC_B97_3               , "gga-xc-b97-3",
	XC_GGA_XC_PBE1W               , "gga-xc-pbe1w",
	XC_GGA_XC_MPWLYP1W            , "gga-xc-mpwlyp1w",
	XC_GGA_XC_PBELYP1W            , "gga-xc-pbelyp1w",
	XC_GGA_XC_SB98_1a             , "gga-xc-sb98-1a",
	XC_GGA_XC_SB98_1b             , "gga-xc-sb98-1b",
	XC_GGA_XC_SB98_1c             , "gga-xc-sb98-1c",
	XC_GGA_XC_SB98_2a             , "gga-xc-sb98-2a",
	XC_GGA_XC_SB98_2b             , "gga-xc-sb98-2b",
	XC_GGA_XC_SB98_2c             , "gga-xc-sb98-2c",
	XC_GGA_XC_MOHLYP              , "gga-xc-mohlyp",
	XC_GGA_XC_MOHLYP2             , "gga-xc-mohlyp2",
	XC_GGA_XC_TH_FL               , "gga-xc-th-fl",
	XC_GGA_XC_TH_FC               , "gga-xc-th-fc",
	XC_GGA_XC_TH_FCFO             , "gga-xc-th-fcfo",
	XC_GGA_XC_TH_FCO              , "gga-xc-th-fco",

	XC_HYB_GGA_XC_B3PW91          , "hyb-gga-xc-b3pw91",
	XC_HYB_GGA_XC_B3LYP           , "hyb-gga-xc-b3lyp",
	XC_HYB_GGA_XC_B3P86           , "hyb-gga-xc-b3p86",
	XC_HYB_GGA_XC_O3LYP           , "hyb-gga-xc-o3lyp",
	XC_HYB_GGA_XC_mPW1K           , "hyb-gga-xc-mpw1k",
	XC_HYB_GGA_XC_PBEH            , "hyb-gga-xc-pbeh",
	XC_HYB_GGA_XC_B97             , "hyb-gga-xc-b97",
	XC_HYB_GGA_XC_B97_1           , "hyb-gga-xc-b97-1",
	XC_HYB_GGA_XC_B97_2           , "hyb-gga-xc-b97-2",
	XC_HYB_GGA_XC_X3LYP           , "hyb-gga-xc-x3lyp",
	XC_HYB_GGA_XC_B1WC            , "hyb-gga-xc-b1wc",
	XC_HYB_GGA_XC_B97_K           , "hyb-gga-xc-b97-k",
	XC_HYB_GGA_XC_B97_3           , "hyb-gga-xc-b97-3",
	XC_HYB_GGA_XC_MPW3PW          , "hyb-gga-xc-mpw3pw",
	XC_HYB_GGA_XC_B1LYP           , "hyb-gga-xc-b1lyp",
	XC_HYB_GGA_XC_B1PW91          , "hyb-gga-xc-b1pw91",
	XC_HYB_GGA_XC_mPW1PW          , "hyb-gga-xc-mpw1pw",
	XC_HYB_GGA_XC_MPW3LYP         , "hyb-gga-xc-mpw3lyp",
	XC_HYB_GGA_XC_SB98_1a         , "hyb-gga-xc-sb98-1a",
	XC_HYB_GGA_XC_SB98_1b         , "hyb-gga-xc-sb98-1b",
	XC_HYB_GGA_XC_SB98_1c         , "hyb-gga-xc-sb98-1c",
	XC_HYB_GGA_XC_SB98_2a         , "hyb-gga-xc-sb98-2a",
	XC_HYB_GGA_XC_SB98_2b         , "hyb-gga-xc-sb98-2b",
	XC_HYB_GGA_XC_SB98_2c         , "hyb-gga-xc-sb98-2c",
	XC_HYB_GGA_XC_HSE03           , "hyb-gga-xc-hse03",
	XC_HYB_GGA_XC_HSE06           , "hyb-gga-xc-hse06",
	XC_HYB_GGA_XC_HJS_PBE         , "hyb-gga-xc-hjs-pbe",
	XC_HYB_GGA_XC_HJS_PBE_SOL     , "hyb-gga-xc-hjs-pbe-sol",
	XC_HYB_GGA_XC_HJS_B88         , "hyb-gga-xc-hjs-b88",
	XC_HYB_GGA_XC_HJS_B97X        , "hyb-gga-xc-hjs-b97x",
	XC_HYB_GGA_XC_CAM_B3LYP       , "hyb-gga-xc-cam-b3lyp",
	XC_HYB_GGA_XC_TUNED_CAM_B3LYP , "hyb-gga-xc-tuned-cam-b3lyp",
	XC_HYB_GGA_XC_BHANDH          , "hyb-gga-xc-bhandh",
	XC_HYB_GGA_XC_BHANDHLYP       , "hyb-gga-xc-bhandhlyp",
	XC_HYB_GGA_XC_MB3LYP_RC04     , "hyb-gga-xc-mb3lyp-rc04",

	XC_HYB_MGGA_XC_M05            , "hyb-mgga-xc-m05",
	XC_HYB_MGGA_XC_M05_2X         , "hyb-mgga-xc-m05-2x",
	XC_HYB_MGGA_XC_B88B95         , "hyb-mgga-xc-b88b95",
	XC_HYB_MGGA_XC_B86B95         , "hyb-mgga-xc-b86b95",
	XC_HYB_MGGA_XC_PW86B95        , "hyb-mgga-xc-pw86b95",
	XC_HYB_MGGA_XC_BB1K           , "hyb-mgga-xc-bb1k",
	XC_HYB_MGGA_XC_M06_HF         , "hyb-mgga-xc-m06-hf",
	XC_HYB_MGGA_XC_MPW1B95        , "hyb-mgga-xc-mpw1b95",
	XC_HYB_MGGA_XC_MPWB1K         , "hyb-mgga-xc-mpwb1k",
	XC_HYB_MGGA_XC_X1B95          , "hyb-mgga-xc-x1b95",
	XC_HYB_MGGA_XC_XB1K           , "hyb-mgga-xc-xb1k",
	XC_HYB_MGGA_XC_M06            , "hyb-mgga-xc-m06",
	XC_HYB_MGGA_XC_M06_2X         , "hyb-mgga-xc-m06-2x"
);

//! Kinetic energy functionals
EnumStringMap<int> xcMap_K
(	XC_LDA_K_TF                   , "lda-k-tf",
	XC_LDA_K_LP                   , "lda-k-lp",

	XC_GGA_K_APBE                 , "gga-k-apbe",
	XC_GGA_K_TW1                  , "gga-k-tw1",
	XC_GGA_K_TW2                  , "gga-k-tw2",
	XC_GGA_K_TW3                  , "gga-k-tw3",
	XC_GGA_K_TW4                  , "gga-k-tw4",
	XC_GGA_K_VW                   , "gga-k-vw",
	XC_GGA_K_GE2                  , "gga-k-ge2",
	XC_GGA_K_GOLDEN               , "gga-k-golden",
	XC_GGA_K_YT65                 , "gga-k-yt65",
	XC_GGA_K_BALTIN               , "gga-k-baltin",
	XC_GGA_K_LIEB                 , "gga-k-lieb",
	XC_GGA_K_ABSR1                , "gga-k-absr1",
	XC_GGA_K_ABSR2                , "gga-k-absr2",
	XC_GGA_K_GR                   , "gga-k-gr",
	XC_GGA_K_LUDENA               , "gga-k-ludena",
	XC_GGA_K_GP85                 , "gga-k-gp85",
	XC_GGA_K_PEARSON              , "gga-k-pearson",
	XC_GGA_K_OL1                  , "gga-k-ol1",
	XC_GGA_K_OL2                  , "gga-k-ol2",
	XC_GGA_K_FR_B88               , "gga-k-fr-b88",
	XC_GGA_K_FR_PW86              , "gga-k-fr-pw86",
	XC_GGA_K_DK                   , "gga-k-dk",
	XC_GGA_K_PERDEW               , "gga-k-perdew",
	XC_GGA_K_VSK                  , "gga-k-vsk",
	XC_GGA_K_VJKS                 , "gga-k-vjks",
	XC_GGA_K_ERNZERHOF            , "gga-k-ernzerhof",
	XC_GGA_K_LC94                 , "gga-k-lc94",
	XC_GGA_K_LLP                  , "gga-k-llp",
	XC_GGA_K_THAKKAR              , "gga-k-thakkar"
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
string getLibXCdescription_K(const string& name) { return getLibXCdescription(name, xcMap_K); }
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
		if(key.length()) //Otherwise default functional set by ExCorr constructor
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
		if((exCorr.xcExchange || exCorr.xcExcorr) && !exCorr.xcCorr)
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

	void printStatus(const ExCorr& exCorr)
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
		forbid("fix-electron-potential");
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
			#ifdef LIBXC_ENABLED
			addDescriptions(xcMap_K.optionList(), getLibXCdescription_K) +
			#endif
			".\nThe available options for <exchange-correlation> are identical to elec-ex-corr\n"
			"and defaults to lda-pz.";
		hasDefault = true;
		emptyParamError = "   A kinetic energy functional must be specified.";
		require("elec-ex-corr");
	}
	
	void process(ParamList& pl, Everything& e)
	{	ExCorr& fluidExCorr = e.eVars.fluidParams.exCorr;
		//Get kinetic energy functional:
		string key; pl.get(key, string(), "kinetic");
		if(!key.length()) fluidExCorr.kineticType = KineticTF; //default: Thomas-Fermi
		else
		{	if(kineticTypeMap.getEnum(key.c_str(), fluidExCorr.kineticType)) {} //Found internal kinetic functional
			#ifdef LIBXC_ENABLED
			else if(xcMap_K.getEnum(key.c_str(), fluidExCorr.xcKinetic)) { fluidExCorr.kineticType = KineticLibXC; } //Found LibXC kinetic functional
			#endif
			else throw key + " is not a recognized kinetic energy functional";
		}
		//Set default exchange-correlation to be LDA
		fluidExCorr.exCorrType =  ExCorrLDA_PZ;
		fluidExCorr.xcName = exCorrTypeMap.getString(ExCorrLDA_PZ);
		  

		//old code which set default exchange correlation to be that of electronic system
		/*	fluidExCorr.exCorrType = e.exCorr.exCorrType;
			fluidExCorr.xcName = e.exCorr.xcName;		
		#ifdef LIBXC_ENABLED
		fluidExCorr.xcExchange = e.exCorr.xcExchange;
		fluidExCorr.xcCorr = e.exCorr.xcCorr;
		fluidExCorr.xcExcorr = e.exCorr.xcExcorr;
		#endif*/
		CommandElecExCorr::process(pl, fluidExCorr);
	}
	
	void printStatus(Everything& e, int iRep)
	{	const ExCorr& fluidExCorr = e.eVars.fluidParams.exCorr;
		switch(fluidExCorr.exCorrType)
		{
			#ifdef LIBXC_ENABLED
			case KineticLibXC: logPrintf("%s ", xcMap_K.getString(fluidExCorr.xcKinetic)); break;
			#endif
			default: logPrintf("%s ", kineticTypeMap.getString(fluidExCorr.kineticType));
		}
		CommandElecExCorr::printStatus(fluidExCorr);
	}
}
commandFluidExCorr;
