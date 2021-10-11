#include <electronic/VanDerWaalsD3.h>
#include <electronic/Everything.h>

//! DFT-D3 parameters for all elements.
// Note that we use the original DFT-D3 parameters and damping scheme here.
// Per-element arrays are initialized at the bottom of this file for clarity.
// All parameters are extracted from the dft-d3 reference code by Grimme et al.
namespace D3
{	const int Zmax = 94; //!< maximum parameterized atomic number
	const int alpha6 = 14; //!< exponent in damping function of r^-6 term
	const int alpha8 = 16; //!< exponent in damping function of r^-8 term
	const double k1 = 16.; //!< damping factor in coordination number calculation
	const double k2 = 4./3.; //!< covalent radius scale factor in coordination number calculation
	const double k3 = -4.; //!< exponential factor in coordination number interpolation of C6
	extern const double sqrtQ[Zmax]; //!< sqrt(Q) for each element, involved in computing C8 from C6
	extern const double k2Rcov[Zmax]; //!< covalent radius in bohrs for each element, scaled by k2 = 4/3

	//! Set scale parameters s6, s8 and damping parameters sr6, sr8 for r^-6, r^-8 terms by name of XC functional xcName.
	//! Note that xcName is canonicalized upon return, especially for LibXC functionals, to unify internal and LibXC names.
	void setXCscale(string& xcName, double& s6, double& sr6, double& s8, double& sr8);
};


VanDerWaalsD3::VanDerWaalsD3(const Everything& e) : VanDerWaals(e)
{
	logPrintf("\nInitializing DFT-D3 calculator:\n");
	
	//Get parameters for exchange-correlation functional
	string xcName = e.exCorr.getName();
	D3::setXCscale(xcName, s6, sr6, s8, sr8);
	logPrintf("\tParameters set for %s functional\n", xcName.c_str());
	logPrintf("\ts6: %6.3lf  s_r6: %6.3lf\n", s6, sr6);
	logPrintf("\ts8: %6.3lf  s_r8: %6.3lf\n", s8, sr8);

	Citations::add("DFT-D3 dispersion correction", "S. Grimme, J. Antony, S. Ehrlich and H. Krieg, J. Chem. Phys. 132, 154104 (2010)");
}


VanDerWaalsD3::~VanDerWaalsD3()
{
}


double VanDerWaalsD3::getScaleFactor(string exCorrName, double scaleOverride) const
{	return 0.;  //global scale factor not used in D3
}


double VanDerWaalsD3::energyAndGrad(std::vector<Atom>& atoms, const double scaleFac, matrix3<>* E_RRT) const
{	die("Not yet implemented.\n");
	return 0.;
}


//-------- D3 element-dependent parameter declarations --------

namespace D3
{
	const double sqrtQ[Zmax] = {
		2.00734898,   1.56637132,  5.01986934,  3.85379032,  3.64446594,
		3.10492822,   2.71175247,  2.59361680,  2.38825250,  2.21522516,
		6.58585536,   5.46295967,  5.65216669,  4.88284902,  4.29727576,
		4.04108902,   3.72932356,  3.44677275,  7.97762753,  7.07623947,
		6.60844053,   6.28791364,  6.07728703,  5.54643096,  5.80491167,
		5.58415602,   5.41374528,  5.28497229,  5.22592821,  5.09817141,
		6.12149689,   5.54083734,  5.06696878,  4.87005108,  4.59089647,
		4.31176304,   9.55461698,  8.67396077,  7.97210197,  7.43439917,
		6.58711862,   6.19536215,  6.01517290,  5.81623410,  5.65710424,
		5.52640661,   5.44263305,  5.58285373,  7.02081898,  6.46815523,
		5.98089120,   5.81686657,  5.53321815,  5.25477007, 11.02204549,
		10.15679528,  9.35167836,  9.06926079,  8.97241155,  8.90092807,
		8.85984840,   8.81736827,  8.79317710,  7.89969626,  8.80588454,
		8.42439218,   8.54289262,  8.47583370,  8.45090888,  8.47339339,
		7.83525634,   8.20702843,  7.70559063,  7.32755997,  7.03887381,
		6.68978720,   6.05450052,  5.88752022,  5.70661499,  5.78450695,
		7.79780729,   7.26443867,  6.78151984,  6.67883169,  6.39024318,
		6.09527958,  11.79156076, 11.10997644,  9.51377795,  8.67197068,
		8.77140725,   8.65402716,  8.53923501,  8.85024712};

	const double k2Rcov[Zmax] = {
		0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865,
		1.88972601, 1.78894056, 1.58736983, 1.61256616, 1.68815527,
		3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820,
		2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730,
		3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923,
		2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188,
		2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349,
		2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216,
		3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717,
		2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967,
		3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625,
		4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657,
		3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833,
		3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098,
		3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878,
		2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790,
		3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584,
		3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289,
		3.82984466, 3.85504098, 3.88023730, 3.90543362};

	//List of supported functionals
	enum XC { XC_LDA, XC_PBE, XC_PBESOL, XC_RPBE, XC_SSB, XC_HCTH_120,
		XC_TPSS, XC_M06_L, XC_HF, XC_PBE0, XC_PBE38, XC_HSE06,
		XC_B3PW91, XC_B3LYP, XC_CAM_B3LYP, XC_PW6B95, XC_TPSS0,
		XC_TPSSH, XC_PWB6K, XC_MPW1B95, XC_MPWB1K, XC_BMK, XC_LC_WPBE,
		XC_M05, XC_M05_2X, XC_M06, XC_M06_2X, XC_M06_HF};
		
	//Map onto shortened names:
	EnumStringMap<XC> xcMap(
		XC_LDA, "lda", XC_PBE, "gga-PBE", XC_PBESOL, "gga-PBEsol",
		XC_RPBE, "gga-RPBE", XC_SSB, "gga-SSB", XC_HCTH_120, "gga-HCTH-120",
		XC_TPSS, "mgga-TPSS", XC_M06_L, "mgga-M06-L",
		XC_HF, "hartree-fock", XC_PBE0, "hyb-PBE0", XC_PBE38, "hyb-PBE38",
		XC_HSE06, "hyb-HSE06", XC_B3PW91, "hyb-B3PW91", XC_B3LYP, "hyb-B3LYP",
		XC_CAM_B3LYP, "hyb-CAM-B3LYP", XC_PW6B95, "hyb-PW6B95",
		XC_TPSS0, "hyb-TPSS0", XC_TPSSH, "hyb-TPSSH", XC_PWB6K, "hyb-PWB6K",
		XC_MPW1B95, "hyb-MPW1B95", XC_MPWB1K, "hyb-MPW1BK", XC_BMK, "hyb-BMK",
		XC_LC_WPBE, "hyb-LC-PBE", XC_M05, "hyb-M05", XC_M05_2X, "hyb-M05-2X",
		XC_M06, "hyb-M06", XC_M06_2X, "hyb-M06-2X", XC_M06_HF, "hyb-M06-HF");

	//Replace first occurence of target in s with replacement
	inline void string_replace(string& s, string target, string replacement)
	{	size_t pos = s.find(target);
		if(pos != string::npos)
			s.replace(pos, target.length(), replacement);
	}
	
	void setXCscale(string& xcName, double& s6, double& sr6, double& s8, double& sr8)
	{	//Canonicalize name:
		if(xcName.substr(0, 3) == "lda") xcName = "lda"; //remove LDA suffixes
		#ifdef LIBXC_ENABLED
		//--- simplify libxc names to match internal format
		//--- remove separate correlation functional name after ':', if any
		size_t sepPos = xcName.find(":");
		if(sepPos != string::npos)
			xcName = xcName.substr(0, sepPos);
		//--- remove xc or x from name
		string_replace(xcName, "-xc-", "-");
		string_replace(xcName, "-x-", "-");
		//--- remove gga or mgga for hybrid functionals
		if(xcName.substr(0, 3) == "hyb")
		{	string_replace(xcName, "-gga-", "-");
			string_replace(xcName, "-mgga-", "-");
		}
		#endif
		//Find XC functional in supported list:
		XC xc;
		if(not xcMap.getEnum(xcName.c_str(), xc))
			die("\nDFT-D3 parameterization not available for %s functional.\n\n", xcName.c_str());
		//Set parameters:
		s6 = 1.;
		sr8 = 1.;
		switch(xc)
		{	case XC_LDA: { sr6 =0.999; s8 =-1.957; sr8=0.697; break; }
			//GGAs:
			case XC_PBE: { sr6=1.217; s8=0.722; break; }
			case XC_PBESOL: { sr6=1.345; s8=0.612; break; }
			case XC_RPBE: { sr6=0.872; s8=0.514; break; }
			case XC_SSB: { sr6=1.215; s8=0.663; break; }
			case XC_HCTH_120: { sr6=1.221; s8=1.206; break; }
			//mGGAs:
			case XC_TPSS: { sr6=1.166; s8=1.105; break; }
			//Hybrids:
			case XC_HF: { sr6=1.158; s8=1.746; break; }
			case XC_PBE0: { sr6=1.287; s8=0.928; break; }
			case XC_PBE38: { sr6=1.333; s8=0.998; break; }
			case XC_HSE06: { sr6=1.129; s8=0.109; break; }
			case XC_B3PW91: { sr6=1.176; s8=1.775; break; }
			case XC_B3LYP: { sr6=1.261; s8=1.703; break; }
			case XC_PW6B95: { sr6=1.532; s8=0.862; break; }
			case XC_TPSS0: { sr6=1.252; s8=1.242; break; }
			case XC_TPSSH: { sr6=1.223; s8=1.219; break; }
			case XC_PWB6K: { sr6=1.660; s8=0.550; break; }
			case XC_MPW1B95: { sr6=1.605; s8=1.118; break; }
			case XC_MPWB1K: { sr6=1.671; s8=1.061; break; }
			case XC_BMK: { sr6=1.931; s8=2.168; break; }
			case XC_CAM_B3LYP: { sr6=1.378; s8=1.217; break; }
			case XC_LC_WPBE: { sr6=1.355; s8=1.279; break; }
			case XC_M05: { sr6=1.373; s8=0.595; break; }
			case XC_M05_2X: { sr6=1.417; s8=0.000; break; }
			case XC_M06_L: { sr6=1.581; s8=0.000; break; }
			case XC_M06: { sr6=1.325; s8=0.000; break; }
			case XC_M06_2X: { sr6=1.619; s8=0.000; break; }
			case XC_M06_HF: { sr6=1.446; s8=0.000; break; }
		}
	}
}
