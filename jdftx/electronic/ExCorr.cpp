/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

#include <electronic/ExCorr.h>
#include <electronic/Basis.h>
#include <electronic/operators.h>
#include <electronic/Everything.h>
#include <electronic/ExactExchange.h>
#include <electronic/ExCorr_internal.h>
#include <electronic/ExCorr_internal_LDA.h>
#include <electronic/ExCorr_internal_GGA.h>
#include <electronic/ExCorr_internal_mGGA.h>
#include <core/Thread.h>
#include <core/GpuUtil.h>
#include <core/DataMultiplet.h>

//---------------- Subset wrapper for MPI parallelization --------------------

void Functional::evaluateSub(int iStart, int iStop,
	std::vector<const double*> n, std::vector<const double*> sigma,
	std::vector<const double*> lap, std::vector<const double*> tau,
	double* E, std::vector<double*> E_n, std::vector<double*> E_sigma,
	std::vector<double*> E_lap, std::vector<double*> E_tau) const
{
	struct Offset
	{	int iStart;
		std::vector<const double*> operator()(std::vector<const double*> v) { auto vOut=v; for(auto& p: vOut) if(p) p+=iStart; return vOut; }
		std::vector<double*> operator()(std::vector<double*> v) { auto vOut=v; for(auto& p: vOut) if(p) p+=iStart; return vOut; }
	}
	offset;
	offset.iStart = iStart;
	evaluate(iStop-iStart, offset(n), offset(sigma), offset(lap), offset(tau), E ? E+iStart : 0, offset(E_n), offset(E_sigma), offset(E_lap), offset(E_tau));
}

//---------------- LDA thread launcher / gpu switch --------------------

FunctionalLDA::FunctionalLDA(LDA_Variant variant, double scaleFac) : Functional(scaleFac), variant(variant)
{	switch(variant)
	{	case LDA_X_Slater:  logPrintf("Initalized Slater LDA exchange.\n"); break;
		case LDA_C_PZ:      logPrintf("Initalized Perdew-Zunger LDA correlation.\n"); break;
		case LDA_C_PW:      logPrintf("Initalized Perdew-Wang LDA correlation.\n"); break;
		case LDA_C_PW_prec: logPrintf("Initalized Perdew-Wang LDA correlation (extended precision).\n"); break;
		case LDA_C_VWN:     logPrintf("Initalized Vosko-Wilk-Nusair LDA correlation.\n"); break;
		case LDA_XC_Teter:  logPrintf("Initalized Teter93 LSD exchange+correlation.\n"); break;
		case LDA_KE_TF:     logPrintf("Initalized Thomas-Fermi LDA kinetic energy.\n"); break;
	}
}

template<LDA_Variant variant, int nCount>
void LDA(int N, array<const double*,nCount> n, double* E, array<double*,nCount> E_n, double scaleFac)
{	threadedLoop(LDA_calc<variant,nCount>::compute, N, n, E, E_n, scaleFac);
}
void LDA(LDA_Variant variant, int N, std::vector<const double*> n, double* E, std::vector<double*> E_n, double scaleFac)
{	SwitchTemplate_spin(SwitchTemplate_LDA, variant, n.size(), LDA, (N, n, E, E_n, scaleFac) )
}
#ifdef GPU_ENABLED
void LDA_gpu(LDA_Variant variant, int N, std::vector<const double*> n, double* E, std::vector<double*> E_n, double scaleFac);
#endif

void FunctionalLDA::evaluate(int N, std::vector<const double*> n, std::vector<const double*> sigma,
	std::vector<const double*> lap, std::vector<const double*> tau,
	double* E, std::vector<double*> E_n, std::vector<double*> E_sigma,
	std::vector<double*> E_lap, std::vector<double*> E_tau) const
{	//LDA: so ignore sigma, lap, tau and call the above functions (CPU/GPU as appropriate)
	assert(n.size()==1 || n.size()==2);
	callPref(LDA)(variant, N, n, E, E_n, scaleFac);
}

//---------------- GGA thread launcher / gpu switch --------------------

FunctionalGGA::FunctionalGGA(GGA_Variant variant, double scaleFac) : Functional(scaleFac), variant(variant)
{	switch(variant)
	{	case GGA_X_PBE: logPrintf("Initalized PBE GGA exchange.\n"); break;
		case GGA_C_PBE: logPrintf("Initalized PBE GGA correlation.\n"); break;
		case GGA_X_PBEsol: logPrintf("Initalized PBEsol GGA exchange.\n"); break;
		case GGA_C_PBEsol: logPrintf("Initalized PBEsol GGA correlation.\n"); break;
		case GGA_X_PW91: logPrintf("Initalized PW91 GGA exchange.\n"); break;
		case GGA_C_PW91: logPrintf("Initalized PW91 GGA correlation.\n"); break;
		case GGA_X_wPBE_SR: logPrintf("Initalized omega-PBE short-ranged GGA exchange.\n"); break;
		case GGA_KE_VW: logPrintf("Initialized von Weisacker kinetic energy gradient correction.\n"); break; 
		case GGA_KE_PW91: logPrintf("Initialized PW91 GGA kinetic energy.\n"); break; 
	}
}

template<GGA_Variant variant, bool spinScaling, int nCount>
void GGA(int N, array<const double*,nCount> n, array<const double*,2*nCount-1> sigma,
	double* E, array<double*,nCount> E_n, array<double*,2*nCount-1> E_sigma, double scaleFac)
{	threadedLoop(GGA_calc<variant,spinScaling,nCount>::compute, N, n, sigma, E, E_n, E_sigma, scaleFac);
}
void GGA(GGA_Variant variant, int N, std::vector<const double*> n, std::vector<const double*> sigma,
	double* E, std::vector<double*> E_n, std::vector<double*> E_sigma, double scaleFac)
{	SwitchTemplate_spin(SwitchTemplate_GGA, variant, n.size(), GGA, (N, n, sigma, E, E_n, E_sigma, scaleFac) )
}
#ifdef GPU_ENABLED
void GGA_gpu(GGA_Variant variant, int N, std::vector<const double*> n, std::vector<const double*> sigma,
	double* E, std::vector<double*> E_n, std::vector<double*> E_sigma, double scaleFac);
#endif

void FunctionalGGA::evaluate(int N, std::vector<const double*> n, std::vector<const double*> sigma,
	std::vector<const double*> lap, std::vector<const double*> tau,
	double* E, std::vector<double*> E_n, std::vector<double*> E_sigma,
	std::vector<double*> E_lap, std::vector<double*> E_tau) const
{	//GGA: so ignore lap, tau and call the above functions (CPU/GPU as appropriate)
	assert(n.size()==1 || n.size()==2);
	callPref(GGA)(variant, N, n, sigma, E, E_n, E_sigma, scaleFac);
}


//---------------- metaGGA thread launcher / gpu switch --------------------

FunctionalMGGA::FunctionalMGGA(mGGA_Variant variant, double scaleFac) : Functional(scaleFac), variant(variant)
{	switch(variant)
	{	case mGGA_X_TPSS: logPrintf("Initalized TPSS mGGA exchange.\n"); break;
		case mGGA_C_TPSS: logPrintf("Initalized TPSS mGGA correlation.\n"); break;
		case mGGA_X_revTPSS: logPrintf("Initalized revTPSS mGGA exchange.\n"); break;
		case mGGA_C_revTPSS: logPrintf("Initalized revTPSS mGGA correlation.\n"); break;
	}
}

template<mGGA_Variant variant, bool spinScaling, int nCount>
void mGGA(int N, array<const double*,nCount> n, array<const double*,2*nCount-1> sigma,
	array<const double*,nCount> lap, array<const double*,nCount> tau,
	double* E, array<double*,nCount> E_n, array<double*,2*nCount-1> E_sigma,
	array<double*,nCount> E_lap, array<double*,nCount> E_tau, double scaleFac)
{	threadedLoop(mGGA_calc<variant,spinScaling,nCount>::compute, N,
		n, sigma, lap, tau, E, E_n, E_sigma, E_lap, E_tau, scaleFac);
}
void mGGA(mGGA_Variant variant, int N, std::vector<const double*> n, std::vector<const double*> sigma,
	std::vector<const double*> lap, std::vector<const double*> tau,
	double* E, std::vector<double*> E_n, std::vector<double*> E_sigma,
	std::vector<double*> E_lap, std::vector<double*> E_tau, double scaleFac)
{	SwitchTemplate_spin(SwitchTemplate_mGGA, variant, n.size(), mGGA, (N,
		n, sigma, lap, tau, E, E_n, E_sigma, E_lap, E_tau, scaleFac) )
}
#ifdef GPU_ENABLED
void mGGA_gpu(mGGA_Variant variant, int N, std::vector<const double*> n, std::vector<const double*> sigma,
	std::vector<const double*> lap, std::vector<const double*> tau,
	double* E, std::vector<double*> E_n, std::vector<double*> E_sigma,
	std::vector<double*> E_lap, std::vector<double*> E_tau, double scaleFac);
#endif

void FunctionalMGGA::evaluate(int N, std::vector<const double*> n, std::vector<const double*> sigma,
	std::vector<const double*> lap, std::vector<const double*> tau,
	double* E, std::vector<double*> E_n, std::vector<double*> E_sigma,
	std::vector<double*> E_lap, std::vector<double*> E_tau) const
{	//mGGA: so ignore lap, tau and call the above functions (CPU/GPU as appropriate)
	assert(n.size()==1 || n.size()==2);
	callPref(mGGA)(variant, N, n, sigma, lap, tau, E, E_n, E_sigma, E_lap, E_tau, scaleFac);
}



//---------------------------- LibXC wrapper -----------------------------
#ifdef LIBXC_ENABLED
#include <xc.h>

//! LibXC wrapper that provides an interface similar to that of the internal functionals
//! with minor differences to handle the different data order and conventions used by LibXC
class FunctionalLibXC
{	xc_func_type funcUnpolarized, funcPolarized;
public:
	bool needsSigma() const { return funcUnpolarized.info->family != XC_FAMILY_LDA; } //!< gradients needed for GGA and HYB_GGA
	bool needsLap() const { return funcUnpolarized.info->family == XC_FAMILY_MGGA; } //!< MGGAs may need laplacian
	bool needsTau() const { return funcUnpolarized.info->family == XC_FAMILY_MGGA; } //!< MGGAs need KE density
	bool isKinetic() const { return funcUnpolarized.info->kind == XC_KINETIC; }
	double exxScale() const { return funcUnpolarized.cam_alpha; }
	double exxOmega() const { return funcUnpolarized.cam_omega; }
	
	FunctionalLibXC(int xcCode, const char* typeName)
	{	if(xc_func_init(&funcUnpolarized, xcCode, XC_UNPOLARIZED) != 0)
			die("Error initializing LibXC unpolarized %s functional\n", typeName);
		if(xc_func_init(&funcPolarized, xcCode, XC_POLARIZED) != 0)
			die("Error initializing LibXC polarized %s functional\n", typeName);
		
		logPrintf("Initialized LibXC %s functional '%s'\n", typeName, funcUnpolarized.info->name);
		
		Citations::add("LibXC library of exchange-correlation functions",
			"M. A. L. Marques, M. J. T. Oliveira and T. Burnus, Comput. Phys. Commun. 183, 2272 (2012)");
		Citations::add(
			funcUnpolarized.info->name + string(" ") + typeName + string(" functional"),
			funcUnpolarized.info->refs);
	}
	~FunctionalLibXC()
	{	xc_func_end(&funcUnpolarized);
		xc_func_end(&funcPolarized);
	}
	
	//! Like Functional::evaluate, except different spin components are stored together
	//! and the computed energy is per-particle (e) instead of per volume (E).
	void evaluate(int nCount, int N,
		const double* n, const double* sigma, const double* lap, const double* tau,
		double* e, double* E_n, double* E_sigma, double* E_lap, double* E_tau) const
	{
		assert(nCount==1 || nCount==2);
		const xc_func_type& func = (nCount==1) ? funcUnpolarized : funcPolarized;
		int sigmaCount = 2*nCount-1; //1 for unpolarized, 3 for polarized
		int Nn = N * nCount;
		int Nsigma = N * sigmaCount;
		//Alloctae temporaries:
		std::vector<double> eTemp(N);
		std::vector<double> E_nTemp(E_n ? Nn : 0);
		std::vector<double> E_sigmaTemp(E_n && needsSigma() ? Nsigma : 0);
		std::vector<double> E_lapTemp(E_n && needsLap() ? Nn : 0);
		std::vector<double> E_tauTemp(E_n && needsTau() ? Nn : 0);
		//Invoke appropriate LibXC function in scratch space:
		if(needsTau())
		{	//HACK: LibXC (v1.0) metaGGAs mess up if gradients are not computed, so always compute them:
			E_nTemp.resize(Nn); E_sigmaTemp.resize(Nsigma); E_lapTemp.resize(Nn); E_tauTemp.resize(Nn);
			xc_mgga_exc_vxc(&func, N, n, sigma, lap, tau,
				&eTemp[0], &E_nTemp[0], &E_sigmaTemp[0], &E_lapTemp[0], &E_tauTemp[0]);
			if(!E_n) { E_nTemp.clear(); E_sigmaTemp.clear(); E_lapTemp.clear(); E_tauTemp.clear(); }
		}
		else if(needsSigma())
		{	if(E_n) xc_gga_exc_vxc(&func, N, n, sigma, &eTemp[0], &E_nTemp[0], &E_sigmaTemp[0]); //need gradient
			else xc_gga_exc(&func, N, n, sigma, &eTemp[0]);
		}
		else
		{	if(E_n) xc_lda_exc_vxc(&func, N, n, &eTemp[0], &E_nTemp[0]); //need gradient
			else xc_lda_exc(&func, N, n, &eTemp[0]);
		}
		//Accumulate onto final results
		eblas_daxpy(N, 1., &eTemp[0], 1, e, 1);
		if(E_nTemp.size()) eblas_daxpy(Nn, 1., &E_nTemp[0], 1, E_n, 1);
		if(E_sigmaTemp.size()) eblas_daxpy(Nsigma, 1., &E_sigmaTemp[0], 1, E_sigma, 1);
		if(E_lapTemp.size()) eblas_daxpy(Nn, 1., &E_lapTemp[0], 1, E_lap, 1);
		if(E_tauTemp.size()) eblas_daxpy(Nn, 1., &E_tauTemp[0], 1, E_tau, 1);
	}
};

//! Convert a collection of scalar fields into an interleaved vector field.
//! result can be freed using delete[]
template<unsigned M> double* transpose(const DataRptrCollection& inVec)
{	assert(inVec.size()==M);
	const unsigned N = inVec[0]->nElem;
	const double* in[M]; for(unsigned m=0; m<M; m++) in[m] = inVec[m]->data();
	double *out = new double[M*N], *outPtr = out;
	for(unsigned n=0; n<N; n++)
		for(unsigned m=0; m<M; m++)
			*(outPtr++) = in[m][n];
	return out;
}

//! Convert an interleaved vector field to a collection of scalar fields
template<unsigned M> void transpose(double* in, DataRptrCollection& outVec)
{	assert(outVec.size()==M);
	const unsigned N = outVec[0]->nElem;
	double* out[M]; for(unsigned m=0; m<M; m++) out[m] = outVec[m]->data();
	double *inPtr = in;
	for(unsigned n=0; n<N; n++)
		for(unsigned m=0; m<M; m++)
			out[m][n] = *(inPtr++);
}

#endif //LIBXC_ENABLED

//----------------------- Functional List --------------------

struct FunctionalList
{	std::vector<std::shared_ptr<Functional> > internal; //!<Functionals with an internal implementation
	void add(LDA_Variant variant, double scaleFac=1.0)
	{	internal.push_back(std::make_shared<FunctionalLDA>(variant, scaleFac));
	}
	void add(GGA_Variant variant, double scaleFac=1.0)
	{	internal.push_back(std::make_shared<FunctionalGGA>(variant, scaleFac));
	}
	void add(mGGA_Variant variant, double scaleFac=1.0)
	{	internal.push_back(std::make_shared<FunctionalMGGA>(variant, scaleFac));
	}
	
	#ifdef LIBXC_ENABLED
	std::vector<std::shared_ptr<FunctionalLibXC> > libXC; //!<Functionals which use LibXC for evaluation
	void add(int xcCode, const char* name)
	{	libXC.push_back(std::make_shared<FunctionalLibXC>(xcCode, name));
	}
	#endif
};


//-------------- ExCorr members ----------------

ExCorr::ExCorr() : exCorrType(ExCorrGGA_PBE), kineticType(KineticNone), xcName("gga-PBE"),
exxScale(0.), exxOmega(0.),
functionals(std::make_shared<FunctionalList>())
#ifdef LIBXC_ENABLED
, xcExchange(0), xcCorr(0), xcExcorr(0), xcKinetic(0)
#endif
{
}

void ExCorr::setup(const Everything& everything)
{	e = &everything;
	
	string citeReason = xcName + " exchange-correlation functional";

	switch(exCorrType)
	{
		#ifdef LIBXC_ENABLED
		case ExCorrLibXC:
		{	if(xcExcorr)
			{	functionals->add(xcExcorr, "exchange-correlation");
				exxScale = functionals->libXC.back()->exxScale(); //set exact exchange factor
				exxOmega = functionals->libXC.back()->exxOmega(); //set exact exchange factor
			}
			else
			{	functionals->add(xcExchange, "exchange");
				functionals->add(xcCorr, "correlation");
				//exact exchange factor 0 by default (all hybrids are listed as combined XC functionals)
			}
			break;
		}
		#endif //LIBXC_ENABLED
		
		case ExCorrLDA_PZ:
			functionals->add(LDA_X_Slater);
			functionals->add(LDA_C_PZ);
			Citations::add(citeReason, "J.P. Perdew and A. Zunger, Phys. Rev. B 23, 5048 (1981)");
			break;
		case ExCorrLDA_PW:
			functionals->add(LDA_X_Slater);
			functionals->add(LDA_C_PW);
			Citations::add(citeReason, "J.P. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992)");
			break;
		case ExCorrLDA_PW_prec:
			functionals->add(LDA_X_Slater);
			functionals->add(LDA_C_PW_prec);
			Citations::add(citeReason, "J.P. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992)");
			break;
		case ExCorrLDA_VWN:
			functionals->add(LDA_X_Slater);
			functionals->add(LDA_C_VWN);
			Citations::add(citeReason, "S.H. Vosko, L. Wilk and M. Nusair, Can. J. Phys. 58, 1200 (1980)");
			break;
		case ExCorrLDA_Teter:
			functionals->add(LDA_XC_Teter);
			Citations::add(citeReason, "S. Goedecker, M. Teter and J. Hutter, Phys. Rev. B 54, 1703 (1996)");
			break;
		case ExCorrGGA_PBE:
			functionals->add(GGA_X_PBE);
			functionals->add(GGA_C_PBE);
			Citations::add(citeReason, "J.P. Perdew, K. Burke and M. Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)");
			break;
		case ExCorrGGA_PBEsol:
			functionals->add(GGA_X_PBEsol);
			functionals->add(GGA_C_PBEsol);
			Citations::add(citeReason, "J.P. Perdew et al., Phys. Rev. Lett. 100, 136406 (2008)");
			break;
		case ExCorrGGA_PW91:
			functionals->add(GGA_X_PW91);
			functionals->add(GGA_C_PW91);
			Citations::add(citeReason, "J.P. Perdew et al., Phys. Rev. B 46, 6671 (1992)");
			break;
		case ExCorrMGGA_TPSS:
			functionals->add(mGGA_X_TPSS);
			functionals->add(mGGA_C_TPSS);
			Citations::add(citeReason, "J. Tao, J.P. Perdew, V.N. Staroverov and G. Scuseria, Phys. Rev. Lett. 91, 146401 (2003)");
			break;
		case ExCorrMGGA_revTPSS:
			functionals->add(mGGA_X_revTPSS);
			functionals->add(mGGA_C_revTPSS);
			Citations::add(citeReason, "J.P. Perdew et al., Phys. Rev. Lett. 103, 026403 (2009)");
			break;
		case ExCorrHYB_PBE0:
			exxScale = 1./4;
			functionals->add(GGA_X_PBE, 3./4);
			functionals->add(GGA_C_PBE);
			Citations::add(citeReason, "M. Ernzerhof and G. E. Scuseria, J. Chem. Phys. 110, 5029 (1999)");
			break;
		case ExCorrHYB_HSE06:
			exxOmega = 0.11;
			exxScale = 1./4;
			functionals->add(GGA_X_wPBE_SR, -1./4);
			functionals->add(GGA_X_PBE);
			functionals->add(GGA_C_PBE);
			Citations::add(citeReason, "A.V. Krukau, O.A. Vydrov, A.F. Izmaylov and G.E. Scuseria, J. Chem. Phys. 125, 224106 (2006)");
			break;
		case ExCorrHYB_HSE12:
			exxOmega = 0.185;
			exxScale = 0.313;
			functionals->add(GGA_X_wPBE_SR, -exxScale);
			functionals->add(GGA_X_PBE);
			functionals->add(GGA_C_PBE);
			Citations::add(citeReason, "J.E. Moussa, P.A. Schultz and J.R. Chelikowsky, J. Chem. Phys. 136, 204117 (2012)");
			break;
		case ExCorrHYB_HSE12s:
			exxOmega = 0.408;
			exxScale = 0.425;
			functionals->add(GGA_X_wPBE_SR, -exxScale);
			functionals->add(GGA_X_PBE);
			functionals->add(GGA_C_PBE);
			Citations::add(citeReason, "J.E. Moussa, P.A. Schultz and J.R. Chelikowsky, J. Chem. Phys. 136, 204117 (2012)");
			break;
		case ExCorrHF:
			exxScale = 1.;
			break;
	}
	
	if(exxScale)
	{	logPrintf("Will include %lg x ", exxScale);
		if(!exxOmega) logPrintf("exact exchange.\n");
		else logPrintf("screened exact exchange with range-parameter %lg\n", exxOmega);
	}
	
	switch(kineticType)
	{	case KineticNone:
			break;
		case KineticTF:
			functionals->add(LDA_KE_TF);
			Citations::add("Thomas-Fermi kinetic energy functional",
				"L.H. Thomas, Proc. Cambridge Phil. Soc. 23, 542 (1927)\n"
				"E. Fermi, Rend. Accad. Naz. Lincei 6, 602 (1927)");
			break;
		case KineticVW:
			functionals->add(LDA_KE_TF);
			functionals->add(GGA_KE_VW,1.);
			Citations::add("Thomas-Fermi-von-Weisacker kinetic energy functional",
				"C.F.v. Weizsacker, Z. Phys. 96, 431 (1935)");
			break;
		case KineticPW91:
			functionals->add(GGA_KE_PW91);
			Citations::add("PW91K kinetic energy functional",
				"A. Lembarki and H. Chermette, Phys. Rev. A 50, 5328 (1994)");
			break;
		#ifdef LIBXC_ENABLED
		case KineticLibXC:
			functionals->add(xcKinetic, "kinetic");
			break;
		#endif //LIBXC_ENABLED
	}
}


string ExCorr::getName() const
{	return xcName;
}

double ExCorr::exxFactor() const
{	return exxScale;
}

double ExCorr::exxRange() const
{	return exxOmega;
}

bool ExCorr::needsKEdensity() const
{	for(auto func: functionals->internal)
		if(func->needsTau())
			return true;
	#ifdef LIBXC_ENABLED
	for(auto func: functionals->libXC)
		if(func->needsTau())
			return true;
	#endif
	return false;
}

double ExCorr::operator()(const DataRptrCollection& n, DataRptrCollection* Vxc, bool includeKinetic,
		const DataRptrCollection* tauPtr, DataRptrCollection* Vtau) const
{
	static StopWatch watch("ExCorrTotal"), watchComm("ExCorrCommunication"), watchFunc("ExCorrFunctional");
	watch.start();
	
	const int nCount = n.size();
	assert(nCount==1 || nCount==2);
	const int sigmaCount = 2*nCount-1;
	const GridInfo& gInfo = e->gInfo;
	
	//------- Prepare inputs, allocate outputs -------
	
	//Energy density per volume:
	DataRptr E; nullToZero(E, gInfo);
	
	//Gradient w.r.t spin densities:
	DataRptrCollection E_n(nCount);
	if(Vxc)
	{	Vxc->clear();
		nullToZero(E_n, gInfo, nCount);
	}
	
	//Check for GGAs and meta GGAs:
	bool needsSigma = false, needsTau=false, needsLap=false;
	for(auto func: functionals->internal)
		if(includeKinetic || !func->isKinetic())
		{	needsSigma |= func->needsSigma();
			needsLap |= func->needsLap();
			needsTau |= func->needsTau();
		}
	#ifdef LIBXC_ENABLED
	for(auto func: functionals->libXC)
		if(includeKinetic || !func->isKinetic())
		{	needsSigma |= func->needsSigma();
			needsLap |= func->needsLap();
			needsTau |= func->needsTau();
		}
	#endif
	
	//Additional inputs/outputs for GGAs (gradient contractions and gradients w.r.t those)
	DataRptrCollection sigma(sigmaCount), E_sigma(sigmaCount);
	std::vector<DataRptrVec> Dn(nCount);
	int iDirStart = (3*mpiUtil->iProcess())/mpiUtil->nProcesses();
	int iDirStop = (3*(mpiUtil->iProcess()+1))/mpiUtil->nProcesses();
	if(needsSigma)
	{	//Compute the gradients of the (spin-)densities:
		for(int s=0; s<nCount; s++)
		{	const DataGptr Jn = J(n[s]);
			for(int i=iDirStart; i<iDirStop; i++)
				Dn[s][i] = I(D(Jn,i),true);
		}
		//Compute the required contractions:
		for(int s1=0; s1<nCount; s1++)
			for(int s2=s1; s2<nCount; s2++)
			{	for(int i=iDirStart; i<iDirStop; i++)
					sigma[s1+s2] += Dn[s1][i] * Dn[s2][i];
				watchComm.start();
				nullToZero(sigma[s1+s2], gInfo);
				sigma[s1+s2]->allReduce(MPIUtil::ReduceSum);
				watchComm.stop();
			}
		//Allocate gradient if required:
		if(Vxc) nullToZero(E_sigma, gInfo, sigmaCount);
	}
	
	//Additional inputs/outputs for MGGAs (Laplacian, orbital KE and gradients w.r.t those)
	DataRptrCollection lap(nCount), E_lap(nCount);
	if(needsLap)
	{	//Compute laplacian
		for(int s=0; s<nCount; s++)
			lap[s] = (1./gInfo.detR) * I(L(J(n[s])));
		//Allocate gradient w.r.t laplacian if required
		if(Vxc) nullToZero(E_lap, gInfo, nCount);
	}
	DataRptrCollection tau(nCount), E_tau(nCount);
	if(needsTau)
	{	//make sure orbital KE density has been provided 
		assert(tauPtr);
		tau = *tauPtr;
		//allocate gradients w.r.t KE density if required
		if(Vxc)
		{	assert(Vtau); //if computing gradients, all gradients must be computed
			Vtau->clear();
			nullToZero(E_tau, gInfo, nCount);
		}
	}
	
	//Cap negative densities to 0:
	DataRptrCollection nCapped = clone(n);
	double nMin, nMax;
	for(int s=0; s<nCount; s++)
		callPref(eblas_capMinMax)(gInfo.nr, nCapped[s]->dataPref(), nMin, nMax, 0.);
	
	#ifdef LIBXC_ENABLED
	//------------------ Evaluate LibXC functionals ---------------
	if(functionals->libXC.size())
	{	//Prepare input/output data on the CPU in transposed order (spins contiguous)
		double *eData = E->data(), *nData=0, *sigmaData=0, *lapData=0, *tauData=0;
		double *E_nData=0, *E_sigmaData=0, *E_lapData=0, *E_tauData=0;
		if(nCount == 1)
		{	nData = nCapped[0]->data();
			if(needsSigma) sigmaData = sigma[0]->data();
			if(needsLap) lapData = lap[0]->data();
			if(needsTau) tauData = tau[0]->data();
			if(Vxc)
			{	E_nData = E_n[0]->data();
				if(needsSigma) E_sigmaData = E_sigma[0]->data();
				if(needsLap) E_lapData = E_lap[0]->data();
				if(needsTau) E_tauData = E_tau[0]->data();
			}
		}
		else //Need to interleave input spin-vector fields:
		{	nData = transpose<2>(nCapped);
			if(needsSigma) sigmaData = transpose<3>(sigma);
			if(needsLap) lapData = transpose<2>(lap);
			if(needsTau) tauData = transpose<2>(tau);
			if(Vxc)
			{	E_nData = new double[2*gInfo.nr]; eblas_zero(2*gInfo.nr, E_nData);
				if(needsSigma) { E_sigmaData = new double[3*gInfo.nr]; eblas_zero(3*gInfo.nr, E_sigmaData); }
				if(needsLap) { E_lapData = new double[2*gInfo.nr]; eblas_zero(2*gInfo.nr, E_lapData); }
				if(needsTau) { E_tauData = new double[2*gInfo.nr]; eblas_zero(2*gInfo.nr, E_tauData); }
			}
		}
		
		//Calculate all the required functionals:
		watchFunc.start();
		for(auto func: functionals->libXC)
			if(includeKinetic || !func->isKinetic())
				func->evaluate(nCount, gInfo.nr, nData, sigmaData, lapData, tauData,
					eData, E_nData, E_sigmaData, E_lapData, E_tauData);
		watchFunc.stop();
		
		//Uninterleave spin-vector field results:
		if(nCount != 1)
		{	delete[] nData;
			if(needsSigma) delete[] sigmaData;
			if(needsLap) delete[] lapData;
			if(needsTau) delete[] tauData;
			if(Vxc)
			{	transpose<2>(E_nData, E_n); delete[] E_nData;
				if(needsSigma) { transpose<3>(E_sigmaData, E_sigma); delete[] E_sigmaData; }
				if(needsLap) { transpose<2>(E_lapData, E_lap); delete[] E_lapData; }
				if(needsTau) { transpose<2>(E_tauData, E_tau); delete[] E_tauData; }
			}
		}
	
		//Convert per-particle energy to energy density per volume
		E = E * (nCount==1 ? nCapped[0] : nCapped[0]+nCapped[1]);
	}
	#endif //LIBXC_ENABLED
	
	//---------------- Compute internal functionals ----------------
	watchFunc.start();
	for(auto func: functionals->internal)
		if(includeKinetic || !func->isKinetic())
			func->evaluateSub(gInfo.irStart, gInfo.irStop,
				constDataPref(nCapped), constDataPref(sigma), constDataPref(lap), constDataPref(tau),
				E->dataPref(), dataPref(E_n), dataPref(E_sigma), dataPref(E_lap), dataPref(E_tau));
	watchFunc.stop();
	
	//---------------- Collect resb results over processes ----------------
	watchComm.start();
	E->allReduce(MPIUtil::ReduceSum);
	for(DataRptr& x: E_n) if(x) x->allReduce(MPIUtil::ReduceSum);
	for(DataRptr& x: E_sigma) if(x) x->allReduce(MPIUtil::ReduceSum);
	for(DataRptr& x: E_lap) if(x) x->allReduce(MPIUtil::ReduceSum);
	for(DataRptr& x: E_tau) if(x) x->allReduce(MPIUtil::ReduceSum);
	watchComm.stop();

	//--------------- Gradient propagation ---------------------
	
	//Propagate gradient w.r.t sigma and lap to grad_n (for GGA/mGGAs)
	if(Vxc && (needsSigma || needsLap))
	{	DataGptr E_nTilde[2]; //contribution to the potential in fourier space
		if(needsSigma)
		{	for(int s1=0; s1<nCount; s1++)
				for(int s2=s1; s2<nCount; s2++)
					for(int i=iDirStart; i<iDirStop; i++)
					{	if(s1==s2) E_nTilde[s1] -= D(Idag(2*(E_sigma[2*s1] * Dn[s1][i])), i);
						else
						{	E_nTilde[s1] -= D(Idag(E_sigma[s1+s2] * Dn[s2][i]), i);
							E_nTilde[s2] -= D(Idag(E_sigma[s1+s2] * Dn[s1][i]), i);
						}
					}
		}
		if(needsLap && mpiUtil->isHead()) //perform only on one process, so that we can allReduce below
		{	for(int s=0; s<nCount; s++)
				E_nTilde[s] += (1./gInfo.detR) * L(Idag(E_lap[s]));
		}
		for(int s=0; s<nCount; s++)
		{	watchComm.start();
			nullToZero(E_nTilde[s], gInfo);
			E_nTilde[s]->allReduce(MPIUtil::ReduceSum);
			watchComm.stop();
			E_n[s] += Jdag(E_nTilde[s],true);
		}
	}
	
	if(Vxc) *Vxc = E_n;
	if(Vtau) *Vtau = E_tau;
	watch.stop();
	return integral(E);
}

//Unpolarized wrapper to above function:
double ExCorr::operator()(const DataRptr& n, DataRptr* Vxc, bool includeKinetic,
		const DataRptr* tau, DataRptr* Vtau) const
{	DataRptrCollection VxcArr(1), tauArr(1), VtauArr(1);
	if(tau) tauArr[0] = *tau;
	double Exc =  (*this)(DataRptrCollection(1, n), Vxc ? &VxcArr : 0, includeKinetic,
		tau ? &tauArr :0, Vtau ? &VtauArr : 0);
	if(Vxc) *Vxc = VxcArr[0];
	if(Vtau) *Vtau = VtauArr[0];
	return Exc;
}

inline void setMask(size_t iStart, size_t iStop, const double* n, double* mask, double nCut)
{	for(size_t i=iStart; i<iStop; i++) mask[i] = (n[i]<nCut ? 0. : 1.);
}

void ExCorr::getSecondDerivatives(const DataRptr& n, DataRptr& e_nn, DataRptr& e_sigma, DataRptr& e_nsigma, DataRptr& e_sigmasigma, double nCut) const
{
	//Check for GGAs and meta GGAs:
	bool needsSigma = false, needsTau=false, needsLap=false;
	for(auto func: functionals->internal)
		if(!func->isKinetic())
		{	needsSigma |= func->needsSigma();
			needsLap |= func->needsLap();
			needsTau |= func->needsTau();
		}
	#ifdef LIBXC_ENABLED
	for(auto func: functionals->libXC)
		if(!func->isKinetic())
		{	needsSigma |= func->needsSigma();
			needsLap |= func->needsLap();
			needsTau |= func->needsTau();
		}
	#endif
	
	if(needsTau || needsLap) die("Second derivatives implemented only for LDAs and GGAs.\n");
	
	const double eps = 1e-7; //Order sqrt(double-precision epsilon)
	const double scalePlus = 1.+eps;
	const double scaleMinus = 1.-eps;
	const DataRptr nPlus = scalePlus * n;
	const DataRptr nMinus = scaleMinus * n;
	
	//Compute mask to zero out low density regions
	DataRptr mask(DataR::alloc(e->gInfo));
	threadLaunch(setMask, e->gInfo.nr, n->data(), mask->data(), nCut);
	
	//Compute gradient-squared for GGA
	DataRptr sigma, sigmaPlus, sigmaMinus;
	if(needsSigma)
	{	sigma = lengthSquared(gradient(n));
		sigmaPlus = scalePlus * sigma;
		sigmaMinus = scaleMinus * sigma;
	}
	
	//Configurations of n and sigma, and the gradients w.r.t them:
	struct Config
	{	const DataRptr *n, *sigma;
		DataRptr e_n, e_sigma;
	};
	std::vector<Config> configs(5);
	configs[0].n = &n;      configs[0].sigma = &sigma;      //original point
	configs[1].n = &nPlus;  configs[1].sigma = &sigma;      // + dn
	configs[2].n = &nMinus; configs[2].sigma = &sigma;      // - dn
	configs[3].n = &n;      configs[3].sigma = &sigmaPlus;  // + dsigma
	configs[4].n = &n;      configs[4].sigma = &sigmaMinus; // - dsigma
	
	//Compute the gradients at all the configurations:
	DataRptr eTmp; nullToZero(eTmp, e->gInfo); //temporary energy return value (ignored)
	for(int i=0; i<(needsSigma ? 5 : 3); i++)
	{	Config& c = configs[i];
		std::vector<const double*> nData(1), sigmaData(1), lapData(1), tauData(1);
		std::vector<double*> e_nData(1), e_sigmaData(1), e_lapData(1), e_tauData(1);
		double* eData = eTmp->dataPref();
		nData[0] = (*c.n)->dataPref();
		nullToZero(c.e_n, e->gInfo);
		e_nData[0] = c.e_n->dataPref();
		if(needsSigma)
		{	sigmaData[0] = (*c.sigma)->dataPref();
			nullToZero(c.e_sigma, e->gInfo);
			e_sigmaData[0] = c.e_sigma->dataPref();
		}
		#ifdef LIBXC_ENABLED
		//Compute LibXC functionals:
		for(auto func: functionals->libXC)
			if(!func->isKinetic())
				func->evaluate(1, e->gInfo.nr, nData[0], sigmaData[0], lapData[0], tauData[0],
					eData, e_nData[0], e_sigmaData[0], e_lapData[0], e_tauData[0]);
		#endif
		//Compute internal functionals:
		for(auto func: functionals->internal)
			if(!func->isKinetic())
				func->evaluate(e->gInfo.nr, nData, sigmaData, lapData, tauData,
					eData, e_nData, e_sigmaData, e_lapData, e_tauData);
	}
	
	//Compute finite difference derivatives:
	DataRptr nDen = (0.5/eps) * inv(n) * mask;
	e_nn = nDen * (configs[1].e_n - configs[2].e_n);
	if(needsSigma)
	{	DataRptr sigmaDen = (0.5/eps) * inv(sigma) * mask;
		e_sigma = configs[0].e_sigma*mask; //First derivative available analytically
		e_nsigma = 0.5*(nDen * (configs[1].e_sigma - configs[2].e_sigma) + sigmaDen * (configs[3].e_n - configs[4].e_n));
		e_sigmasigma = sigmaDen * (configs[3].e_sigma - configs[4].e_sigma);
	}
	else
	{	e_sigma = 0;
		e_nsigma = 0;
		e_sigmasigma = 0;
	}
}
