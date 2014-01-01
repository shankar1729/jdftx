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

#ifndef JDFTX_ELECTRONIC_EXCORR_INTERNAL_H
#define JDFTX_ELECTRONIC_EXCORR_INTERNAL_H

//! @file ExCorr_internal.h
//! Internal abstractions of, and helper routines for the internal exchange and correlation routines

#include <electronic/operators_internal.h>

static const double nCutoff = 1e-16; //!< ignore densities below this value

//! Switch a function fTemplate templated over a functional variant and spin count
//! SwitchTemplate_functional is a macro such as SwitchTemplate_LDA
//! (This is needed to switch from a run-time nCount to a compile-time template argument)
#define SwitchTemplate_spin(SwitchTemplate_functional,variant,nCount, fTemplate,argList) \
	switch(nCount) \
	{	case 1: { SwitchTemplate_functional(variant,1, fTemplate,argList) break; } \
		case 2: { SwitchTemplate_functional(variant,2, fTemplate,argList) break; } \
		default: break; \
	}


//! Abstract base class for functionals
class Functional
{
protected:
	double scaleFac; //!< scale factor (to support mixing for hybrid functionals)
public:
	Functional(double scaleFac=1.0) : scaleFac(scaleFac) {}
	virtual bool needsSigma() const=0; //!< return true if density gradients are used
	virtual bool needsLap() const=0; //!< return true if laplacian of density is used (MGGA)
	virtual bool needsTau() const=0; //!< return true if orbital kinetic energy density is used (MGGA)
	virtual bool isKinetic() const=0; //!< whether this is a kinetic energy functional
	
	//! Compute exchange-correlation energy densities and optionally gradients if E_n[0] is non-null.
	//! Note that if E_n[0] is non-null, then so must all components of E_n and E_sigma (if a GGA).
	//! All the energies and gradients must be scaled by scaleFac.
	//! @param N number of points
	//! @param n (spin-)densities: 1 or 2 vectors for unpolarized / polarized
	//! @param sigma contracted (spin-)density gradient contractions: 1 or 3 vectors for unpolarized / polarized
	//! @param lap laplacian of (spin-)denisties
	//! @param tau kinetic energy density per psin channel (form orbitals)
	//! @param E accumulate energy density per volume
	//! @param E_n accumulate gradient w.r.t n's (if E_n[0] is non-null)
	//! @param E_sigma accumulate gradient w.r.t sigma's (if E_n[0] is non-null)
	//! @param E_lap accumulate gradient w.r.t lap's (if E_n[0] is non-null)
	//! @param E_tau accumulate gradient w.r.t tau's (if E_n[0] is non-null)
	virtual void evaluate(int N, std::vector<const double*> n, std::vector<const double*> sigma,
		std::vector<const double*> lap, std::vector<const double*> tau,
		double* E, std::vector<double*> E_n, std::vector<double*> E_sigma,
		std::vector<double*> E_lap, std::vector<double*> E_tau) const=0;
	
	//!Call evaluate for a subset of data points:
	void evaluateSub(int iStart, int iStop,
		std::vector<const double*> n, std::vector<const double*> sigma,
		std::vector<const double*> lap, std::vector<const double*> tau,
		double* E, std::vector<double*> E_n, std::vector<double*> E_sigma,
		std::vector<double*> E_lap, std::vector<double*> E_tau) const;
};


//! LDA spin interpolation function f(zeta) and its derivative
__hostanddev__ double spinInterpolation(double zeta, double& f_zeta)
{	const double scale = 1./(pow(2.,4./3) - 2);
	double zetaPlusCbrt = pow(1+zeta, 1./3);
	double zetaMinusCbrt = pow(1-zeta, 1./3);
	f_zeta = scale*(zetaPlusCbrt - zetaMinusCbrt)*(4./3);
	return scale*((1+zeta)*zetaPlusCbrt + (1-zeta)*zetaMinusCbrt - 2.);
}

//! Spin-interpolate an LDA functional given its paramagnetic and ferromagnetic functors
template<typename Para, typename Ferro> __hostanddev__ 
double spinInterpolate(double rs, double zeta, double& e_rs, double& e_zeta, const Para& para, const Ferro& ferro)
{	double ePara_rs, ePara = para(rs, ePara_rs);
	if(!zeta)
	{	//Return paramagnetic result:
		e_rs = ePara_rs;
		e_zeta = 0.;
		return ePara;
	}
	else //Mix in ferromagnetic result:
	{	double eFerro_rs, eFerro = ferro(rs, eFerro_rs);
		double f_zeta, f = spinInterpolation(zeta, f_zeta); //spin interpolation function
		e_rs = ePara_rs + f*(eFerro_rs - ePara_rs);
		e_zeta = f_zeta*(eFerro - ePara);
		return ePara + f*(eFerro - ePara);
	}
}

//! Spin-interpolate an LDA functional given its paramagnetic, ferromagnetic and spin-stiffness functors
//! (This is the spin-inteprolation technique used in the VWN and PW correlation functionals)
//! (For numerical compatibility with the original PW routine, the f"(0) scale factor may be over-ridden)
template<typename Para, typename Ferro, typename Stiff> __hostanddev__ 
double spinInterpolate(double rs, double zeta, double& e_rs, double& e_zeta,
	const Para& para, const Ferro& ferro, const Stiff& stiff,
	const double fDblPrime0 = 4./(9*(pow(2., 1./3)-1)))
{
	double ePara_rs, ePara = para(rs, ePara_rs); //Paramagentic
	if(!zeta) //return paramagentic result
	{	e_rs = ePara_rs;
		e_zeta = 0.;
		return ePara;
	}
	else //Mix in ferromagnetic and zeta-derivative results:
	{	double eFerro_rs, eFerro = ferro(rs, eFerro_rs); //Ferromagnetic
		double eStiff_rs, eStiff = stiff(rs, eStiff_rs); //Spin-derivative
		//Compute mix factors:
		double f_zeta, f = spinInterpolation(zeta, f_zeta); //spin interpolation function
		double zeta2=zeta*zeta, zeta3=zeta2*zeta, zeta4=zeta2*zeta2; //powers of zeta
		const double scale = -1./fDblPrime0;
		double w1 = zeta4*f,             w1_zeta = 4*zeta3*f + zeta4*f_zeta;
		double w2 = scale*((1-zeta4)*f), w2_zeta = scale*(-4*zeta3*f + (1-zeta4)*f_zeta);
		//Mix:
		e_rs = ePara_rs + w1*(eFerro_rs-ePara_rs) + w2*eStiff_rs;
		e_zeta = w1_zeta*(eFerro-ePara) + w2_zeta*eStiff;
		return ePara + w1*(eFerro-ePara) + w2*eStiff;
	}
}

#endif // JDFTX_ELECTRONIC_EXCORR_INTERNAL_H

