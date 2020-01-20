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

#ifndef JDFTX_CORE_COULOMB_INTERNAL_H
#define JDFTX_CORE_COULOMB_INTERNAL_H

//! @addtogroup LongRange
//! @{

//! @file Coulomb_internal.h Shared inline functions / internal declarations for Coulomb framework

#include <core/matrix3.h>
#include <core/Spline.h>
#include <gsl/gsl_integration.h>

//Common citation for Coulomb truncation
#define wsTruncationPaper "R. Sundararaman and T.A. Arias, Phys. Rev. B 87, 165122 (2013)"

//Ion-margin error message
#define ionMarginMessage "Expand unit cell, or if absolutely sure, reduce coulomb-truncation-ion-margin.\n"

//! Periodic coulomb interaction (4 pi/G^2)
struct CoulombPeriodic_calc
{	__hostanddev__ double operator()(const vector3<int>& iG, const matrix3<>& GGT) const
	{	double Gsq = GGT.metric_length_squared(iG);
		return Gsq ? (4*M_PI)/Gsq : 0.;
	}
	__hostanddev__ symmetricMatrix3<> latticeGradient(const vector3<int>& iG, const matrix3<>& GGT) const
	{	double Gsq = GGT.metric_length_squared(iG);
		return (Gsq ? (8*M_PI)/(Gsq*Gsq) : 0.) * outer(vector3<>(iG));
	}
};

//! Slab-truncated coulomb interaction
struct CoulombSlab_calc
{	int iDir; double hlfL;
	CoulombSlab_calc(int iDir, double hlfL) : iDir(iDir), hlfL(hlfL) {}
	__hostanddev__ double operator()(const vector3<int>& iG, const matrix3<>& GGT) const
	{	double Gsq = GGT.metric_length_squared(iG);
		double Gplane = Gsq - GGT(iDir,iDir) * iG[iDir]*iG[iDir]; //G along the non-truncated directions
		Gplane = Gplane>0. ? sqrt(Gplane) : 0.; //safe sqrt to prevent NaN from roundoff errors
		return (4*M_PI) * (Gsq ? (1. - exp(-Gplane*hlfL) * cos(M_PI*iG[iDir]))/Gsq : -0.5*hlfL*hlfL);
	}
	__hostanddev__ symmetricMatrix3<> latticeGradient(const vector3<int>& iG, const matrix3<>& GGT) const
	{	double Gsq = GGT.metric_length_squared(iG);
		symmetricMatrix3<> result;
		if(Gsq)
		{	double GsqInv = 1./Gsq;
			vector3<int> iGplane = iG; iGplane[iDir] = 0.; //iG along non-truncated directions
			double Gplane = sqrt(GGT.metric_length_squared(iGplane));
			double GplaneInv = Gplane ? 1./Gplane : 0.;
			double expCosTerm = exp(-Gplane*hlfL) * cos(M_PI*iG[iDir]);
			//Compute planar part:
			double prefac1 = 2.*(1.-expCosTerm)*GsqInv;
			symmetricMatrix3<> result = (prefac1 - expCosTerm*hlfL*GplaneInv) * outer(vector3<>(iGplane));
			//Fix normal part:
			((double*)&result)[iDir] = prefac1*iG[iDir]*iG[iDir] + expCosTerm*hlfL*Gplane/GGT(iDir,iDir);
			return (4*M_PI*GsqInv) * result;
		}
		else
		{	symmetricMatrix3<> result;
			((double*)&result)[iDir] = (-4*M_PI)*hlfL*hlfL/GGT(iDir,iDir);
			return result;
		}
	}
};

//! Sphere-truncated coulomb interaction
struct CoulombSpherical_calc
{	double Rc;
	CoulombSpherical_calc(double Rc) : Rc(Rc) {}
	__hostanddev__ double operator()(const vector3<int>& iG, const matrix3<>& GGT) const
	{	double Gsq = GGT.metric_length_squared(iG);
		return Gsq ? (4*M_PI) * (1. - cos(Rc*sqrt(Gsq)))/Gsq : (2*M_PI)*Rc*Rc;
	}
	__hostanddev__ symmetricMatrix3<> latticeGradient(const vector3<int>& iG, const matrix3<>& GGT) const
	{	double Gsq = GGT.metric_length_squared(iG);
		double GRc = sqrt(Gsq)*Rc, cosGRc, sinGRc;
		sincos(GRc, &sinGRc, &cosGRc);
		return (Gsq ? (4*M_PI)*(2*(1.-cosGRc) - GRc*sinGRc)/(Gsq*Gsq) : 0.) * outer(vector3<>(iG));
	}
};

//! Lattice derivative calculation for ionKernel used in embedded mode
struct CoulombIonKernel_calc
{	double expFac;
	CoulombIonKernel_calc(double ionWidth) : expFac(0.5*ionWidth*ionWidth) {}
	__hostanddev__ symmetricMatrix3<> latticeGradient(const vector3<int>& iG, const matrix3<>& GGT) const
	{	double Gsq = GGT.metric_length_squared(iG);
		return (Gsq ? (8*M_PI)*(1.-exp(-expFac*Gsq)*(1.+expFac*Gsq))/(Gsq*Gsq) : 0.) * outer(vector3<>(iG));
	}
};


#ifdef GPU_ENABLED
void coulombAnalytic_gpu(vector3<int> S, const matrix3<>& GGT, const CoulombPeriodic_calc& calc, complex* data);
void coulombAnalytic_gpu(vector3<int> S, const matrix3<>& GGT, const CoulombSlab_calc& calc, complex* data);
void coulombAnalytic_gpu(vector3<int> S, const matrix3<>& GGT, const CoulombSpherical_calc& calc, complex* data);
void coulombAnalyticStress_gpu(vector3<int> S, const matrix3<>& GGT, const CoulombPeriodic_calc& calc, const complex* X, const complex* Y, symmetricMatrix3<>* grad_RRT);
void coulombAnalyticStress_gpu(vector3<int> S, const matrix3<>& GGT, const CoulombSlab_calc& calc, const complex* X, const complex* Y, symmetricMatrix3<>* grad_RRT);
void coulombAnalyticStress_gpu(vector3<int> S, const matrix3<>& GGT, const CoulombSpherical_calc& calc, const complex* X, const complex* Y, symmetricMatrix3<>* grad_RRT);
void coulombAnalyticStress_gpu(vector3<int> S, const matrix3<>& GGT, const CoulombIonKernel_calc& calc, const complex* X, const complex* Y, symmetricMatrix3<>* grad_RRT);
void coulombNumericalStress_gpu(vector3<int> S, const matrix3<>& GGT, const symmetricMatrix3<>* Vc_RRT, const complex* X, const complex* Y, symmetricMatrix3<>* grad_RRT);
#endif
void coulombAnalytic(vector3<int> S, const matrix3<>& GGT, const CoulombPeriodic_calc& calc, complex* data);
void coulombAnalytic(vector3<int> S, const matrix3<>& GGT, const CoulombSlab_calc& calc, complex* data);
void coulombAnalytic(vector3<int> S, const matrix3<>& GGT, const CoulombSpherical_calc& calc, complex* data);
void coulombAnalyticStress(vector3<int> S, const matrix3<>& GGT, const CoulombPeriodic_calc& calc, const complex* X, const complex* Y, symmetricMatrix3<>* grad_RRT);
void coulombAnalyticStress(vector3<int> S, const matrix3<>& GGT, const CoulombSlab_calc& calc, const complex* X, const complex* Y, symmetricMatrix3<>* grad_RRT);
void coulombAnalyticStress(vector3<int> S, const matrix3<>& GGT, const CoulombSpherical_calc& calc, const complex* X, const complex* Y, symmetricMatrix3<>* grad_RRT);
void coulombAnalyticStress(vector3<int> S, const matrix3<>& GGT, const CoulombIonKernel_calc& calc, const complex* X, const complex* Y, symmetricMatrix3<>* grad_RRT);
void coulombNumericalStress(vector3<int> S, const matrix3<>& GGT, const symmetricMatrix3<>* Vc_RRT, const complex* X, const complex* Y, symmetricMatrix3<>* grad_RRT);

//! Compute erf(x)/x (with x~0 handled properly)
__hostanddev__ double erf_by_x(double x)
{	double xSq = x*x;
	if(xSq<1e-6) return (1./sqrt(M_PI))*(2. + xSq*(-2./3 + 0.2*xSq));
	else return erf(x)/x;
}

//! Compute (1/x) d(erf(x)/x)/dx (with x~0 handled properly)
__hostanddev__ double erf_by_xPrime_by_x(double x)
{	double xSq = x*x;
	if(xSq<1e-6) return (2./sqrt(M_PI))*(-2./3 + 0.4*xSq);
	else return ((2./sqrt(M_PI))*x*exp(-xSq) - erf(x))/(x*xSq);
}



//--------------- Special function for cylinder/wire modes ------------
//                 (implemented in CoulombWire.cpp)

//! Compute Cbar_k^sigma - the gaussian convolved cylindrical coulomb kernel - by numerical quadrature
struct Cbar
{	Cbar();
	~Cbar();
	double operator()(double k, double sigma, double rho, double rho0=1.); //!< Compute Cbar_k^sigma(rho)
private:
	static const size_t maxIntervals = 1000; //!< Size of integration workspace
	gsl_integration_workspace* iWS; //!< Integration workspace
	static double integrandSmallRho(double t, void* params); //!< Integrand for rho < sigma
	static double integrandLargeRho(double t, void* params); //!< Integrand for rho > sigma
};

//! Look-up table for Cbar_k^sigma(rho) for specific values of k and sigma
//! If prime=true, construct lookup table for -d(Cbar_k^sigma(rho))/dk instead
struct Cbar_k_sigma
{	Cbar_k_sigma(double k, double sigma, double rhoMax, double rho0=1., bool prime=false);
	//! Get value:
	inline double value(double rho) const
	{	double f = QuinticSpline::value(coeff.data(), drhoInv * rho);
		return isLog ? exp(f) : f;
	}
	//! Get derivative:
	inline double deriv(double rho) const
	{	double fp = QuinticSpline::deriv(coeff.data(), drhoInv * rho) * drhoInv;
		return isLog ? fp * value(rho) : fp;
	}
private:
	double drhoInv; bool isLog;
	std::vector<double> coeff;
};


//---------------------- Exchange Kernels --------------------

//In each of the following functions, kSq is the square of the appropriate
//wave vector (includes reciprocal lattice vector and k-point difference),
//and will not be zero (the G=0 term is handled in the calling routine)

//! Radial fourier transform of erfc(omega r)/r (not valid at G=0)
__hostanddev__ double erfcTilde(double Gsq, double omegaSq)
{	return (4*M_PI) * (omegaSq ? (1.-exp(-0.25*Gsq/omegaSq)) : 1.) / Gsq;
}


//! Periodic exchange
struct ExchangePeriodic_calc
{	__hostanddev__ double operator()(double kSq) const
	{	return (4*M_PI) / kSq;
	}
	__hostanddev__ double latticeGradientPrefac(double kSq, double& grad_lnDetR) const
	{	return (8*M_PI) / (kSq*kSq);
	}
};

//! Erfc-screened Periodic exchange
struct ExchangePeriodicScreened_calc
{	double inv4omegaSq; //!< 1/(4 omega^2)
	ExchangePeriodicScreened_calc(double omega) : inv4omegaSq(0.25/(omega*omega)) {}
	
	__hostanddev__ double operator()(double kSq) const
	{	return (4*M_PI) * (1.-exp(-inv4omegaSq*kSq)) / kSq;
	}
	__hostanddev__ double latticeGradientPrefac(double kSq, double& grad_lnDetR) const
	{	return (8*M_PI) * (1.-exp(-inv4omegaSq*kSq)*(1.+inv4omegaSq*kSq)) / (kSq*kSq);
	}
};

//! Spherical-truncated exchange
struct ExchangeSpherical_calc
{	double Rc;
	ExchangeSpherical_calc(double Rc) : Rc(Rc) {}
	
	__hostanddev__ double operator()(double kSq) const
	{	return (4*M_PI) * (1. - cos(Rc * sqrt(kSq))) / kSq;
	}
	__hostanddev__ double latticeGradientPrefac(double kSq, double& grad_lnDetR) const
	{	double kR = sqrt(kSq)*Rc, coskR, sinkR;
		sincos(kR, &sinkR, &coskR);
		grad_lnDetR = (4*M_PI) * (kR*sinkR) / (3*kSq);
		return (4*M_PI) * (2*(1.-coskR) - kR*sinkR)/(kSq*kSq);
	}
};

//! Erfc-screened Spherical-truncated exchange
struct ExchangeSphericalScreened_calc
{	double* coeff; //!< quintic spline coefficients
	double dGinv; //!< inverse of coefficient spacing
	size_t nSamples; //!< number of coefficients
	ExchangeSphericalScreened_calc() : coeff(0) {}
	
	__hostanddev__ double operator()(double kSq) const
	{	double t = dGinv * sqrt(kSq);
		if(t >= nSamples) return 0.;
		else return QuinticSpline::value(coeff, t);
	}

	//Additional members set for stress calculation:
	double Rc; //!< truncation radius
	double erfcOmegaRc_4piBy3; //!< (4 pi/3) erfc(omega*Rc), where omega is the screening parameter
	
	__hostanddev__ double latticeGradientPrefac(double kSq, double& grad_lnDetR) const
	{	double k = sqrt(kSq);
		grad_lnDetR = erfcOmegaRc_4piBy3 * sin(k*Rc) * Rc/k;
		double t = dGinv * k;
		if(t >= nSamples) return 0.;
		else return QuinticSpline::deriv(coeff, t) * (-t/kSq);
	}
};

//! Slab-truncated exchange
struct ExchangeSlab_calc
{	int iDir; double hlfL;
	double* coeff; double dGinv; size_t nSamples, nCoeff; //quintic-spline coefficients for screened mode
	ExchangeSlab_calc() : coeff(0) {}
	__hostanddev__ double operator()(const vector3<int>& iG, const matrix3<>& GGT, const vector3<>& kDiff, double Vzero, double thresholdSq) const
	{	vector3<> g = iG + kDiff; //net G-vector in reciprocal lattice coordinates including k-point
		double Gsq = GGT.metric_length_squared(g);
		if(Gsq < thresholdSq)
			return Vzero;
		double Gplane = Gsq - GGT(iDir,iDir) * iG[iDir]*iG[iDir]; //G along the non-truncated directions (note kDiff[iDir]=0)
		Gplane = Gplane>0. ? sqrt(Gplane) : 0.; //safe sqrt to prevent NaN from roundoff errors
		double Vc = (4*M_PI) * (1. - exp(-Gplane*hlfL) * cos(M_PI*iG[iDir]))/Gsq; //Unscreened exchange (calculate analytically)
		if(coeff)
		{	//Correct for Screened exchange using lookup table:
			const double* coeffPlane = coeff + abs(iG[iDir]) * nCoeff; //coefficients for this plane
			double t = dGinv * Gplane;
			double prefac = iG[iDir] ? 1. : 1./Gplane;
			if(t<nSamples) Vc += prefac * QuinticSpline::value(coeffPlane, t);
		}
		return Vc;
	}
};


void exchangeAnalytic(vector3<int> S, const matrix3<>& GGT, const ExchangePeriodic_calc& calc, complex* data, const vector3<>& kDiff, double Vzero, double thresholdSq);
void exchangeAnalytic(vector3<int> S, const matrix3<>& GGT, const ExchangePeriodicScreened_calc& calc, complex* data, const vector3<>& kDiff, double Vzero, double thresholdSq);
void exchangeAnalytic(vector3<int> S, const matrix3<>& GGT, const ExchangeSpherical_calc& calc, complex* data, const vector3<>& kDiff, double Vzero, double thresholdSq);
void exchangeAnalytic(vector3<int> S, const matrix3<>& GGT, const ExchangeSphericalScreened_calc& calc, complex* data, const vector3<>& kDiff, double Vzero, double thresholdSq);
void exchangeAnalytic(vector3<int> S, const matrix3<>& GGT, const ExchangeSlab_calc& calc, complex* data, const vector3<>& kDiff, double Vzero, double thresholdSq);
#ifdef GPU_ENABLED
void exchangeAnalytic_gpu(vector3<int> S, const matrix3<>& GGT, const ExchangePeriodic_calc& calc, complex* data, const vector3<>& kDiff, double Vzero, double thresholdSq);
void exchangeAnalytic_gpu(vector3<int> S, const matrix3<>& GGT, const ExchangePeriodicScreened_calc& calc, complex* data, const vector3<>& kDiff, double Vzero, double thresholdSq);
void exchangeAnalytic_gpu(vector3<int> S, const matrix3<>& GGT, const ExchangeSpherical_calc& calc, complex* data, const vector3<>& kDiff, double Vzero, double thresholdSq);
void exchangeAnalytic_gpu(vector3<int> S, const matrix3<>& GGT, const ExchangeSphericalScreened_calc& calc, complex* data, const vector3<>& kDiff, double Vzero, double thresholdSq);
void exchangeAnalytic_gpu(vector3<int> S, const matrix3<>& GGT, const ExchangeSlab_calc& calc, complex* data, const vector3<>& kDiff, double Vzero, double thresholdSq);
#endif

void exchangeAnalyticStress(vector3<int> S, const matrix3<>& G, const ExchangePeriodic_calc& calc, const complex* X, symmetricMatrix3<>* grad_RRT, const vector3<>& kDiff, double thresholdSq);
void exchangeAnalyticStress(vector3<int> S, const matrix3<>& G, const ExchangePeriodicScreened_calc& calc, const complex* X, symmetricMatrix3<>* grad_RRT, const vector3<>& kDiff, double thresholdSq);
void exchangeAnalyticStress(vector3<int> S, const matrix3<>& G, const ExchangeSpherical_calc& calc, const complex* X, symmetricMatrix3<>* grad_RRT, const vector3<>& kDiff, double thresholdSq);
void exchangeAnalyticStress(vector3<int> S, const matrix3<>& G, const ExchangeSphericalScreened_calc& calc, const complex* X, symmetricMatrix3<>* grad_RRT, const vector3<>& kDiff, double thresholdSq);
//void exchangeAnalyticStress(vector3<int> S, const matrix3<>& G, const ExchangeSlab_calc& calc, const complex* X, symmetricMatrix3<>* grad_RRT, const vector3<>& kDiff, double thresholdSq);
#ifdef GPU_ENABLED
void exchangeAnalyticStress_gpu(vector3<int> S, const matrix3<>& G, const ExchangePeriodic_calc& calc, const complex* X, symmetricMatrix3<>* grad_RRT, const vector3<>& kDiff, double thresholdSq);
void exchangeAnalyticStress_gpu(vector3<int> S, const matrix3<>& G, const ExchangePeriodicScreened_calc& calc, const complex* X, symmetricMatrix3<>* grad_RRT, const vector3<>& kDiff, double thresholdSq);
void exchangeAnalyticStress_gpu(vector3<int> S, const matrix3<>& G, const ExchangeSpherical_calc& calc, const complex* X, symmetricMatrix3<>* grad_RRT, const vector3<>& kDiff, double thresholdSq);
void exchangeAnalyticStress_gpu(vector3<int> S, const matrix3<>& G, const ExchangeSphericalScreened_calc& calc, const complex* X, symmetricMatrix3<>* grad_RRT, const vector3<>& kDiff, double thresholdSq);
//void exchangeAnalyticStress(vector3<int> S, const matrix3<>& G, const ExchangeSlab_calc& calc, const complex* X, symmetricMatrix3<>* grad_RRT, const vector3<>& kDiff, double thresholdSq);
#endif


//! Multiply a complexScalarFieldTilde's data by a RealKernel (real-symmetry reduced)
__hostanddev__ void multRealKernel_calc(size_t i, const vector3<int>& iG,
	const vector3<int>& S, const double* kernel, complex* data)
{	//Compute index on the real kernel:
	vector3<int> iGreal = iG;
	if(iGreal[2]<0) iGreal = -iGreal; //inversion symmetry in G-space for real-kernels
	if(iGreal[1]<0) iGreal[1] += S[1];
	if(iGreal[0]<0) iGreal[0] += S[0];
	size_t iReal = iGreal[2] + size_t(1+S[2]/2) * (iGreal[1] + S[1]*iGreal[0]);
	//Multiply:
	data[i] *= kernel[iReal];
}
void multRealKernel(vector3<int> S, const double* kernel, complex* data);
#ifdef GPU_ENABLED
void multRealKernel_gpu(vector3<int> S, const double* kernel, complex* data);
#endif

//! Multiply a complexScalarFieldTilde's data by a kernel sampled with offset and rotation by rot
__hostanddev__ void multTransformedKernel_calc(size_t i, const vector3<int>& iG,
	const vector3<int>& S, const double* kernel, complex* data, const vector3<int>& offset)
{	vector3<int> iGkernel = (iG - offset); //Compute index on the real kernel
	for(int k=0; k<3; k++) if(iGkernel[k]<0) iGkernel[k] += S[k]; //Reduce to [0,S-1) in each dimension
	size_t iReal = iGkernel[2] + S[2]*size_t(iGkernel[1] + S[1]*iGkernel[0]); //net index into kernel
	data[i] *= kernel[iReal];
}
void multTransformedKernel(vector3<int> S, const double* kernel, complex* data, const vector3<int>& offset);
#ifdef GPU_ENABLED
void multTransformedKernel_gpu(vector3<int> S, const double* kernel, complex* data, const vector3<int>& offset);
#endif


//! Compute stress corresponding to multRealKernel()
__hostanddev__ void realKernelStress_calc(size_t i, const vector3<int>& iG,
	const vector3<int>& S, const symmetricMatrix3<>* kernel_RRT, const complex* X, symmetricMatrix3<>* grad_RRT)
{	//Compute index on the real kernel:
	vector3<int> iGreal = iG;
	if(iGreal[2]<0) iGreal = -iGreal; //inversion symmetry in G-space for real-kernels
	if(iGreal[1]<0) iGreal[1] += S[1];
	if(iGreal[0]<0) iGreal[0] += S[0];
	size_t iReal = iGreal[2] + size_t(1+S[2]/2) * (iGreal[1] + S[1]*iGreal[0]);
	//Multiply:
	grad_RRT[i] = kernel_RRT[iReal] * X[i].norm();
}
void realKernelStress(vector3<int> S, const symmetricMatrix3<>* kernel_RRT, const complex* X, symmetricMatrix3<>* grad_RRT);
#ifdef GPU_ENABLED
void realKernelStress_gpu(vector3<int> S, const symmetricMatrix3<>* kernel_RRT, const complex* X, symmetricMatrix3<>* grad_RRT);
#endif

//! @}
#endif // JDFTX_CORE_COULOMB_INTERNAL_H
