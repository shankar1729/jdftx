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

#ifndef JDFTX_ELECTRONIC_RADIALFUNCTION_H
#define JDFTX_ELECTRONIC_RADIALFUNCTION_H

#include <core/Spline.h>

struct RadialFunctionR;

//! G-space radial function stored on a uniform grid (of |G|)
struct RadialFunctionG
{
	double dGinv; //!< inverse sample spacing
	int nCoeff; //!< number of coefficients
	std::vector<double> coeff; //!< coefficients on cpu
	#ifdef GPU_ENABLED
	double* coeffGpu; //!< coefficients on gpu
	const double* coeffPref() const { return coeffGpu; }
	#else
	const double* coeffPref() const { return coeff.data(); }
	#endif

	//Access appropriate coefficients depending on whether called from CPU/GPU
	__hostanddev__ const double* getCoeff() const
	{
		#ifdef __CUDA_ARCH__
		return coeffGpu;
		#else
		return coeff.data();
		#endif
	}

	RadialFunctionG();
	void init(int l, int nSamples, double dG, const char* filename, double scale=1.0); //!< read and initialize from an ascii file (DFT PSP format)
	void init(int l, const std::vector<double>& samples, double dG); //!< initialize from an array of samples in memory
	void set(const std::vector<double>& coeff, double dGInv); //!< set the coefficients (and update the GPU versions etc.)
	void updateGmax(int l, int nSamples); //!< if created from a RadialFunctionR, increase nCoeff if necessary (call when lattice is modified)
	void free(bool rFuncDelete=true);
	
	//! Blip (quintic spline evaluation)
	__hostanddev__ double operator()(double G) const
	{	double Gindex = G * dGinv;
		if(Gindex >= nCoeff-5) return 0.;
		else return QuinticSpline::value(getCoeff(), Gindex);
	}
	
	RadialFunctionR* rFunc; //!< copy of the real-space radial version (if created from one)
	
	#ifndef __in_a_cu_file__
	//! Helper functional for initializing using a function
	template<typename Func, typename... Args> void init(int l, double dG, double Gmax, const Func& func, Args... args)
	{	std::vector<double> samples(unsigned(ceil(Gmax/dG))+5);
		for(unsigned i=0; i<samples.size(); i++) samples[i] = func(i*dG, args...);
		init(l, samples, dG);
	}
	
	//Fourier transform of cuspless exponential (for use with init)
	static inline double cusplessExpTilde(double G, double norm, double a)
	{	double aG = a*G;
		double den = 1./(1.+aG*aG);
		return norm * den*den*den;
	}

       //Fourier transform of exponential (for use with init)
       static inline double exponentialTilde(double G, double norm, double a)
       {	
	        double aG = a*G;
	        double den = 1./(1.+aG*aG);
	        return norm * den*den;
       }

	//Fourier transform of gaussian (for use with init)
	static inline double gaussTilde(double G, double norm, double sigma)
	{	double sigmaG = sigma*G;
		return norm * exp(-0.5*sigmaG*sigmaG);
	}
	
	explicit operator bool() const { return nCoeff; } //!< test null-ness
	#endif
};


//! A function on a non-uniform real-space radial grid
struct RadialFunctionR
{	std::vector<double> r; //!< radial location
	std::vector<double> dr; //!< radial weight
	std::vector<double> f; //!< sample value
	
	RadialFunctionR(int nSamples=0); //!< allocate for a specified sample count
	RadialFunctionR(const std::vector<double>& r, double dlogr); //!< initialize logarithmic grid
	void set(std::vector<double> r, std::vector<double> dr); //!< update the sample locations, weights
	
	//! Compute optimum dr given r (must be in ascending order) which will be exact
	//! for the integration of cubic splines with natural boundary conditions.
	//! It is recommended that r[0]=0 or r[0]<<r[1], but it'll work fine even otherwise.
	void initWeights();
	
	//! Evaluate spherical bessel transform of order l at wavevector G:
	//! @$ func(G) = \int dr 4\pi r^2 j_l(G r) f(r) @$
	double transform(int l, double G) const;
	
	//! Initialize a uniform G radial function from the logPrintf grid function according to
	//! @$ func(G) = \int dr 4\pi r^2 j_l(G r) f(r) @$
	void transform(int l, double dG, int nGrid, RadialFunctionG& func) const;
};

#endif // JDFTX_ELECTRONIC_RADIALFUNCTION_H
