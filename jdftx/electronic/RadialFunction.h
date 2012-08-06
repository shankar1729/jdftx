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

#include <core/scalar.h>
#include <vector>

//! G-space radial function stored on a uniform grid (of |G|)
class RadialFunctionG
{
	double dGinv; //!< inverse sample spacing
	int nCoeff; //!< number of coefficients
	double* coeff; //!< coefficients (either on cpu or gpu depending on isGpuEnabled())
public:
	RadialFunctionG();
	operator bool() const; //!< test null-ness
	void init(int l, int nSamples, double dG, const char* filename, double scale=1.0); //!< read and initialize from an ascii file (DFT PSP format)
	void init(int l, const std::vector<double>& samples, double dG); //!< initialize from an array of samples in memory
	void free();
	
	//! Blip (cubic spline evaluation)
	__hostanddev__ double operator()(double G) const
	{	double Gindex = G * dGinv;
		int j = (int)Gindex;
		double tR = Gindex-j; //right weight for interval
		double tL = 1.-tR; //left weight for interval
		//Load blip coefficients:
		if(j+4 > nCoeff) return 0.;
		double c0 = coeff[j], c1 = coeff[j+1], c2 = coeff[j+2], c3 = coeff[j+3];
		//Convert to bernstein polynomial coefficients:
		double b0 = c0 + 4*c1 + c2;
		double b1 = 4*c1 + 2*c2;
		double b2 = 2*c1 + 4*c2;
		double b3 = c1 + 4*c2 + c3;
		//Evaluate bernstein polynomial by de Casteljau's reduction (for best numerical stability)
		//3->2
		double d0 = tL*b0 + tR*b1;
		double d1 = tL*b1 + tR*b2;
		double d2 = tL*b2 + tR*b3;
		//2->1
		double e0 = tL*d0 + tR*d1;
		double e1 = tL*d1 + tR*d2;
		//1->0 (result with scale factor)
		return 0.25*(tL*e0 + tR*e1);
	}
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
