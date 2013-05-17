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

#ifndef JDFTX_FLUID_ERFFMTWEIGHT_H
#define JDFTX_FLUID_ERFFMTWEIGHT_H

#include <electronic/SphericalHarmonics.h>
#include <cmath>

//! Utility for creating soft FMT weight functions
class ErfFMTweight
{
public:
	ErfFMTweight(double R, double sigma) : R(R),sigma(sigma)
	{	//root and curvature of the saddle point integrand used in w0, w1v and w2m:
		r0 = sqrt(pow(0.5*R,2) - pow(sigma,2)) + 0.5*R;
		c0 = 1 - pow(sigma/r0,2);
	}

	//! set the weights at a given G
	void operator()(double G, double& w0, double& w1, double& w2, double& w3, double& w1v, double& w2m) const
	{	//Weights without saddle point approx:
		double prefac = exp(-0.5*pow(G*sigma,2));
		double j0 = bessel_jl(0, G*R);
		double j2 = bessel_jl(2, G*R);
		double softness = pow(sigma/R,2);
		//---
		w1 = R * prefac * j0;
		w2 = 4*M_PI*pow(R,2) * prefac * (j0 + softness*cos(G*R));
		w3 = (4*M_PI*pow(R,3)/3) * prefac * (j0*(1+3*softness) + j2);

		//Weights with saddle point approx:
		prefac = exp(-0.5*pow(sigma ? (r0-R)/sigma : 0.0, 2))/sqrt(c0) * exp(-(0.5/c0)*pow(G*sigma,2));
		j0 = bessel_jl(0, G*r0);
		j2 = bessel_jl(2, G*r0);
		double j4 = bessel_jl(4, G*r0);
		softness = pow(sigma/r0,2)/c0;
		//---
		w0  = prefac * j0;
		w1v = -pow(r0,2)/3 * prefac * (j0*(1+3*softness) + j2);
		w2m = (4*M_PI*pow(r0,4)/15) * prefac * (j0*(1+softness*(10+softness*15)) + 10*j2*(1.0/7+softness) + (3.0/7)*j4);
	}

private:
	double R, sigma, r0, c0;
};

#endif // JDFTX_FLUID_ERFFMTWEIGHT_H
