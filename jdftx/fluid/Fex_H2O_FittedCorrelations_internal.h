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

#ifndef JDFTX_FLUID_FEX_H2O_FITTEDCORRELATIONS_INTERNAL_H
#define JDFTX_FLUID_FEX_H2O_FITTEDCORRELATIONS_INTERNAL_H

#include <core/scalar.h>

namespace Fex_H2O_FittedCorrelations_internal
{
	//Coefficients from Table II of [FittedCorrelations]:
	static const double f0 = 3.83929979127276e-18;
	static const double f1 = -1.01124264551199e-10;
	static const double f2 = 8.87857550440045e-04;
	static const double f3 = -2.59856289573731e+03;
	static const double f4 = 9.17183837107696e+05;
	static const double f5 = -1.27079357385654e+08;
	static const double f6 = 6.58720860464796e+09;
	static const double pObar = -19.0/20;
	static const double pHbar = +12.0/20;
	static const double pMean = +27.0/20;

	__hostanddev__ double fex(double N)
	{	return N<0. ? 0. : f0 + N*(f1 + N*(f2 + N*(f3 + N*(f4 + N*(f5 + N*f6)))));
	}
	__hostanddev__ double fexDot(double N)
	{	return N<0. ? 0. : f1 + N*(2*f2 + N*(3*f3 + N*(4*f4 + N*(5*f5 + N*(6*f6)))));
	}
}
//Compute the gaussian weighted density energy and gradients
__hostanddev__
double Fex_H2O_FittedCorrelations_calc(int i, const double* NObar, const double* NHbar, double* Phi_NObar, double* Phi_NHbar)
{	using namespace Fex_H2O_FittedCorrelations_internal;
	double Nmean = (1.0/3)*(NObar[i] + NHbar[i]);
	double fDotMean = fexDot(Nmean);
	Phi_NObar[i] = ((1.0/3)*pMean*fDotMean + pObar*fexDot(NObar[i]));
	Phi_NHbar[i] = ((1.0/3)*pMean*fDotMean + pHbar*fexDot(NHbar[i]*0.5)*0.5);
	return pMean*fex(Nmean) + pObar*fex(NObar[i]) + pHbar*fex(NHbar[i]*0.5);
}

#endif // JDFTX_FLUID_FEX_H2O_FITTEDCORRELATIONS_INTERNAL_H
