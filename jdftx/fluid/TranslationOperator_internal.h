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

#ifndef JDFTX_FLUID_TRANSLATIONOPERATORINTERNAL_H
#define JDFTX_FLUID_TRANSLATIONOPERATORINTERNAL_H

#include <core/vector3.h>

//! Get the index into data arrays for a vector3<int> which could be in the 2nd zone due to a translation
__hostanddev__ int wrappedIndex(const vector3<int> i, const vector3<int> S)
{
	vector3<int> iWrapped = i;
	for(int k=0; k<3; k++) if(iWrapped[k]>=S[k]) iWrapped[k]-=S[k];
	return iWrapped[2] + S[2]*(iWrapped[1] + S[1]*iWrapped[0]);
}

__hostanddev__
void constantSplineTaxpy_calc(int yIndex, const vector3<int> iy, const vector3<int> S,
	double alpha, const double* x, double* y, const vector3<int> Tint)
{
	y[yIndex] += alpha * x[wrappedIndex(iy+Tint,S)];
}

__hostanddev__
void linearSplineTaxpy_calc(int yIndex, const vector3<int> iy, const vector3<int> S,
	double alpha, const double* x, double* y, const vector3<int> Tint, const vector3<> Tfrac)
{
	//Weights for linear interpolation:
	double w0[] = {1-Tfrac[0], Tfrac[0]};
	double w1[] = {1-Tfrac[1], Tfrac[1]};
	double w2[] = {1-Tfrac[2], Tfrac[2]};
	vector3<int> ix = iy + Tint;
	//Interpolate:
	double temp0 = 0.0;
	for(int i0=0; i0<2; i0++) //loop unrolled by default
	{	double temp1 = 0.0;
		for(int i1=0; i1<2; i1++) //loop unrolled by default
		{	double temp2 = 0.0;
			for(int i2=0; i2<2; i2++) //loop unrolled by default
			{	temp2 += w2[i2] * x[wrappedIndex(ix+vector3<int>(i0,i1,i2),S)];
			}
			temp1 += w1[i1] * temp2;
		}
		temp0 += w0[i0] * temp1;
	}
	y[yIndex] += alpha * temp0;
}

__hostanddev__
void fourierTranslate_calc(int i, const vector3<int> iG, const vector3<int> S, const vector3<> Gt, complex* xTilde)
{	xTilde[i] *= cis(-dot(iG,Gt));
}


#endif // JDFTX_FLUID_TRANSLATIONOPERATORINTERNAL_H
