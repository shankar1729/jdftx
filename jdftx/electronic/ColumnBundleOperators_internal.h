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

#ifndef JDFTX_ELECTRONIC_COLUMNBUNDLEOPERATORS_INTERNAL_H
#define JDFTX_ELECTRONIC_COLUMNBUNDLEOPERATORS_INTERNAL_H

#include <core/matrix3.h>

//! @cond

__hostanddev__ void reducedL_calc(int j, int nbasis, int ncols, const complex* Y, complex* LY,
	const matrix3<> GGT, const vector3<int>* iGarr, const vector3<> k, double detR)
{	for (int i=0; i < ncols; i++)
		LY[nbasis*i+j] = (-detR * GGT.metric_length_squared(iGarr[j]+k)) * Y[nbasis*i+j];
}
__hostanddev__ void reducedLinv_calc(int j, int nbasis, int ncols, const complex* Y, complex* LinvY,
	const matrix3<> GGT, const vector3<int>* iGarr, const vector3<> k, double detR)
{	for (int i=0; i < ncols; i++)
	{	double G2 = GGT.metric_length_squared(iGarr[j]+k);
		LinvY[nbasis*i+j] = (G2 ? -1./(detR*G2) : 0.) * Y[nbasis*i+j];
	}
}

__hostanddev__ void precond_inv_kinetic_calc(int j, int nbasis, int ncols, complex* Ydata,
	double KErollover, const matrix3<> GGT, const vector3<int>* iGarr, const vector3<> k, double invdetR)
{
	double x = 0.5*GGT.metric_length_squared(iGarr[j]+k)/KErollover;
	double precondFactor = 1.+x*(1.+x*(1.+x*(1.+x*(1.+x*(1.+x*(1.+x*(1.+x)))))));
	precondFactor = precondFactor*invdetR/(1.+x*precondFactor);
	for(int i=0; i < ncols; i++)
		Ydata[nbasis*i+j] *= precondFactor;
}

__hostanddev__ void precond_inv_kinetic_band_calc(int j, int nbasis, int ncols, complex* Ydata, const double* KEref,
	const matrix3<>& GGT, const vector3<int>* iGarr, const vector3<>& k)
{
	double KE = 0.5*GGT.metric_length_squared(iGarr[j]+k);
	for(int i=0; i<ncols; i++)
	{	double x = KE/KEref[i];
		Ydata[nbasis*i+j] *= (27.+x*(18.+x*(12.+x*8.))) / (27.+x*(18.+x*(12.+x*(8.+x*16))));
	}
}

__hostanddev__
void translate_calc(int j, int nbasis, int ncols, complex* Y, const vector3<int>* iGarr, const vector3<>& k, const vector3<>& dr)
{	complex tFactor = cis(-2*M_PI*dot(iGarr[j]+k,dr));
	for(int i=0; i<ncols; i++)
		Y[nbasis*i+j] *= tFactor;
}

__hostanddev__
void translateColumns_calc(int j, int nbasis, int ncols, complex* Y, const vector3<int>* iGarr, const vector3<>& k, const vector3<>* dr)
{	for(int i=0; i<ncols; i++)
		Y[nbasis*i+j] *= cis(-2*M_PI*dot(iGarr[j]+k,dr[i]));
}

__hostanddev__ void reducedD_calc(int j, int nbasis, int ncols, const complex* Ydata, complex* DYdata,
	const vector3<int>* iGarr, double kdotGe, const vector3<> Ge)
{
	complex Di(0, kdotGe+dot(iGarr[j],Ge)); // i (k+G).e where e is the cartesian direction of interest
	for(int i=0; i < ncols; i++)
		DYdata[nbasis*i+j] = Di*Ydata[nbasis*i+j];
}

__hostanddev__ void reducedDD_calc(int j, int nbasis, int ncols, const complex* Ydata, complex* DYdata,
	const vector3<int>* iGarr, double kdotGe1, double kdotGe2, const vector3<> Ge1, const vector3<> Ge2)
{
	complex Di(0, kdotGe1+dot(iGarr[j],Ge1)); // i (k+G).ei
	complex Dj(0, kdotGe2+dot(iGarr[j],Ge2)); // i (k+G).ej
	for(int i=0; i < ncols; i++)
		DYdata[nbasis*i+j] = Di*Dj*Ydata[nbasis*i+j];
}

//! @endcond
#endif // JDFTX_ELECTRONIC_COLUMNBUNDLEOPERATORS_INTERNAL_H
