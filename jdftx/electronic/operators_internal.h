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

#ifndef JDFTX_ELECTRONIC_OPERATORS_INTERNAL_H
#define JDFTX_ELECTRONIC_OPERATORS_INTERNAL_H

#include <core/matrix3.h>
#include <electronic/RadialFunction.h>

__hostanddev__ void D_calc(int i, const vector3<int>& iG, const complex* in, complex* out,
	const vector3<>& Ge)
{	out[i] = in[i] * complex(0, dot(iG,Ge));
}
__hostanddev__ void DD_calc(int i, const vector3<int>& iG, const complex* in, complex* out,
	const vector3<>& Ge1, const vector3<>& Ge2)
{	out[i] =  in[i] * (-dot(iG,Ge1) * dot(iG,Ge2));
}

__hostanddev__ complex blochPhase_calc(const vector3<int>& iv, const vector3<>& invS, const vector3<>& k)
{	return cis(2*M_PI*dot(k, vector3<>(iv[0]*invS[0], iv[1]*invS[1], iv[2]*invS[2])));
}

template<typename scalar> __hostanddev__
void pointGroupGather_calc(int i, const vector3<int>& iv, const vector3<int>& S,
	const scalar* in, scalar* out, const matrix3<int>& mMesh)
{	vector3<int> ivRot = mMesh * iv; //rotated index vector
	//Project back into range:
	for(int j=0; j<3; j++)
	{	ivRot[j] = ivRot[j] % S[j];
		if(ivRot[j]<0) ivRot[j] += S[j];
	}
	int iRot = ivRot[2]+S[2]*(ivRot[1]+S[1]*ivRot[0]); //rotated index
	out[i] = in[iRot];
}

__hostanddev__ complex radialFunction_calc(const vector3<int>& iG, const matrix3<>& GGT, const RadialFunctionG& f, const vector3<>& r0)
{	return f(sqrt(GGT.metric_length_squared(iG))) * cis(-2*M_PI*dot(iG,r0));
}

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

__hostanddev__ void precond_inv_kinetic_calc(int j, int nbasis, int ncols, const complex* Ydata, complex* KYdata,
	double KErollover, const matrix3<> GGT, const vector3<int>* iGarr, const vector3<> k, double invdetR)
{
	double x = 0.5*GGT.metric_length_squared(iGarr[j]+k)/KErollover;
	double precondFactor = 1.0+x*(1.0+x*(1.0+x*(1.0+x*(1.0+x*(1.0+x*(1.0+x*(1.0+x)))))));
	precondFactor = precondFactor*invdetR/(1.0+x*precondFactor);
	for(int i=0; i < ncols; i++)
		KYdata[nbasis*i+j] = precondFactor * Ydata[nbasis*i+j];
}

__hostanddev__
void translate_calc(int j, int nbasis, int ncols, complex* Y, const vector3<int>* iGarr, const vector3<>& k, const vector3<>& dr)
{	complex tFactor = cis(-2*M_PI*dot(iGarr[j]+k,dr));
	for(int i=0; i<ncols; i++)
		Y[nbasis*i+j] *= tFactor;
}

__hostanddev__ void reducedD_calc(int j, int nbasis, int ncols, const complex* Ydata, complex* DYdata,
	const vector3<int>* iGarr, double kdotGe, const vector3<> Ge)
{
	complex Di(0, kdotGe+dot(iGarr[j],Ge)); // i (k+G).e where e is the cartesian direction of interest
	for(int i=0; i < ncols; i++)
		DYdata[nbasis*i+j] = Di*Ydata[nbasis*i+j];
}

#endif // JDFTX_ELECTRONIC_OPERATORS_INTERNAL_H
