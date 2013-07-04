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
#include <electronic/SphericalHarmonics.h>
#include <vector>

//! Struct to wrap a fixed size array for passing to templated functions
//! (Pretty much std::array, but that is not yet supported in CUDA)
template<typename T, int N>
struct array
{	T arr[N];
	array(const std::vector<T>& vec) { for(int s=0; s<N; s++) arr[s]=vec[s]; }
	__hostanddev__ array(T t=0) { for(int s=0; s<N; s++) arr[s]=t; }
	__hostanddev__ T& operator[](int i) { return arr[i]; }
	__hostanddev__ const T& operator[](int i) const { return arr[i]; }
};


__hostanddev__ void D_calc(int i, const vector3<int>& iG, const complex* in, complex* out,
	const vector3<>& Ge)
{	out[i] = in[i] * complex(0, dot(iG,Ge));
}
__hostanddev__ void DD_calc(int i, const vector3<int>& iG, const complex* in, complex* out,
	const vector3<>& Ge1, const vector3<>& Ge2)
{	out[i] =  in[i] * (-dot(iG,Ge1) * dot(iG,Ge2));
}

//! Switch a function templated over l for all supported l with parenthesis enclosed argument list argList
#define SwitchTemplate_l(l,fTemplate,argList) \
	switch(l) \
	{	case 0: fTemplate<0> argList; break; \
		case 1: fTemplate<1> argList; break; \
		case 2: fTemplate<2> argList; break; \
		case 3: fTemplate<3> argList; break; \
	}

template<int l, int lpm> struct lGradient_staticLoop
{	static __hostanddev__ void set(int i, const vector3<>& g, const complex& phasedIn, const array<complex*,2*l+1>& out)
	{	out[lpm][i] = phasedIn * Ylm<l,lpm-l>(g);
		lGradient_staticLoop<l,lpm-1>::set(i, g, phasedIn, out);
	}
};
template<int l> struct lGradient_staticLoop<l,-1> { static __hostanddev__ void set(int i, const vector3<>& g, const complex& phasedIn, const array<complex*,2*l+1>& out) {} }; //end recursion

template<int l> __hostanddev__ void lGradient_calc(int i, const vector3<int>& iG, bool isNyq, const complex* in, const array<complex*,2*l+1>& out, const matrix3<>& G)
{	const complex phase = cis(l*0.5*M_PI); // iota^l (computable at compile time)
	lGradient_staticLoop<l,l+l>::set(i, iG*G, isNyq ? 0. : phase * in[i], out);
}

template<int l, int lpm> struct lDivergence_staticLoop
{	static __hostanddev__ complex get(int i, const vector3<>& g, const array<const complex*,2*l+1>& in)
	{	return in[lpm][i] * Ylm<l,lpm-l>(g) +  lDivergence_staticLoop<l,lpm-1>::get(i, g, in);
	}
};
template<int l> struct lDivergence_staticLoop<l,-1> { static __hostanddev__ complex get(int i, const vector3<>& g, const array<const complex*,2*l+1>& in) { return complex(); } }; //end recursion

template<int l> __hostanddev__ void lDivergence_calc(int i, const vector3<int>& iG, bool isNyq, const array<const complex*,2*l+1>& in, complex* out, const matrix3<>& G)
{	const complex phase = cis(l*0.5*M_PI); // iota^l (computable at compile time)
	out[i] = (isNyq ? 0. : phase) * lDivergence_staticLoop<l,l+l>::get(i, iG*G, in);
}


__hostanddev__ complex blochPhase_calc(const vector3<int>& iv, const vector3<>& invS, const vector3<>& k)
{	return cis(2*M_PI*dot(k, vector3<>(iv[0]*invS[0], iv[1]*invS[1], iv[2]*invS[2])));
}

template<typename scalar> __hostanddev__
void pointGroupScatter_calc(int i, const vector3<int>& iv, const vector3<int>& S,
	const scalar* in, scalar* out, const matrix3<int>& mMesh)
{	vector3<int> ivRot = mMesh * iv; //rotated index vector
	//Project back into range:
	for(int j=0; j<3; j++)
	{	ivRot[j] = ivRot[j] % S[j];
		if(ivRot[j]<0) ivRot[j] += S[j];
	}
	int iRot = ivRot[2]+S[2]*(ivRot[1]+S[1]*ivRot[0]); //rotated index
	out[iRot] = in[i];
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
