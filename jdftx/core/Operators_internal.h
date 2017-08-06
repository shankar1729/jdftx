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

#ifndef JDFTX_CORE_OPERATORS_INTERNAL_H
#define JDFTX_CORE_OPERATORS_INTERNAL_H

#include <core/matrix3.h>
#include <core/tensor3.h>
#include <core/LoopMacros.h>
#include <core/SphericalHarmonics.h>
#include <core/RadialFunction.h>

//! @cond

//! Compute the full G-space indices corresponding to a half Gspace index and its mirror image:
#define COMPUTE_fullRefIndices \
	vector3<int> iv(iG), ivRef; \
	for(int k=0; k<3; k++) \
	{	if(iv[k]<0) iv[k] += S[k]; \
		ivRef[k] = iv[k] ? S[k] - iv[k] : 0; \
	} \
	int iFull = iv[2] + S[2]*(iv[1] + S[1]*iv[0]); \
	int iFullRef = ivRef[2] + S[2]*(ivRef[1] + S[1]*ivRef[0]);


__hostanddev__
void RealG_calc(int iHalf, const vector3<int> iG, const vector3<int> S,
	const complex* vFull, complex* vHalf, double scaleFac)
{
	COMPUTE_fullRefIndices
	vHalf[iHalf] = scaleFac*0.5*(vFull[iFull] + vFull[iFullRef].conj());
}

__hostanddev__
void ImagG_calc(int iHalf, const vector3<int> iG, const vector3<int> S,
	const complex* vFull, complex* vHalf, double scaleFac)
{
	COMPUTE_fullRefIndices
	vHalf[iHalf] = complex(0,-scaleFac*0.5)*(vFull[iFull] - vFull[iFullRef].conj());
}

__hostanddev__
void ComplexG_calc(int iHalf, const vector3<int> iG, const vector3<int> S,
	const complex* vHalf, complex *vFull, double scaleFac)
{
	COMPUTE_fullRefIndices
	complex temp = scaleFac*vHalf[iHalf];
	vFull[iFull] = temp; //Copy value into corresponding location in full-space
	vFull[iFullRef] = temp.conj(); //Also set complex conjugate into the mirror location
}

__hostanddev__
void changeGrid_calc(const vector3<int>& iG, const vector3<int>& Sin, const vector3<int>& Sout, const complex* in, complex* out)
{	//Compute index:
	#define COMPUTE_index(suffix) \
		int i##suffix = 0; \
		for(int k=0; k<2; k++) \
		{	if(2*iG[k]<1-S##suffix[k] || 2*iG[k]>S##suffix[k]) return; \
			i##suffix = i##suffix * S##suffix[k] + (iG[k]<0 ? (iG[k]+S##suffix[k]) : iG[k]); \
		} \
		if(2*iG[2]>S##suffix[2]) return; \
		else i##suffix = i##suffix*(1+S##suffix[2]/2) + iG[2];
	COMPUTE_index(in)
	COMPUTE_index(out)
	#undef COMPUTE_index
	out[iout] = in[iin];
}
__hostanddev__
void changeGridFull_calc(const vector3<int>& iG, const vector3<int>& Sin, const vector3<int>& Sout, const complex* in, complex* out)
{	//Compute index:
	#define COMPUTE_index(suffix) \
		int i##suffix = 0; \
		for(int k=0; k<3; k++) \
		{	if(2*iG[k]<1-S##suffix[k] || 2*iG[k]>S##suffix[k]) return; \
			i##suffix = i##suffix * S##suffix[k] + (iG[k]<0 ? (iG[k]+S##suffix[k]) : iG[k]); \
		}
	COMPUTE_index(in)
	COMPUTE_index(out)
	#undef COMPUTE_index
	out[iout] = in[iin];
}

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
		case 4: fTemplate<4> argList; break; \
		case 5: fTemplate<5> argList; break; \
		case 6: fTemplate<6> argList; break; \
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

__hostanddev__ complex radialFunction_calc(const vector3<int>& iG, const matrix3<>& GGT, const RadialFunctionG& f, const vector3<>& r0)
{	return f(sqrt(GGT.metric_length_squared(iG))) * cis(-2*M_PI*dot(iG,r0));
}


//----------- Functions from VectorField.h -------

__hostanddev__ void gradient_calc(int i, const vector3<int> iG, bool nyq, const matrix3<> G,
	const complex* Xtilde, vector3<complex*>& gradTilde)
{	complex iota(0.0, nyq ? 0.0 : 1.0); //zero nyquist frequencies
	storeVector((iG*G) * (iota*Xtilde[i]), gradTilde, i);
}

__hostanddev__ void divergence_calc(int i, const vector3<int> iG, bool nyq, const matrix3<> G,
	vector3<const complex*>& Vtilde, complex* divTilde)
{	complex iota(0.0, nyq ? 0.0 : 1.0); //zero nyquist frequencies
	divTilde[i] = iota * dot(iG*G, loadVector(Vtilde,i));
}


__hostanddev__ void tensorGradient_calc(int i, const vector3<int> iG, bool nyq, const matrix3<> G,
	const complex* Xtilde, tensor3<complex*>& gradTilde)
{
	complex minus_Xtilde = nyq ? complex(0,0) : -Xtilde[i]; //zero nyquist frequencies
	vector3<> Gvec = iG*G;
	double Gsq = Gvec.length_squared();
	gradTilde.xy()[i] = minus_Xtilde*Gvec.x()*Gvec.y();
	gradTilde.yz()[i] = minus_Xtilde*Gvec.y()*Gvec.z();
	gradTilde.zx()[i] = minus_Xtilde*Gvec.z()*Gvec.x();
	gradTilde.xxr()[i] = minus_Xtilde*(Gvec.x()*Gvec.x() - (1.0/3)*Gsq);
	gradTilde.yyr()[i] = minus_Xtilde*(Gvec.y()*Gvec.y() - (1.0/3)*Gsq);
}

__hostanddev__ void tensorDivergence_calc(int i, const vector3<int> iG, bool nyq, const matrix3<> G,
	tensor3<const complex*>& Vtilde, complex* divTilde)
{
	complex temp = complex(0,0);
	if(!nyq)
	{	vector3<> Gvec = iG*G;
		temp += Vtilde.xy()[i]*( 2*Gvec.x()*Gvec.y() );
		temp += Vtilde.yz()[i]*( 2*Gvec.y()*Gvec.z() );
		temp += Vtilde.zx()[i]*( 2*Gvec.z()*Gvec.x() );
		temp += Vtilde.xxr()[i]*( Gvec.x()*Gvec.x() - Gvec.z()*Gvec.z() );
		temp += Vtilde.yyr()[i]*( Gvec.y()*Gvec.y() - Gvec.z()*Gvec.z() );
	}
	divTilde[i] = -temp;
}

//! @endcond
#endif //JDFTX_CORE_OPERATORS_INTERNAL_H
