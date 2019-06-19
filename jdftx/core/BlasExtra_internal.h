/*-------------------------------------------------------------------
Copyright 2016 Ravishankar Sundararaman

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

#ifndef JDFTX_CORE_BLASEXTRA_INTERNAL_H
#define JDFTX_CORE_BLASEXTRA_INTERNAL_H

#include <core/matrix3.h>

//! @cond

//------ Helpers for scatter / axpy -----
template<typename scalar, bool conjx, bool havew, bool conjw> struct Conjugator {};
template<> struct Conjugator<double,false,false,false> { __hostanddev__ double operator()(const double* x, int ix, const double* w, int iw) const { return x[ix]; } };
template<> struct Conjugator<double,false,true,false> { __hostanddev__ double operator()(const double* x, int ix, const double* w, int iw) const { return x[ix]*w[iw]; } };
template<> struct Conjugator<complex,false,false,false> { __hostanddev__ complex operator()(const complex* x, int ix, const complex* w, int iw) const { return x[ix]; } };
template<> struct Conjugator<complex,false,true,false> { __hostanddev__ complex operator()(const complex* x, int ix, const complex* w, int iw) const { return x[ix]*w[iw]; } };
template<> struct Conjugator<complex,false,true,true> { __hostanddev__ complex operator()(const complex* x, int ix, const complex* w, int iw) const { return x[ix]*conj(w[iw]); } };
template<> struct Conjugator<complex,true,false,false> { __hostanddev__ complex operator()(const complex* x, int ix, const complex* w, int iw) const { return conj(x[ix]); } };
template<> struct Conjugator<complex,true,true,false> { __hostanddev__ complex operator()(const complex* x, int ix, const complex* w, int iw) const { return conj(x[ix])*w[iw]; } };
template<> struct Conjugator<complex,true,true,true> { __hostanddev__ complex operator()(const complex* x, int ix, const complex* w, int iw) const { return conj(x[ix]*w[iw]); } };

#define DEFINE_SPARSE_AXPY(type, suffix) \
	template<typename scalar2> void eblas_##type##_axpy##suffix(const int Nindex, scalar2 a, const int* index, const complex* x, complex* y, bool conjx, const complex* w, bool conjw) \
	{	if(conjx) \
		{	if(w) \
			{	if(conjw) eblas_##type##_axpy##suffix(Nindex, a, index, x, y, w, Conjugator<complex,true,true,true>()); \
				else eblas_##type##_axpy##suffix(Nindex, a, index, x, y, w, Conjugator<complex,true,true,false>()); \
			} \
			else eblas_##type##_axpy##suffix(Nindex, a, index, x, y, w, Conjugator<complex,true,false,false>()); \
		} \
		else \
		{	if(w) \
			{	if(conjw) eblas_##type##_axpy##suffix(Nindex, a, index, x, y, w, Conjugator<complex,false,true,true>()); \
				else eblas_##type##_axpy##suffix(Nindex, a, index, x, y, w, Conjugator<complex,false,true,false>()); \
			} \
			else eblas_##type##_axpy##suffix(Nindex, a, index, x, y, w, Conjugator<complex,false,false,false>()); \
		} \
	} \
	void eblas_##type##_zdaxpy##suffix(const int Nindex, double a, const int* index, const complex* x, complex* y, bool conjx, const complex* w, bool conjw) \
	{	eblas_##type##_axpy##suffix(Nindex, a, index, x, y, conjx, w, conjw); \
	} \
	void eblas_##type##_zaxpy##suffix(const int Nindex, complex a, const int* index, const complex* x, complex* y, bool conjx, const complex* w, bool conjw) \
	{	eblas_##type##_axpy##suffix(Nindex, a, index, x, y, conjx, w, conjw); \
	} \
	void eblas_##type##_daxpy##suffix(const int Nindex, double a, const int* index, const double* x, double* y, const double* w) \
	{	if(w) eblas_##type##_axpy##suffix(Nindex, a, index, x, y, w, Conjugator<double,false,true,false>()); \
		else eblas_##type##_axpy##suffix(Nindex, a, index, x, y, w, Conjugator<double,false,false,false>()); \
	}

	
//---- symmetrize helper routines ----

template<typename scalar> __hostanddev__ void eblas_symmetrize_calc(size_t i, int n, const int* symmIndex, scalar* x, double nInv)
{	scalar xSum=0.0;
	for(int j=0; j<n; j++) xSum += x[symmIndex[n*i+j]];
	xSum *= nInv; //average n in the equivalence class
	for(int j=0; j<n; j++) x[symmIndex[n*i+j]] = xSum;
}

__hostanddev__ void eblas_symmetrize_phase_calc(size_t i, int n, const int* symmIndex, const int* symmMult, const complex* phase, complex* x)
{	complex xSum = 0.;
	for(int j=0; j<n; j++)
		xSum += x[symmIndex[n*i+j]] * phase[n*i+j];
	xSum *= 1./(n*symmMult[i]); //average n in the equivalence class, with weight for accumulation below accounted)
	for(int j=0; j<n; j++)
		x[symmIndex[n*i+j]] = 0.;
	for(int j=0; j<n; j++)
		x[symmIndex[n*i+j]] += xSum * phase[n*i+j].conj();
}

//! Quadruplet of complex arrays corresponding to spin density matrix channels
class complexPtr4
{	complex *up, *dn, *re, *im; //!< UpUp, DnDn, Re(UpDn), Im(UpDn) components respectively
public:
	complexPtr4(const std::vector<complex*>& v) : up(v[0]), dn(v[1]), re(v[2]), im(v[3]) {}
	__hostanddev__ void get(int i, const complex& alpha, complex& s, vector3<complex>& v) const //get scalar and vector parts scaled by alpha
	{	s = alpha*(up[i]+dn[i]);
		v = alpha*vector3<complex>(2.*re[i], -2.*im[i], up[i]-dn[i]);
	}
	__hostanddev__ void accum(int i, const complex& alpha, const complex& s, const vector3<complex>& v) //set scalar and vector parts scaled by alpha
	{	complex alphaHlf = 0.5*alpha;
		up[i] += alphaHlf*(s+v[2]);
		dn[i] += alphaHlf*(s-v[2]);
		re[i] += alphaHlf*v[0];
		im[i] -= alphaHlf*v[1];
	}
	__hostanddev__ void zero(int i) //set to zero
	{	up[i] = 0.;
		dn[i] = 0.;
		re[i] = 0.;
		im[i] = 0.;
	}
};

__hostanddev__ void eblas_symmetrize_phase_rot_calc(size_t i, int n, const int* symmIndex, const int* symmMult, const complex* phase, const matrix3<>* rotSpin, complexPtr4 x)
{	//Gather sums:
	complex sSum = 0.;
	vector3<complex> vSum;
	for(int j=0; j<n; j++)
	{	complex s; vector3<complex> v;
		x.get(symmIndex[n*i+j], phase[n*i+j], s, v);
		sSum += s;
		vSum += rotSpin[j] * v;
	}
	//Normalize to average:
	double scaleFac = 1./(n*symmMult[i]); //account for multiplicity in accumulation below
	sSum *= scaleFac;
	vSum *= scaleFac;
	//Zero target:
	for(int j=0; j<n; j++)
		x.zero(symmIndex[n*i+j]);
	//Scatter averages:
	for(int j=0; j<n; j++)
		x.accum(symmIndex[n*i+j], phase[n*i+j].conj(), sSum, vSum * rotSpin[j]); //rotation conjugate to above
}

//! @endcond

#endif // JDFTX_CORE_BLASEXTRA_INTERNAL_H
