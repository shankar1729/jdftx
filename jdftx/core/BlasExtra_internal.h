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

#include <core/scalar.h>

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

//! @endcond

#endif // JDFTX_CORE_BLASEXTRA_INTERNAL_H
