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

#ifndef JDFTX_CORE_BLASEXTRA_H
#define JDFTX_CORE_BLASEXTRA_H

/** @file BlasExtra.h
@brief Commonly used BLAS-like routines
*/

#include <core/cblas_wrapper.h>
#include <cstdlib>
#include <cstdio>
#include <cfloat>
#include <core/scalar.h>
#include <core/Thread.h>

#ifdef GPU_ENABLED
#include <cublas.h>
#include <cuda_runtime.h>
#endif

/** @brief Templated elementwise multiply Y *= X for arrays X, Y
@tparam Ty primiitive data-type for array Y
@tparam Tx primiitive data-type for array X
@param N number of elements in X and Y
@param X pointer to the first element of array X
@param incX stride along the X array
@param Y pointer to the first element of array Y
@param incY stride along the Y array (must be non-zero)
*/
template<typename Ty, typename Tx> void eblas_mul(const int N, const Tx* X, const int incX, Ty* Y, const int incY);

//!@brief Specialization of #eblas_mul for double[] *= double[]
inline void eblas_dmul(const int N, const double* X, const int incX, double* Y, const int incY) { eblas_mul(N,X,incX,Y,incY); }
//!@brief Specialization of #eblas_mul for complex[] *= complex[]
inline void eblas_zmul(const int N, const complex* X, const int incX, complex* Y, const int incY) { eblas_mul(N,X,incX,Y,incY); }
//!@brief Specialization of #eblas_mul for complex[] *= double[]
inline void eblas_zmuld(const int N, const double* X, const int incX, complex* Y, const int incY) { eblas_mul(N,X,incX,Y,incY); }
#ifdef GPU_ENABLED
//GPU versions of the above functions implemented in BlasExtra.cu
//!@brief Equivalent of #eblas_dmul for GPU data pointers
void eblas_dmul_gpu(const int N, const double* X, const int incX, double* Y, const int incY);
//!@brief Equivalent of #eblas_zmul for GPU data pointers
void eblas_zmul_gpu(const int N, const complex* X, const int incX, complex* Y, const int incY);
//!@brief Equivalent of #eblas_zmuld for GPU data pointers
void eblas_zmuld_gpu(const int N, const double* X, const int incX, complex* Y, const int incY);
#endif

/** @brief Templated elementwise divide Y /= X for arrays X, Y
@tparam Ty primiitive data-type for array Y
@tparam Tx primiitive data-type for array X
@param N number of elements in X and Y
@param X pointer to the first element of array X
@param incX stride along the X array
@param Y pointer to the first element of array Y
@param incY stride along the Y array (must be non-zero)
*/
template<typename Ty, typename Tx> void eblas_div(const int N, const Tx* X, const int incX, Ty* Y, const int incY);

//!@brief Specialization of #eblas_div for double[] /= double[]
inline void eblas_ddiv(const int N, const double* X, const int incX, double* Y, const int incY) { eblas_div(N,X,incX,Y,incY); }
//!@brief Specialization of #eblas_div for #complex[] /= #complex[]
inline void eblas_zdiv(const int N, const complex* X, const int incX, complex* Y, const int incY) { eblas_div(N,X,incX,Y,incY); }
//!@brief Specialization of #eblas_div for #complex[] /= double[]
inline void eblas_zdivd(const int N, const double* X, const int incX, complex* Y, const int incY) { eblas_div(N,X,incX,Y,incY); }
#ifdef GPU_ENABLED
//GPU versions of the above functions implemented in BlasExtra.cu
//!@brief Equivalent of #eblas_ddiv for GPU data pointers
void eblas_ddiv_gpu(const int N, const double* X, const int incX, double* Y, const int incY);
//!@brief Equivalent of #eblas_zdiv for GPU data pointers
void eblas_zdiv_gpu(const int N, const complex* X, const int incX, complex* Y, const int incY);
//!@brief Equivalent of #eblas_zdivd for GPU data pointers
void eblas_zdivd_gpu(const int N, const double* X, const int incX, complex* Y, const int incY);
#endif

//! Elementwise linear combination Z = sX * X + sY * Y
void eblas_lincomb(const int N,
	const complex& sX, const complex* X, const int incX,
	const complex& sY, const complex* Y, const int incY,
	complex* Z, const int incZ);

#ifdef GPU_ENABLED
//! Elementwise linear combination Z = sX * X + sY * Y (gpu version)
void eblas_lincomb_gpu(const int N,
	const complex& sX, const complex* X, const int incX,
	const complex& sY, const complex* Y, const int incY,
	complex* Z, const int incZ);
#endif

/** @brief Threaded complex matrix multiply (threaded wrapper around zgemm)
All the parameters have the same meaning as in zgemm, except element order is always Column Major (FORTRAN order!)
*/
void eblas_zgemm(CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, int M, int N, int K,
	const complex& alpha, const complex *A, const int lda, const complex *B, const int ldb,
	const complex& beta, complex *C, const int ldc);
#ifdef GPU_ENABLED
//! cublasZgemm wrapper to overload eblas_zgemm
void eblas_zgemm_gpu(CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, int M, int N, int K,
	const complex& alpha, const complex *A, const int lda, const complex *B, const int ldb,
	const complex& beta, complex *C, const int ldc);
#endif

//Sparse-dense vector operations:
//! Perform a scatter axpy operation: y(index) += a*x (in Octave notation)
void eblas_scatter_zdaxpy(const int Nindex, double a, const int* index, const complex* x, complex* y);
//! Perform a gather axpy operation: y += a*x(index) (in Octave notation)
void eblas_gather_zdaxpy(const int Nindex, double a, const int* index, const complex* x, complex* y);
#ifdef GPU_ENABLED
//! Perform a GPU scatter axpy operation: y(index) += a*x (in Octave notation)
void eblas_scatter_zdaxpy_gpu(const int Nindex, double a, const int* index, const complex* x, complex* y);
//! Perform a GPU gather axpy operation: y += a*x(index) (in Octave notation)
void eblas_gather_zdaxpy_gpu(const int Nindex, double a, const int* index, const complex* x, complex* y);
#endif


//! Elementwise y += a|x|^2
void eblas_accumNorm(int N, const double& a, const complex* x, double* y);
#ifdef GPU_ENABLED
//! Elementwise y += a|x|^2 on the GPU
void eblas_accumNorm_gpu(int N, const double& a, const complex* x, double* y);
#endif


//Threaded-wrappers for BLAS1 functions (Cblas)
template<typename T> void eblas_copy(T* dest, const T* src, int N) { memcpy(dest, src, N*sizeof(T)); }
void eblas_zero(int N, double* x); //!< wraps memset
void eblas_zero(int N, complex* x); //!< wraps memset
void eblas_dscal(int N, double a, double* x, int incx); //!< wraps cblas_dscal
void eblas_zdscal(int N, double a, complex* x, int incx); //!< wraps cblas_zdscal
void eblas_zscal(int N, const complex& a, complex* x, int incx); //!< wraps cblas_zscal
void eblas_daxpy(int N, double a, const double* x, int incx, double* y, int incy); //!< wraps cblas_daxpy
void eblas_zaxpy(int N, const complex& a, const complex* x, int incx, complex* y, int incy); //!< wraps cblas_zaxpy
complex eblas_zdotc(int N, const complex* x, int incx, const complex* y, int incy); //!< wraps cblas_zdotc_sub
double eblas_ddot(int N, const double* x, int incx, const double* y, int incy); //!< wraps cblas_ddot
double eblas_dznrm2(int N, const complex* x, int incx); //!< wraps cblas_dznrm2
double eblas_dnrm2(int N, const double* x, int incx); //!< wraps cblas_dnrm2

//Wrappers/alternate names for CUBLAS functions (auto-selectable with Cblas ones above using callPref)
#ifdef GPU_ENABLED
template<typename T> void eblas_copy_gpu(T* dest, const T* src, int N) { cudaMemcpy(dest, src, N*sizeof(T), cudaMemcpyDeviceToDevice); }
void eblas_zero_gpu(int N, double* x); //!< wraps cudaMemset
void eblas_zero_gpu(int N, complex* x); //!< wraps cudaMemset
#define eblas_dscal_gpu cublasDscal
void eblas_zdscal_gpu(int N, double a, complex* x, int incx); //!< wraps cublasZdscal
void eblas_zscal_gpu(int N, const complex& a, complex* x, int incx); //!< wraps cublasZscal
#define eblas_daxpy_gpu cublasDaxpy
void eblas_zaxpy_gpu(int N, const complex& a, const complex* x, int incx, complex* y, int incy); //!< wraps cublasZaxpy
complex eblas_zdotc_gpu(int N, const complex* x, int incx, const complex* y, int incy); //!< wraps cublasZdotc
#define eblas_ddot_gpu cublasDdot
double eblas_dznrm2_gpu(int N, const complex* x, int incx); //!< wraps cublasDznrm2
#define eblas_dnrm2_gpu cublasDnrm2

#endif

#ifdef GPU_ENABLED
	#define callPref(functionName) functionName##_gpu //gpu versions available, so use them preferentially
#else
	#define callPref(functionName) functionName //only cpu version available
#endif


//! Find the minimum, xMin, and maximum, xMax, of an array x, and optionally restrict its range to [capLo,capHi]
void eblas_capMinMax(const int N, double* x, double& xMin, double& xMax, double capLo=-DBL_MAX, double capHi=+DBL_MAX);
#ifdef GPU_ENABLED
//! Equivalent of eblas_capMinMax for the gpu
void eblas_capMinMax_gpu(const int N, double* x, double& xMin, double& xMax, double capLo=-DBL_MAX, double capHi=+DBL_MAX);
#endif


//###################################################################################################
//####  Implementation  ####
//##########################

//!@cond

//Elementwise multiply implementation:
template<typename Ty, typename Tx>
void eblas_mul_sub(size_t iMin, size_t iMax, const Tx* X, const int incX, Ty* Y, const int incY)
{	for(size_t i=iMin; i<iMax; i++) Y[incY*i] *= X[incX*i];
}

template<typename Ty, typename Tx>
void eblas_mul(const int N, const Tx* X, const int incX, Ty* Y, const int incY)
{	if(incY==0) die("incY cannot be = 0")
	threadLaunch((N<100000 || (!threadOperators)) ? 1 : 0, //force single threaded for small problem sizes
		eblas_mul_sub<Ty,Tx>, N, X, incX, Y, incY);
}

//Elementwise divide implementation:
template<typename Ty, typename Tx>
void eblas_div_sub(size_t iMin, size_t iMax, const Tx* X, const int incX, Ty* Y, const int incY)
{	for(size_t i=iMin; i<iMax; i++) Y[incY*i] /= X[incX*i];
}
template<typename Ty, typename Tx>
void eblas_div(const int N, const Tx* X, const int incX, Ty* Y, const int incY)
{	if(incY==0) die("incY cannot be = 0")
	threadLaunch((N<100000 || (!threadOperators)) ? 1 : 0, //force single threaded for small problem sizes
		eblas_div_sub<Ty,Tx>, N, X, incX, Y, incY);
}
//!@endcond

#endif // JDFTX_CORE_BLASEXTRA_H
