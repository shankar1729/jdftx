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

#include <gsl/gsl_cblas.h>
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

//!@brief Specialization of eblas_mul() for double[] *= double[]
inline void eblas_dmul(const int N, const double* X, const int incX, double* Y, const int incY) { eblas_mul(N,X,incX,Y,incY); }
//!@brief Specialization of eblas_mul() for complex[] *= complex[]
inline void eblas_zmul(const int N, const complex* X, const int incX, complex* Y, const int incY) { eblas_mul(N,X,incX,Y,incY); }
//!@brief Specialization of eblas_mul() for complex[] *= double[]
inline void eblas_zmuld(const int N, const double* X, const int incX, complex* Y, const int incY) { eblas_mul(N,X,incX,Y,incY); }
#ifdef GPU_ENABLED
//GPU versions of the above functions implemented in BlasExtra.cu
//!@brief Equivalent of eblas_dmul() for GPU data pointers
void eblas_dmul_gpu(const int N, const double* X, const int incX, double* Y, const int incY);
//!@brief Equivalent of eblas_zmul() for GPU data pointers
void eblas_zmul_gpu(const int N, const complex* X, const int incX, complex* Y, const int incY);
//!@brief Equivalent of eblas_zmuld() for GPU data pointers
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

//!@brief Specialization of eblas_div() for double[] /= double[]
inline void eblas_ddiv(const int N, const double* X, const int incX, double* Y, const int incY) { eblas_div(N,X,incX,Y,incY); }
//!@brief Specialization of eblas_div() for #complex[] /= #complex[]
inline void eblas_zdiv(const int N, const complex* X, const int incX, complex* Y, const int incY) { eblas_div(N,X,incX,Y,incY); }
//!@brief Specialization of eblas_div() for #complex[] /= double[]
inline void eblas_zdivd(const int N, const double* X, const int incX, complex* Y, const int incY) { eblas_div(N,X,incX,Y,incY); }
#ifdef GPU_ENABLED
//GPU versions of the above functions implemented in BlasExtra.cu
//!@brief Equivalent of eblas_ddiv() for GPU data pointers
void eblas_ddiv_gpu(const int N, const double* X, const int incX, double* Y, const int incY);
//!@brief Equivalent of eblas_zdiv() for GPU data pointers
void eblas_zdiv_gpu(const int N, const complex* X, const int incX, complex* Y, const int incY);
//!@brief Equivalent of eblas_zdivd() for GPU data pointers
void eblas_zdivd_gpu(const int N, const double* X, const int incX, complex* Y, const int incY);
#endif

//! @brief Elementwise linear combination Z = sX * X + sY * Y
//! @param N Number of elements in X, Y and X
//! @param sX Scale factor for input X
//! @param X Data pointer for input X
//! @param incX Pointer increment for input X
//! @param sY Scale factor for input Y
//! @param Y Data pointer for input Y
//! @param incY Pointer increment for input Y
//! @param Z Data pointer for input Z
//! @param incZ Pointer increment for input Z
void eblas_lincomb(const int N,
	const complex& sX, const complex* X, const int incX,
	const complex& sY, const complex* Y, const int incY,
	complex* Z, const int incZ);

#ifdef GPU_ENABLED
//! @brief Equivalent of eblas_lincomb() for GPU data pointers
void eblas_lincomb_gpu(const int N,
	const complex& sX, const complex* X, const int incX,
	const complex& sY, const complex* Y, const int incY,
	complex* Z, const int incZ);
#endif

//! @brief Threaded complex matrix multiply (threaded wrapper around zgemm)
//! All the parameters have the same meaning as in cblas_zgemm, except element order is always Column Major (FORTRAN order!)
void eblas_zgemm(CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, int M, int N, int K,
	const complex& alpha, const complex *A, const int lda, const complex *B, const int ldb,
	const complex& beta, complex *C, const int ldc);
#ifdef GPU_ENABLED
//! @brief Wrap cublasZgemm to provide the same interface as eblas_zgemm()
void eblas_zgemm_gpu(CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, int M, int N, int K,
	const complex& alpha, const complex *A, const int lda, const complex *B, const int ldb,
	const complex& beta, complex *C, const int ldc);
#endif

//Sparse<->dense vector operations:
//! @brief Scatter y(index) += a * x
//! @param Nindex Length of index array
//! @param a Scale factor (real)
//! @param index 0-based index array (with length at least Nindex)
//! @param x Input array that is sampled consecutively (with length at least Nindex)
//! @param y Output array that is sampled with the index array (with length at least max(index)+1)
//! @param conjugate If true, accumulate from conj(x) instead of x
void eblas_scatter_zdaxpy(const int Nindex, double a, const int* index, const complex* x, complex* y, bool conjugate=false);
//! @brief Equivalent of eblas_scatter_zdaxpy() with a complex scale factor
void eblas_scatter_zaxpy(const int Nindex, complex a, const int* index, const complex* x, complex* y, bool conjugate=false);
//! @brief Equivalent of eblas_scatter_zdaxpy() for real data arrays
void eblas_scatter_daxpy(const int Nindex, double a, const int* index, const double* x, double* y);

//! @brief Gather y += a * x(index)
//! @param Nindex Length of index array
//! @param a Scale factor (real)
//! @param index 0-based index array (with length at least Nindex)
//! @param x Input array that is sampled with the index array (with length at least max(index)+1)
//! @param y Output array that is sampled consecutively (with length at least Nindex)
//! @param conjugate If true, accumulate from conj(x) instead of x
void eblas_gather_zdaxpy(const int Nindex, double a, const int* index, const complex* x, complex* y, bool conjugate=false);
//! @brief Equivalent of eblas_gather_zdaxpy() with a complex scale factor
void eblas_gather_zaxpy(const int Nindex, complex a, const int* index, const complex* x, complex* y, bool conjugate=false);
//! @brief Equivalent of eblas_scatter_zdaxpy() for real data arrays
void eblas_gather_daxpy(const int Nindex, double a, const int* index, const double* x, double* y);

#ifdef GPU_ENABLED
//! @brief Equivalent of eblas_scatter_zdaxpy() for GPU data pointers
void eblas_scatter_zdaxpy_gpu(const int Nindex, double a, const int* index, const complex* x, complex* y, bool conjugate=false);
//! @brief Equivalent of eblas_scatter_zaxpy() for GPU data pointers
void eblas_scatter_zaxpy_gpu(const int Nindex, complex a, const int* index, const complex* x, complex* y, bool conjugate=false);
//! @brief Equivalent of eblas_scatter_daxpy() for GPU data pointers
void eblas_scatter_daxpy_gpu(const int Nindex, double a, const int* index, const double* x, double* y);
//! @brief Equivalent of eblas_gather_zdaxpy() for GPU data pointers
void eblas_gather_zdaxpy_gpu(const int Nindex, double a, const int* index, const complex* x, complex* y, bool conjugate=false);
//! @brief Equivalent of eblas_gather_zaxpy() for GPU data pointers
void eblas_gather_zaxpy_gpu(const int Nindex, complex a, const int* index, const complex* x, complex* y, bool conjugate=false);
//! @brief Equivalent of eblas_gather_daxpy() for GPU data pointers
void eblas_gather_daxpy_gpu(const int Nindex, double a, const int* index, const double* x, double* y);
#endif

//! @brief Accumulate elementwise norm of a complex array x into y i.e. y += a x conj(x) 
//! @param N Length of array
//! @param a scale factor
//! @param x Input complex data array
//! @param y Ouput real data array
void eblas_accumNorm(int N, const double& a, const complex* x, double* y);
//! @brief Accumulate elementwise product of two complex arrays xU and xC into real and imaginary parts yRe and yIm i.e. (yRe + i yIm) += a xU conj(xC) 
//! @param N Length of array
//! @param a scale factor
//! @param xU Unconjugated input complex data array
//! @param xC Conjugated input complex data array
//! @param yRe Ouput real-part data array
//! @param yIm Ouput imaginary-part data array
void eblas_accumProd(int N, const double& a, const complex* xU, const complex* xC, double* yRe, double* yIm);
#ifdef GPU_ENABLED
//! @brief Equivalent of eblas_accumNorm() for GPU data pointers
void eblas_accumNorm_gpu(int N, const double& a, const complex* x, double* y);
//! @brief Equivalent of eblas_accumProd() for GPU data pointers
void eblas_accumProd_gpu(int N, const double& a, const complex* xU, const complex* xC, double* yRe, double* yIm);
#endif

//! @brief Symmetrize an array x, using N n-fold equivalence classes in symmIndex
//! @param N Length of array x
//! @param n Length of symmetry equivalence classes
//! @param symmIndex Every consecutive set of n indices in this array forms an equivalence class
//! @param x Data array to be symmetrized in place
void eblas_symmetrize(int N, int n, const int* symmIndex, double* x);
//! @brief Equivalent of eblas_symmetrize() for complex data pointers
void eblas_symmetrize(int N, int n, const int* symmIndex, complex* x);
#ifdef GPU_ENABLED
//! @brief Equivalent of eblas_symmetrize() for real GPU data pointers
void eblas_symmetrize_gpu(int N, int n, const int* symmIndex, double* x);
//! @brief Equivalent of eblas_symmetrize() for complex GPU data pointers
void eblas_symmetrize_gpu(int N, int n, const int* symmIndex, complex* x);
#endif

//Threaded-wrappers for BLAS1 functions (Cblas)
//! @brief Copy a data array
//! @tparam T Data type of the input and output arrays
//! @param dest Destination data pointer
//! @param src Source data pointer
//! @param N Number of elements of the data type T to copy
template<typename T> void eblas_copy(T* dest, const T* src, int N) { memcpy(dest, src, N*sizeof(T)); }
//! @brief Zero a data array
//! @param N Number of elements to zero
//! @param x Data pointer
void eblas_zero(int N, double* x);
//! @brief Equivalent of eblas_zero() for complex data arrays
void eblas_zero(int N, complex* x);
//! @brief Scale a real array: threaded wrapper to the cblas_dscal BLAS1 function
void eblas_dscal(int N, double a, double* x, int incx);
//! @brief Scale a complex array by a real scale factor: threaded wrapper to the cblas_zdscal BLAS1 function
void eblas_zdscal(int N, double a, complex* x, int incx);
//! @brief Scale a complex array by a complex scale factor: threaded wrapper to the cblas_zscal BLAS1 function
void eblas_zscal(int N, const complex& a, complex* x, int incx);
//! @brief Scaled-accumulate on real arrays: threaded wrapper to the cblas_daxpy BLAS1 function
void eblas_daxpy(int N, double a, const double* x, int incx, double* y, int incy);
//! @brief Scaled-accumulate on complex arrays: threaded wrapper to the cblas_zaxpy BLAS1 function
void eblas_zaxpy(int N, const complex& a, const complex* x, int incx, complex* y, int incy);
//! @brief Dot product of complex arrays: threaded wrapper to the cblas_zdotc BLAS1 function
complex eblas_zdotc(int N, const complex* x, int incx, const complex* y, int incy);
//! @brief Dot product of real arrays: threaded wrapper to the cblas_ddot BLAS1 function
double eblas_ddot(int N, const double* x, int incx, const double* y, int ncy);
//! @brief 2-norm of a complex array: threaded wrapper to the cblas_dznrm2 BLAS1 function
double eblas_dznrm2(int N, const complex* x, int incx);
//! @brief 2-norm of a real array: threaded wrapper to the cblas_dnrm2 BLAS1 function
double eblas_dnrm2(int N, const double* x, int incx);

//Wrappers/alternate names for CUBLAS functions (auto-selectable with Cblas ones above using callPref)
#ifdef GPU_ENABLED
//! @brief Equivalent of eblas_copy() for GPU data pointers
template<typename T> void eblas_copy_gpu(T* dest, const T* src, int N) { cudaMemcpy(dest, src, N*sizeof(T), cudaMemcpyDeviceToDevice); }
//! @brief Equivalent of eblas_zero() for GPU data pointers
void eblas_zero_gpu(int N, double* x);
//! @brief Equivalent of eblas_zero() for GPU data pointers
void eblas_zero_gpu(int N, complex* x);
//! @brief Equivalent of eblas_dscal() for GPU data pointers
#define eblas_dscal_gpu cublasDscal
//! @brief Equivalent of eblas_zdscal() for GPU data pointers
void eblas_zdscal_gpu(int N, double a, complex* x, int incx);
//! @brief Equivalent of eblas_zscal for GPU data pointers
void eblas_zscal_gpu(int N, const complex& a, complex* x, int incx);
//! @brief Equivalent of eblas_daxpy() for GPU data pointers
#define eblas_daxpy_gpu cublasDaxpy
//! @brief Equivalent of eblas_zaxpy() for GPU data pointers
void eblas_zaxpy_gpu(int N, const complex& a, const complex* x, int incx, complex* y, int incy);
//! @brief Equivalent of eblas_zdotc() for GPU data pointers
complex eblas_zdotc_gpu(int N, const complex* x, int incx, const complex* y, int incy);
//! @brief Equivalent of eblas_ddot() for GPU data pointers
#define eblas_ddot_gpu cublasDdot
//! @brief Equivalent of eblas_dznrm2() for GPU data pointers
double eblas_dznrm2_gpu(int N, const complex* x, int incx);
//! @brief Equivalent of eblas_dnrm2() for GPU data pointers
#define eblas_dnrm2_gpu cublasDnrm2

#endif

//! @brief Select between functionName and functionName_gpu for the CPU and GPU executables respectively
#ifdef GPU_ENABLED
	#define callPref(functionName) functionName##_gpu //gpu versions available, so use them preferentially
#else
	#define callPref(functionName) functionName //only cpu version available
#endif


//! @brief Find the minimum and maximum of a data array and optionally cap it from above and/or below
//! @param N Length of data array
//! @param x Input data array, that might be modified if capLo or capHi are finite
//! @param xMin On output, the minimum of the input data array
//! @param xMax On output, the maximum of the input data array
//! @param capLo If finite, cap the data array on output from below at this value (no capping for the default value)
//! @param capHi If finite, cap the data array on output from above at this value (no capping for the default value)
void eblas_capMinMax(const int N, double* x, double& xMin, double& xMax, double capLo=-DBL_MAX, double capHi=+DBL_MAX);
#ifdef GPU_ENABLED
//! @brief Equivalent of eblas_capMinMax() for GPU data pointers
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
	threadLaunch((N<100000) ? 1 : 0, //force single threaded for small problem sizes
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
	threadLaunch((N<100000) ? 1 : 0, //force single threaded for small problem sizes
		eblas_div_sub<Ty,Tx>, N, X, incX, Y, incY);
}
//!@endcond

#endif // JDFTX_CORE_BLASEXTRA_H
