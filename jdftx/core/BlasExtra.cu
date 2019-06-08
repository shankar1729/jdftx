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

#include <core/GpuKernelUtils.h>
#include <core/BlasExtra_internal.h>
#include <algorithm>
#include <cublas_v2.h>
#include <cfloat>
#include <gsl/gsl_cblas.h>

template<typename Tx, typename Ty> __global__
void eblas_mul_kernel(const int N, const Tx* X, const int incX, Ty* Y, const int incY)
{	int i = kernelIndex1D();
	if(i<N) Y[i*incY] *= X[i*incX];
}
void eblas_dmul_gpu(const int N, const double* X, const int incX, double* Y, const int incY)
{	GpuLaunchConfig1D glc(eblas_mul_kernel<double,double>, N);
	eblas_mul_kernel<double,double><<<glc.nBlocks,glc.nPerBlock>>>(N, X,incX, Y,incY);
	gpuErrorCheck();
}
void eblas_zmul_gpu(const int N, const complex* X, const int incX, complex* Y, const int incY)
{	GpuLaunchConfig1D glc(eblas_mul_kernel<complex,complex>, N);
	eblas_mul_kernel<complex,complex><<<glc.nBlocks,glc.nPerBlock>>>(N, X,incX, Y,incY);
	gpuErrorCheck();
}
void eblas_zmuld_gpu(const int N, const double* X, const int incX, complex* Y, const int incY)
{	GpuLaunchConfig1D glc(eblas_mul_kernel<double,complex>, N);
	eblas_mul_kernel<double,complex><<<glc.nBlocks,glc.nPerBlock>>>(N, X,incX, Y,incY);
	gpuErrorCheck();
}


__global__
void eblas_lincomb_kernel(const int N,
	complex sX, const complex* X, const int incX,
	complex sY, const complex* Y, const int incY,
	complex* Z, const int incZ)
{	int i = kernelIndex1D();
	if(i<N) Z[i*incZ] = sX*X[i*incX] + sY*Y[i*incY];
}
void eblas_lincomb_gpu(const int N,
	const complex& sX, const complex* X, const int incX,
	const complex& sY, const complex* Y, const int incY,
	complex* Z, const int incZ)
{	GpuLaunchConfig1D glc(eblas_lincomb_kernel, N);
	eblas_lincomb_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, sX,X,incX, sY,Y,incY, Z,incZ);
	gpuErrorCheck();
}

inline cublasOperation_t cublasTranspose(CBLAS_TRANSPOSE trans)
{	switch(trans)
	{	case CblasNoTrans: return CUBLAS_OP_N;
		case CblasTrans: return CUBLAS_OP_T;
		case CblasConjTrans: return CUBLAS_OP_C;
	}
	return CUBLAS_OP_N;
}
void eblas_zgemm_gpu(CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, int M, int N, int K,
	const complex& alpha, const complex *A, const int lda, const complex *B, const int ldb,
	const complex& beta, complex *C, const int ldc)
{	cublasZgemm(cublasHandle, cublasTranspose(TransA), cublasTranspose(TransB), M, N, K,
		(const double2*)&alpha, (const double2*)A, lda, (const double2*)B, ldb,
		(const double2*)&beta, (double2*)C, ldc);
}

template<typename scalar, typename scalar2, typename Conjugator> __global__ 
void eblas_scatter_axpy_kernel(const int N, scalar2 a, const int* index, const scalar* x, scalar* y, const scalar* w, const Conjugator& conjugator)
{	int i = kernelIndex1D();
	if(i<N) y[index[i]] += a * conjugator(x,i, w,i);
}
template<typename scalar, typename scalar2, typename Conjugator> __global__
void eblas_gather_axpy_kernel(const int N, scalar2 a, const int* index, const scalar* x, scalar* y, const scalar* w, const Conjugator& conjugator)
{	int i = kernelIndex1D();
	if(i<N) y[i] += a * conjugator(x,index[i], w,i);
}
#define DEFINE_SPARSE_AXPY_GPU_LAUNCHER(type) \
	template<typename scalar, typename scalar2, typename Conjugator> \
	void eblas_##type##_axpy_gpu(const int N, scalar2 a, const int* index, const scalar* x, scalar* y, const scalar* w, const Conjugator& conjugator) \
	{	GpuLaunchConfig1D glc(eblas_##type##_axpy_kernel<scalar,scalar2,Conjugator>, N); \
		eblas_##type##_axpy_kernel<scalar,scalar,Conjugator><<<glc.nBlocks,glc.nPerBlock>>>(N, a, index, x, y, w, conjugator); \
		gpuErrorCheck(); \
	}
DEFINE_SPARSE_AXPY_GPU_LAUNCHER(scatter)
DEFINE_SPARSE_AXPY_GPU_LAUNCHER(gather)
DEFINE_SPARSE_AXPY(scatter, _gpu)
DEFINE_SPARSE_AXPY(gather, _gpu)


__global__
void eblas_accumNorm_kernel(int N, double a, const complex* x, double* y)
{	int i = kernelIndex1D();
	if(i<N) y[i] += a * norm(x[i]);
}
void eblas_accumNorm_gpu(int N, const double& a, const complex* x, double* y)
{	GpuLaunchConfig1D glc(eblas_accumNorm_kernel, N);
	eblas_accumNorm_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, a, x, y);
	gpuErrorCheck();
}

__global__
void eblas_accumProd_kernel(int N, double a, const complex* xU, const complex* xC, double* yRe, double* yIm)
{	int i = kernelIndex1D();
	if(i<N)
	{	complex z = a * xU[i] * xC[i].conj();
		yRe[i] += z.real();
		yIm[i] += z.imag();
	}
}
void eblas_accumProd_gpu(int N, const double& a, const complex* xU, const complex* xC, double* yRe, double* yIm)
{	GpuLaunchConfig1D glc(eblas_accumProd_kernel, N);
	eblas_accumProd_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, a, xU, xC, yRe, yIm);
	gpuErrorCheck();
}

template<typename scalar> __global__
void eblas_symmetrize_kernel(int N, int n, const int* symmIndex, scalar* x, double nInv)
{	int i=kernelIndex1D();
	if(i<N) eblas_symmetrize_calc(i, n, symmIndex, x, nInv);
}
template<typename scalar> void eblas_symmetrize_gpu(int N, int n, const int* symmIndex, scalar* x)
{	GpuLaunchConfig1D glc(eblas_symmetrize_kernel<scalar>, N);
	eblas_symmetrize_kernel<scalar><<<glc.nBlocks,glc.nPerBlock>>>(N, n, symmIndex, x, 1./n);
	gpuErrorCheck();
}
void eblas_symmetrize_gpu(int N, int n, const int* symmIndex, double* x) { eblas_symmetrize_gpu<double>(N, n, symmIndex, x); }
void eblas_symmetrize_gpu(int N, int n, const int* symmIndex, complex* x) { eblas_symmetrize_gpu<complex>(N, n, symmIndex, x); }

__global__
void eblas_symmetrize_phase_kernel(int N, int n, const int* symmIndex, const int* symmMult, const complex* phase, complex* x)
{	int i=kernelIndex1D();
	if(i<N) eblas_symmetrize_phase_calc(i, n, symmIndex, symmMult, phase, x);
}
void eblas_symmetrize_gpu(int N, int n, const int* symmIndex, const int* symmMult, const complex* phase, complex* x)
{	GpuLaunchConfig1D glc(eblas_symmetrize_phase_kernel, N);
	eblas_symmetrize_phase_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, n, symmIndex, symmMult, phase, x);
	gpuErrorCheck();
}

__global__
void eblas_symmetrize_phase_rot_kernel(int N, int n, const int* symmIndex, const int* symmMult, const complex* phase, const matrix3<>* rotSpin, complexPtr4 x)
{	int i=kernelIndex1D();
	if(i<N) eblas_symmetrize_phase_rot_calc(i, n, symmIndex, symmMult, phase, rotSpin, x);
}
void eblas_symmetrize_gpu(int N, int n, const int* symmIndex, const int* symmMult, const complex* phase, const matrix3<>* rotSpin, std::vector<complex*> x)
{	GpuLaunchConfig1D glc(eblas_symmetrize_phase_rot_kernel, N);
	eblas_symmetrize_phase_rot_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, n, symmIndex, symmMult, phase, rotSpin, complexPtr4(x));
	gpuErrorCheck();
}

//BLAS-1 wrappers:
void eblas_dscal_gpu(int N, double a, double* x, int incx)
{	cublasDscal_v2(cublasHandle, N, &a, x, incx);
}
void eblas_zdscal_gpu(int N, double a, complex* x, int incx)
{	cublasZdscal_v2(cublasHandle, N, &a, (double2*)x, incx);
}
void eblas_zscal_gpu(int N, const complex& a, complex* x, int incx)
{	cublasZscal_v2(cublasHandle, N, (const double2*)&a, (double2*)x, incx);
}
void eblas_daxpy_gpu(int N, double a, const double* x, int incx, double* y, int incy)
{	cublasDaxpy_v2(cublasHandle, N, &a, x, incx, y, incy);
}
void eblas_zaxpy_gpu(int N, const complex& a, const complex* x, int incx, complex* y, int incy)
{	cublasZaxpy_v2(cublasHandle, N, (const double2*)&a, (const double2*)x, incx, (double2*)y, incy);
}
complex eblas_zdotc_gpu(int N, const complex* x, int incx, const complex* y, int incy)
{	complex result;
	cublasZdotc_v2(cublasHandle, N, (const double2*)x, incx, (const double2*)y, incy, (double2*)&result);
	return result;
}
double eblas_ddot_gpu(int N, const double* x, int incx, const double* y, int incy)
{	double result;
	cublasDdot_v2(cublasHandle, N, x, incx, y, incy, &result);
	return result;
}
double eblas_dznrm2_gpu(int N, const complex* x, int incx)
{	double result;
	cublasDznrm2_v2(cublasHandle, N, (const double2*)x, incx, &result);
	return result;
}
double eblas_dnrm2_gpu(int N, const double* x, int incx)
{	double result;
	cublasDnrm2_v2(cublasHandle, N, x, incx, &result);
	return result;
}


//Min-max:

//forward declare the cpu version (used at the end for colletcing results):
void eblas_capMinMax(const int N, double* x, double& xMin, double& xMax, double capLo=-DBL_MAX, double capHi=+DBL_MAX);

extern __shared__ double xMinLoc[];
__global__
void eblas_capMinMax_kernel(int N, double* x, double* xMinBlk, double* xMaxBlk, double capLo, double capHi)
{	int i=kernelIndex1D();
	int iThread = threadIdx.x;
	//Store the original value as the local min and max:
	double* xMaxLoc = xMinLoc + blockDim.x;
	if(i<N)
	{	xMinLoc[iThread] = x[i];
		xMaxLoc[iThread] = x[i];
		//Cap the value in the array:
		if(x[i]<capLo) x[i]=capLo;
		if(x[i]>capHi) x[i]=capHi;
	}
	else
	{	xMinLoc[iThread] = +DBL_MAX;
		xMaxLoc[iThread] = -DBL_MAX;
	}
	//Min-max within block:
	int extent = blockDim.x/2;
	int stride = (blockDim.x+1)/2;
	while(extent)
	{	__syncthreads();
		if(iThread<extent)
		{	if(xMinLoc[iThread+stride]<xMinLoc[iThread]) xMinLoc[iThread]=xMinLoc[iThread+stride];
			if(xMaxLoc[iThread+stride]>xMaxLoc[iThread]) xMaxLoc[iThread]=xMaxLoc[iThread+stride];
		}
		extent = stride/2;
		stride = (stride+1)/2;
	}
	__syncthreads();
	//Save to the global min-max:
	if(iThread==0)
	{	int iBlock = blockIdx.y*gridDim.x+blockIdx.x;
		xMinBlk[iBlock] = xMinLoc[0];
		xMaxBlk[iBlock] = xMaxLoc[0];
	}
}
void eblas_capMinMax_gpu(const int N, double* x, double& xMin, double& xMax, double capLo, double capHi)
{	GpuLaunchConfig1D glc(eblas_capMinMax_kernel, N);
	int nBlocksTot = glc.nBlocks.x * glc.nBlocks.y;
	//First perform the capping and obtain the max and min per block:
	double* xMinBlk; cudaMalloc(&xMinBlk, sizeof(double)*nBlocksTot*2);
	double* xMaxBlk = xMinBlk + nBlocksTot;
	int sharedMemBytes = 2*glc.nPerBlock.x*sizeof(double);
	eblas_capMinMax_kernel<<<glc.nBlocks,glc.nPerBlock,sharedMemBytes>>>(N,x,xMinBlk,xMaxBlk,capLo,capHi);
	//Finish on the CPU:
	double* xMinCpu = new double[2*nBlocksTot];
	double* xMaxCpu = xMinCpu + nBlocksTot;
	cudaMemcpy(xMinCpu, xMinBlk, sizeof(double)*nBlocksTot*2, cudaMemcpyDeviceToHost);
	cudaFree(xMinBlk);
	xMin = *std::min_element(xMinCpu, xMinCpu+nBlocksTot);
	xMax = *std::max_element(xMaxCpu, xMaxCpu+nBlocksTot);
	delete[] xMinCpu;
}
