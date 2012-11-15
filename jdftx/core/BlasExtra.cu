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
#include <core/scalar.h>
#include <algorithm>
#include <cublas.h>
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

inline char cublasTranspose(CBLAS_TRANSPOSE trans)
{	switch(trans)
	{	case CblasNoTrans: return 'N';
		case CblasTrans: return 'T';
		case CblasConjTrans: return 'C';
	}
	return 0;
}
void eblas_zgemm_gpu(CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, int M, int N, int K,
	const complex& alpha, const complex *A, const int lda, const complex *B, const int ldb,
	const complex& beta, complex *C, const int ldc)
{	cublasZgemm(cublasTranspose(TransA), cublasTranspose(TransB), M, N, K,
		alpha, (const double2*)A, lda, (const double2*)B, ldb, beta, (double2*)C, ldc);
}


__global__
void eblas_scatter_zdaxpy_kernel(const int N, double a, const int* index, const complex* x, complex* y)
{	int i = kernelIndex1D();
	if(i<N) y[index[i]] += x[i] * a;
}
void eblas_scatter_zdaxpy_gpu(const int N, double a, const int* index, const complex* x, complex* y)
{	GpuLaunchConfig1D glc(eblas_scatter_zdaxpy_kernel, N);
	eblas_scatter_zdaxpy_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, a, index, x, y);
	gpuErrorCheck();
}

__global__
void eblas_gather_zdaxpy_kernel(const int N, double a, const int* index, const complex* x, complex* y)
{	int i = kernelIndex1D();
	if(i<N) y[i] += x[index[i]] * a;
}
void eblas_gather_zdaxpy_gpu(const int N, double a, const int* index, const complex* x, complex* y)
{	GpuLaunchConfig1D glc(eblas_gather_zdaxpy_kernel, N);
	eblas_gather_zdaxpy_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, a, index, x, y);
	gpuErrorCheck();
}

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

//BLAS-1 wrappers:
void eblas_zero_gpu(int N, complex* x)
{	cudaMemset(x, 0, N*sizeof(complex));
}
void eblas_zero_gpu(int N, double* x)
{	cudaMemset(x, 0, N*sizeof(double));
}
void eblas_zdscal_gpu(int N, double a, complex* x, int incx)
{	cublasZdscal(N, a, (double2*)x, incx);
}
void eblas_zscal_gpu(int N, const complex& a, complex* x, int incx)
{	cublasZscal(N, a, (double2*)x, incx);
}
void eblas_zaxpy_gpu(int N, const complex& a, const complex* x, int incx, complex* y, int incy)
{	cublasZaxpy(N, a, (const double2*)x, incx, (double2*)y, incy);
}
complex eblas_zdotc_gpu(int N, const complex* x, int incx, const complex* y, int incy)
{	return cublasZdotc(N, (const double2*)x, incx, (const double2*)y, incy);
}
double eblas_dznrm2_gpu(int N, const complex* x, int incx)
{	return cublasDznrm2(N, (const double2*)x, incx);
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
