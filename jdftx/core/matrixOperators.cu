/*-------------------------------------------------------------------
Copyright 2018 Ravishankar Sundararaman

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
#include <cublas_v2.h>

__global__
void matrixSubGet_kernel(int nr, int iStart, int iStep, int iDelta, int jStart, int jStep, int jDelta, const complex* in, complex* out)
{	int i = kernelIndex(x); if(i>=iDelta) return;
	int j = kernelIndex(y); if(j>=jDelta) return;
	out[iDelta*j+i] = in[nr*(j*jStep+jStart) + (i*iStep+iStart)];
	
}
void matrixSubGet_gpu(int nr, int iStart, int iStep, int iDelta, int jStart, int jStep, int jDelta, const complex* in, complex* out)
{	GpuLaunchConfig3D glc(matrixSubGet_kernel, vector3<int>(1,jDelta,iDelta));
	matrixSubGet_kernel<<<glc.nBlocks,glc.nPerBlock>>>(nr, iStart,iStep,iDelta, jStart,jStep,jDelta, in, out);
	gpuErrorCheck();
}

#define DECLARE_matrixSubSetAccum(NAME, OP) \
	__global__ \
	void matrixSub ##NAME## _kernel(int nr, int iStart, int iStep, int iDelta, int jStart, int jStep, int jDelta, const complex* in, complex* out) \
	{	int i = kernelIndex(x); if(i>=iDelta) return; \
		int j = kernelIndex(y); if(j>=jDelta) return; \
		out[nr*(j*jStep+jStart) + (i*iStep+iStart)] OP in[iDelta*j+i]; \
	} \
	void matrixSub ##NAME## _gpu(int nr, int iStart, int iStep, int iDelta, int jStart, int jStep, int jDelta, const complex* in, complex* out) \
	{	GpuLaunchConfig3D glc(matrixSub ##NAME## _kernel, vector3<int>(1,jDelta,iDelta)); \
		matrixSub ##NAME## _kernel<<<glc.nBlocks,glc.nPerBlock>>>(nr, iStart,iStep,iDelta, jStart,jStep,jDelta, in, out); \
		gpuErrorCheck(); \
	}
DECLARE_matrixSubSetAccum(Set, =)
DECLARE_matrixSubSetAccum(Accum, +=)
#undef DECLARE_matrixSubSetAccum


//----- diagonal matrix multiplies
template<typename scalar> __global__ void mulMD_kernel(int nRows, int nCols, const complex* M, const scalar* D, complex* out)
{	int i = kernelIndex(x); if(i>=nRows) return;
	int j = kernelIndex(y); if(j>=nCols) return;
	int index = j*nRows + i;
	out[index] = M[index] * D[j];
}
template<typename scalar> __global__ void mulDM_kernel(int nRows, int nCols, const scalar* D, const complex* M, complex* out)
{	int i = kernelIndex(x); if(i>=nRows) return;
	int j = kernelIndex(y); if(j>=nCols) return;
	int index = j*nRows + i;
	out[index] = D[i] * M[index];
}
template<typename scalar> void mulMD_gpu(int nRows, int nCols, const complex* M, const scalar* D, complex* out)
{	GpuLaunchConfig3D glc(mulMD_kernel<scalar>, vector3<int>(1,nCols,nRows));
	mulMD_kernel<scalar><<<glc.nBlocks,glc.nPerBlock>>>(nRows, nCols, M, D, out);
	gpuErrorCheck();
}
template<typename scalar> void mulDM_gpu(int nRows, int nCols, const scalar* D, const complex* M, complex* out)
{	GpuLaunchConfig3D glc(mulDM_kernel<scalar>, vector3<int>(1,nCols,nRows));
	mulDM_kernel<scalar><<<glc.nBlocks,glc.nPerBlock>>>(nRows, nCols, D, M, out);
	gpuErrorCheck();
}
void mulMDdouble_gpu(int nRows, int nCols, const complex* M, const double* D, complex* out) { mulMD_gpu<double>(nRows, nCols, M, D, out); }
void mulDMdouble_gpu(int nRows, int nCols, const double* D, const complex* M, complex* out) { mulDM_gpu<double>(nRows, nCols, D, M, out); }
void mulMDcomplex_gpu(int nRows, int nCols, const complex* M, const complex* D, complex* out) { mulMD_gpu<complex>(nRows, nCols, M, D, out); }
void mulDMcomplex_gpu(int nRows, int nCols, const complex* D, const complex* M, complex* out) { mulDM_gpu<complex>(nRows, nCols, D, M, out); }


__global__
void relativeHermiticityError_kernel(int N, const complex* data, double* buf)
{	int i = kernelIndex1D();
	if(i<N)
	{	double errNum = 0., errDen = 0.;
		for(int j=0; j<N; j++)
		{	int index = N*i + j;
			int indexT = N*j + i;
			errNum += norm(data[index]-data[indexT].conj());
			errDen += norm(data[index]);
		}
		buf[i] = errNum;
		buf[i+N] = errDen;
	}
}
double relativeHermiticityError_gpu(int N, const complex* data)
{	GpuLaunchConfig1D glc(relativeHermiticityError_kernel, N);
	double* buf; cudaMalloc(&buf, sizeof(double)*(2*N)); //buffer to store results per row
	relativeHermiticityError_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, data, buf);
	gpuErrorCheck();
	double errNum = 0., errDen = 0.;
	cublasDasum(cublasHandle, N, buf, 1, &errNum);
	cublasDasum(cublasHandle, N, buf+N, 1, &errDen);
	cudaFree(buf);
	return sqrt(errNum / (errDen*N));
}
