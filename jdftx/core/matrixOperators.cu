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

__global__
void matrixSubSet_kernel(int nr, int iStart, int iStep, int iDelta, int jStart, int jStep, int jDelta, const complex* in, complex* out)
{	int i = kernelIndex(x); if(i>=iDelta) return;
	int j = kernelIndex(y); if(j>=jDelta) return;
	out[nr*(j*jStep+jStart) + (i*iStep+iStart)] = in[iDelta*j+i];
	
}
void matrixSubSet_gpu(int nr, int iStart, int iStep, int iDelta, int jStart, int jStep, int jDelta, const complex* in, complex* out)
{	GpuLaunchConfig3D glc(matrixSubSet_kernel, vector3<int>(1,jDelta,iDelta));
	matrixSubSet_kernel<<<glc.nBlocks,glc.nPerBlock>>>(nr, iStart,iStep,iDelta, jStart,jStep,jDelta, in, out);
	gpuErrorCheck();
}

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
