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
#include <electronic/JDFT1_shapeFunc.h>

__global__
void JDFT1_shapeFunc_kernel(int N, const double* nCavity, double* shape, const double nc, const double sigma)
{	int i = kernelIndex1D(); if(i<N) JDFT1_shapeFunc_sub(i, nCavity, shape, nc, sigma);
}
void JDFT1_shapeFunc_gpu(int N, const double* nCavity, double* shape, const double nc, const double sigma)
{	GpuLaunchConfig1D glc(JDFT1_shapeFunc_kernel, N);
	JDFT1_shapeFunc_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, nCavity, shape, nc, sigma);
	gpuErrorCheck();
}

__global__
void JDFT1_shapeFunc_grad_kernel(int N, const double* nCavity, const double* grad_shape, double* grad_nCavity, const double nc, const double sigma)
{	int i = kernelIndex1D(); if(i<N) JDFT1_shapeFunc_grad_sub(i, nCavity, grad_shape, grad_nCavity, nc, sigma);
}
void JDFT1_shapeFunc_grad_gpu(int N, const double* nCavity, const double* grad_shape, double* grad_nCavity, const double nc, const double sigma)
{	GpuLaunchConfig1D glc(JDFT1_shapeFunc_grad_kernel, N);
	JDFT1_shapeFunc_grad_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, nCavity, grad_shape, grad_nCavity, nc, sigma);
	gpuErrorCheck();
}

