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
#include <electronic/PCM_internal.h>

__global__
void pcmShapeFunc_kernel(int N, const double* nCavity, double* shape, const double nc, const double sigma)
{	int i = kernelIndex1D(); if(i<N) pcmShapeFunc_calc(i, nCavity, shape, nc, sigma);
}
void pcmShapeFunc_gpu(int N, const double* nCavity, double* shape, const double nc, const double sigma)
{	GpuLaunchConfig1D glc(pcmShapeFunc_kernel, N);
	pcmShapeFunc_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, nCavity, shape, nc, sigma);
	gpuErrorCheck();
}

__global__
void pcmShapeFunc_grad_kernel(int N, const double* nCavity, const double* grad_shape, double* grad_nCavity, const double nc, const double sigma)
{	int i = kernelIndex1D(); if(i<N) pcmShapeFunc_grad_calc(i, nCavity, grad_shape, grad_nCavity, nc, sigma);
}
void pcmShapeFunc_grad_gpu(int N, const double* nCavity, const double* grad_shape, double* grad_nCavity, const double nc, const double sigma)
{	GpuLaunchConfig1D glc(pcmShapeFunc_grad_kernel, N);
	pcmShapeFunc_grad_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, nCavity, grad_shape, grad_nCavity, nc, sigma);
	gpuErrorCheck();
}
