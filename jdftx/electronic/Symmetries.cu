/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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
#include <cstdio>

__global__
void symmetrize_kernel(int N, int nRot, double* x, int* symmIndex, double nRotInv)
{	int i=kernelIndex1D();
	if(i<N)
	{	double xSum = 0.0;
		for(int j=0; j<nRot; j++) xSum += x[symmIndex[nRot*i+j]];
		xSum *= nRotInv; //average n in the equivalence class
		for(int j=0; j<nRot; j++) x[symmIndex[nRot*i+j]] = xSum;
	}
}
void symmetrize_gpu(int N, int nRot, double* x, int* symmIndex)
{	GpuLaunchConfig1D glc(symmetrize_kernel, N);
	symmetrize_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, nRot, x, symmIndex, 1.0/nRot);
	gpuErrorCheck();
}
