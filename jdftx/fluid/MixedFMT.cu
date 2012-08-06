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

#include <fluid/MixedFMT_internal.h>
#include <core/GpuKernelUtils.h>
#include <core/LoopMacros.h>

__global__
void tensorKernel_kernel(int zBlock, const vector3<int> S, const matrix3<> G, const complex* nTilde, tensor3<complex*> mTilde)
{	COMPUTE_halfGindices
	tensorKernel_calc(i, iG, IS_NYQUIST, G, nTilde, mTilde);
}
void tensorKernel_gpu(const vector3<int> S, const matrix3<> G, const complex* nTilde, tensor3<complex*> mTilde)
{	GpuLaunchConfigHalf3D glc(tensorKernel_kernel, S);
	for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
		tensorKernel_kernel<<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, G, nTilde, mTilde);
	gpuErrorCheck();
}

__global__
void tensorKernel_grad_kernel(int zBlock, const vector3<int> S, const matrix3<> G, tensor3<const complex*> grad_mTilde, complex* grad_nTilde)
{	COMPUTE_halfGindices
	tensorKernel_grad_calc(i, iG, IS_NYQUIST, G, grad_mTilde, grad_nTilde);
}
void tensorKernel_grad_gpu(const vector3<int> S, const matrix3<> G, tensor3<const complex*> grad_mTilde, complex* grad_nTilde)
{	GpuLaunchConfigHalf3D glc(tensorKernel_grad_kernel, S);
	for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
		tensorKernel_grad_kernel<<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, G, grad_mTilde, grad_nTilde);
	gpuErrorCheck();
}

__global__
void phiFMT_kernel(int N, double* phiArr,
	const double *n0arr, const double *n1arr, const double *n2arr, const double *n3arr,
	vector3<const double*> n1vArr, vector3<const double*> n2vArr, tensor3<const double*> n2mArr,
	double *grad_n0arr, double *grad_n1arr, double *grad_n2arr, double *grad_n3arr,
	vector3<double*> grad_n1vArr, vector3<double*> grad_n2vArr, tensor3<double*> grad_n2mArr)
{	int i = kernelIndex1D();
	if(i<N)
		phiArr[i] = phiFMT_calc(i, n0arr, n1arr, n2arr, n3arr, n1vArr, n2vArr, n2mArr,
			grad_n0arr, grad_n1arr, grad_n2arr, grad_n3arr, grad_n1vArr, grad_n2vArr, grad_n2mArr);
}
void phiFMT_gpu(int N, double* phiArr,
	const double *n0arr, const double *n1arr, const double *n2arr, const double *n3arr,
	vector3<const double*> n1vArr, vector3<const double*> n2vArr, tensor3<const double*> n2mArr,
	double *grad_n0arr, double *grad_n1arr, double *grad_n2arr, double *grad_n3arr,
	vector3<double*> grad_n1vArr, vector3<double*> grad_n2vArr, tensor3<double*> grad_n2mArr)
{	GpuLaunchConfig1D glc(phiFMT_kernel, N);
	phiFMT_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, phiArr, n0arr, n1arr, n2arr, n3arr, n1vArr, n2vArr, n2mArr,
			grad_n0arr, grad_n1arr, grad_n2arr, grad_n3arr, grad_n1vArr, grad_n2vArr, grad_n2mArr);
}

__global__
void phiBond_kernel(int N, double Rhm, double scale, double* phiArr,
	const double *n0arr, const double *n2arr, const double *n3arr, vector3<const double*> n2vArr,
	double *grad_n0arr, double *grad_n2arr, double *grad_n3arr, vector3<double*> grad_n2vArr)
{	int i = kernelIndex1D();
	if(i<N)
		phiArr[i] = phiBond_calc(i, Rhm, scale,
			n0arr, n2arr, n3arr, n2vArr, grad_n0arr, grad_n2arr, grad_n3arr, grad_n2vArr);
}
void phiBond_gpu(int N, double Rhm, double scale, double* phiArr,
	const double *n0arr, const double *n2arr, const double *n3arr, vector3<const double*> n2vArr,
	double *grad_n0arr, double *grad_n2arr, double *grad_n3arr, vector3<double*> grad_n2vArr)
{	GpuLaunchConfig1D glc(phiBond_kernel, N);
	phiBond_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, Rhm, scale, phiArr,
		n0arr, n2arr, n3arr, n2vArr, grad_n0arr, grad_n2arr, grad_n3arr, grad_n2vArr);
}
