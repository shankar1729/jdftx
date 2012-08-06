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

#include <core/LoopMacros.h>
#include <core/GpuKernelUtils.h>
#include <fluid/TranslationOperator_internal.h>

__global__
void constantSplineTaxpy_kernel(int zBlock, const vector3<int> S,
	double alpha, const double* x, double* y, const vector3<int> Tint)
{	COMPUTE_rIndices
	constantSplineTaxpy_calc(i, iv, S, alpha, x, y, Tint);
}
void constantSplineTaxpy_gpu(const vector3<int> S,
	double alpha, const double* x, double* y, const vector3<int> Tint)
{	GpuLaunchConfig3D glc(constantSplineTaxpy_kernel, S);
	for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
		constantSplineTaxpy_kernel<<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, alpha, x, y, Tint);
	gpuErrorCheck();
}

__global__
void linearSplineTaxpy_kernel(int zBlock, const vector3<int> S,
	double alpha, const double* x, double* y, const vector3<int> Tint, const vector3<> Tfrac)
{	COMPUTE_rIndices
	linearSplineTaxpy_calc(i, iv, S, alpha, x, y, Tint, Tfrac);
}
void linearSplineTaxpy_gpu(const vector3<int> S,
	double alpha, const double* x, double* y, const vector3<int> Tint, const vector3<> Tfrac)
{	GpuLaunchConfig3D glc(linearSplineTaxpy_kernel, S);
	for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
		linearSplineTaxpy_kernel<<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, alpha, x, y, Tint, Tfrac);
	gpuErrorCheck();
}



__global__
void fourierTranslate_kernel(int zBlock, const vector3<int> S, const vector3<> Gt, complex* xTilde)
{	COMPUTE_halfGindices
	fourierTranslate_calc(i, iG, S, Gt, xTilde);
}
void fourierTranslate_gpu(const vector3<int> S, const vector3<> Gt, complex* xTilde)
{	GpuLaunchConfigHalf3D glc(fourierTranslate_kernel, S);
	for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
		fourierTranslate_kernel<<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, Gt, xTilde);
}
