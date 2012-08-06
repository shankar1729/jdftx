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

#include <electronic/ExactExchange_internal.h>
#include <core/GpuKernelUtils.h>
#include <core/LoopMacros.h>

__global__
void screenedCoulombK_kernel(int zBlock, const vector3<int> S, const matrix3<> GGT,
	const complex* in, complex* out, vector3<> kDiff, double weightedVzero, double omegaSq)
{	COMPUTE_fullGindices
	out[i] = in[i] * screenedCoulombK_calc(iG, GGT, kDiff, weightedVzero, omegaSq); 
}
void screenedCoulombK_gpu(const vector3<int>& S, const matrix3<>& GGT,
	const complex* in, complex* out, const vector3<>& kDiff, double weightedVzero, double omegaSq)
{	GpuLaunchConfig3D glc(screenedCoulombK_kernel, S);
	for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
		screenedCoulombK_kernel<<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, GGT, in, out, kDiff, weightedVzero, omegaSq);
	gpuErrorCheck();
}
