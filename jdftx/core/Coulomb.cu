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

#include <core/Coulomb_internal.h>
#include <core/GpuKernelUtils.h>
#include <core/LoopMacros.h>

template<typename Coulomb_calc> __global__
void coulombAnalytic_kernel(int zBlock, vector3<int> S, const matrix3<> GGT, const Coulomb_calc calc, complex* data)
{	COMPUTE_halfGindices
	data[i] *= calc(iG, GGT);
}
#define DECLARE_coulombAnalytic_gpu(Type) \
	void coulombAnalytic_gpu(vector3<int> S, const matrix3<>& GGT, const Coulomb##Type##_calc& calc, complex* data) \
	{	GpuLaunchConfigHalf3D glc(coulombAnalytic_kernel<Coulomb##Type##_calc>, S); \
		for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++) \
			coulombAnalytic_kernel<Coulomb##Type##_calc><<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, GGT, calc, data); \
	}
DECLARE_coulombAnalytic_gpu(Periodic)
DECLARE_coulombAnalytic_gpu(Slab)
DECLARE_coulombAnalytic_gpu(Spherical)
#undef DECLARE_coulombAnalytic_gpu
