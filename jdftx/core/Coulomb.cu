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


template<typename Exchange_calc> __global__
void exchangeAnalytic_kernel(int zBlock, vector3<int> S, const matrix3<> GGT, const Exchange_calc calc,
	complex* data, const vector3<> kDiff, double Vzero, double thresholdSq)
{
	COMPUTE_fullGindices
	double kplusGsq = GGT.metric_length_squared(iG + kDiff);
	data[i] *= kplusGsq<thresholdSq ? Vzero : calc(kplusGsq);
}
//Specialization for slab mode:
template<> void exchangeAnalytic_kernel<ExchangeSlab_calc>(int zBlock,
	vector3<int> S, const matrix3<> GGT, const ExchangeSlab_calc calc,
	complex* data, const vector3<> kDiff, double Vzero, double thresholdSq)
{
	COMPUTE_fullGindices
	data[i] *= calc(iG, GGT, kDiff, Vzero, thresholdSq);
}
#define DECLARE_exchangeAnalytic_gpu(Type) \
	void exchangeAnalytic_gpu(vector3<int> S, const matrix3<>& GGT, const Exchange##Type##_calc& calc, \
		complex* data, const vector3<>& kDiff, double Vzero, double thresholdSq) \
	{	\
		GpuLaunchConfig3D glc(exchangeAnalytic_kernel<Exchange##Type##_calc>, S); \
		for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++) \
			exchangeAnalytic_kernel<Exchange##Type##_calc><<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, GGT, calc, \
				data, kDiff, Vzero, thresholdSq); \
	}
DECLARE_exchangeAnalytic_gpu(Periodic)
DECLARE_exchangeAnalytic_gpu(PeriodicScreened)
DECLARE_exchangeAnalytic_gpu(Spherical)
DECLARE_exchangeAnalytic_gpu(SphericalScreened)
DECLARE_exchangeAnalytic_gpu(Slab)
#undef DECLARE_exchangeAnalytic_gpu

__global__
void multRealKernel_kernel(int zBlock, vector3<int> S, const double* kernel, complex* data)
{	COMPUTE_fullGindices
	multRealKernel_calc(i, iG, S, kernel, data);
}
void multRealKernel_gpu(vector3<int> S, const double* kernel, complex* data)
{	GpuLaunchConfig3D glc(multRealKernel_kernel, S);
	for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
		multRealKernel_kernel<<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, kernel, data);
}

__global__
void multTransformedKernel_kernel(int zBlock, vector3<int> S, const double* kernel, complex* data, const vector3<int> offset)
{	COMPUTE_fullGindices
	multTransformedKernel_calc(i, iG, S, kernel, data, offset);
}
void multTransformedKernel_gpu(vector3<int> S, const double* kernel, complex* data, const vector3<int>& offset)
{	GpuLaunchConfig3D glc(multTransformedKernel_kernel, S);
	for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
		multTransformedKernel_kernel<<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, kernel, data, offset);
}
