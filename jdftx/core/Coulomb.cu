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


template<typename Coulomb_calc> __global__
void coulombAnalyticStress_kernel(int zBlock, vector3<int> S, const matrix3<> GGT, const Coulomb_calc calc,
	const complex* X, const complex* Y, symmetricMatrix3<>* grad_RRT)
{
	COMPUTE_halfGindices
	double weight = ((iG[2]==0) or (2*iG[2]==S[2])) ? 1 : 2; //weight factor for points in reduced reciprocal space of real scalar fields
	grad_RRT[i] = (weight * real(X[i].conj() * Y[i])) * calc.latticeGradient(iG, GGT);
}
#define DECLARE_coulombAnalyticStress_gpu(Type) \
	void coulombAnalyticStress_gpu(vector3<int> S, const matrix3<>& GGT, const Coulomb##Type##_calc& calc, \
		const complex* X, const complex* Y, symmetricMatrix3<>* grad_RRT) \
	{	\
		GpuLaunchConfigHalf3D glc(coulombAnalyticStress_kernel<Coulomb##Type##_calc>, S); \
		for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++) \
			coulombAnalyticStress_kernel<Coulomb##Type##_calc><<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, GGT, calc, X, Y, grad_RRT); \
	}
DECLARE_coulombAnalyticStress_gpu(Periodic)
DECLARE_coulombAnalyticStress_gpu(Slab)
DECLARE_coulombAnalyticStress_gpu(Spherical)
DECLARE_coulombAnalyticStress_gpu(IonKernel)
#undef DECLARE_coulombAnalytic_gpu


__global__ void coulombNumericalStress_kernel(int zBlock, vector3<int> S, const matrix3<> GGT, const symmetricMatrix3<>* Vc_RRT,
	const complex* X, const complex* Y, symmetricMatrix3<>* grad_RRT)
{
	COMPUTE_halfGindices
	double weight = ((iG[2]==0) or (2*iG[2]==S[2])) ? 1 : 2; //weight factor for points in reduced reciprocal space of real scalar fields
	grad_RRT[i] = (weight * real(X[i].conj() * Y[i])) * Vc_RRT[i];
}
void coulombNumericalStress_gpu(vector3<int> S, const matrix3<>& GGT, const symmetricMatrix3<>* Vc_RRT,
	const complex* X, const complex* Y, symmetricMatrix3<>* grad_RRT)
{
	GpuLaunchConfigHalf3D glc(coulombNumericalStress_kernel, S);
	for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
		coulombNumericalStress_kernel<<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, GGT, Vc_RRT, X, Y, grad_RRT);
}


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


template<typename Exchange_calc> __global__
void exchangeAnalyticStress_kernel(int zBlock, vector3<int> S, const matrix3<> G, const Exchange_calc calc,
	const complex* X, symmetricMatrix3<>* grad_RRT, const vector3<> kDiff, double thresholdSq)
{
	COMPUTE_fullGindices
	vector3<> iGkCart = (iG + kDiff) * G;
	double kplusGsq = iGkCart.length_squared();
	double XiNorm = X[i].norm();
	double grad_lnDetR = 0.;
	grad_RRT[i] = outer(iGkCart) * (kplusGsq<thresholdSq ? 0. : calc.latticeGradientPrefac(kplusGsq, grad_lnDetR)*XiNorm);
	if(grad_lnDetR) *((vector3<>*)(grad_RRT+i)) += grad_lnDetR * XiNorm; //diagonal contribution due to detR gradient
}
template<> __global__
void exchangeAnalyticStress_kernel<ExchangeSlab_calc>(int zBlock, vector3<int> S, const matrix3<> G, const ExchangeSlab_calc calc,
	const complex* X, symmetricMatrix3<>* grad_RRT, const vector3<> kDiff, double thresholdSq)
{
	COMPUTE_fullGindices
	grad_RRT[i] = calc.latticeGradient(iG, G, kDiff, thresholdSq) * X[i].norm();
}
#define DECLARE_exchangeAnalyticStress_gpu(Type) \
	void exchangeAnalyticStress_gpu(vector3<int> S, const matrix3<>& G, const Exchange##Type##_calc& calc, \
		const complex* X, symmetricMatrix3<>* grad_RRT, const vector3<>& kDiff, double thresholdSq) \
	{	\
		GpuLaunchConfig3D glc(exchangeAnalyticStress_kernel<Exchange##Type##_calc>, S); \
		for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++) \
			exchangeAnalyticStress_kernel<Exchange##Type##_calc><<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, G, calc, \
				X, grad_RRT, kDiff, thresholdSq); \
	}
DECLARE_exchangeAnalyticStress_gpu(Periodic)
DECLARE_exchangeAnalyticStress_gpu(PeriodicScreened)
DECLARE_exchangeAnalyticStress_gpu(Spherical)
DECLARE_exchangeAnalyticStress_gpu(SphericalScreened)
DECLARE_exchangeAnalyticStress_gpu(Slab)
#undef DECLARE_exchangeAnalyticStress_gpu


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


__global__
void realKernelStress_kernel(int zBlock, vector3<int> S, const symmetricMatrix3<>* kernel_RRT, const complex* X, symmetricMatrix3<>* grad_RRT)
{	COMPUTE_fullGindices
	realKernelStress_calc(i, iG, S, kernel_RRT, X, grad_RRT);
}
void realKernelStress_gpu(vector3<int> S, const symmetricMatrix3<>* kernel_RRT, const complex* X, symmetricMatrix3<>* grad_RRT)
{	GpuLaunchConfig3D glc(realKernelStress_kernel, S);
	for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
		realKernelStress_kernel<<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, kernel_RRT, X, grad_RRT);
}


__global__
void transformedKernelStress_kernel(int zBlock, vector3<int> S, const symmetricMatrix3<>* kernel_RRT, const complex* X, symmetricMatrix3<>* grad_RRT, const vector3<int> offset)
{	COMPUTE_fullGindices
	transformedKernelStress_calc(i, iG, S, kernel_RRT, X, grad_RRT, offset);
}
void transformedKernelStress_gpu(vector3<int> S, const symmetricMatrix3<>* kernel_RRT, const complex* X, symmetricMatrix3<>* grad_RRT, const vector3<int>& offset)
{	GpuLaunchConfig3D glc(transformedKernelStress_kernel, S);
	for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
		transformedKernelStress_kernel<<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, kernel_RRT, X, grad_RRT, offset);
}