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

namespace ShapeFunction
{
	__global__
	void compute_kernel(int N, const double* n, double* shape, const double nc, const double sigma)
	{	int i = kernelIndex1D(); if(i<N) compute_calc(i, n, shape, nc, sigma);
	}
	void compute_gpu(int N, const double* n, double* shape, const double nc, const double sigma)
	{	GpuLaunchConfig1D glc(compute_kernel, N);
		compute_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, n, shape, nc, sigma);
		gpuErrorCheck();
	}

	__global__
	void propagateGradient_kernel(int N, const double* n, const double* grad_shape, double* grad_n, const double nc, const double sigma)
	{	int i = kernelIndex1D(); if(i<N) propagateGradient_calc(i, n, grad_shape, grad_n, nc, sigma);
	}
	void propagateGradient_gpu(int N, const double* n, const double* grad_shape, double* grad_n, const double nc, const double sigma)
	{	GpuLaunchConfig1D glc(propagateGradient_kernel, N);
		propagateGradient_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, n, grad_shape, grad_n, nc, sigma);
		gpuErrorCheck();
	}

	__global__
	void expandDensityHelper_kernel(int N, double alpha, const double* nBar, const double* DnBarSq, double* nEx, double* nEx_nBar, double* nEx_DnBarSq)
	{	int i = kernelIndex1D(); if(i<N) expandDensity_calc(i, alpha, nBar, DnBarSq, nEx, nEx_nBar, nEx_DnBarSq);
	}
	void expandDensityHelper_gpu(int N, double alpha, const double* nBar, const double* DnBarSq, double* nEx, double* nEx_nBar, double* nEx_DnBarSq)
	{	GpuLaunchConfig1D glc(expandDensityHelper_kernel, N);
		expandDensityHelper_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, alpha, nBar, DnBarSq, nEx, nEx_nBar, nEx_DnBarSq);
		gpuErrorCheck();
	}
}

//------------- Helper classes for NonlinearPCM  -------------
namespace NonlinearPCMeval
{
	__global__
	void ScreeningFreeEnergy_kernel(size_t N, double mu0, const double* muPlus, const double* muMinus, const double* s, double* rho, double* A, double* A_muPlus, double* A_muMinus, double* A_s, const Screening eval)
	{	int i = kernelIndex1D(); if(i<N) eval.freeEnergy_calc(i, mu0, muPlus, muMinus, s, rho, A, A_muPlus, A_muMinus, A_s);
	}
	void Screening::freeEnergy_gpu(size_t N, double mu0, const double* muPlus, const double* muMinus, const double* s, double* rho, double* A, double* A_muPlus, double* A_muMinus, double* A_s) const
	{	GpuLaunchConfig1D glc(ScreeningFreeEnergy_kernel, N);
		ScreeningFreeEnergy_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, mu0, muPlus, muMinus, s, rho, A, A_muPlus, A_muMinus, A_s, *this);
		gpuErrorCheck();
	}
	
	__global__
	void ScreeningConvertDerivative_kernel(size_t N, double mu0, const double* muPlus, const double* muMinus, const double* s, const double* A_rho, double* A_muPlus, double* A_muMinus, double* A_s, const Screening eval)
	{	int i = kernelIndex1D(); if(i<N) eval.convertDerivative_calc(i, mu0, muPlus, muMinus, s, A_rho, A_muPlus, A_muMinus, A_s);
	}
	void Screening::convertDerivative_gpu(size_t N, double mu0, const double* muPlus, const double* muMinus, const double* s, const double* A_rho, double* A_muPlus, double* A_muMinus, double* A_s) const
	{	GpuLaunchConfig1D glc(ScreeningConvertDerivative_kernel, N);
		ScreeningConvertDerivative_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, mu0, muPlus, muMinus, s, A_rho, A_muPlus, A_muMinus, A_s, *this);
		gpuErrorCheck();
	}
	
	__global__
	void DielectricFreeEnergy_kernel(size_t N, vector3<const double*> eps, const double* s, vector3<double*> p, double* A, vector3<double*> A_eps, double* A_s, const Dielectric eval)
	{	int i = kernelIndex1D(); if(i<N) eval.freeEnergy_calc(i, eps, s, p, A, A_eps, A_s);
	}
	void Dielectric::freeEnergy_gpu(size_t N, vector3<const double*> eps, const double* s, vector3<double*> p, double* A, vector3<double*> A_eps, double* A_s) const
	{	GpuLaunchConfig1D glc(DielectricFreeEnergy_kernel, N);
		DielectricFreeEnergy_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, eps, s, p, A, A_eps, A_s, *this);
		gpuErrorCheck();
	}
	
	__global__
	void DielectricConvertDerivative_kernel(size_t N, vector3<const double*> eps, const double* s, vector3<const double*> A_p, vector3<double*> A_eps, double* A_s, const Dielectric eval)
	{	int i = kernelIndex1D(); if(i<N) eval.convertDerivative_calc(i, eps, s, A_p, A_eps, A_s);
	}
	void Dielectric::convertDerivative_gpu(size_t N, vector3<const double*> eps, const double* s, vector3<const double*> A_p, vector3<double*> A_eps, double* A_s) const
	{	GpuLaunchConfig1D glc(DielectricConvertDerivative_kernel, N);
		DielectricConvertDerivative_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, eps, s, A_p, A_eps, A_s, *this);
		gpuErrorCheck();
	}
}
