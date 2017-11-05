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
#include <core/LoopMacros.h>
#include <fluid/PCM_internal.h>

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
}

namespace ShapeFunctionCANDLE
{
	__global__
	void compute_or_grad_kernel(int N, bool grad,
		const double* n, vector3<const double*> Dn, vector3<const double*> Dphi, double* shape,
		const double* A_shape, double* A_n, vector3<double*> A_Dn, vector3<double*> A_Dphi, double* A_pCavity,
		const double nc, const double invSigmaSqrt2, const double pCavity)
	{	int i = kernelIndex1D();
		if(i<N) compute_or_grad_calc(i, grad, n, Dn, Dphi, shape, A_shape, A_n, A_Dn, A_Dphi, A_pCavity, nc, invSigmaSqrt2, pCavity);
	}
	void compute_or_grad_gpu(int N, bool grad,
		const double* n, vector3<const double*> Dn, vector3<const double*> Dphi, double* shape,
		const double* A_shape, double* A_n, vector3<double*> A_Dn, vector3<double*> A_Dphi, double* A_pCavity,
		const double nc, const double invSigmaSqrt2, const double pCavity)
	{	GpuLaunchConfig1D glc(compute_or_grad_kernel, N);
		compute_or_grad_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, grad, n, Dn, Dphi, shape, A_shape, A_n, A_Dn, A_Dphi, A_pCavity, nc, invSigmaSqrt2, pCavity);
		gpuErrorCheck();
	}
}

namespace ShapeFunctionSGA13
{
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

namespace ShapeFunctionSoftSphere
{
	__global__
	void compute_kernel(int zBlock, const vector3<int> S, const vector3<> Sinv, const matrix3<> RTR,
		int nAtoms, const vector3<>* x, int nReps, const vector3<int>* reps, const double* radius, double* shape, double sigmaInv)
	{	COMPUTE_rIndices
		compute_calc(i, iv, Sinv, RTR, nAtoms, x, nReps, reps, radius, shape, sigmaInv);
	}
	void compute_gpu(const vector3<int>& S, const matrix3<>& RTR,
		int nAtoms, const vector3<>* x, int nReps, const vector3<int>* reps, const double* radius, double* shape, double sigmaInv)
	{	GpuLaunchConfig3D glc(compute_kernel, S);
		vector3<> Sinv(1./S[0], 1./S[1], 1./S[2]);
		for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
			compute_kernel<<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, Sinv, RTR, nAtoms, x, nReps, reps, radius, shape, sigmaInv);
		gpuErrorCheck();
	}
	
	__global__
	void propagateGradient_kernel(int zBlock, const vector3<int> S, const vector3<> Sinv, const matrix3<> RTR,
		const vector3<> x, int nReps, const vector3<int>* reps, double radius,
		const double* shape, const double* E_shape, vector3<double*> E_x, double* E_radius, double sigmaInv)
	{	COMPUTE_rIndices
		propagateGradient_calc(i, iv, Sinv, RTR, x, nReps, reps, radius, shape, E_shape, E_x, E_radius, sigmaInv);
	}
	void propagateGradient_gpu(const vector3<int>& S, const matrix3<>& RTR,
		const vector3<>& x, int nReps, const vector3<int>* reps, double radius,
		const double* shape, const double* E_shape, vector3<double*> E_x, double* E_radius, double sigmaInv)
	{	GpuLaunchConfig3D glc(propagateGradient_kernel, S);
		vector3<> Sinv(1./S[0], 1./S[1], 1./S[2]);
		for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
			propagateGradient_kernel<<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, Sinv, RTR, x, nReps, reps, radius, shape, E_shape, E_x, E_radius, sigmaInv);
		gpuErrorCheck();
	}
}

namespace ShapeFunctionSCCS
{
	__global__
	void compute_kernel(int N, const double* n, double* shape, const double rhoMin, const double rhoMax, const double epsBulk)
	{	int i = kernelIndex1D(); if(i<N) compute_calc(i, n, shape, rhoMin, rhoMax, epsBulk);
	}
	void compute_gpu(int N, const double* n, double* shape, const double rhoMin, const double rhoMax, const double epsBulk)
	{	GpuLaunchConfig1D glc(compute_kernel, N);
		compute_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, n, shape, rhoMin, rhoMax, epsBulk);
		gpuErrorCheck();
	}

	__global__
	void propagateGradient_kernel(int N, const double* n, const double* grad_shape, double* grad_n, const double rhoMin, const double rhoMax, const double epsBulk)
	{	int i = kernelIndex1D(); if(i<N) propagateGradient_calc(i, n, grad_shape, grad_n, rhoMin, rhoMax, epsBulk);
	}
	void propagateGradient_gpu(int N, const double* n, const double* grad_shape, double* grad_n, const double rhoMin, const double rhoMax, const double epsBulk)
	{	GpuLaunchConfig1D glc(propagateGradient_kernel, N);
		propagateGradient_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, n, grad_shape, grad_n, rhoMin, rhoMax, epsBulk);
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
	void ScreeningPhiToState_kernel(size_t N, const double* phi, const double* s, const RadialFunctionG xLookup, bool setState, double* muPlus, double* muMinus, double* kappaSq, const Screening eval)
	{	int i = kernelIndex1D(); if(i<N) eval.phiToState_calc(i, phi, s, xLookup, setState, muPlus, muMinus, kappaSq);
	}
	void Screening::phiToState_gpu(size_t N, const double* phi, const double* s, const RadialFunctionG& xLookup, bool setState, double* muPlus, double* muMinus, double* kappaSq) const
	{	GpuLaunchConfig1D glc(ScreeningPhiToState_kernel, N);
		ScreeningPhiToState_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, phi, s, xLookup, setState, muPlus, muMinus, kappaSq, *this);
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

	__global__
	void DielectricPhiToState_kernel(size_t N, vector3<const double*> Dphi, const double* s, const RadialFunctionG gLookup, bool setState, vector3<double*> eps, double* epsilon, const Dielectric eval)
	{	int i = kernelIndex1D(); if(i<N) eval.phiToState_calc(i, Dphi, s, gLookup, setState, eps, epsilon);
	}
	void Dielectric::phiToState_gpu(size_t N, vector3<const double*> Dphi, const double* s, const RadialFunctionG& gLookup, bool setState, vector3<double*> eps, double* epsilon) const
	{	GpuLaunchConfig1D glc(DielectricPhiToState_kernel, N);
		DielectricPhiToState_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, Dphi, s, gLookup, setState, eps, epsilon, *this);
		gpuErrorCheck();
	}
}
