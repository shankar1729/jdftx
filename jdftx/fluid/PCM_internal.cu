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
	void ScreeningApply_kernel(size_t N, const RadialFunctionG ionEnergyLookup,
			const double* s, const double* phi, double* A, double* A_phi, double* A_s, const Screening eval)
	{	int i = kernelIndex1D(); if(i<N) eval.apply_calc(i, ionEnergyLookup, s, phi, A, A_phi, A_s);
	}
	void Screening::apply_gpu(size_t N, const RadialFunctionG& ionEnergyLookup,
			const double* s, const double* phi, double* A, double* A_phi, double* A_s) const
	{	GpuLaunchConfig1D glc(ScreeningApply_kernel, N);
		ScreeningApply_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, ionEnergyLookup, s, phi, A, A_phi, A_s, *this);
		gpuErrorCheck();
	}

	__global__
	void DielectricApply_kernel(size_t N, const RadialFunctionG dielEnergyLookup,
			const double* s, vector3<const double*> Dphi, double* A, vector3<double*> A_Dphi, double* A_s,
			const Dielectric eval)
	{	int i = kernelIndex1D(); if(i<N) eval.apply_calc(i, dielEnergyLookup, s, Dphi, A, A_Dphi, A_s);
	}
	void Dielectric::apply_gpu(size_t N, const RadialFunctionG& dielEnergyLookup,
			const double* s, vector3<const double*> Dphi, double* A, vector3<double*> A_Dphi, double* A_s) const
	{	GpuLaunchConfig1D glc(DielectricApply_kernel, N);
		DielectricApply_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, dielEnergyLookup, s, Dphi, A, A_Dphi, A_s, *this);
		gpuErrorCheck();
	}
}
