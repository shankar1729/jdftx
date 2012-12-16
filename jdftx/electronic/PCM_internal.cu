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

__global__
void pcmShapeFunc_kernel(int N, const double* nCavity, double* shape, const double nc, const double sigma)
{	int i = kernelIndex1D(); if(i<N) pcmShapeFunc_calc(i, nCavity, shape, nc, sigma);
}
void pcmShapeFunc_gpu(int N, const double* nCavity, double* shape, const double nc, const double sigma)
{	GpuLaunchConfig1D glc(pcmShapeFunc_kernel, N);
	pcmShapeFunc_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, nCavity, shape, nc, sigma);
	gpuErrorCheck();
}

__global__
void pcmShapeFunc_grad_kernel(int N, const double* nCavity, const double* grad_shape, double* grad_nCavity, const double nc, const double sigma)
{	int i = kernelIndex1D(); if(i<N) pcmShapeFunc_grad_calc(i, nCavity, grad_shape, grad_nCavity, nc, sigma);
}
void pcmShapeFunc_grad_gpu(int N, const double* nCavity, const double* grad_shape, double* grad_nCavity, const double nc, const double sigma)
{	GpuLaunchConfig1D glc(pcmShapeFunc_grad_kernel, N);
	pcmShapeFunc_grad_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, nCavity, grad_shape, grad_nCavity, nc, sigma);
	gpuErrorCheck();
}

//------------- Helper classes for NonlinearPCM  -------------
namespace NonlinearPCMeval
{
	__global__
	void ScreeningFreeEnergy_kernel(size_t N, double mu0, const double* mu, const double* s, double* rho, double* A, double* A_mu, double* A_s, const Screening eval)
	{	int i = kernelIndex1D(); if(i<N) eval.freeEnergy_calc(i, mu0, mu, s, rho, A, A_mu, A_s);
	}
	void Screening::freeEnergy_gpu(size_t N, double mu0, const double* mu, const double* s, double* rho, double* A, double* A_mu, double* A_s) const
	{	GpuLaunchConfig1D glc(ScreeningFreeEnergy_kernel, N);
		ScreeningFreeEnergy_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, mu0, mu, s, rho, A, A_mu, A_s, *this);
		gpuErrorCheck();
	}
	
	__global__
	void ScreeningConvertDerivative_kernel(size_t N, double mu0, const double* mu, const double* s, const double* A_rho, double* A_mu, double* A_s, const Screening eval)
	{	int i = kernelIndex1D(); if(i<N) eval.convertDerivative_calc(i, mu0, mu, s, A_rho, A_mu, A_s);
	}
	void Screening::convertDerivative_gpu(size_t N, double mu0, const double* mu, const double* s, const double* A_rho, double* A_mu, double* A_s) const
	{	GpuLaunchConfig1D glc(ScreeningConvertDerivative_kernel, N);
		ScreeningConvertDerivative_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, mu0, mu, s, A_rho, A_mu, A_s, *this);
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
