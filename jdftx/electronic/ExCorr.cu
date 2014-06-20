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

#include <core/GpuKernelUtils.h>
#include <electronic/ExCorr_internal.h>
#include <electronic/ExCorr_internal_LDA.h>
#include <electronic/ExCorr_internal_GGA.h>
#include <electronic/ExCorr_internal_mGGA.h>

//---------------- Spin-density-matrix transformations for noncollinear magentism --------------------

__global__
void spinDiagonalize_kernel(int N, array<const double*,4> n, array<const double*,4> x, array<double*,2> xDiag)
{	int i = kernelIndex1D();
	if(i<N) spinDiagonalize_calc(i, n, x, xDiag);
}
void spinDiagonalize_gpu(int N, std::vector<const double*> n, std::vector<const double*> x, std::vector<double*> xDiag)
{	GpuLaunchConfig1D glc(spinDiagonalize_kernel, N);
	spinDiagonalize_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, n, x, xDiag);
	gpuErrorCheck();
}

__global__
void spinDiagonalizeGrad_kernel(int N, array<const double*,4> n, array<const double*,4> x, array<const double*,2> E_xDiag, array<double*,4> E_n, array<double*,4> E_x)
{	int i = kernelIndex1D();
	if(i<N) spinDiagonalizeGrad_calc(i, n, x, E_xDiag, E_n, E_x);
}
void spinDiagonalizeGrad_gpu(int N, std::vector<const double*> n, std::vector<const double*> x, std::vector<const double*> E_xDiag, std::vector<double*> E_n, std::vector<double*> E_x)
{	GpuLaunchConfig1D glc(spinDiagonalizeGrad_kernel, N);
	spinDiagonalizeGrad_kernel<<<glc.nBlocks,glc.nPerBlock>>>(N, n, x, E_xDiag, E_n, E_x);
	gpuErrorCheck();
}

//-------------------------- LDA GPU launch mechanism ----------------------------

template<LDA_Variant variant, int nCount> __global__
void LDA_kernel(int N, array<const double*,nCount> n, double* E, array<double*,nCount> E_n, double scaleFac)
{	int i = kernelIndex1D();
	if(i<N) LDA_calc<variant,nCount>::compute(i, n, E, E_n, scaleFac);
}
template<LDA_Variant variant, int nCount>
void LDA_gpu(int N, array<const double*,nCount> n, double* E, array<double*,nCount> E_n, double scaleFac)
{	GpuLaunchConfig1D glc(LDA_kernel<variant,nCount>, N);
	LDA_kernel<variant,nCount><<<glc.nBlocks,glc.nPerBlock>>>(N, n, E, E_n, scaleFac);
	gpuErrorCheck();
}
void LDA_gpu(LDA_Variant variant, int N, std::vector<const double*> n, double* E, std::vector<double*> E_n, double scaleFac)
{	SwitchTemplate_spin(SwitchTemplate_LDA, variant, n.size(), LDA_gpu, (N, n, E, E_n, scaleFac) )
}

//-------------------------- GGA GPU launch mechanism ----------------------------

template<GGA_Variant variant, bool spinScaling, int nCount> __global__
void GGA_kernel(int N, array<const double*,nCount> n, array<const double*,2*nCount-1> sigma,
	double* E, array<double*,nCount> E_n, array<double*,2*nCount-1> E_sigma, double scaleFac)
{	int i = kernelIndex1D();
	if(i<N) GGA_calc<variant,spinScaling,nCount>::compute(i, n, sigma, E, E_n, E_sigma, scaleFac);
}
template<GGA_Variant variant, bool spinScaling, int nCount>
void GGA_gpu(int N, array<const double*,nCount> n, array<const double*,2*nCount-1> sigma,
	double* E, array<double*,nCount> E_n, array<double*,2*nCount-1> E_sigma, double scaleFac)
{	GpuLaunchConfig1D glc(GGA_kernel<variant,spinScaling,nCount>, N);
	GGA_kernel<variant,spinScaling,nCount><<<glc.nBlocks,glc.nPerBlock>>>(N, n, sigma, E, E_n, E_sigma, scaleFac);
	gpuErrorCheck();
}
void GGA_gpu(GGA_Variant variant, int N, std::vector<const double*> n, std::vector<const double*> sigma,
	double* E, std::vector<double*> E_n, std::vector<double*> E_sigma, double scaleFac)
{	SwitchTemplate_spin(SwitchTemplate_GGA, variant, n.size(), GGA_gpu, (N, n, sigma, E, E_n, E_sigma, scaleFac) )
}

//-------------------------- metaGGA GPU launch mechanism ----------------------------

template<mGGA_Variant variant, bool spinScaling, int nCount> __global__
void mGGA_kernel(int N, array<const double*,nCount> n, array<const double*,2*nCount-1> sigma,
	array<const double*,nCount> lap, array<const double*,nCount> tau,
	double* E, array<double*,nCount> E_n, array<double*,2*nCount-1> E_sigma,
	array<double*,nCount> E_lap, array<double*,nCount> E_tau, double scaleFac)
{	int i = kernelIndex1D();
	if(i<N) mGGA_calc<variant,spinScaling,nCount>::compute(i,
		n, sigma, lap, tau, E, E_n, E_sigma, E_lap, E_tau, scaleFac);
}
template<mGGA_Variant variant, bool spinScaling, int nCount>
void mGGA_gpu(int N, array<const double*,nCount> n, array<const double*,2*nCount-1> sigma,
	array<const double*,nCount> lap, array<const double*,nCount> tau,
	double* E, array<double*,nCount> E_n, array<double*,2*nCount-1> E_sigma,
	array<double*,nCount> E_lap, array<double*,nCount> E_tau, double scaleFac)
{	GpuLaunchConfig1D glc(mGGA_kernel<variant,spinScaling,nCount>, N);
	mGGA_kernel<variant,spinScaling,nCount><<<glc.nBlocks,glc.nPerBlock>>>(N,
		n, sigma, lap, tau, E, E_n, E_sigma, E_lap, E_tau, scaleFac);
	gpuErrorCheck();
}
void mGGA_gpu(mGGA_Variant variant, int N, std::vector<const double*> n, std::vector<const double*> sigma,
	std::vector<const double*> lap, std::vector<const double*> tau,
	double* E, std::vector<double*> E_n, std::vector<double*> E_sigma,
	std::vector<double*> E_lap, std::vector<double*> E_tau, double scaleFac)
{	SwitchTemplate_spin(SwitchTemplate_mGGA, variant, n.size(), mGGA_gpu, 
		(N, n, sigma, lap, tau, E, E_n, E_sigma, E_lap, E_tau, scaleFac) )
}

