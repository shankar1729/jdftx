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

#include <electronic/ColumnBundleOperators_internal.h>
#include <core/GpuKernelUtils.h>
#include <core/LoopMacros.h>

__global__
void reducedL_kernel(int nbasis, int ncols, const complex* Y, complex* LY,
	const matrix3<> GGT, const vector3<int>* iGarr, const vector3<> k, double detR)
{	int j = kernelIndex1D();
	if(j<nbasis) reducedL_calc(j, nbasis, ncols, Y, LY, GGT, iGarr, k, detR);
}
void reducedL_gpu(int nbasis, int ncols, const complex* Y, complex* LY,
	const matrix3<> GGT, const vector3<int>* iGarr, const vector3<> k, double detR)
{	GpuLaunchConfig1D glc(reducedL_kernel, nbasis);
	reducedL_kernel<<<glc.nBlocks,glc.nPerBlock>>>(nbasis, ncols, Y, LY, GGT, iGarr, k, detR);
	gpuErrorCheck();
}

__global__
void reducedLinv_kernel(int nbasis, int ncols, const complex* Y, complex* LinvY,
	const matrix3<> GGT, const vector3<int>* iGarr, const vector3<> k, double detR)
{	int j = kernelIndex1D();
	if(j<nbasis) reducedLinv_calc(j, nbasis, ncols, Y, LinvY, GGT, iGarr, k, detR);
}
void reducedLinv_gpu(int nbasis, int ncols, const complex* Y, complex* LinvY,
	const matrix3<> GGT, const vector3<int>* iGarr, const vector3<> k, double detR)
{	GpuLaunchConfig1D glc(reducedLinv_kernel, nbasis);
	reducedLinv_kernel<<<glc.nBlocks,glc.nPerBlock>>>(nbasis, ncols, Y, LinvY, GGT, iGarr, k, detR);
	gpuErrorCheck();
}


__global__
void precond_inv_kinetic_kernel(int nbasis, int ncols, complex* Y, 
	double KErollover, const matrix3<> GGT, const vector3<int>* iGarr, const vector3<> k, double invdetR)
{	int j = kernelIndex1D();
	if(j<nbasis) precond_inv_kinetic_calc(j, nbasis, ncols, Y, KErollover, GGT, iGarr, k, invdetR);
}
void precond_inv_kinetic_gpu(int nbasis, int ncols, complex* Y,
	double KErollover, const matrix3<> GGT, const vector3<int>* iGarr, const vector3<> k, double invdetR)
{	GpuLaunchConfig1D glc(precond_inv_kinetic_kernel, nbasis);
	precond_inv_kinetic_kernel<<<glc.nBlocks,glc.nPerBlock>>>(nbasis, ncols, Y, KErollover, GGT, iGarr, k, invdetR);
	gpuErrorCheck();
}

__global__
void precond_inv_kinetic_band_kernel(int nbasis, int ncols, complex* Y, const double* KEref,
	const matrix3<> GGT, const vector3<int>* iGarr, const vector3<> k)
{	int j = kernelIndex1D();
	if(j<nbasis) precond_inv_kinetic_band_calc(j, nbasis, ncols, Y, KEref, GGT, iGarr, k);
}
void precond_inv_kinetic_band_gpu(int nbasis, int ncols, complex* Y, const double* KEref,
	const matrix3<>& GGT, const vector3<int>* iGarr, const vector3<>& k)
{	GpuLaunchConfig1D glc(precond_inv_kinetic_band_kernel, nbasis);
	precond_inv_kinetic_band_kernel<<<glc.nBlocks,glc.nPerBlock>>>(nbasis, ncols, Y, KEref, GGT, iGarr, k);
	gpuErrorCheck();
}

__global__
void translate_kernel(int nbasis, int ncols, complex* Y, const vector3<int>* iGarr, const vector3<> k, const vector3<> dr)
{	int j = kernelIndex1D();
	if(j<nbasis) translate_calc(j, nbasis, ncols, Y, iGarr, k, dr);
}
void translate_gpu(int nbasis, int ncols, complex* Y, const vector3<int>* iGarr, const vector3<>& k, const vector3<>& dr)
{	GpuLaunchConfig1D glc(translate_kernel, nbasis);
	translate_kernel<<<glc.nBlocks,glc.nPerBlock>>>(nbasis, ncols, Y, iGarr, k, dr);
	gpuErrorCheck();
}

__global__
void translateColumns_kernel(int nbasis, int ncols, complex* Y, const vector3<int>* iGarr, const vector3<> k, const vector3<>* dr)
{	int j = kernelIndex1D();
	if(j<nbasis) translateColumns_calc(j, nbasis, ncols, Y, iGarr, k, dr);
}
void translateColumns_gpu(int nbasis, int ncols, complex* Y, const vector3<int>* iGarr, const vector3<>& k, const vector3<>* dr)
{	GpuLaunchConfig1D glc(translateColumns_kernel, nbasis);
	translateColumns_kernel<<<glc.nBlocks,glc.nPerBlock>>>(nbasis, ncols, Y, iGarr, k, dr);
	gpuErrorCheck();
}


__global__
void reducedD_kernel(int nbasis, int ncols, const complex* Y, complex* DY, 
	const vector3<int>* iGarr, double kdotGe, const vector3<> Ge)
{	int j = kernelIndex1D();
	if(j<nbasis) reducedD_calc(j, nbasis, ncols, Y, DY, iGarr, kdotGe, Ge);
}
void reducedD_gpu(int nbasis, int ncols, const complex* Y, complex* DY, 
	const vector3<int>* iGarr, double kdotGe, const vector3<> Ge)
{	GpuLaunchConfig1D glc(reducedD_kernel, nbasis);
	reducedD_kernel<<<glc.nBlocks,glc.nPerBlock>>>(nbasis, ncols, Y, DY, iGarr, kdotGe, Ge);
	gpuErrorCheck();
}

__global__
void reducedDD_kernel(int nbasis, int ncols, const complex* Y, complex* DDY, 
	const vector3<int>* iGarr, double kdotGe1, double kdotGe2, const vector3<> Ge1, const vector3<> Ge2)
{	int j = kernelIndex1D();
	if(j<nbasis) reducedDD_calc(j, nbasis, ncols, Y, DDY, iGarr, kdotGe1, kdotGe2, Ge1, Ge2);
}
void reducedDD_gpu(int nbasis, int ncols, const complex* Y, complex* DDY, 
	const vector3<int>* iGarr, double kdotGe1, double kdotGe2, const vector3<> Ge1, const vector3<> Ge2)
{	GpuLaunchConfig1D glc(reducedDD_kernel, nbasis);
	reducedDD_kernel<<<glc.nBlocks,glc.nPerBlock>>>(nbasis, ncols, Y, DDY, iGarr, kdotGe1, kdotGe2, Ge1, Ge2);
	gpuErrorCheck();
}
