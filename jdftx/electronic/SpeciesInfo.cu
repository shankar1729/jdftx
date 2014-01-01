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

#include <electronic/SpeciesInfo_internal.h>
#include <core/GpuKernelUtils.h>
#include <core/LoopMacros.h>

#undef _GLIBCXX_ATOMIC_BUILTINS
#undef _GLIBCXX_USE_INT128
#include <thrust/device_ptr.h>
#include <thrust/sort.h>

//Calculate non-local pseudopotential projector
template<int l, int m> __global__
void Vnl_kernel(int nbasis, int atomStride, int nAtoms, vector3<> k, const vector3<int>* iGarr,
	const matrix3<> G, const vector3<>* pos, const RadialFunctionG VnlRadial, complex* V)
{	int n = kernelIndex1D();
	if(n<nbasis) Vnl_calc<l,m>(n, atomStride, nAtoms, k, iGarr, G, pos, VnlRadial, V);
}
template<int l, int m>
void Vnl_gpu(int nbasis, int atomStride, int nAtoms, vector3<> k, const vector3<int>* iGarr,
	const matrix3<> G, const vector3<>* pos, const RadialFunctionG& VnlRadial, complex* V)
{	GpuLaunchConfig1D glc(Vnl_kernel<l,m>, nbasis);
	Vnl_kernel<l,m><<<glc.nBlocks,glc.nPerBlock>>>(nbasis, atomStride, nAtoms, k, iGarr, G, pos, VnlRadial, V);
	gpuErrorCheck();
}
void Vnl_gpu(int nbasis, int atomStride, int nAtoms, int l, int m, vector3<> k, const vector3<int>* iGarr,
	const matrix3<> G, const vector3<>* pos, const RadialFunctionG& VnlRadial, complex* V)
{
	SwitchTemplate_lm(l,m, Vnl_gpu, (nbasis, atomStride, nAtoms, k, iGarr, G, pos, VnlRadial, V) )
}


//Augment electron density by spherical functions
template<int Nlm> __global__ void nAugment_kernel(int zBlock, const vector3<int> S, const matrix3<> G, int iGstart, int iGstop,
	int nCoeff, double dGinv, const double* nRadial, const vector3<> atpos, complex* n)
{	COMPUTE_halfGindices
	if(i<iGstart || i>=iGstop) return;
	nAugment_calc<Nlm>(i, iG, G, nCoeff, dGinv, nRadial, atpos, n);
}
template<int Nlm> void nAugment_gpu(const vector3<int> S, const matrix3<>& G, int iGstart, int iGstop,
	int nCoeff, double dGinv, const double* nRadial, const vector3<>& atpos, complex* n)
{	GpuLaunchConfigHalf3D glc(nAugment_kernel<Nlm>, S);
	for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
		nAugment_kernel<Nlm><<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, G, iGstart, iGstop, nCoeff, dGinv, nRadial, atpos, n);
	gpuErrorCheck();
}
void nAugment_gpu(int Nlm, const vector3<int> S, const matrix3<>& G, int iGstart, int iGstop,
	int nCoeff, double dGinv, const double* nRadial, const vector3<>& atpos, complex* n)
{
	SwitchTemplate_Nlm(Nlm, nAugment_gpu, (S, G, iGstart, iGstop, nCoeff, dGinv, nRadial, atpos, n) )
}


//Function for initializing the index arrays used by nAugmentGrad
__global__ void setNagIndex_kernel(int zBlock, const vector3<int> S, const matrix3<> G, int iGstart, int iGstop, double dGinv, uint64_t* nagIndex)
{	COMPUTE_halfGindices
	if(i<iGstart || i>=iGstop) return;
	nagIndex[i-iGstart] = setNagIndex_calc(iG, S, G, dGinv);
}
__global__ void setNagIndexPtr_kernel(int iMax, int nCoeff, const uint64_t* nagIndex, size_t* nagIndexPtr)
{	int i = kernelIndex1D();
	if(i>=iMax) return;
	setNagIndexPtr_calc(i, iMax, nCoeff, nagIndex, nagIndexPtr);
}
void setNagIndex_gpu(const vector3<int>& S, const matrix3<>& G, int iGstart, int iGstop, int nCoeff, double dGinv, uint64_t*& nagIndex, size_t*& nagIndexPtr)
{	//First initialize the indices:
	size_t nGsub = iGstop-iGstart;
	{	if(!nagIndex) cudaMalloc(&nagIndex, nGsub*sizeof(uint64_t));
		GpuLaunchConfigHalf3D glc(setNagIndex_kernel, S);
		for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
			setNagIndex_kernel<<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, G, iGstart, iGstop, dGinv, nagIndex);
		gpuErrorCheck();
	}
	//Now sort them to be ordered by Gindex
	thrust::sort(thrust::device_ptr<uint64_t>(nagIndex), thrust::device_ptr<uint64_t>(nagIndex+nGsub));
	gpuErrorCheck();
	//Finally initialize the pointers to boundaries between different Gindices:
	{	if(!nagIndexPtr) cudaMalloc(&nagIndexPtr, (nCoeff+1)*sizeof(size_t));
		GpuLaunchConfig1D glc(setNagIndexPtr_kernel, nGsub);
		setNagIndexPtr_kernel<<<glc.nBlocks,glc.nPerBlock>>>(nGsub, nCoeff, nagIndex, nagIndexPtr);
		gpuErrorCheck();
	}
}


//Propagate gradients corresponding to above electron density augmentation
template<int Nlm> __global__ void nAugmentGrad_kernel(const vector3<int> S, const matrix3<> G,
	int nCoeff, double dGinv, const double* nRadial, const vector3<> atpos,
	const complex* ccE_n, double* E_nRadialTemp, vector3<complex*> E_atpos,
	const uint64_t* nagIndex, const size_t* nagIndexPtr)
{
	const int& iCoeff = blockIdx.x;
	double* E_nRadialThread = E_nRadialTemp + (Nlm*nCoeff) * (iCoeff & 0x7); //Write to (iCoeff % 8)th copy (breaks write dependency (of range 6) between nearby iCoeffs)
	const size_t& ptrStart = nagIndexPtr[iCoeff];
	const size_t& ptrStop = nagIndexPtr[iCoeff+1];
	int nIter = (ptrStop-ptrStart + blockDim.x-1) / blockDim.x; //Note (A+B-1)/B for positive integers equals ceil(A/B)
	for(int iIter=0; iIter<nIter; iIter++)
	{	size_t ptr = ptrStart + threadIdx.x + iIter*blockDim.x;
		bool dummy = (ptr >= ptrStop); //can't quit here, because reduction in QuinticSpline::valueGrad needs all threads to get there
		if(dummy) ptr = ptrStop-1; //to avoid unnecssary branches later (but threads with dummy=true should never write to memory)
		nAugmentGrad_calc<Nlm>(nagIndex[ptr], S, G, nCoeff, dGinv, nRadial, atpos, ccE_n, E_nRadialThread, E_atpos, dummy);
	}
}
__global__ void nAugmentGrad_collectKernel(int nData, const double* in, double* out)
{	int i = kernelIndex1D();
	if(i<nData)
	{	double sum = 0.;
		for(int iCopy=0; iCopy<8; iCopy++)
			sum += in[i + iCopy*nData];
		out[i] += sum;
	}
}
template<int Nlm> void nAugmentGrad_gpu(const vector3<int> S, const matrix3<>& G,
	int nCoeff, double dGinv, const double* nRadial, const vector3<>& atpos,
	const complex* ccE_n, double* E_nRadial, vector3<complex*> E_atpos, const uint64_t* nagIndex, const size_t* nagIndexPtr)
{
	//Allocate temporary memory:
	int iDevice; cudaGetDevice(&iDevice);
	cudaDeviceProp prop; cudaGetDeviceProperties(&prop, iDevice);
	cudaFuncAttributes attr; cudaFuncGetAttributes(&attr, nAugmentGrad_kernel<Nlm>);
	int sharedMemPerThread = 6 * sizeof(double);
	int nPerBlock = std::min(prop.warpSize, std::min(attr.maxThreadsPerBlock, int(prop.sharedMemPerBlock/sharedMemPerThread)));
	int nBlocks = nCoeff;
	double* E_nRadialTemp; cudaMalloc(&E_nRadialTemp, sizeof(double)*nCoeff*Nlm*8);
	cudaMemset(E_nRadialTemp, 0, sizeof(double)*nCoeff*Nlm*8);
	gpuErrorCheck();
	//Stage 1: calculate with the scattered accumulate to E_nRadial
	nAugmentGrad_kernel<Nlm><<<nBlocks,nPerBlock,sharedMemPerThread*nPerBlock>>>(S, G, nCoeff, dGinv, nRadial, atpos, ccE_n, E_nRadialTemp, E_atpos, nagIndex, nagIndexPtr);
	gpuErrorCheck();
	//Stage 2: collect from E_nRadialTemp to E_nRadial
	GpuLaunchConfig1D glc(nAugmentGrad_collectKernel, nCoeff*Nlm);
	nAugmentGrad_collectKernel<<<glc.nBlocks,glc.nPerBlock>>>(nCoeff*Nlm, E_nRadialTemp, E_nRadial);
	gpuErrorCheck();
	//Cleanup:
	cudaFree(E_nRadialTemp);
}
void nAugmentGrad_gpu(int Nlm, const vector3<int> S, const matrix3<>& G,
	int nCoeff, double dGinv, const double* nRadial, const vector3<>& atpos,
	const complex* ccE_n, double* E_nRadial, vector3<complex*> E_atpos, const uint64_t* nagIndex, const size_t* nagIndexPtr)
{	
	SwitchTemplate_Nlm(Nlm, nAugmentGrad_gpu, (S, G, nCoeff, dGinv, nRadial, atpos, ccE_n, E_nRadial, E_atpos, nagIndex, nagIndexPtr) )
}


//Structure factor
__global__
void getSG_kernel(int zBlock, const vector3<int> S, int nAtoms, const vector3<>* atpos, double invVol, complex* SG)
{	COMPUTE_halfGindices
	SG[i] = invVol * getSG_calc(iG, nAtoms, atpos);
}
void getSG_gpu(const vector3<int> S, int nAtoms, const vector3<>* atpos, double invVol, complex* SG)
{	GpuLaunchConfigHalf3D glc(getSG_kernel, S);
	for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
		getSG_kernel<<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, nAtoms, atpos, invVol, SG);
	gpuErrorCheck();
}

//Calculate local pseudopotetial, ionic charge, charge ball and partial core density
__global__
void updateLocal_kernel(int zBlock, const vector3<int> S, const matrix3<> GGT,
	complex *Vlocps,  complex *rhoIon, complex *nChargeball, complex *nCore, complex* tauCore,
	int nAtoms, const vector3<>* atpos, double invVol, const RadialFunctionG VlocRadial,
	double Z, const RadialFunctionG nCoreRadial, const RadialFunctionG tauCoreRadial,
	double Zchargeball, double wChargeball)
{
	COMPUTE_halfGindices
	updateLocal_calc(i, iG, GGT, Vlocps, rhoIon, nChargeball,
		nCore, tauCore, nAtoms, atpos, invVol, VlocRadial,
		Z, nCoreRadial, tauCoreRadial, Zchargeball, wChargeball);
}
void updateLocal_gpu(const vector3<int> S, const matrix3<> GGT,
	complex *Vlocps,  complex *rhoIon, complex *nChargeball, complex *nCore, complex* tauCore,
	int nAtoms, const vector3<>* atpos, double invVol, const RadialFunctionG& VlocRadial,
	double Z, const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball)
{	GpuLaunchConfigHalf3D glc(updateLocal_kernel, S);
	for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
		updateLocal_kernel<<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, GGT, Vlocps, rhoIon, nChargeball,
			nCore, tauCore, nAtoms, atpos, invVol, VlocRadial,
			Z, nCoreRadial, tauCoreRadial, Zchargeball, wChargeball);
	gpuErrorCheck();
}

//Propagate gradients w.r.t pseudopotetial, ionic charge, chargeball and partial cores to structure factor
__global__
void gradLocalToSG_kernel(int zBlock, const vector3<int> S, const matrix3<> GGT,
	const complex* ccgrad_Vlocps, const complex* ccgrad_rhoIon, const complex* ccgrad_nChargeball,
	const complex* ccgrad_nCore, const complex* ccgrad_tauCore, complex* ccgrad_SG,
	const RadialFunctionG VlocRadial, double Z,
	const RadialFunctionG nCoreRadial, const RadialFunctionG tauCoreRadial,
	double Zchargeball, double wChargeball)
{
	COMPUTE_halfGindices
	gradLocalToSG_calc(i, iG, GGT, ccgrad_Vlocps, ccgrad_rhoIon, ccgrad_nChargeball,
		ccgrad_nCore, ccgrad_tauCore, ccgrad_SG, VlocRadial, Z,
		nCoreRadial, tauCoreRadial, Zchargeball, wChargeball);
}
void gradLocalToSG_gpu(const vector3<int> S, const matrix3<> GGT,
	const complex* ccgrad_Vlocps, const complex* ccgrad_rhoIon, const complex* ccgrad_nChargeball,
	const complex* ccgrad_nCore, const complex* ccgrad_tauCore, complex* ccgrad_SG,
	const RadialFunctionG& VlocRadial, double Z,
	const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball)
{	GpuLaunchConfigHalf3D glc(gradLocalToSG_kernel, S);
	for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
		gradLocalToSG_kernel<<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, GGT,
			ccgrad_Vlocps, ccgrad_rhoIon, ccgrad_nChargeball,
			ccgrad_nCore, ccgrad_tauCore, ccgrad_SG, VlocRadial, Z,
			nCoreRadial, tauCoreRadial, Zchargeball, wChargeball);
	gpuErrorCheck();
}

//Calculate forces from gradient w.r.t structure factor
__global__
void gradSGtoAtpos_kernel(int zBlock, const vector3<int> S, const vector3<> atpos,
	const complex* ccgrad_SG, vector3<complex*> grad_atpos)
{
	COMPUTE_halfGindices
	gradSGtoAtpos_calc(i, iG, atpos, ccgrad_SG, grad_atpos);
}
void gradSGtoAtpos_gpu(const vector3<int> S, const vector3<> atpos,
	const complex* ccgrad_SG, vector3<complex*> grad_atpos)
{	GpuLaunchConfigHalf3D glc(gradSGtoAtpos_kernel, S);
	for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
		gradSGtoAtpos_kernel<<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, atpos, ccgrad_SG, grad_atpos);
	gpuErrorCheck();
}
