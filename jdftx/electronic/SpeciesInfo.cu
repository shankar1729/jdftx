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

//Calculate non-local pseudopotential projector
template<int l, int m> __global__
void Vnl_kernel(int nbasis, int atomStride, int nAtoms, vector3<> k, const vector3<int>* iGarr, const matrix3<> G,
	const vector3<>* pos, const RadialFunctionG VnlRadial, complex* V, bool computeGrad, vector3<complex*> dV)
{	int n = kernelIndex1D();
	if(n<nbasis) Vnl_calc<l,m>(n, atomStride, nAtoms, k, iGarr, G, pos, VnlRadial, V, computeGrad, dV);
}
template<int l, int m>
void Vnl_gpu(int nbasis, int atomStride, int nAtoms, vector3<> k, const vector3<int>* iGarr, const matrix3<> G,
	const vector3<>* pos, const RadialFunctionG& VnlRadial, complex* V, bool computeGrad, vector3<complex*> dV)
{	GpuLaunchConfig1D glc(Vnl_kernel<l,m>, nbasis);
	Vnl_kernel<l,m><<<glc.nBlocks,glc.nPerBlock>>>(nbasis, atomStride, nAtoms, k, iGarr, G, pos, VnlRadial, V, computeGrad, dV);
	gpuErrorCheck();
}
void Vnl_gpu(int nbasis, int atomStride, int nAtoms, int l, int m, vector3<> k, const vector3<int>* iGarr, const matrix3<> G,
	const vector3<>* pos, const RadialFunctionG& VnlRadial, complex* V, bool computeGrad, vector3<complex*> dV)
{
	SwitchTemplate_lm(l,m, Vnl_gpu, (nbasis, atomStride, nAtoms, k, iGarr, G, pos, VnlRadial, V, computeGrad, dV) )
}


//Augment electron density by spherical functions
template<int Nlm> __global__ void nAugment_kernel(int zBlock, const vector3<int> S, const matrix3<> G,
	int nGloc, double dGinv, const double* nRadial, const vector3<> atpos, complex* n)
{	COMPUTE_halfGindices
	nAugment_calc<Nlm>(i, iG, G, nGloc, dGinv, nRadial, atpos, n);
}
template<int Nlm> void nAugment_gpu(const vector3<int> S, const matrix3<>& G,
	int nGloc, double dGinv, const double* nRadial, const vector3<>& atpos, complex* n)
{	GpuLaunchConfigHalf3D glc(nAugment_kernel<Nlm>, S);
	for(int zBlock=0; zBlock<glc.zBlockMax; zBlock++)
		nAugment_kernel<Nlm><<<glc.nBlocks,glc.nPerBlock>>>(zBlock, S, G, nGloc, dGinv, nRadial, atpos, n);
	gpuErrorCheck();
}
void nAugment_gpu(int Nlm, const vector3<int> S, const matrix3<>& G,
	int nGloc, double dGinv, const double* nRadial, const vector3<>& atpos, complex* n)
{
	SwitchTemplate_Nlm(Nlm, nAugment_gpu, (S, G, nGloc, dGinv, nRadial, atpos, n) )
}


//Propagate gradients corresponding to above electron density augmentation
template<int Nlm> __global__ void nAugmentGrad_kernel(const vector3<int> S, const matrix3<> G,
	int nGloc, double dGinv, const double* nRadial, const vector3<> atpos,
	const complex* ccE_n, double* E_nRadialTemp, vector3<complex*> E_atpos)
{
	int iThread = kernelIndex1D();
	int nThreads = blockDim.x * gridDim.x;
	double* E_nRadialThread = E_nRadialTemp + (Nlm*nGloc)*iThread;
	size_t nG = S[0]*S[1]*(S[2]/2+1);
	size_t iStart = ((iThread) * nG) / nThreads;
	size_t iStop = ((iThread+1) * nG) / nThreads;
	THREAD_halfGspaceLoop( (nAugmentGrad_calc<Nlm>)(i, iG, G, nGloc, dGinv, nRadial, atpos, ccE_n, E_nRadialThread, E_atpos, iG[2]==0||2*iG[2]==S[2] ? 1 : 2); )
}
__global__ void nAugmentGrad_collectKernel(int nData, int nCopies, const double* in, double* out)
{	int i = kernelIndex1D();
	if(i<nData)
	{	double sum = 0.;
		for(int iCopy=0; iCopy<nCopies; iCopy++)
			sum += in[i + iCopy*nData];
		out[i] += sum;
	}
}
template<int Nlm> void nAugmentGrad_gpu(const vector3<int> S, const matrix3<>& G,
	int nGloc, double dGinv, const double* nRadial, const vector3<>& atpos,
	const complex* ccE_n, double* E_nRadial, vector3<complex*> E_atpos)
{
	//Allocate temporary memory:
	cudaDeviceProp prop;
	int iDevice; cudaGetDevice(&iDevice);
	cudaGetDeviceProperties(&prop, iDevice);
	int nPerBlock = 8;
	int nBlocks = prop.multiProcessorCount * 4;
	int nThreads = nBlocks * nPerBlock;
	double* E_nRadialTemp; cudaMalloc(&E_nRadialTemp, sizeof(double)*nGloc*Nlm*nThreads);
	cudaMemset(E_nRadialTemp, 0, sizeof(double)*nGloc*Nlm*nThreads);
	gpuErrorCheck();
	//Stage 1: calculate with the scattered accumulate to per-thread temporary space
	nAugmentGrad_kernel<Nlm><<<nBlocks,nPerBlock>>>(S, G, nGloc, dGinv, nRadial, atpos, ccE_n, E_nRadialTemp, E_atpos);
	gpuErrorCheck();
	//Stage 2: collect from E_nRadialTemp to E_nRadial
	GpuLaunchConfig1D glc(nAugmentGrad_collectKernel, nGloc*Nlm);
	nAugmentGrad_collectKernel<<<glc.nBlocks,glc.nPerBlock>>>(nGloc*Nlm, nThreads, E_nRadialTemp, E_nRadial);
	gpuErrorCheck();
	//Cleanup:
	cudaFree(E_nRadialTemp);
}
void nAugmentGrad_gpu(int Nlm, const vector3<int> S, const matrix3<>& G,
	int nGloc, double dGinv, const double* nRadial, const vector3<>& atpos,
	const complex* ccE_n, double* E_nRadial, vector3<complex*> E_atpos)
{	
	SwitchTemplate_Nlm(Nlm, nAugmentGrad_gpu, (S, G, nGloc, dGinv, nRadial, atpos, ccE_n, E_nRadial, E_atpos) )
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
