/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman
Copyright 1996-2003 Sohrab Ismail-Beigi

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
#include <core/LoopMacros.h>
#include <core/Thread.h>
#include <core/BlasExtra.h>
#include <algorithm>

//Initialize non-local projector from a radial function at a particular l,m
template<int l, int m>
void Vnl(int nbasis, int atomStride, int nAtoms, const vector3<> k, const vector3<int>* iGarr, const matrix3<> G,
	const vector3<>* pos, const RadialFunctionG& VnlRadial, complex* V, bool computeGrad, vector3<complex*> dV)
{	threadedLoop(Vnl_calc<l,m>, nbasis, atomStride, nAtoms, k, iGarr, G, pos, VnlRadial, V, computeGrad, dV);
}
void Vnl(int nbasis, int atomStride, int nAtoms, int l, int m, const vector3<> k, const vector3<int>* iGarr, const matrix3<> G,
	const vector3<>* pos, const RadialFunctionG& VnlRadial, complex* V, bool computeGrad, vector3<complex*> dV)
{	SwitchTemplate_lm(l,m, Vnl, (nbasis, atomStride, nAtoms, k, iGarr, G, pos, VnlRadial, V, computeGrad, dV) )
}

//Augment electron density by spherical functions
template<int Nlm> void nAugment_sub(size_t iStart, size_t iStop, const vector3<int> S, const matrix3<>& G,
	int nCoeff, double dGinv, const double* nRadial, const vector3<>& atpos, complex* n)
{
	THREAD_halfGspaceLoop( (nAugment_calc<Nlm>)(i, iG, G, nCoeff, dGinv, nRadial, atpos, n); )
}
template<int Nlm> void nAugment(const vector3<int> S, const matrix3<>& G,
	int nCoeff, double dGinv, const double* nRadial, const vector3<>& atpos, complex* n)
{
	threadLaunch(nAugment_sub<Nlm>, S[0]*S[1]*(S[2]/2+1), S, G, nCoeff, dGinv, nRadial, atpos, n);
}
void nAugment(int Nlm, const vector3<int> S, const matrix3<>& G,
	int nCoeff, double dGinv, const double* nRadial, const vector3<>& atpos, complex* n)
{	
	SwitchTemplate_Nlm(Nlm, nAugment, (S, G, nCoeff, dGinv, nRadial, atpos, n) )
}

//Function for initializing the index arrays used by nAugmentGrad
void setNagIndex_sub(size_t iStart, size_t iStop, const vector3<int> S, const matrix3<> G, double dGinv, uint64_t* nagIndex)
{	THREAD_halfGspaceLoop( setNagIndex_calc(i, iG, S, G, dGinv, nagIndex); )
}
void setNagIndexPtr_sub(int iStart, int iStop, size_t nG, int nCoeff, const uint64_t* nagIndex, size_t* nagIndexPtr)
{	for(int i=iStart; i<iStop; i++)
		setNagIndexPtr_calc(i, nG, nCoeff, nagIndex, nagIndexPtr);
}
void setNagIndex(const vector3<int>& S, const matrix3<>& G, int nCoeff, double dGinv, uint64_t*& nagIndex, size_t*& nagIndexPtr)
{	//First initialize the indices:
	size_t nG = S[0]*S[1]*(S[2]/2+1);
	{	if(!nagIndex) nagIndex = new uint64_t[nG];
		threadLaunch(setNagIndex_sub, nG, S, G, dGinv, nagIndex);
	}
	//Now sort them to be ordered by Gindex
	std::sort(nagIndex, nagIndex+nG);
	//Finally initialize the pointers to boundaries between different Gindices:
	{	if(!nagIndexPtr) nagIndexPtr = new size_t[nCoeff+1];
		threadLaunch(setNagIndexPtr_sub, nG, nG, nCoeff, nagIndex, nagIndexPtr);
	}
}


//Propagate gradients corresponding to above electron density augmentation
template<int Nlm> void nAugmentGrad_sub(int jStart, int jStop, const vector3<int> S, const matrix3<>& G,
	int nCoeff, double dGinv, const double* nRadial, const vector3<>& atpos,
	const complex* ccE_n, double* E_nRadial, vector3<complex*> E_atpos,
	const uint64_t* nagIndex, const size_t* nagIndexPtr, int pass)
{
	for(int j=jStart; j<jStop; j++)
	{	const int iCoeff = pass+6*j;
		for(size_t ptr=nagIndexPtr[iCoeff]; ptr<nagIndexPtr[iCoeff+1]; ptr++)
			nAugmentGrad_calc<Nlm>(nagIndex[ptr], S, G, nCoeff, dGinv, nRadial, atpos, ccE_n, E_nRadial, E_atpos, false);
	}
}
template<int Nlm> void nAugmentGrad(const vector3<int> S, const matrix3<>& G,
	int nCoeff, double dGinv, const double* nRadial, const vector3<>& atpos,
	const complex* ccE_n, double* E_nRadial, vector3<complex*> E_atpos,
	const uint64_t* nagIndex, const size_t* nagIndexPtr)
{	//Accumulate in non-overlapping passes:
	for(int pass=0; pass<6; pass++)
	{	int nBlocks = (nCoeff-pass+5)/6; // ceil((nCoeff-pass)/6)
		threadLaunch(nAugmentGrad_sub<Nlm>, nBlocks, S, G, nCoeff, dGinv, nRadial, atpos, ccE_n, E_nRadial, E_atpos, nagIndex, nagIndexPtr, pass);
	}
}
void nAugmentGrad(int Nlm, const vector3<int> S, const matrix3<>& G,
	int nCoeff, double dGinv, const double* nRadial, const vector3<>& atpos,
	const complex* ccE_n, double* E_nRadial, vector3<complex*> E_atpos, const uint64_t* nagIndex, const size_t* nagIndexPtr)
{	
	SwitchTemplate_Nlm(Nlm, nAugmentGrad, (S, G, nCoeff, dGinv, nRadial, atpos, ccE_n, E_nRadial, E_atpos, nagIndex, nagIndexPtr) )
}


//Structure factor
void getSG_sub(size_t iStart, size_t iStop, const vector3<int> S,
	int nAtoms, const vector3<>* atpos, double invVol, complex* SG)
{	THREAD_halfGspaceLoop( SG[i] = invVol * getSG_calc(iG, nAtoms, atpos); )
}
void getSG(const vector3<int> S, int nAtoms, const vector3<>* atpos, double invVol, complex* SG)
{	threadLaunch(getSG_sub, S[0]*S[1]*(S[2]/2+1), S, nAtoms, atpos, invVol, SG);
}

//Local pseudopotential, ionic charge, chargeball and partial cores (CPU thread and launcher)
void updateLocal_sub(size_t iStart, size_t iStop, const vector3<int> S, const matrix3<> GGT,
	complex *Vlocps,  complex *rhoIon, complex *nChargeball, complex *nCore, complex* tauCore,
	int nAtoms, const vector3<>* atpos, double invVol, const RadialFunctionG& VlocRadial,
	double Z, const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball)
{	THREAD_halfGspaceLoop(
		updateLocal_calc(i, iG, GGT,
			Vlocps, rhoIon, nChargeball, nCore, tauCore,
			nAtoms, atpos, invVol, VlocRadial,
			Z, nCoreRadial, tauCoreRadial, Zchargeball, wChargeball); )
}
void updateLocal(const vector3<int> S, const matrix3<> GGT,
	complex *Vlocps,  complex *rhoIon, complex *nChargeball, complex *nCore, complex* tauCore,
	int nAtoms, const vector3<>* atpos, double invVol, const RadialFunctionG& VlocRadial,
	double Z, const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball)
{	threadLaunch(updateLocal_sub, S[0]*S[1]*(S[2]/2+1), S, GGT,
		Vlocps, rhoIon, nChargeball, nCore, tauCore,
		nAtoms, atpos, invVol, VlocRadial,
		Z, nCoreRadial, tauCoreRadial, Zchargeball, wChargeball);
}

//Forces due to local pseudopotential, ionic charge, chargeball and partial cores (to gradient w.r.t structure factor)
void gradLocalToSG_sub(size_t iStart, size_t iStop, const vector3<int> S, const matrix3<> GGT,
	const complex* ccgrad_Vlocps, const complex* ccgrad_rhoIon, const complex* ccgrad_nChargeball,
	const complex* ccgrad_nCore, const complex* ccgrad_tauCore, complex* ccgrad_SG,
	const RadialFunctionG& VlocRadial, double Z,
	const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball)
{	THREAD_halfGspaceLoop(
		gradLocalToSG_calc(i, iG, GGT,
		ccgrad_Vlocps, ccgrad_rhoIon, ccgrad_nChargeball,
		ccgrad_nCore, ccgrad_tauCore, ccgrad_SG, VlocRadial,
		Z, nCoreRadial, tauCoreRadial, Zchargeball, wChargeball); )
}
void gradLocalToSG(const vector3<int> S, const matrix3<> GGT,
	const complex* ccgrad_Vlocps, const complex* ccgrad_rhoIon, const complex* ccgrad_nChargeball,
	const complex* ccgrad_nCore, const complex* ccgrad_tauCore, complex* ccgrad_SG,
	const RadialFunctionG& VlocRadial, double Z,
	const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball)
{	threadLaunch(gradLocalToSG_sub, S[0]*S[1]*(S[2]/2+1), S, GGT,
		ccgrad_Vlocps, ccgrad_rhoIon, ccgrad_nChargeball,
		ccgrad_nCore, ccgrad_tauCore, ccgrad_SG, VlocRadial,
		Z, nCoreRadial, tauCoreRadial, Zchargeball, wChargeball);
}

//Gradient w.r.t structure factor -> gradient w.r.t atom positions
void gradSGtoAtpos_sub(size_t iStart, size_t iStop, const vector3<int> S, const vector3<> atpos,
	const complex* ccgrad_SG, vector3<complex*> grad_atpos)
{	THREAD_halfGspaceLoop( gradSGtoAtpos_calc(i, iG, atpos, ccgrad_SG, grad_atpos); )
}
void gradSGtoAtpos(const vector3<int> S, const vector3<> atpos,
	const complex* ccgrad_SG, vector3<complex*> grad_atpos)
{	threadLaunch(gradSGtoAtpos_sub, S[0]*S[1]*(S[2]/2+1), S, atpos, ccgrad_SG, grad_atpos);
}
