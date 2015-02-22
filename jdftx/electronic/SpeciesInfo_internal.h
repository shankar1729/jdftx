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

#ifndef JDFTX_ELECTRONIC_SPECIESINFO_INTERNAL_H
#define JDFTX_ELECTRONIC_SPECIESINFO_INTERNAL_H

//!@file SpeciesInfo_internal.h Shared GPU/CPU code for ion/pseudopotential related calculations

#include <core/matrix3.h>
#include <electronic/RadialFunction.h>
#include <electronic/SphericalHarmonics.h>
#include <stdint.h>

//! Compute Vnl and optionally its gradients for a subset of the basis space, and for multiple atomic positions
template<int l, int m> __hostanddev__
void Vnl_calc(int n, int atomStride, int nAtoms, const vector3<>& k, const vector3<int>* iGarr,
	const matrix3<>& G, const vector3<>* pos, const RadialFunctionG& VnlRadial, complex* Vnl)
{
	vector3<> kpG = k + iGarr[n]; //k+G in reciprocal lattice coordinates:
	vector3<> qvec = kpG * G; //k+G in cartesian coordinates
	double q = qvec.length();
	vector3<> qhat = qvec * (q ? 1.0/q : 0.0); //the unit vector along qvec (set qhat to 0 for q=0 (doesn't matter))
	double prefac = Ylm<l,m>(qhat) * VnlRadial(q); //prefactor to structure factor
	//Loop over columns (multiple atoms at same l,m):
	for(int atom=0; atom<nAtoms; atom++)
		Vnl[atom*atomStride+n] = prefac * cis((-2*M_PI)*dot(pos[atom],kpG));
}
void Vnl(int nbasis, int atomStride, int nAtoms, int l, int m, const vector3<> k, const vector3<int>* iGarr,
	const matrix3<> G, const vector3<>* pos, const RadialFunctionG& VnlRadial, complex* Vnl);
#ifdef GPU_ENABLED
void Vnl_gpu(int nbasis, int atomStride, int nAtoms, int l, int m, const vector3<> k, const vector3<int>* iGarr,
	const matrix3<> G, const vector3<>* pos, const RadialFunctionG& VnlRadial, complex* Vnl);
#endif


//! Perform the loop:
//!   for(lm=0; lm < Nlm; lm++) (*f)(tag< lm >);
//! at compile time using templates. Note with lm := l*(l+1)+m to loop over spherical harmonics, Nlm = (lMax+1)^2.
//! The functor f is templated over lm by using the tag class StaticLoopYlmTag< lm >.
template<int lm> struct StaticLoopYlmTag{};
template<int Nlm, typename Functor, int lmInv> struct StaticLoopYlm
{	__hostanddev__ static void exec(Functor* f)
	{	(*f)(StaticLoopYlmTag<Nlm-lmInv>());
		StaticLoopYlm< Nlm, Functor, lmInv-1 >::exec(f);
	}
};
template<int Nlm, typename Functor> struct StaticLoopYlm< Nlm, Functor, 0>
{	__hostanddev__ static void exec(Functor* f) { return; }
};
template<int Nlm, typename Functor> __hostanddev__ void staticLoopYlm(Functor* f)
{	StaticLoopYlm< Nlm, Functor, Nlm >::exec(f);
}

#define SwitchTemplate_Nlm(Nlm, func, args) \
	switch(Nlm) \
	{	case 1:  func<1>  args; break; \
		case 4:  func<4>  args; break; \
		case 9:  func<9>  args; break; \
		case 16: func<16> args; break; \
		case 25: func<25> args; break; \
		case 49: func<49> args; break; \
		default: fprintf(stderr, "Invalid Nlm in SwitchTemplate_Nlm"); exit(1); \
	}


//! Augment electron density by spherical functions (radial functions multiplied by spherical harmonics)
//! and propagate gradient w.r.t to it to that w.r.t the atom position (accumulate)
//! (In MPI mode, each process only collects contributions for a subset of G-vectors)
struct nAugmentFunctor
{	vector3<> qhat; double q;
	int nCoeff; double dGinv; const double* nRadial;
	complex n;
	
	__hostanddev__ nAugmentFunctor(const vector3<>& qvec, int nCoeff, double dGinv, const double* nRadial)
	: nCoeff(nCoeff), dGinv(dGinv), nRadial(nRadial)
	{	q = qvec.length();
		qhat = qvec * (q ? 1.0/q : 0.0); //the unit vector along qvec (set qhat to 0 for q=0 (doesn't matter))
	}
	
	template<int lm> __hostanddev__ void operator()(const StaticLoopYlmTag<lm>&)
	{	//Compute phase (-i)^l:
		complex mIota(0,-1), phase(1,0);
		for(int l=0; l*(l+2) < lm; l++) phase *= mIota;
		//Accumulate result:
		double Gindex = q * dGinv;
		if(Gindex < nCoeff-5)
			n += phase * Ylm<lm>(qhat) * QuinticSpline::value(nRadial+lm*nCoeff, Gindex); 
	}
};
template<int Nlm> __hostanddev__
void nAugment_calc(int i, const vector3<int>& iG, const matrix3<>& G,
	int nCoeff, double dGinv, const double* nRadial, const vector3<>& atpos, complex* n)
{
	nAugmentFunctor functor(iG*G, nCoeff, dGinv, nRadial);
	staticLoopYlm<Nlm>(&functor);
	n[i] += functor.n * cis((-2*M_PI)*dot(atpos,iG));
}
void nAugment(int Nlm,
	const vector3<int> S, const matrix3<>& G, int iGstart, int iGstop,
	int nCoeff, double dGinv, const double* nRadial, const vector3<>& atpos, complex* n);
#ifdef GPU_ENABLED
void nAugment_gpu(int Nlm,
	const vector3<int> S, const matrix3<>& G, int iGstart, int iGstop,
	int nCoeff, double dGinv, const double* nRadial, const vector3<>& atpos, complex* n);
#endif

//Function for initializing the index arrays used by nAugmentGrad
//! (In MPI mode, only a subset of G-vectors are indexed on each process (to correspond to nAUgment))
__hostanddev__ uint64_t setNagIndex_calc(const vector3<int>& iG, const vector3<int>& S, const matrix3<>& G, double dGinv)
{	uint64_t Gindex = uint64_t((iG*G).length() * dGinv);
	vector3<int> iv = iG; for(int k=0; k<3; k++) if(iv[k]<0) iv[k] += S[k];
	return (Gindex << 48) //Putting Gindex in the higher word allows sorting by it first, and then by grid point index
		+ (uint64_t(iv[0]) << 32) + (uint64_t(iv[1]) << 16) + uint64_t(iv[2]);
}
__hostanddev__ void setNagIndexPtr_calc(int i, int iMax, int nCoeff, const uint64_t* nagIndex, size_t* nagIndexPtr)
{	int Gindex = int(nagIndex[i] >> 48);
	int GindexNext = (i+1<iMax) ? int(nagIndex[i+1] >> 48) : nCoeff;
	if(i==0) for(int j=0; j<=Gindex; j++) nagIndexPtr[j] = 0;
	for(int j=Gindex; j<GindexNext; j++)
		nagIndexPtr[j+1] = i+1;
}
void setNagIndex(const vector3<int>& S, const matrix3<>& G, int iGstart, int iGstop, int nCoeff, double dGinv, uint64_t*& nagIndex, size_t*& nagIndexPtr);
#ifdef GPU_ENABLED
void setNagIndex_gpu(const vector3<int>& S, const matrix3<>& G, int iGstart, int iGstop, int nCoeff, double dGinv, uint64_t*& nagIndex, size_t*& nagIndexPtr);
#endif

//Gradient propragation corresponding to nAugment:
//(The MPI division happens implicitly here, because nagIndex is limited to each process's share (see above))
struct nAugmentGradFunctor
{	vector3<> qhat; double q;
	int nCoeff; double dGinv; const double* nRadial;
	complex E_n, nE_n;
	double* E_nRadial;
	int dotPrefac; //prefactor in dot-product (1 or 2 for each reciprocal space point, because of real symmetry)

	__hostanddev__ nAugmentGradFunctor(const vector3<>& qvec, int nCoeff, double dGinv, const double* nRadial, const complex& E_n, double* E_nRadial, int dotPrefac)
	: nCoeff(nCoeff), dGinv(dGinv), nRadial(nRadial), E_n(E_n), E_nRadial(E_nRadial), dotPrefac(dotPrefac)
	{	q = qvec.length();
		qhat = qvec * (q ? 1.0/q : 0.0); //the unit vector along qvec (set qhat to 0 for q=0 (doesn't matter))
	}
	
	template<int lm> __hostanddev__ void operator()(const StaticLoopYlmTag<lm>&)
	{	//Compute phase (-i)^l:
		complex mIota(0,-1), phase(1,0);
		for(int l=0; l*(l+2) < lm; l++) phase *= mIota;
		//Accumulate result:
		double Gindex = q * dGinv;
		if(Gindex < nCoeff-5)
		{	complex term = phase * Ylm<lm>(qhat) * E_n;
			QuinticSpline::valueGrad(dotPrefac * term.real(), E_nRadial+lm*nCoeff, Gindex);
			if(nRadial) nE_n += term * QuinticSpline::value(nRadial+lm*nCoeff, Gindex); //needed again only when computing forces
		}
	}
};
template<int Nlm> __hostanddev__
void nAugmentGrad_calc(uint64_t key, const vector3<int>& S, const matrix3<>& G,
	int nCoeff, double dGinv, const double* nRadial, const vector3<>& atpos,
	const complex* ccE_n, double* E_nRadial, vector3<complex*> E_atpos, bool dummyGpuThread=false)
{
	//Obtain 3D index iG and array offset i for this point (similar to COMPUTE_halfGindices)
	vector3<int> iG;
	iG[2] = int(0xFFFF & key); key >>= 16;
	iG[1] = int(0xFFFF & key); key >>= 16;
	iG[0] = int(0xFFFF & key);
	size_t i = iG[2] + (S[2]/2+1)*size_t(iG[1] + S[1]*iG[0]);
	for(int j=0; j<3; j++) if(2*iG[j]>S[j]) iG[j]-=S[j];
	int dotPrefac = (iG[2]==0||2*iG[2]==S[2]) ? 1 : 2;
	
	nAugmentGradFunctor functor(iG*G, nCoeff, dGinv, nRadial, dummyGpuThread ? complex() : ccE_n[i].conj() * cis((-2*M_PI)*dot(atpos,iG)), E_nRadial, dotPrefac);
	staticLoopYlm<Nlm>(&functor);
	if(nRadial && !dummyGpuThread) accumVector((functor.nE_n * complex(0,-2*M_PI)) * iG, E_atpos, i);
}
void nAugmentGrad(int Nlm, const vector3<int> S, const matrix3<>& G,
	int nCoeff, double dGinv, const double* nRadial, const vector3<>& atpos,
	const complex* ccE_n, double* E_nRadial, vector3<complex*> E_atpos,
	const uint64_t* nagIndex, const size_t* nagIndexPtr);
#ifdef GPU_ENABLED
void nAugmentGrad_gpu(int Nlm, const vector3<int> S, const matrix3<>& G,
	int nCoeff, double dGinv, const double* nRadial, const vector3<>& atpos,
	const complex* ccE_n, double* E_nRadial, vector3<complex*> E_atpos,
	const uint64_t* nagIndex, const size_t* nagIndexPtr);
#endif


//!Get structure factor for a specific iG, given a list of atoms
__hostanddev__ complex getSG_calc(const vector3<int>& iG, const int& nAtoms, const vector3<>* atpos)
{	complex SG = complex(0,0);
	for(int atom=0; atom<nAtoms; atom++)
		SG += cis(-2*M_PI*dot(iG,atpos[atom]));
	return SG;
}
//!Get structure factor in a ScalarFieldTilde's data/dataGpu (with 1/vol normalization factor)
void getSG(const vector3<int> S, int nAtoms, const vector3<>* atpos, double invVol, complex* SG);
#ifdef GPU_ENABLED
void getSG_gpu(const vector3<int> S, int nAtoms, const vector3<>* atpos, double invVol, complex* SG);
#endif

//! Calculate local pseudopotential, ionic density and chargeball due to one species at a given G-vector
__hostanddev__ void updateLocal_calc(int i, const vector3<int>& iG, const matrix3<>& GGT,
	complex *Vlocps, complex *rhoIon, complex *nChargeball, complex* nCore, complex* tauCore,
	int nAtoms, const vector3<>* atpos, double invVol, const RadialFunctionG& VlocRadial,
	double Z, const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball)
{
	double Gsq = GGT.metric_length_squared(iG);

	//Compute structure factor (scaled by 1/detR):
	complex SGinvVol = getSG_calc(iG, nAtoms, atpos) * invVol;

	//Short-ranged part of Local potential (long-ranged part added on later in IonInfo.cpp):
	Vlocps[i] += SGinvVol * VlocRadial(sqrt(Gsq));

	//Nuclear charge (optionally widened to a gaussian later in IonInfo.cpp):
	rhoIon[i] += SGinvVol * (-Z);

	//Chargeball:
	if(nChargeball)
		nChargeball[i] += SGinvVol * Zchargeball * exp(-0.5*Gsq*pow(wChargeball,2));

	//Partial core:
	if(nCore) nCore[i] += SGinvVol * nCoreRadial(sqrt(Gsq));
	if(tauCore) tauCore[i] += SGinvVol * tauCoreRadial(sqrt(Gsq));
}
void updateLocal(const vector3<int> S, const matrix3<> GGT,
	complex *Vlocps,  complex *rhoIon, complex *n_chargeball, complex* n_core, complex* tauCore,
	int nAtoms, const vector3<>* atpos, double invVol, const RadialFunctionG& VlocRadial,
	double Z, const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball);
#ifdef GPU_ENABLED
void updateLocal_gpu(const vector3<int> S, const matrix3<> GGT,
	complex *Vlocps,  complex *rhoIon, complex *n_chargeball, complex* n_core, complex* tauCore,
	int nAtoms, const vector3<>* atpos, double invVol, const RadialFunctionG& VlocRadial,
	double Z, const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball);
#endif


//! Propagate (complex conjugates of) gradients w.r.t Vlocps, rhoIon etc to complex conjugate gradient w.r.t SG (the strutcure factor)
__hostanddev__ void gradLocalToSG_calc(int i, const vector3<int> iG, const matrix3<> GGT,
	const complex* ccgrad_Vlocps, const complex* ccgrad_rhoIon, const complex* ccgrad_nChargeball,
	const complex* ccgrad_nCore, const complex* ccgrad_tauCore, complex* ccgrad_SG,
	const RadialFunctionG& VlocRadial, double Z,
	const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball)
{
	double Gsq = GGT.metric_length_squared(iG);
	complex ccgrad_SGinvVol(0,0); //result for this G value (gradient w.r.t structure factor/volume)

	//Local potential (short ranged part in the radial function - Z/r):
	ccgrad_SGinvVol += ccgrad_Vlocps[i] * VlocRadial(sqrt(Gsq));

	//Nuclear charge
	if(ccgrad_rhoIon)
		ccgrad_SGinvVol += ccgrad_rhoIon[i] * (-Z);

	//Chargeball:
	if(ccgrad_nChargeball)
		ccgrad_SGinvVol += ccgrad_nChargeball[i] * Zchargeball * exp(-0.5*Gsq*pow(wChargeball,2));

	//Partial core:
	if(ccgrad_nCore) ccgrad_SGinvVol += ccgrad_nCore[i] * nCoreRadial(sqrt(Gsq));
	if(ccgrad_tauCore) ccgrad_SGinvVol += ccgrad_tauCore[i] * tauCoreRadial(sqrt(Gsq));
	
	//Store result:
	ccgrad_SG[i] = ccgrad_SGinvVol;
}
void gradLocalToSG(const vector3<int> S, const matrix3<> GGT,
	const complex* ccgrad_Vlocps, const complex* ccgrad_rhoIon, const complex* ccgrad_nChargeball,
	const complex* ccgrad_nCore, const complex* ccgrad_tauCore, complex* ccgrad_SG,
	const RadialFunctionG& VlocRadial, double Z,
	const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball);
#ifdef GPU_ENABLED
void gradLocalToSG_gpu(const vector3<int> S, const matrix3<> GGT,
	const complex* ccgrad_Vlocps, const complex* ccgrad_rhoIon, const complex* ccgrad_nChargeball,
	const complex* ccgrad_nCore, const complex* ccgrad_tauCore, complex* ccgrad_SG,
	const RadialFunctionG& VlocRadial, double Z,
	const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball);
#endif


//! Propagate the complex conjugate gradient w.r.t the structure factor to the given atomic position
//! this is still per G-vector, need to sum grad_atpos over G to get the force on that atom
__hostanddev__ void gradSGtoAtpos_calc(int i, const vector3<int> iG, const vector3<> atpos,
	const complex* ccgrad_SG, vector3<complex*> grad_atpos)
{
	complex term = complex(0,-2*M_PI) * cis(-2*M_PI*dot(iG,atpos)) * ccgrad_SG[i].conj();
	storeVector(iG * term, grad_atpos, i);
}
void gradSGtoAtpos(const vector3<int> S, const vector3<> atpos,
	const complex* ccgrad_SG, vector3<complex*> grad_atpos);
#ifdef GPU_ENABLED
void gradSGtoAtpos_gpu(const vector3<int> S, const vector3<> atpos,
	const complex* ccgrad_SG, vector3<complex*> grad_atpos);
#endif


#endif // JDFTX_ELECTRONIC_SPECIESINFO_INTERNAL_H
