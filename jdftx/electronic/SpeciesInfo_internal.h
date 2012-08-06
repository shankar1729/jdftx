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

//! Compute Vnl and optionally its gradients for a subset of the basis space, and for multiple atomic positions
template<int l, int m> __hostanddev__
void Vnl_calc(int n, int atomStride, int nAtoms, const vector3<>& k, const vector3<int>* iGarr, const matrix3<>& G,
	const vector3<>* pos, const RadialFunctionG& VnlRadial, complex* Vnl, bool computeGrad, vector3<complex*> dV)
{
	vector3<> kpG = k + iGarr[n]; //k+G in reciprocal lattice coordinates:
	vector3<> qvec = kpG * G; //k+G in cartesian coordinates
	double q = qvec.length();
	vector3<> qhat = qvec * (q ? 1.0/q : 0.0); //the unit vector along qvec (set qhat to 0 for q=0 (doesn't matter))
	//Compute the prefactor to the sturtcure factor in Vnl:
	double prefac = Ylm<l,m>(qhat) * VnlRadial(q);
	//Loop over columns (multiple atoms at same l,m):
	for(int atom=0; atom<nAtoms; atom++)
	{	//Multiply above prefactor by the structure factor for current atom
		complex temp = prefac * cis((-2*M_PI)*dot(pos[atom],kpG));
		Vnl[atom*atomStride+n] = temp;
		//Also set the gradients if requested
		if(computeGrad)
		{	temp *= complex(0,-2*M_PI);
			storeVector(temp*kpG, dV, atom*atomStride+n);
		}
	}
}
void Vnl(int nbasis, int atomStride, int nAtoms, int l, int m, const vector3<> k, const vector3<int>* iGarr, const matrix3<> G,
	const vector3<>* pos, const RadialFunctionG& VnlRadial, complex* Vnl, bool computeGrad, vector3<complex*> dV);
#ifdef GPU_ENABLED
void Vnl_gpu(int nbasis, int atomStride, int nAtoms, int l, int m, const vector3<> k, const vector3<int>* iGarr, const matrix3<> G,
	const vector3<>* pos, const RadialFunctionG& VnlRadial, complex* Vnl, bool computeGrad, vector3<complex*> dV);
#endif


//! Compute the Q radial function (product of two spherical basis functions projected to a given total l)
//! and propagate gradient w.r.t to it to that w.r.t the atom position (accumulate)
template<int l1,int m1, int l2,int m2, int l> __hostanddev__
void Qr_calc(int i, const vector3<int>& iG, const matrix3<>& G,  const RadialFunctionG& Qradial,
	const vector3<>& atpos, complex* Q, const complex* ccgrad_Q, vector3<complex*> grad_atpos, double gradScale)
{
	vector3<> qvec = iG * G; //to cartesian coordinates
	double q = qvec.length();
	vector3<> qhat = qvec * (q ? 1.0/q : 0.0); //the unit vector along qvec (set qhat to 0 for q=0 (doesn't matter))
	//Contribution at a given G-vector:
	complex contrib = YlmProd<l1,m1,l2,m2,l>(qhat) * Qradial(q) * cis((-2*M_PI)*dot(atpos,iG) - 0.5*M_PI*l);
	Q[i] = contrib;
	//Optional gradient accumulation:
	if(ccgrad_Q)
		accumVector((gradScale * ccgrad_Q[i].conj() * contrib * complex(0,-2*M_PI)) * iG, grad_atpos, i);
}
bool Qr(int l1, int m1, int l2, int m2, int l,
	const vector3<int> S, const matrix3<>& G, const RadialFunctionG& Qradial,
	const vector3<>& atpos, complex* Q, const complex* ccgrad_Q, vector3<complex*> grad_atpos, double gradScale=1.);
#ifdef GPU_ENABLED
bool Qr_gpu(int l1, int m1, int l2, int m2, int l,
	const vector3<int> S, const matrix3<>& G, const RadialFunctionG& Qradial,
	const vector3<>& atpos, complex* Q, const complex* ccgrad_Q, vector3<complex*> grad_atpos, double gradScale=1.);
#endif


//! Calculate local pseudopotential, ionic density and chargeball due to one species at a given G-vector
__hostanddev__ void updateLocal_calc(int i, const vector3<int>& iG, const matrix3<>& GGT,
	complex *Vlocps, complex *rhoIon, complex *nChargeball, complex* nCore, complex* tauCore,
	int nAtoms, const vector3<>* atpos, double invVol, const RadialFunctionG& VlocRadial,
	double Z, double ionWidth, const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball)
{
	double Gsq = GGT.metric_length_squared(iG);

	//Compute structure factor (scaled by 1/detR):
	complex SGinvVol = complex(0,0);
	for(int atom=0; atom<nAtoms; atom++)
		SGinvVol += cis(-2*M_PI*dot(iG,atpos[atom]));
	SGinvVol *= invVol;

	//Local potential (short ranged part in the radial function - Z/r):
	Vlocps[i] += SGinvVol * ( VlocRadial(sqrt(Gsq)) - Z * (i ? 4*M_PI/Gsq : 2*M_PI*pow(ionWidth,2)) );

	//Nuclear charge (optionally widened to a gaussian (if ionWidth!=0)):
	rhoIon[i] += SGinvVol * (-Z) * exp(-0.5*Gsq*pow(ionWidth,2));

	//Chargeball:
	if(nChargeball)
		nChargeball[i] += SGinvVol * Zchargeball
			* pow(sqrt(2*M_PI)*wChargeball,3)*exp(-0.5*Gsq*pow(wChargeball,2));

	//Partial core:
	if(nCore) nCore[i] += SGinvVol * nCoreRadial(sqrt(Gsq));
	if(tauCore) tauCore[i] += SGinvVol * tauCoreRadial(sqrt(Gsq));
}
void updateLocal(const vector3<int> S, const matrix3<> GGT,
	complex *Vlocps,  complex *rhoIon, complex *n_chargeball, complex* n_core, complex* tauCore,
	int nAtoms, const vector3<>* atpos, double invVol, const RadialFunctionG& VlocRadial,
	double Z, double ionWidth, const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball);
#ifdef GPU_ENABLED
void updateLocal_gpu(const vector3<int> S, const matrix3<> GGT,
	complex *Vlocps,  complex *rhoIon, complex *n_chargeball, complex* n_core, complex* tauCore,
	int nAtoms, const vector3<>* atpos, double invVol, const RadialFunctionG& VlocRadial,
	double Z, double ionWidth, const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball);
#endif


//! Propagate (complex conjugates of) gradients w.r.t Vlocps, rhoIon etc to gradient w.r.t SG (the strutcure factor)
__hostanddev__ void gradLocalToSG_calc(int i, const vector3<int> iG, const matrix3<> GGT,
	const complex* ccgrad_Vlocps, const complex* ccgrad_rhoIon, const complex* ccgrad_nChargeball,
	const complex* ccgrad_nCore, const complex* ccgrad_tauCore, complex* grad_SG,
	const RadialFunctionG& VlocRadial, double Z, double ionWidth,
	const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball)
{
	double Gsq = GGT.metric_length_squared(iG);
	complex ccgrad_SGinvVol(0,0); //result for this G value (gradient w.r.t structure factor/volume)

	//Local potential (short ranged part in the radial function - Z/r):
	ccgrad_SGinvVol += ccgrad_Vlocps[i] * ( VlocRadial(sqrt(Gsq)) - Z * (i ? 4*M_PI/Gsq : 0) );

	//Nuclear charge (optionally widened to a gaussian (if ionWidth!=0)):
	if(ccgrad_rhoIon)
		ccgrad_SGinvVol += ccgrad_rhoIon[i] * (-Z) * exp(-0.5*Gsq*pow(ionWidth,2));

	//Chargeball:
	if(ccgrad_nChargeball)
		ccgrad_SGinvVol += ccgrad_nChargeball[i] * Zchargeball
			* pow(sqrt(2*M_PI)*wChargeball,3)*exp(-0.5*Gsq*pow(wChargeball,2));

	//Partial core:
	if(ccgrad_nCore) ccgrad_SGinvVol += ccgrad_nCore[i] * nCoreRadial(sqrt(Gsq));
	if(ccgrad_tauCore) ccgrad_SGinvVol += ccgrad_tauCore[i] * tauCoreRadial(sqrt(Gsq));
	
	//Store result:
	grad_SG[i] = ccgrad_SGinvVol.conj();
}
void gradLocalToSG(const vector3<int> S, const matrix3<> GGT,
	const complex* ccgrad_Vlocps, const complex* ccgrad_rhoIon, const complex* ccgrad_nChargeball,
	const complex* ccgrad_nCore, const complex* ccgrad_tauCore, complex* grad_SG,
	const RadialFunctionG& VlocRadial, double Z, double ionWidth,
	const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball);
#ifdef GPU_ENABLED
void gradLocalToSG_gpu(const vector3<int> S, const matrix3<> GGT,
	const complex* ccgrad_Vlocps, const complex* ccgrad_rhoIon, const complex* ccgrad_nChargeball,
	const complex* ccgrad_nCore, const complex* ccgrad_tauCore, complex* grad_SG,
	const RadialFunctionG& VlocRadial, double Z, double ionWidth,
	const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball);
#endif


//! Propagate the gradient form the structure factor to the given atomic position
//! this is still per G-vector, need to sum grad_atpos over G to get the force on that atom
__hostanddev__ void gradSGtoAtpos_calc(int i, const vector3<int> iG, const vector3<> atpos,
	const complex* grad_SG, vector3<complex*> grad_atpos)
{
	complex term = complex(0,-2*M_PI) * cis(-2*M_PI*dot(iG,atpos)) * grad_SG[i];
	storeVector(iG * term, grad_atpos, i);
}
void gradSGtoAtpos(const vector3<int> S, const vector3<> atpos,
	const complex* grad_SG, vector3<complex*> grad_atpos);
#ifdef GPU_ENABLED
void gradSGtoAtpos_gpu(const vector3<int> S, const vector3<> atpos,
	const complex* grad_SG, vector3<complex*> grad_atpos);
#endif


#endif // JDFTX_ELECTRONIC_SPECIESINFO_INTERNAL_H
