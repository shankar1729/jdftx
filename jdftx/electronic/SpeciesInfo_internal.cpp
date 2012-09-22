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


//Initialize radial augmentation function for a pair product of 2 Ylm's, projected to a given total l
template<int l1,int m1, int l2,int m2, int l>
void Qr_sub(size_t iStart, size_t iStop, const vector3<int> S, const matrix3<>& G,  const RadialFunctionG& Qradial,
	const vector3<>& atpos, complex* Q, const complex* ccgrad_Q, vector3<complex*> grad_atpos, double gradScale)
{
	THREAD_halfGspaceLoop( (Qr_calc<l1,m1,l2,m2,l>)(i, iG, G, Qradial, atpos, Q, ccgrad_Q, grad_atpos, gradScale); )
}
template<int l1,int m1, int l2,int m2, int l>
void Qr(const vector3<int> S, const matrix3<>& G,  const RadialFunctionG& Qradial,
	const vector3<>& atpos, complex* Q, const complex* ccgrad_Q, vector3<complex*> grad_atpos, double gradScale, bool& nonZero)
{
	threadLaunch(Qr_sub<l1,m1,l2,m2,l>, S[0]*S[1]*(S[2]/2+1), S, G, Qradial, atpos, Q, ccgrad_Q, grad_atpos, gradScale);
	nonZero = true;
}
bool Qr(int l1, int m1, int l2, int m2, int l,
	const vector3<int> S, const matrix3<>& G, const RadialFunctionG& Qradial,
	const vector3<>& atpos, complex* Q, const complex* ccgrad_Q, vector3<complex*> grad_atpos, double gradScale)
{	
	bool nonZero = false;
	SwitchTemplate_lmPair(l1,m1, l2,m2, l, Qr, (S, G, Qradial, atpos, Q, ccgrad_Q, grad_atpos, gradScale, nonZero) )
	return nonZero;
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
	const complex* ccgrad_nCore, const complex* ccgrad_tauCore, complex* grad_SG,
	const RadialFunctionG& VlocRadial, double Z,
	const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball)
{	THREAD_halfGspaceLoop(
		gradLocalToSG_calc(i, iG, GGT,
		ccgrad_Vlocps, ccgrad_rhoIon, ccgrad_nChargeball,
		ccgrad_nCore, ccgrad_tauCore, grad_SG, VlocRadial,
		Z, nCoreRadial, tauCoreRadial, Zchargeball, wChargeball); )
}
void gradLocalToSG(const vector3<int> S, const matrix3<> GGT,
	const complex* ccgrad_Vlocps, const complex* ccgrad_rhoIon, const complex* ccgrad_nChargeball,
	const complex* ccgrad_nCore, const complex* ccgrad_tauCore, complex* grad_SG,
	const RadialFunctionG& VlocRadial, double Z,
	const RadialFunctionG& nCoreRadial, const RadialFunctionG& tauCoreRadial,
	double Zchargeball, double wChargeball)
{	threadLaunch(gradLocalToSG_sub, S[0]*S[1]*(S[2]/2+1), S, GGT,
		ccgrad_Vlocps, ccgrad_rhoIon, ccgrad_nChargeball,
		ccgrad_nCore, ccgrad_tauCore, grad_SG, VlocRadial,
		Z, nCoreRadial, tauCoreRadial, Zchargeball, wChargeball);
}

//Gradient w.r.t structure factor -> gradient w.r.t atom positions
void gradSGtoAtpos_sub(size_t iStart, size_t iStop, const vector3<int> S, const vector3<> atpos,
	const complex* grad_SG, vector3<complex*> grad_atpos)
{	THREAD_halfGspaceLoop( gradSGtoAtpos_calc(i, iG, atpos, grad_SG, grad_atpos); )
}
void gradSGtoAtpos(const vector3<int> S, const vector3<> atpos,
	const complex* grad_SG, vector3<complex*> grad_atpos)
{	threadLaunch(gradSGtoAtpos_sub, S[0]*S[1]*(S[2]/2+1), S, atpos, grad_SG, grad_atpos);
}
