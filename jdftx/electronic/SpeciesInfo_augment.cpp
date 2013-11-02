/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman

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

#include <electronic/SpeciesInfo.h>
#include <electronic/SpeciesInfo_internal.h>
#include <electronic/Everything.h>
#include <electronic/matrix.h>
#include <electronic/operators.h>
#include <electronic/ColumnBundle.h>

//------- additional SpeciesInfo functions for ultrasoft pseudopotentials (density and overlap augmentation) -------


void SpeciesInfo::augmentOverlap(const ColumnBundle& Cq, ColumnBundle& OCq) const
{	static StopWatch watch("augmentOverlap"); watch.start();
	if(!atpos.size()) return; //unused species
	if(!Qint.size()) return; //no overlap augmentation
	const GridInfo &gInfo = e->gInfo;
	const Basis& basis = *Cq.basis;
	
	for(int l=0; l<int(VnlRadial.size()); l++)
	{	unsigned nProj = VnlRadial[l].size(); if(!nProj) continue; //skip l if no projectors
		//Copy Qint into a block-diagonal form for all atoms:
		tiledBlockMatrix Q(this->Qint[l], atpos.size());
		//Allocate temporaries:
		ColumnBundle V = Cq.similar(nProj * atpos.size());
		for(int m=-l; m<=l; m++)
		{	// Calculate the nonlocal projectors:
			for(unsigned p=0; p<nProj; p++)
			{	size_t offs = p * basis.nbasis;
				size_t atomStride = nProj * basis.nbasis;
				callPref(Vnl)(basis.nbasis, atomStride, atpos.size(), l, m, V.qnum->k, basis.iGarrPref, gInfo.G,
					atposPref, VnlRadial[l][p], V.dataPref()+offs, false, vector3<complex*>());
			}
			// Augment overlap:
			OCq += V * (Q * (V^Cq));
		}
	}
	watch.stop();
}

#define augmentDensity_COMMON_INIT \
	if(!atpos.size()) return; /*unused species*/ \
	if(!Qint.size()) return; /*no overlap augmentation*/ \
	/*Determine dimensions:*/ \
	int nProj = 0, lMax = 0; \
	for(unsigned l=0; l<VnlRadial.size(); l++) \
	{	nProj += (2*l+1)*VnlRadial[l].size(); \
		if(VnlRadial[l].size()) lMax=l; \
	} \
	int Nlm = (2*lMax+1)*(2*lMax+1); \
	int nCoeff = Qradial.cbegin()->second.nCoeff; \
	int nCoeffAtom = nCoeff * Nlm; \
	int nCoeffSpin = nCoeffAtom * atpos.size();

void SpeciesInfo::augmentDensityInit()
{	augmentDensity_COMMON_INIT
	size_t nCoeffTot = nCoeffSpin * (e->eInfo.spinType==SpinNone ? 1 : 2);
	if(!nAug)
	{
		#ifdef GPU_ENABLED
		cudaMalloc(&nAug, nCoeffTot*sizeof(double));
		cudaMalloc(&E_nAug, nCoeffTot*sizeof(double));
		#else
		nAug = new double[nCoeffTot];
		E_nAug = new double[nCoeffTot];
		#endif
	}
	callPref(eblas_zero)(nCoeffTot, nAug);
}

void SpeciesInfo::augmentDensityCleanup()
{	if(nAug)
	{
		#ifdef GPU_ENABLED
		cudaFree(nAug); nAug=0;
		cudaFree(E_nAug); E_nAug=0;
		if(nagIndex) cudaFree(nagIndex); nagIndex=0;
		if(nagIndexPtr) cudaFree(nagIndexPtr); nagIndexPtr=0;
		#else
		delete[] nAug; nAug=0;
		delete[] E_nAug; E_nAug=0;
		if(nagIndex) delete[] nagIndex; nagIndex=0;
		if(nagIndexPtr) delete[] nagIndexPtr; nagIndexPtr=0;
		#endif
	}
}

void SpeciesInfo::augmentDensitySpherical(const diagMatrix& Fq, const ColumnBundle& Cq)
{	static StopWatch watch("augmentDensitySpherical"), watchMakeProj("adsMakeProj"), watchMatMul("adsMatMul"); watch.start(); 
	augmentDensity_COMMON_INIT
	const GridInfo &gInfo = e->gInfo;
	const Basis& basis = *Cq.basis;
	//Allocate temporaries:
	ColumnBundle V = Cq.similar(nProj);
	//Loop over atoms:
	for(unsigned atom=0; atom<atpos.size(); atom++)
	{	//Initialize all projectors at this atom:
		watchMakeProj.start();
		int iProj = 0;
		for(int l=0; l<int(VnlRadial.size()); l++)
			for(unsigned p=0; p<VnlRadial[l].size(); p++)
				for(int m=-l; m<=l; m++)
				{	size_t offs = iProj * basis.nbasis;
					size_t atomStride = nProj * basis.nbasis;
					callPref(Vnl)(basis.nbasis, atomStride, 1, l, m, V.qnum->k, basis.iGarrPref, gInfo.G,
						atposPref+atom, VnlRadial[l][p], V.dataPref()+offs, false, vector3<complex*>());
					iProj++;
				}
		watchMakeProj.stop();
		watchMatMul.start();
		matrix VdagC = V ^ Cq;
		matrix Rho = VdagC * Fq * dagger(VdagC); //density matrix in projector basis
		watchMatMul.stop();
		
		//Collect contributions as spherical functions:
		double* nAugCur = nAug + Cq.qnum->index()*nCoeffSpin + atom*nCoeffAtom;
		//Triple loop over first projector:
		int i1 = 0;
		for(int l1=0; l1<int(VnlRadial.size()); l1++)
		for(int p1=0; p1<int(VnlRadial[l1].size()); p1++)
		for(int m1=-l1; m1<=l1; m1++)
		{	//Triple loop over second projector:
			int i2 = 0;
			for(int l2=0; l2<int(VnlRadial.size()); l2++)
			for(int p2=0; p2<int(VnlRadial[l2].size()); p2++)
			for(int m2=-l2; m2<=l2; m2++)
			{	if(i2<=i1) //rest handled by i1<->i2 symmetry
				{	std::vector<YlmProdTerm> terms = expandYlmProd(l1,m1, l2,m2);
					double prefac = Cq.qnum->weight * ((i1==i2 ? 1 : 2)/gInfo.detR)
								* (Rho.data()[Rho.index(i2,i1)] * cis(0.5*M_PI*(l2-l1))).real();
					for(const YlmProdTerm& term: terms)
					{	QijIndex qIndex = { l1, p1, l2, p2, term.l };
						auto Qijl = Qradial.find(qIndex);
						if(Qijl==Qradial.end()) continue; //no entry at this l
						callPref(eblas_daxpy)(nCoeff, term.coeff * prefac, Qijl->second.coeffPref(),1, nAugCur+nCoeff*(term.l*(term.l+1) + term.m),1);
					}
				}
				i2++;
			}
			i1++;
		}
	}
	watch.stop();
}

void SpeciesInfo::augmentDensityGrid(DataRptrCollection& n) const
{	static StopWatch watch("augmentDensityGrid"); watch.start(); 
	augmentDensity_COMMON_INIT
	const GridInfo &gInfo = e->gInfo;
	double dGinv = 1./dGloc;
	for(unsigned s=0; s<n.size(); s++)
	{	DataGptr nAugTilde; nullToZero(nAugTilde, gInfo);
		for(unsigned atom=0; atom<atpos.size(); atom++)
		{	const double* nAugCur = nAug + s*nCoeffSpin + atom*nCoeffAtom;
			callPref(nAugment)(Nlm, gInfo.S, gInfo.G, nCoeff, dGinv, nAugCur, atpos[atom], nAugTilde->dataPref());
		}
		n[s] += I(nAugTilde,true);
	}
	watch.stop();
}

void SpeciesInfo::augmentDensityGridGrad(const DataRptrCollection& E_n, std::vector<vector3<> >* forces)
{	static StopWatch watch("augmentDensityGridGrad"); watch.start();
	augmentDensity_COMMON_INIT
	if(!nAug) augmentDensityInit();
	const GridInfo &gInfo = e->gInfo;
	double dGinv = 1./dGloc;
	callPref(eblas_zero)(nCoeffSpin * (e->eInfo.spinType==SpinNone ? 1 : 2), E_nAug);
	DataGptrVec E_atpos; if(forces) nullToZero(E_atpos, gInfo);
	for(unsigned s=0; s<E_n.size(); s++)
	{	DataGptr ccE_n = Idag(E_n[s]);
		for(unsigned atom=0; atom<atpos.size(); atom++)
		{	const double* nAugCur = forces ? nAug + s*nCoeffSpin + atom*nCoeffAtom : 0;
			double* E_nAugCur = E_nAug + s*nCoeffSpin + atom*nCoeffAtom;
			if(forces) initZero(E_atpos);
			callPref(nAugmentGrad)(Nlm, gInfo.S, gInfo.G, nCoeff, dGinv, nAugCur, atpos[atom],
				ccE_n->dataPref(), E_nAugCur, forces ? E_atpos.dataPref() : vector3<complex*>(), nagIndex, nagIndexPtr);
			if(forces) for(int k=0; k<3; k++) (*forces)[atom][k] -= sum(E_atpos[k]);
		}
	}
	watch.stop();
}

void SpeciesInfo::augmentDensitySphericalGrad(const diagMatrix& Fq, const ColumnBundle& Cq, ColumnBundle& HCq, std::vector<vector3<> >* forces, const matrix& gradCdagOCq) const
{	static StopWatch watch("augmentDensitySphericalGrad"); watch.start();
	augmentDensity_COMMON_INIT
	const GridInfo &gInfo = e->gInfo;
	const Basis& basis = *Cq.basis;
	ColumnBundle V = Cq.similar(nProj);
	ColumnBundle dV[3]; if(forces) for(int k=0; k<3; k++) dV[k] = V.similar();
	//Loop over atoms:
	for(unsigned atom=0; atom<atpos.size(); atom++)
	{	vector3<complex*> dVdata;
		if(forces) for(int k=0; k<3; k++) dVdata[k] = dV[k].dataPref();
		//Initialize all projectors at this atom:
		int iProj = 0;
		for(int l=0; l<int(VnlRadial.size()); l++)
			for(unsigned p=0; p<VnlRadial[l].size(); p++)
				for(int m=-l; m<=l; m++)
				{	size_t offs = iProj * basis.nbasis;
					size_t atomStride = nProj * basis.nbasis;
					callPref(Vnl)(basis.nbasis, atomStride, 1, l, m, V.qnum->k, basis.iGarrPref, gInfo.G,
						atposPref+atom, VnlRadial[l][p], V.dataPref()+offs, forces, dVdata+offs);
					iProj++;
				}
		matrix VdagC = V ^ Cq;
		matrix Rho = VdagC * Fq * dagger(VdagC); //density matrix in projector basis
		matrix Qint(nProj, nProj); Qint.zero(); //Full nProj x nProj version of this->Qint[l]
		matrix E_Rho(nProj, nProj); E_Rho.zero(); //gradient w.r.t density matrix in projector basis
		
		const double* E_nAugCur = E_nAug + Cq.qnum->index()*nCoeffSpin + atom*nCoeffAtom;
		
		int i1 = 0;
		//Triple loop over first projector:
		for(int l1=0; l1<int(VnlRadial.size()); l1++)
		for(int p1=0; p1<int(VnlRadial[l1].size()); p1++)
		for(int m1=-l1; m1<=l1; m1++)
		{	//Triple loop over second projector:
			int i2 = 0;
			for(int l2=0; l2<int(VnlRadial.size()); l2++)
			for(int p2=0; p2<int(VnlRadial[l2].size()); p2++)
			for(int m2=-l2; m2<=l2; m2++)
			{	if(i2<=i1) //rest handled by i1<->i2 symmetry 
				{	std::vector<YlmProdTerm> terms = expandYlmProd(l1,m1, l2,m2);
					double E_Rho_i1i2sum = 0.;
					for(const YlmProdTerm& term: terms)
					{	QijIndex qIndex = { l1, p1, l2, p2, term.l };
						auto Qijl = Qradial.find(qIndex);
						if(Qijl==Qradial.end()) continue; //no entry at this l
						E_Rho_i1i2sum += term.coeff * callPref(eblas_ddot)(nCoeff, Qijl->second.coeffPref(),1, E_nAugCur+nCoeff*(term.l*(term.l+1) + term.m),1);
					}
					complex E_Rho_i1i2 = E_Rho_i1i2sum * (1./gInfo.detR) * cis(0.5*M_PI*(l2-l1));
					E_Rho.data()[E_Rho.index(i2,i1)] += E_Rho_i1i2.conj();
					if(i1!=i2) E_Rho.data()[E_Rho.index(i1,i2)] += E_Rho_i1i2;
				}
				if(l1==l2 && m1==m2)
				{	Qint.data()[Qint.index(i1,i2)] = this->Qint[l1].data()[this->Qint[l1].index(p1,p2)];
					Qint.data()[Qint.index(i2,i1)] = this->Qint[l1].data()[this->Qint[l1].index(p2,p1)];
				}
				i2++;
			}
			i1++;
		}
		matrix E_RhoVdagC = E_Rho * VdagC;
		if(HCq) HCq += V * E_RhoVdagC;
		if(forces)
		{	for(int k=0; k<3; k++)
			{	matrix dVdagC = dV[k]^Cq;
				(*forces)[atom][k] -= 2.*Cq.qnum->weight *
						( trace(E_RhoVdagC * Fq * dagger(dVdagC)).real() //Contribution via dV
						+ trace(Qint * VdagC * gradCdagOCq * dagger(dVdagC)).real() ); //Contribution via overlap
			}
		}
	}
	watch.stop();
}


bool SpeciesInfo::QijIndex::operator<(const SpeciesInfo::QijIndex& other) const
{	//Bring both indices to the upper triangular part:
	QijIndex q1 = *this; q1.sortIndices();
	QijIndex q2 = other; q2.sortIndices();
	//Compare:
	if(q1.l1 < q2.l1) return true;
	if(q1.l1 == q2.l1)
	{	if(q1.p1 < q2.p1) return true;
		if(q1.p1 == q2.p1)
		{	if(q1.l2 < q2.l2) return true;
			if(q1.l2 == q2.l2)
			{	if(q1.p2 < q2.p2) return true;
				if(q1.p2 == q2.p2)
					return q1.l < q2.l;
			}
		}
	}
	return false;
}
void SpeciesInfo::QijIndex::sortIndices()
{	if(l1>l2)
	{	std::swap(l1,l2);
		std::swap(p1,p2);
	}
	else if(l1==l2 && p1>p2)
	{	std::swap(p1,p2);
	}
}
