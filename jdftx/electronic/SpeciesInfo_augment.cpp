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


void SpeciesInfo::augmentOverlap(const ColumnBundle& Cq, ColumnBundle& OCq, matrix* VdagCqPtr) const
{	static StopWatch watch("augmentOverlap"); watch.start();
	if(!atpos.size()) return; //unused species
	if(!Qint.size()) return; //no overlap augmentation
	std::shared_ptr<ColumnBundle> V = getV(Cq);
	matrix VdagCq = (*V) ^ Cq;
	if(VdagCqPtr) *VdagCqPtr = VdagCq; //cache for later usage
	OCq += (*V) * (tiledBlockMatrix(QintAll,atpos.size()) * VdagCq);
	watch.stop();
}

#define augmentDensity_COMMON_INIT \
	if(!atpos.size()) return; /*unused species*/ \
	if(!Qint.size()) return; /*no overlap augmentation*/ \
	/*Determine dimensions:*/ \
	int lMax = 0; \
	for(unsigned l=0; l<VnlRadial.size(); l++) \
		if(VnlRadial[l].size()) lMax=l; \
	int Nlm = (2*lMax+1)*(2*lMax+1);

#define augmentDensityGrid_COMMON_INIT \
	augmentDensity_COMMON_INIT \
	int nCoeffHlf = (Qradial.cbegin()->second.nCoeff+1)/2; /*pack real radial functions into complex numbers*/ \
	int nCoeff = 2*nCoeffHlf;


void SpeciesInfo::augmentDensityInit()
{	augmentDensity_COMMON_INIT
	size_t nSpinAtomLM = (e->eInfo.spinType==SpinNone ? 1 : 2) * atpos.size() * Nlm;
	if(!nAug)
	{	nAug.init(Qradial.size(), nSpinAtomLM, isGpuEnabled());
		E_nAug.init(Qradial.size(), nSpinAtomLM, isGpuEnabled());
	}
	nAug.zero();
	E_nAug.zero();
}

void SpeciesInfo::augmentDensitySpherical(const QuantumNumber& qnum, const diagMatrix& Fq, const matrix& VdagCq)
{	static StopWatch watch("augmentDensitySpherical"); watch.start(); 
	augmentDensity_COMMON_INIT
	int nProj = MnlAll.nRows();
	const GridInfo &gInfo = e->gInfo;
	complex* nAugData = nAug.data();
	
	//Loop over atoms:
	for(unsigned atom=0; atom<atpos.size(); atom++)
	{	//Get projections at this atom:
		matrix atomVdagC = VdagCq(atom*nProj,(atom+1)*nProj, 0,VdagCq.nCols());
		matrix Rho = atomVdagC * Fq * dagger(atomVdagC); //density matrix in projector basis on this atom
		
		int atomOffs = Nlm*(atom + qnum.index()*atpos.size());
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
					double prefac = qnum.weight * ((i1==i2 ? 1 : 2)/gInfo.detR)
								* (Rho.data()[Rho.index(i2,i1)] * cis(0.5*M_PI*(l2-l1))).real();
					for(const YlmProdTerm& term: terms)
					{	QijIndex qIndex = { l1, p1, l2, p2, term.l };
						auto Qijl = Qradial.find(qIndex);
						if(Qijl==Qradial.end()) continue; //no entry at this l
						nAugData[nAug.index(Qijl->first.index, atomOffs + term.l*(term.l+1) + term.m)] += term.coeff * prefac;
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
	augmentDensityGrid_COMMON_INIT
	const GridInfo &gInfo = e->gInfo;
	double dGinv = 1./dGloc;
	matrix nAugTot = nAug; nAugTot.allReduce(MPIUtil::ReduceSum); //collect radial functions from all processes, and split by G-vectors below
	matrix nAugRadial = QradialMat * nAugTot; //transform from radial functions to spline coeffs
	double* nAugRadialData = (double*)nAugRadial.dataPref();
	for(unsigned s=0; s<n.size(); s++)
	{	DataGptr nAugTilde; nullToZero(nAugTilde, gInfo);
		for(unsigned atom=0; atom<atpos.size(); atom++)
		{	int atomOffs = nCoeff * Nlm * (atom + atpos.size()*s);
			callPref(nAugment)(Nlm, gInfo.S, gInfo.G, gInfo.iGstart, gInfo.iGstop, nCoeff, dGinv, nAugRadialData+atomOffs, atpos[atom], nAugTilde->dataPref());
		}
		n[s] += I(nAugTilde,true);
	}
	watch.stop();
}

void SpeciesInfo::augmentDensityGridGrad(const DataRptrCollection& E_n, std::vector<vector3<> >* forces)
{	static StopWatch watch("augmentDensityGridGrad"); watch.start();
	augmentDensityGrid_COMMON_INIT
	if(!nAug) augmentDensityInit();
	const GridInfo &gInfo = e->gInfo;
	double dGinv = 1./dGloc;
	matrix E_nAugRadial = zeroes(nCoeffHlf, (e->eInfo.spinType==SpinNone ? 1 : 2) * atpos.size() * Nlm);
	double* E_nAugRadialData = (double*)E_nAugRadial.dataPref();
	matrix nAugRadial; const double* nAugRadialData=0;
	if(forces)
	{	matrix nAugTot = nAug; nAugTot.allReduce(MPIUtil::ReduceSum);
		nAugRadial = QradialMat * nAugTot;
		nAugRadialData = (const double*)nAugRadial.dataPref();
	}
	DataGptrVec E_atpos; if(forces) nullToZero(E_atpos, gInfo);
	for(unsigned s=0; s<E_n.size(); s++)
	{	DataGptr ccE_n = Idag(E_n[s]);
		for(unsigned atom=0; atom<atpos.size(); atom++)
		{	int atomOffs = nCoeff * Nlm * (atom + atpos.size()*s);
			if(forces) initZero(E_atpos);
			callPref(nAugmentGrad)(Nlm, gInfo.S, gInfo.G, nCoeff, dGinv, forces? (nAugRadialData+atomOffs) :0, atpos[atom],
				ccE_n->dataPref(), E_nAugRadialData+atomOffs, forces ? E_atpos.dataPref() : vector3<complex*>(), nagIndex, nagIndexPtr);
			if(forces) for(int k=0; k<3; k++) (*forces)[atom][k] -= sum(E_atpos[k]);
		}
	}
	E_nAug = dagger(QradialMat) * E_nAugRadial;  //propagate from spline coeffs to radial functions
	E_nAug.allReduce(MPIUtil::ReduceSum);
	watch.stop();
}

void SpeciesInfo::augmentDensitySphericalGrad(const QuantumNumber& qnum, const diagMatrix& Fq, const matrix& VdagCq, matrix& HVdagCq) const
{	static StopWatch watch("augmentDensitySphericalGrad"); watch.start();
	augmentDensity_COMMON_INIT
	int nProj = MnlAll.nRows();
	const GridInfo &gInfo = e->gInfo;
	const complex* E_nAugData = E_nAug.data();

	matrix E_RhoVdagC(VdagCq.nRows(),VdagCq.nCols(),isGpuEnabled());
	
	//Loop over atoms:
	for(unsigned atom=0; atom<atpos.size(); atom++)
	{	matrix E_Rho = zeroes(nProj, nProj); //gradient w.r.t density matrix in projector basis
		
		int atomOffs = Nlm*(atom + qnum.index()*atpos.size());
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
					double E_Rho_i1i2sum = 0.;
					for(const YlmProdTerm& term: terms)
					{	QijIndex qIndex = { l1, p1, l2, p2, term.l };
						auto Qijl = Qradial.find(qIndex);
						if(Qijl==Qradial.end()) continue; //no entry at this l
						E_Rho_i1i2sum += term.coeff * E_nAugData[E_nAug.index(Qijl->first.index, atomOffs + term.l*(term.l+1) + term.m)].real();
					}
					complex E_Rho_i1i2 = E_Rho_i1i2sum * (1./gInfo.detR) * cis(0.5*M_PI*(l2-l1));
					E_Rho.data()[E_Rho.index(i2,i1)] += E_Rho_i1i2.conj();
					if(i1!=i2) E_Rho.data()[E_Rho.index(i1,i2)] += E_Rho_i1i2;
				}
				i2++;
			}
			i1++;
		}
		
		matrix atomVdagC = VdagCq(atom*nProj,(atom+1)*nProj, 0,VdagCq.nCols());
		matrix E_atomRhoVdagC = E_Rho * atomVdagC;
		E_RhoVdagC.set(atom*nProj,(atom+1)*nProj, 0,VdagCq.nCols(), E_atomRhoVdagC);
	}
	HVdagCq += E_RhoVdagC;
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
