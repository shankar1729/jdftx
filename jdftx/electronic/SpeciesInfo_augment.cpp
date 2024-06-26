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
#include <electronic/ColumnBundle.h>
#include <core/matrix.h>

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

void SpeciesInfo::augmentOverlapDeriv(const ColumnBundle& Cq, ColumnBundle& OCq, ColumnBundle& V, ColumnBundle& dV) const
{	static StopWatch watch("augmentOverlapDeriv"); watch.start();
	if(!atpos.size()) return; //unused species
	if(!Qint.size()) return; //no overlap augmentation
	OCq += V * (QintAll * (dV^Cq));
	OCq += dV * (QintAll * (V^Cq));
	watch.stop();
}

void SpeciesInfo::augmentOverlapSecondDeriv(const ColumnBundle& Cq, ColumnBundle& OCq, ColumnBundle& V, ColumnBundle& dV_A, ColumnBundle& dV_B, ColumnBundle& dsqV) const
{	static StopWatch watch("augmentOverlapDeriv"); watch.start();
	if(!atpos.size()) return; //unused species
	if(!Qint.size()) return; //no overlap augmentation
	OCq += V * (QintAll * (dsqV^Cq));
	OCq += dV_A * (QintAll * (dV_B^Cq));
	OCq += dV_B * (QintAll * (dV_A^Cq));
	OCq += dsqV * (QintAll * (V^Cq));
	watch.stop();
}

void SpeciesInfo::augmentSpinOverlap(const ColumnBundle& Cq, vector3<matrix>& Sq) const
{	if(!atpos.size()) return; //unused species
	if(!Qint.size()) return; //no overlap augmentation
	assert(e->eInfo.isNoncollinear());
	matrix VdagCq = (*getV(Cq)) ^ Cq; //pseudopotential projection
	//Pauli sigma matrices:
	std::vector<matrix> pauli(3, zeroes(2,2));
	pauli[0].set(0,1, +1); pauli[0].set(1,0, +1);
	pauli[1].set(0,1, complex(0,-1.)); pauli[1].set(1,0, complex(0,1.));
	pauli[2].set(0,0, +1); pauli[2].set(1,1, -1);
	//Loop over spin components:
	for(int iDir=0; iDir<3; iDir++)
	{	//Version of QintAll based on Pauli matrix insted of eye(2):
		int nProj = QintAll.nRows();
		matrix Qspin = zeroes(nProj,nProj);
		int lOffset = 0;
		for(unsigned l=0; l<VnlRadial.size(); l++)
		{	unsigned nMS = 2*(2*l+1); //number of m and spins at each l
			int iProj = lOffset;
			for(unsigned ni=0; ni<VnlRadial[l].size(); ni++)
			{	int jProj = lOffset;
				for(unsigned nj=0; nj<VnlRadial[l].size(); nj++)
				{	if(Qint[l])
					{	for(unsigned iMS=0; iMS<nMS; iMS+=2)
							Qspin.set(iProj+iMS,iProj+iMS+2, jProj+iMS,jProj+iMS+2,
								Qint[l].data()[Qint[l].index(ni,nj)] * pauli[iDir]);
					}
					jProj += nMS;
				}
				iProj += nMS;
			}
			lOffset = iProj;
		}
		if(isRelativistic())
			Qspin = fljAll * Qspin * fljAll;
		//Augment:
		Sq[iDir] += dagger(VdagCq) * (tiledBlockMatrix(Qspin,atpos.size()) * VdagCq);
	}
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


void SpeciesInfo::augmentDensityInit(int atom)
{	augmentDensity_COMMON_INIT
	size_t nSpinAtomLM = e->eInfo.nDensities * (atom >= 0? 1 : atpos.size()) * Nlm;
	if(!nAug)
	{	nAug.init(Qradial.size(), nSpinAtomLM, isGpuEnabled());
		E_nAug.init(Qradial.size(), nSpinAtomLM, isGpuEnabled());
	}
	nAug.zero();
	E_nAug.zero();
}

void SpeciesInfo::augmentDensitySpherical(const QuantumNumber& qnum, const diagMatrix& Fq, const matrix& VdagCq, const matrix* VdagdCqL, const matrix* VdagdCqR, int atom)
{	static StopWatch watch("augmentDensitySpherical"); watch.start(); 
	augmentDensity_COMMON_INIT
	int nProj = MnlAll.nRows();
	const GridInfo &gInfo = e->gInfo;
	complex* nAugData = nAug.data();
	unsigned nAtoms = (atom >= 0) ? 1 : atpos.size();
	bool derivativeMode = (VdagdCqL or VdagdCqR);
	if(derivativeMode) assert(VdagdCqL);
	
	//Loop over atoms:
	for(unsigned at=0; at<nAtoms; at++)
	{	//Get projections and calculate density matrix at this atom:
		matrix RhoAll;
		matrix atomVdagC = VdagCq(at*nProj,(at+1)*nProj, 0,VdagCq.nCols());
		if(not derivativeMode)
			RhoAll = atomVdagC * Fq * dagger(atomVdagC); //density matrix in projector basis on this atom
		else
		{	//Compute and store first derivative of RhoAll instead:
			matrix atomVdagdCL = (*VdagdCqL)(at*nProj,(at+1)*nProj, 0,VdagdCqL->nCols());
			const matrix& atomVdagdCR = VdagdCqR
				? matrix((*VdagdCqR)(at*nProj,(at+1)*nProj, 0,VdagdCqR->nCols()))
				: atomVdagdCL;
			RhoAll = atomVdagC * Fq * dagger(atomVdagdCR) + atomVdagdCL * Fq * dagger(atomVdagC); //Compute and store first 		
		}
		if(isRelativistic()) RhoAll = fljAll * RhoAll * fljAll; //transformation for relativistic pseudopotential
		std::vector<matrix> Rho(e->eInfo.nDensities); //RhoAll split by spin(-density-matrix) components
		if(e->eInfo.isNoncollinear())
		{	matrix RhoUp = RhoAll(0,2,nProj, 0,2,nProj);
			matrix RhoDn = RhoAll(1,2,nProj, 1,2,nProj);
			if(Rho.size()==1)
				Rho[0] = RhoUp + RhoDn; //unpolarized noncollinear mode
			else
			{	matrix RhoUpDn = RhoAll(0,2,nProj, 1,2,nProj);
				matrix RhoDnUp = RhoAll(1,2,nProj, 0,2,nProj);
				Rho[0] = RhoUp;
				Rho[1] = RhoDn;
				Rho[2] = (RhoUpDn + RhoDnUp) * 0.5; //'real part' of UpDn
				Rho[3] = (RhoUpDn - RhoDnUp) * complex(0,-0.5); //'imaginary part' of UpDn
			}
		}
		else std::swap(Rho[qnum.index()], RhoAll); //in this case each qnum contributes to a specific spin component
		
		//Calculate spherical function contributions from density matrix:
		for(size_t s=0; s<Rho.size(); s++) if(Rho[s])
		{	int atomOffs = Nlm * (at + s*nAtoms);
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
									* (Rho[s].data()[Rho[s].index(i2,i1)] * cis(0.5*M_PI*(l2-l1))).real();
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
	}
	watch.stop();
}

void SpeciesInfo::augmentDensityGrid(ScalarFieldArray& n, int atom, const vector3<>* atposDeriv) const
{	static StopWatch watch("augmentDensityGrid"); watch.start(); 
	augmentDensityGrid_COMMON_INIT
	const GridInfo &gInfo = e->gInfo;
	double dGinv = 1./gInfo.dGradial;
	matrix nAugTot = nAug; mpiWorld->allReduceData(nAugTot, MPIUtil::ReduceSum); //collect radial functions from all processes, and split by G-vectors below
	matrix nAugRadial = QradialMat * nAugTot; //transform from radial functions to spline coeffs
	double* nAugRadialData = (double*)nAugRadial.dataPref();
	for(unsigned s=0; s<n.size(); s++)
	{	ScalarFieldTilde nAugTilde; nullToZero(nAugTilde, gInfo);
		unsigned atoms = (atom >= 0)? 1 : atpos.size();
		for(unsigned atomindex=0; atomindex<atoms; atomindex++)
		{	int atomOffs = (atom >= 0)? nCoeff * Nlm * s : nCoeff * Nlm * (atomindex + atpos.size()*s);
			callPref(nAugment)(Nlm, gInfo.S, gInfo.G, gInfo.iGstart, gInfo.iGstop, nCoeff, dGinv, nAugRadialData+atomOffs, atpos[(atom >= 0)? atom : atomindex], nAugTilde->dataPref(), atposDeriv);
		}
		n[s] += I(nAugTilde);
	}
	watch.stop();
}

void SpeciesInfo::augmentDensityGridGrad(const ScalarFieldArray& E_n, std::vector<vector3<>>* forces, matrix3<>* Eaug_RRT)
{	static StopWatch watch("augmentDensityGridGrad"); watch.start();
	augmentDensityGrid_COMMON_INIT
	if(!nAug) augmentDensityInit();
	const GridInfo &gInfo = e->gInfo;
	double dGinv = 1./gInfo.dGradial;
	matrix E_nAugRadial = zeroes(nCoeffHlf, e->eInfo.nDensities * atpos.size() * Nlm);
	double* E_nAugRadialData = (double*)E_nAugRadial.dataPref();
	matrix nAugRadial; const double* nAugRadialData=0;
	if(forces or Eaug_RRT)
	{	matrix nAugTot = nAug; mpiWorld->allReduceData(nAugTot, MPIUtil::ReduceSum);
		nAugRadial = QradialMat * nAugTot;
		nAugRadialData = (const double*)nAugRadial.dataPref();
	}
	VectorFieldTilde E_atpos; if(forces) nullToZero(E_atpos, gInfo);
	ScalarFieldTildeArray E_RRT(6); if(Eaug_RRT) nullToZero(E_RRT, gInfo);
	for(unsigned s=0; s<E_n.size(); s++)
	{	ScalarFieldTilde ccE_n = Idag(E_n[s]);
		for(unsigned atom=0; atom<atpos.size(); atom++)
		{	int atomOffs = nCoeff * Nlm * (atom + atpos.size()*s);
			if(forces) initZero(E_atpos);
			callPref(nAugmentGrad)(Nlm, gInfo.S, gInfo.G, nCoeff, dGinv,
				nAugRadialData ? (nAugRadialData+atomOffs) : 0,
				atpos[atom], ccE_n->dataPref(), E_nAugRadialData+atomOffs,
				forces ? E_atpos.dataPref() : vector3<complex*>(),
				Eaug_RRT ? array<complex*,6>(dataPref(E_RRT)) : array<complex*,6>(),
				0, nagIndex.dataPref(), nagIndexPtr.dataPref());
			if(forces) for(int k=0; k<3; k++) (*forces)[atom][k] -= sum(E_atpos[k]);
		}
	}
	if(Eaug_RRT)
	{	symmetricMatrix3<> E_RRTsum;
		double* E_RRTsumData = (double*)&E_RRTsum;
		for(int ij=0; ij<6; ij++)
			E_RRTsumData[ij] = sum(E_RRT[ij]);
		*Eaug_RRT += matrix3<>(E_RRTsum);
	}
	E_nAug = dagger(QradialMat) * E_nAugRadial;  //propagate from spline coeffs to radial functions
	mpiWorld->allReduceData(E_nAug, MPIUtil::ReduceSum);
	watch.stop();
}


void SpeciesInfo::augmentDensityGridGradDeriv(const ScalarFieldArray& E_n, int atom, const vector3<>* atposDeriv)
{	static StopWatch watch("augmentDensityGridGradDeriv"); watch.start();
	augmentDensityGrid_COMMON_INIT
	if(!nAug) augmentDensityInit();
	const GridInfo &gInfo = e->gInfo;
	double dGinv = 1./gInfo.dGradial;
	matrix E_nAugRadial = zeroes(nCoeffHlf, e->eInfo.nDensities * atpos.size() * Nlm);
	double* E_nAugRadialData = (double*)E_nAugRadial.dataPref();
	for(unsigned s=0; s<E_n.size(); s++)
	{	ScalarFieldTilde ccE_n = Idag(E_n[s]);
		int atomOffs = nCoeff * Nlm * (atom + atpos.size()*s);
		callPref(nAugmentGrad)(Nlm, gInfo.S, gInfo.G, nCoeff, dGinv, 0,
			atpos[atom], ccE_n->dataPref(), E_nAugRadialData+atomOffs,
			vector3<complex*>(),
			array<complex*,6>(),
			atposDeriv, nagIndex.dataPref(), nagIndexPtr.dataPref());
	}
	E_nAug = dagger(QradialMat) * E_nAugRadial;  //propagate from spline coeffs to radial functions
	mpiWorld->allReduceData(E_nAug, MPIUtil::ReduceSum);
	watch.stop();
}

void SpeciesInfo::augmentDensitySphericalGrad(const QuantumNumber& qnum, const matrix& VdagCq, matrix& HVdagCq, int atom) const
{	static StopWatch watch("augmentDensitySphericalGrad"); watch.start();
	augmentDensity_COMMON_INIT
	int nProj = MnlAll.nRows();
	const GridInfo &gInfo = e->gInfo;
	const complex* E_nAugData = E_nAug.data();

	matrix E_RhoVdagC(VdagCq.nRows(),VdagCq.nCols(),isGpuEnabled());
	
	//Loop over atoms:
	unsigned atoms = atpos.size();
	if (atom >= 0) atoms = 1;
	
	for(unsigned Vatomindex=0; Vatomindex<atoms; Vatomindex++)
	{
		int atomindex = Vatomindex;
		if (atom >= 0)
			atomindex = atom;
		
		//Prepare gradient w.r.t density matrix in basis of current atom's projectors (split by spinor components, if any)
		std::vector<matrix> E_Rho(e->eInfo.nDensities);
		if(e->eInfo.isNoncollinear()) E_Rho.assign(E_Rho.size(), zeroes(nProj/2, nProj/2));
		else E_Rho[qnum.index()] = zeroes(nProj, nProj);
		
		//Propagate gradients from spherical functions to density matrix:
		for(size_t s=0; s<E_Rho.size(); s++) if(E_Rho[s])
		{	int atomOffs = Nlm*(atomindex + s*atpos.size());
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
						E_Rho[s].data()[E_Rho[s].index(i2,i1)] += E_Rho_i1i2.conj();
						if(i1!=i2) E_Rho[s].data()[E_Rho[s].index(i1,i2)] += E_Rho_i1i2;
					}
					i2++;
				}
				i1++;
			}
		}
		
		//Collate density matrix from spinor components (if necessary)
		matrix E_RhoAll;
		if(e->eInfo.isNoncollinear())
		{	E_RhoAll = zeroes(nProj, nProj);
			E_RhoAll.set(0,2,nProj, 0,2,nProj, E_Rho[0]);
			E_RhoAll.set(1,2,nProj, 1,2,nProj, E_Rho[E_Rho.size()>1 ? 1 : 0]);
			if(E_Rho.size()>1) //full noncollinear mode (with magnetization)
			{	E_RhoAll.set(0,2,nProj, 1,2,nProj, 0.5*E_Rho[2] + complex(0,0.5)*E_Rho[3]);
				E_RhoAll.set(1,2,nProj, 0,2,nProj, 0.5*E_Rho[2] - complex(0,0.5)*E_Rho[3]);
			}
		}
		else std::swap(E_RhoAll, E_Rho[qnum.index()]);
		
		//Propagate gradients from densiy matrix to projections:
		if(isRelativistic()) E_RhoAll = fljAll * E_RhoAll * fljAll; //transformation for relativistic pseudopotential
		matrix atomVdagC = VdagCq(Vatomindex*nProj,(Vatomindex+1)*nProj, 0,VdagCq.nCols());
		matrix E_atomRhoVdagC = E_RhoAll * atomVdagC;
		E_RhoVdagC.set(Vatomindex*nProj,(Vatomindex+1)*nProj, 0,VdagCq.nCols(), E_atomRhoVdagC);
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
