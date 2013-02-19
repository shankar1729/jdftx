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

#include <electronic/Everything.h>
#include <electronic/SpeciesInfo.h>
#include <electronic/Symmetries.h>
#include <electronic/ElecInfo.h>
#include <electronic/IonInfo.h>
#include <electronic/Basis.h>
#include <electronic/operators.h>
#include <core/LatticeUtils.h>
#include <core/GridInfo.h>
#include <core/Thread.h>

//Symmetry detection code based on that of Nikolaj Moll, April 1999 (symm.c in Sohrab's version)

Symmetries::Symmetries()
{	nSymmIndex = 0;
	shouldPrintMatrices = false;
}

Symmetries::~Symmetries()
{	if(nSymmIndex)
	{	
		#ifdef GPU_ENABLED
		cudaFree(symmIndex);
		#else
		delete[] symmIndex;
		#endif
	}
}

void Symmetries::setup(const Everything& everything)
{	e = &everything;
	logPrintf("\n---------- Setting up symmetries ----------\n");

	//Calculate and check symmetries if needed:
	switch(mode)
	{	case SymmetriesAutomatic: //Automatic symmetries
			calcSymmetries();
			break;
		case SymmetriesManual: //Manually specified matrices
			if(!sym.size())
				die("\nManual symmetries specified without specifying any symmetry matrices.\n");
			sortSymmetries(); //make sure first symmetry is identity
			checkSymmetries(); //make sure atoms respect the specified symmetries
			break;
		default: //No symmetry (only matrix is identity)
			sym.assign(1, matrix3<int>(1,1,1)); 
	}
	
	initAtomMaps(); // Map atoms to symmetry related ones
}

void Symmetries::setupMesh()
{	checkFFTbox(); //Check that the FFT box is commensurate with the symmetries and initialize mesh matrices
	checkKmesh(); //Check symmetries of k-point mesh (and warn if lower than basis symmetries)
	initSymmIndex(); //Initialize the equivalence classes for scalar field symmetrization (using mesh matrices)
}


std::list<QuantumNumber> Symmetries::reduceKmesh(const std::vector<QuantumNumber>& qnums) const
{	//Convert to a list for efficient removal:
	std::list<QuantumNumber> qlist(qnums.begin(), qnums.end());
	typedef std::list<QuantumNumber>::iterator Iter;
	bool usedInversion = false; //whether inversion was used in reducing k-point mesh
	std::vector<int> invertList;
	invertList.push_back(+1);
	if(sym.size()>1 && e->eInfo.spinType==SpinNone)
		invertList.push_back(-1);
	for(int invert: invertList) //First try without inversion, and then if required, try again with inversion
	{	//For each state, gobble up the weights of all subsequent equivalent states:
		for(Iter i=qlist.begin(); i!=qlist.end(); i++)
		{	//start at the next element
			Iter j=i; j++;
			while(j!=qlist.end())
			{	bool found = false;
				for(const matrix3<int>& m: sym)
					if(circDistanceSquared((~m)*i->k, invert*j->k)<symmThresholdSq)
					{	found = true;
						if(invert<0) usedInversion = true;
						break;
					}
				if(found)
				{	i->weight += j->weight;
					j = qlist.erase(j); //after this j points to the element just after the one rmeoved
				}
				else j++;
			}
		}
	}
	if(usedInversion) logPrintf("Adding inversion symmetry to k-mesh for non-inversion-symmetric unit cell.\n");
	else invertList.resize(1); //drop explicit inversion if not required
	((Symmetries*)this)->kpointInvertList = invertList; //Set kpointInvertList
	return qlist;
}

//Symmetrize scalar fields:
void symmetrize_sub(size_t iStart, size_t iStop, int nRot, double* x, int* symmIndex)
{	double nrotInv = 1.0/nRot;
	for(size_t i=iStart; i<iStop; i++)
	{	double xSum=0.0;
		for(int j=0; j<nRot; j++) xSum += x[symmIndex[nRot*i+j]];
		xSum *= nrotInv; //average n in the equivalence class
		for(int j=0; j<nRot; j++) x[symmIndex[nRot*i+j]] = xSum;
	}
}
#ifdef GPU_ENABLED //implemented in Symmetries.cu:
extern void symmetrize_gpu(int nSymmClasses, int nRot, double* x, int* symmIndex);
#endif
void Symmetries::symmetrize(DataRptr& x) const
{	if(sym.size()==1) return; // No symmetries, nothing to do
	int nSymmClasses = nSymmIndex / sym.size(); //number of equivalence classes
	#ifdef GPU_ENABLED
	symmetrize_gpu(nSymmClasses, sym.size(), x->dataGpu(), symmIndex);
	#else
	threadLaunch(symmetrize_sub, nSymmClasses, sym.size(), x->data(), symmIndex);
	#endif
}

//Symmetrize forces:
void Symmetries::symmetrize(IonicGradient& f) const
{	if(sym.size() <= 1) return;
	for(unsigned sp=0; sp<f.size(); sp++)
	{	std::vector<vector3<> > tempForces(f[sp].size());
		for(unsigned atom=0; atom<f[sp].size(); atom++)
			for(unsigned iRot=0; iRot<sym.size(); iRot++)
				tempForces[atom] += (~sym[iRot]) * f[sp][atomMap[sp][atom][iRot]];
		for(unsigned atom=0; atom<f[sp].size(); atom++)
			f[sp][atom] = tempForces[atom] / sym.size();
	}
}

const std::vector< matrix3<int> >& Symmetries::getMatrices() const
{	return sym;
}
const std::vector< matrix3<int> >& Symmetries::getMeshMatrices() const
{	return symMesh;
}
const std::vector<int>& Symmetries::getKpointInvertList() const
{	return kpointInvertList;
}


void Symmetries::calcSymmetries()
{
	const IonInfo& iInfo = e->iInfo;
	logPrintf("Searching for point group symmetries:\n");

	//Find symmetries of bravais lattice
	matrix3<> Rreduced; matrix3<int> transmission;
	std::vector<matrix3<int>> symLattice = getSymmetries(e->gInfo.R, e->coulombParams.isTruncated(), &Rreduced, &transmission);
	if(nrm2(Rreduced - e->gInfo.R) > symmThreshold * nrm2(Rreduced)) //i.e. if R != Rreduced
	{	logPrintf("Non-trivial transmission matrix:\n"); transmission.print(globalLog," %2d ");
		logPrintf("with reduced lattice vectors:\n"); Rreduced.print(globalLog," %12.6f ");
	}
	logPrintf("\n%lu symmetries of the bravais lattice\n", symLattice.size()); logFlush();

	//Find symmetries commensurate with atom positions:
	vector3<> rCenter;
	sym = basisReduce(symLattice, rCenter);
	logPrintf("reduced to %lu symmetries with basis\n", sym.size());
	
	//Make sure identity is the first symmetry
	sortSymmetries();

	//Print symmetry matrices
	if(shouldPrintMatrices)
	{	for(const matrix3<int>& m: sym)
		{	m.print(globalLog, " %2d ");
			logPrintf("\n");
		}
	}
	logFlush();
	
	if(shouldMoveAtoms) //Check for better symmetry centers:
	{	std::vector< vector3<> > rCenterCandidates;
		//Check atom positions and midpoints of atom pairs as candidate symmetry centers:
		for(auto sp: iInfo.species)
			for(unsigned n1 = 0; n1 < sp->atpos.size(); n1++)
			{	rCenterCandidates.push_back(sp->atpos[n1]);
				for(unsigned n2 = 0; n2 < n1; n2++)
					rCenterCandidates.push_back(0.5*(sp->atpos[n1] + sp->atpos[n2]));
			}
		//Check if any of the candidates leads to more symmetries than current rCenter:
		size_t origSymSize = sym.size();
		for(vector3<> rProposed: rCenterCandidates)
		{	std::vector< matrix3<int> > symTemp = basisReduce(symLattice, rProposed);
			if(symTemp.size() > sym.size())
			{	rCenter = rProposed;
				sym = symTemp;
			}
		}
		//Print positions and quit if a better symmetry center is found:
		if(sym.size()>origSymSize)
		{	logPrintf(
				"\nTranslating atoms by [ %lg %lg %lg ] (in lattice coordinates) will\n"
				"increase symmetry count from %lu to %lu. Translated atom positions follow:\n",
				-rCenter[0], -rCenter[1], -rCenter[2], origSymSize, sym.size());
			for(auto sp: iInfo.species)
				for(vector3<>& pos: sp->atpos)
					pos -= rCenter;
			iInfo.printPositions(globalLog);
			die("Use the suggested ionic positions, or set <moveAtoms>=no in command symmetry.\n");
		}
	}
}

std::vector< matrix3<int> > Symmetries::basisReduce(const std::vector< matrix3<int> >& symLattice, vector3<> offset) const
{	std::vector< matrix3<int> > symBasis;
	//Loop over lattice symmetries:
	for(const matrix3<int>& m: symLattice)
	{	bool symmetric = true;
		for(auto sp: e->iInfo.species) //For each species
		{	for(unsigned a1=0; a1<sp->atpos.size(); a1++) //For each atom
			{	bool foundImage = false;
				vector3<> mapped_pos1 = offset + m * (sp->atpos[a1]-offset);
				for(unsigned a2=0; a2<sp->atpos.size(); a2++)
					if(circDistanceSquared(mapped_pos1, sp->atpos[a2]) < symmThresholdSq) //Check if the image exists
						if(!sp->initialMagneticMoments.size()
							|| (sp->initialMagneticMoments[a1]==sp->initialMagneticMoments[a2]) ) //Check spin symmetry
						{
							foundImage = true;
							break;
						}
				if(!foundImage) { symmetric = false; break; }
			}
			if(!symmetric) break;
		}
		if(symmetric) //For each species, each atom maps onto another
			symBasis.push_back(m);
	}
	return symBasis;
}

void Symmetries::checkKmesh() const
{	//Find subgroup of sym which leaves k-mesh invariant
	std::vector< matrix3<int> > symKmesh;
	for(const matrix3<int>& m: sym)
	{	bool symmetric = true;
		//For each k-point, search if the image under m belongs to the k-mesh
		for(const QuantumNumber& q1: e->eInfo.qnums)
		{	bool foundImage = false;
			for(const QuantumNumber& q2: e->eInfo.qnums)
				if(circDistanceSquared((~m)*q1.k,q2.k)<symmThresholdSq
					&& fabs(q1.weight-q2.weight)<symmThreshold)
				{	foundImage = true;
					break;
				}
			if(!foundImage) { symmetric = false; break; }
		}
		if(symmetric) symKmesh.push_back(m); //m maps k-mesh onto itself
	}
	
	if(symKmesh.size() < sym.size())
	{	logPrintf("\nWARNING: k-mesh symmetries are a subgroup of size %lu\n", symKmesh.size());
		if(shouldPrintMatrices)
		{	for(const matrix3<int>& m: symKmesh)
			{	m.print(globalLog, " %2d ");
				logPrintf("\n");
			}
		}
		logPrintf("The effectively sampled k-mesh is a superset of the specified one,\n"
			"and the answers need not match those with symmetries turned off.\n");
	}
}

void Symmetries::initSymmIndex()
{	const GridInfo& gInfo = e->gInfo;
	if(sym.size()==1) return;

	std::vector<int> symmIndexVec;
	symmIndexVec.reserve(gInfo.nr);
	vector3<int> r;
	std::vector<bool> done(gInfo.nr, false);
	//Loop over all points not already handled as an image of a previous one:
	for(r[0]=0; r[0]<gInfo.S[0]; r[0]+=1)
		for(r[1]=0; r[1]<gInfo.S[1]; r[1]+=1)
			for(r[2]=0; r[2]<gInfo.S[2]; r[2]+=1)
			{	int index = gInfo.fullRindex(r);
				if(!done[index])
				{	//Loop over symmetry matrices:
					for(const matrix3<int>& m: symMesh)
					{	vector3<int> rNew = m * r;
						//project back into range:
						for(int i=0; i<3; i++)
							rNew[i] = rNew[i] % gInfo.S[i];
						int index2 = gInfo.fullGindex(rNew); //fullGindex handles wrapping negative indices
						symmIndexVec.push_back(index2);
						done[index2] = true;
					}
				}
			}
	//Set the final pointers:
	nSymmIndex = symmIndexVec.size();
	#ifdef GPU_ENABLED
	cudaMalloc(&symmIndex, nSymmIndex*sizeof(int));
	cudaMemcpy(symmIndex, &symmIndexVec[0], nSymmIndex*sizeof(int), cudaMemcpyHostToDevice);
	#else
	symmIndex = new int[nSymmIndex];
	memcpy(symmIndex, &symmIndexVec[0], nSymmIndex*sizeof(int));
	#endif
}

void Symmetries::sortSymmetries()
{	//Ensure first matrix is identity:
	matrix3<int> identity(1,1,1);
	for(unsigned i=1; i<sym.size(); i++)
		if(sym[i]==identity)
			std::swap(sym[0], sym[i]);
}


void Symmetries::checkFFTbox()
{	const vector3<int>& S = e->gInfo.S;
	symMesh.resize(sym.size());
	for(unsigned iRot = 0; iRot<sym.size(); iRot++)
	{	//the mesh coordinate symmetry matrices are Diag(S) * m * Diag(inv(S))
		//and these must be integral for the mesh to be commensurate:
		symMesh[iRot] = Diag(S) * sym[iRot];
		//Right-multiply by Diag(inv(S)) and ensure integer results:
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				if(symMesh[iRot](i,j) % S[j] == 0)
					symMesh[iRot](i,j) /= S[j];
				else
				{	logPrintf("FFT box not commensurate with symmetry matrix:\n");
					sym[iRot].print(globalLog, " %2d ");
					die("FFT box not commensurate with symmetries\n");
				}
	}
}


void Symmetries::checkSymmetries() const
{	logPrintf("Checking manually specified symmetry matrices.\n");
	for(const matrix3<int>& m: sym) //For each symmetry matrix
		for(auto sp: e->iInfo.species) //For each species
			for(auto pos1: sp->atpos) //For each atom
			{	bool foundImage = false;
				vector3<> mapped_pos1 = m * pos1;
				for(auto pos2: sp->atpos) //Check if the image exists:
					if(circDistanceSquared(mapped_pos1, pos2) < symmThresholdSq)
					{	foundImage = true;
						break;
					}
				if(!foundImage) die("Symmetries do not agree with atomic positions!\n");
			}
}

void Symmetries::initAtomMaps()
{	const IonInfo& iInfo = e->iInfo;
	if(sym.size()==1) return;
	if(shouldPrintMatrices) logPrintf("\nMapping of atoms according to symmetries:\n");
	atomMap.resize(iInfo.species.size());
	
	for(unsigned sp = 0; sp < iInfo.species.size(); sp++)
	{
		const SpeciesInfo& spInfo = *(iInfo.species[sp]);
		atomMap[sp].resize(spInfo.atpos.size());
		
		for(unsigned at1=0; at1<spInfo.atpos.size(); at1++)
		{	if(shouldPrintMatrices) logPrintf("%s %3u: ", spInfo.name.c_str(), at1);
			atomMap[sp][at1].resize(sym.size());
			
			for(unsigned iRot = 0; iRot<sym.size(); iRot++)
			{	vector3<> mapped_pos1 = sym[iRot] * spInfo.atpos[at1];
				
				for(unsigned at2=0; at2<spInfo.atpos.size(); at2++)
					if(circDistanceSquared(mapped_pos1, spInfo.atpos[at2]) < symmThresholdSq)
					{	atomMap[sp][at1][iRot] = at2;
				
						if(not spInfo.constraints[at1].isEquivalent(spInfo.constraints[at2], e->gInfo.R*sym[iRot]*e->gInfo.invR))
							die("Species %s atoms %u and %u are related by symmetry "
							"but have different move scale factors or inconsistent move constraints.\n\n",
								spInfo.name.c_str(), at1, at2);
							
					}
				
				if(shouldPrintMatrices) logPrintf(" %3u", atomMap[sp][at1][iRot]);
			}
			if(shouldPrintMatrices) logPrintf("\n");
		}
	}
	logFlush();
}
