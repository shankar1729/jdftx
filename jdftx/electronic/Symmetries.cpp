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
#include <electronic/SphericalHarmonics.h>
#include <core/LatticeUtils.h>
#include <core/GridInfo.h>
#include <core/Thread.h>

static const int lMaxSpherical = 3;

Symmetries::Symmetries() : symSpherical(lMaxSpherical+1)
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
	if(mode != SymmetriesNone) logPrintf("\n---------- Setting up symmetries ----------\n");

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
	initSymmIndex(); //Initialize the equivalence classes for scalar field symmetrization (using mesh matrices)
}

//Pack and unpack kpoint map entry to a single 64-bit integer
inline unsigned long long kmapPack(size_t iSrc, int invert, int iSym) { return (((unsigned long long)iSrc)<<8) | ((invert<0 ? 1 : 0)<<7) | iSym; }
inline void kmapUnpack(unsigned long long kmap, size_t& iSrc, int& invert, int& iSym)
{	iSrc = size_t(kmap >> 8);
	invert = (0x80 & kmap) ? -1 : +1;
	iSym = int(0x7F & kmap);
}

//Conversion functions needed by the PeriodicLookup template:
inline vector3<> getCoord(const QuantumNumber& qnum) { return qnum.k; } //for k-point reduction
inline vector3<> getCoord(const vector3<>& pos) { return pos; } //for atom symmetries

std::vector<QuantumNumber> Symmetries::reduceKmesh(const std::vector<QuantumNumber>& qnums) const
{	static StopWatch watch("reduceKmesh"); watch.start();
	if(mode == SymmetriesNone)
	{	((Symmetries*)this)->kpointInvertList.assign(1, +1);
		watch.stop();
		return qnums;
	}
	//Compile list of inversions to check:
	std::vector<int> invertList;
	invertList.push_back(+1);
	invertList.push_back(-1);
	for(const matrix3<int>& m: sym)
		if(m==matrix3<int>(-1,-1,-1))
		{	invertList.resize(1); //inversion explicitly found in symmetry list, so remove from invertList
			break;
		}
	//Compile kpoint map:
	std::vector<unsigned long long>& kmap = ((Symmetries*)this)->kmap;
	kmap.assign(qnums.size(), ~0ULL); //list of source k-point and symmetry operation (ordered to prefer no inversion, earliest k-point and then earliest symmetry matrix)
	std::vector<int> isSymKmesh(sym.size(), true); //whether each symmetry matrix leaves the k-mesh invariant
	size_t iSrcStart = (mpiUtil->iProcess() * qnums.size()) / mpiUtil->nProcesses();
	size_t iSrcStop = ((mpiUtil->iProcess()+1) * qnums.size()) / mpiUtil->nProcesses();
	PeriodicLookup<QuantumNumber> plook(qnums, e->gInfo.GGT);
	for(size_t iSrc=iSrcStart; iSrc<iSrcStop; iSrc++)
		for(int invert: invertList)
			for(int iSym=0; iSym<int(sym.size()); iSym++)
			{	size_t iDest = plook.find(invert * qnums[iSrc].k * sym[iSym]);
				if(iDest != string::npos)
					kmap[iDest] = std::min(kmap[iDest], kmapPack(iSrc, invert, iSym));
				else
				{	if(invert>0) isSymKmesh[iSym] = false;
				}
			}
	//Sync map across processes
	mpiUtil->allReduce(kmap.data(), kmap.size(), MPIUtil::ReduceMin);
	mpiUtil->allReduce(isSymKmesh.data(), isSymKmesh.size(), MPIUtil::ReduceLAnd);
	//Print symmetry-incommensurate kmesh warning if necessary:
	size_t nSymKmesh = std::count(isSymKmesh.begin(), isSymKmesh.end(), true);
	if(nSymKmesh < sym.size()) //if even one of them is false
	{	logPrintf("\nWARNING: k-mesh symmetries are a subgroup of size %lu\n", nSymKmesh);
		if(shouldPrintMatrices)
		{	for(int iSym=0; iSym<int(sym.size()); iSym++)
				if(isSymKmesh[iSym])
				{	sym[iSym].print(globalLog, " %2d ");
					logPrintf("\n");
				}
		}
		logPrintf("The effectively sampled k-mesh is a superset of the specified one,\n"
			"and the answers need not match those with symmetries turned off.\n");
	}
	//Compile set of source k-points and whether inversion is necessary:
	bool usedInversion = false;
	std::map<size_t,size_t> iSrcMap; //map from original to reduced kpoints
	for(unsigned long long kmapEntry: kmap)
	{	size_t iSrc; int invert, iSym;
		kmapUnpack(kmapEntry, iSrc, invert, iSym);
		iSrcMap[iSrc] = 0;
		if(invert<0) usedInversion = true;
	}
	size_t iReduced=0; for(auto& mapEntry: iSrcMap) mapEntry.second = (iReduced++);
	//Set invertList:
	if(usedInversion) logPrintf("Adding inversion symmetry to k-mesh for non-inversion-symmetric unit cell.\n");
	else invertList.resize(1); //drop explicit inversion if not required
	((Symmetries*)this)->kpointInvertList = invertList; //Set kpointInvertList
	//Compile list of reduced kpoints:
	std::vector<QuantumNumber> qRed(iSrcMap.size());
	for(const auto& mapEntry: iSrcMap) qRed[mapEntry.second] = qnums[mapEntry.first];
	//Update iSrc in map with the reduced value, and accumulate weights:
	for(size_t iDest=0; iDest<qnums.size(); iDest++)
	{	unsigned long long& kmapEntry = kmap[iDest];
		size_t iSrc; int invert, iSym;
		kmapUnpack(kmapEntry, iSrc, invert, iSym);
		size_t iReduced = iSrcMap[iSrc];
		if(iDest != iSrc) qRed[iReduced].weight += qnums[iDest].weight; //collect weight
		kmapEntry = kmapPack(iReduced, invert, iSym);
	}
	watch.stop();
	return qRed;
}

//Symmetrize scalar fields:
void Symmetries::symmetrize(DataRptr& x) const
{	if(sym.size()==1) return; // No symmetries, nothing to do
	int nSymmClasses = nSymmIndex / sym.size(); //number of equivalence classes
	callPref(eblas_symmetrize)(nSymmClasses, sym.size(), symmIndex, x->dataPref());
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

//Symmetrize Ylm-basis matrices:
void Symmetries::symmetrizeSpherical(matrix& X, int sp) const
{	int nAtoms = atomMap[sp].size();
	int l = (X.nRows()/nAtoms-1)/2; //matrix dimension = (2l+1)*nAtoms
	int nm = 2*l+1;
	int nTot = nm*nAtoms;
	assert(X.nCols()==nTot);
	if(!l || sym.size()==1) return; //symmetries do nothing
	const std::vector<matrix>& sym_l = getSphericalMatrices(l);
	matrix result;
	for(unsigned iRot=0; iRot<sym_l.size(); iRot++)
	{	//Construct transformation matrix including atom maps:
		matrix m = zeroes(nTot, nTot);
		for(int atom=0; atom<nAtoms; atom++)
			m.set(atomMap[sp][atom][iRot],nAtoms,nTot, atom,nAtoms,nTot, sym_l[iRot]);
		//Apply
		result += m * X * dagger(m);
	}
	X = (1./sym_l.size()) * result;
}


const std::vector< matrix3<int> >& Symmetries::getMatrices() const
{	return sym;
}
const std::vector< matrix3<int> >& Symmetries::getMeshMatrices() const
{	return symMesh;
}

const std::vector<matrix>& Symmetries::getSphericalMatrices(int l) const
{	if(l>lMaxSpherical) die("l=%d > lMax=%d supported for density matrix symmetrization\n", l, lMaxSpherical);
	if(!symSpherical[l].size()) //Not yet initialized, do so now:
	{	//Create a basis of unit vectors for which Ylm are linearly independent:
		std::vector< vector3<> > nHat(2*l+1);
		nHat[0] = vector3<>(0,0,1);
		for(int m=1; m<=l; m++)
		{	double phi = 2./l; //chosen empirically to get small basis-matrix condition numbers for l<=3
			double theta = m*2./l;
			nHat[2*m-1] = vector3<>(sin(theta), 0, cos(theta));
			nHat[2*m-0] = vector3<>(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
		}
		//Construct basis matrix at nHat:
		matrix bOrig(2*l+1, 2*l+1); complex* bOrigData = bOrig.data();
		for(unsigned nIndex=0; nIndex<nHat.size(); nIndex++)
			for(int m=-l; m<=l; m++)
				bOrigData[bOrig.index(l+m,nIndex)] = Ylm(l, m, nHat[nIndex]);
		matrix bOrigInv = inv(bOrig);
		//For each rotation matrix, construct rotated basis, and deduce rotation matrix in l,m basis:
		std::vector<matrix>& out = (((Symmetries*)this)->symSpherical)[l];
		out.resize(sym.size());
		for(unsigned iRot=0; iRot<sym.size(); iRot++)
		{	matrix3<> rot = e->gInfo.R * sym[iRot] * inv(e->gInfo.R); //cartesian rotation matrix
			matrix bRot(2*l+1, 2*l+1); complex* bRotData = bRot.data();
			for(unsigned nIndex=0; nIndex<nHat.size(); nIndex++)
				for(int m=-l; m<=l; m++)
					bRotData[bRot.index(l+m,nIndex)] = Ylm(l, m, rot * nHat[nIndex]);
			out[iRot] = bRot * bOrigInv;
		}
	}
	return symSpherical[l];
}


const std::vector<int>& Symmetries::getKpointInvertList() const
{	return kpointInvertList;
}

const std::vector< std::vector< std::vector<int> > >& Symmetries::getAtomMap() const
{	return atomMap;
}

void Symmetries::printKmap(FILE* fp) const
{	fprintf(fp, "#ReducedKpt #Symmetry inversion      (0-based indices, unreduced k-points in C array order)\n");
	int width = int(ceil(log(e->eInfo.nStates)/log(10)));
	for(const unsigned long long& kmapEntry: kmap)
	{	size_t iSrc; int invert, iSym;
		kmapUnpack(kmapEntry, iSrc, invert, iSym);
		fprintf(fp, "%*lu %2d %+d\n", width, iSrc, iSym, invert);
	}
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
		{	die("\nTranslating atoms by [ %lg %lg %lg ] (in lattice coordinates) will increase symmetry count\n"
				"from %lu to %lu. Update the ionic positions, or set <moveAtoms>=no in command symmetry.\n",
				-rCenter[0], -rCenter[1], -rCenter[2], origSymSize, sym.size());
		}
	}
}

bool magMomEquivalent(const vector3<>& a, const vector3<>& b)
{	return (a-b).length_squared() < symmThresholdSq;
}

std::vector< matrix3<int> > Symmetries::basisReduce(const std::vector< matrix3<int> >& symLattice, vector3<> offset) const
{	std::vector< matrix3<int> > symBasis;
	//Loop over lattice symmetries:
	for(const matrix3<int>& m: symLattice)
	{	bool symmetric = true;
		for(auto sp: e->iInfo.species) //For each species
		{	PeriodicLookup< vector3<> > plook(sp->atpos, e->gInfo.RTR);
			const std::vector< vector3<> >& M = sp->initialMagneticMoments;
			for(size_t a1=0; a1<sp->atpos.size(); a1++) //For each atom
			{	if(string::npos == plook.find(offset + m*(sp->atpos[a1]-offset), M.size()?M[a1]:vector3<>(), M.size()?&M:0, magMomEquivalent)) //match position and magentic moment
				{	symmetric = false;
					break;
				}
			}
			if(!symmetric) break;
		}
		if(symmetric) //For each species, each atom maps onto another
			symBasis.push_back(m);
	}
	return symBasis;
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
					die("FFT box not commensurate with symmetries.\n");
				}
	}
	
	//Check embedded truncation center:
	if(e->coulombParams.embed)
	{	const vector3<>& c = e->coulombParams.embedCenter;
		for(unsigned iRot = 0; iRot<sym.size(); iRot++)
			if(circDistanceSquared(c, sym[iRot]*c) > symmThresholdSq)
			{	logPrintf("Coulomb truncation embedding center is not invariant under symmetry matrix:\n");
				sym[iRot].print(globalLog, " %2d ");
				die("Coulomb truncation embedding center is not invariant under symmetries.\n");
			}
		
		//Find the nearest grid point to embedCenter that is commensurate with symmetries:
		vector3<int> iv0; for(int k=0; k<3; k++) iv0[k] = round(c[k]*S[k]);
		matrix3<> invDiagS = inv(Diag(vector3<>(S)));
		vector3<int> dv;
		bool done = false;
		for(int d=0; d<=(S[0]+S[1]+S[2])/2+1 && !done; d++) //search outwards (sorted by a Manhattan metric)
			for(dv[0]=-d; dv[0]<=d && !done; dv[0]++)
			{	int d0 = d - abs(dv[0]);
				for(dv[1]=-d0; dv[1]<=d0 && !done; dv[1]++)
				{	int d1 = d0 - abs(dv[1]);
					for(dv[2]=-d1; dv[2]<=d1 && !done; dv[2]+=2*std::max(1,d1)) //only points that satisfy abs(dv[0])+abs(dv[1])+abs(dv[2])==d
					{	vector3<int> iv = iv0 + dv;
						vector3<> x = invDiagS * iv;
						bool valid = true;
						for(const matrix3<int>& m: sym)
							if(circDistanceSquared(x, m*x) > symmThresholdSq)
							{	valid = false;
								break;
							}
						if(valid)
						{	((CoulombParams&)e->coulombParams).embedCenter = x;
							done = true; //terminate the search (4 loops)
						}
					}
				}
			}
		if(!done)
			die("Could not find a (integer) grid point to use as the truncation embedding center that\n"
				"is invariant under symmetries. HINT: center on the origin, or disable symmetries.\n");
	}
}


void Symmetries::checkSymmetries() const
{	logPrintf("Checking manually specified symmetry matrices.\n");
	for(const matrix3<int>& m: sym) //For each symmetry matrix
		for(auto sp: e->iInfo.species) //For each species
		{	PeriodicLookup< vector3<> > plook(sp->atpos, e->gInfo.RTR);
			const std::vector< vector3<> >& M = sp->initialMagneticMoments;
			for(size_t a1=0; a1<sp->atpos.size(); a1++) //For each atom
			{	if(string::npos == plook.find(m * sp->atpos[a1], M.size()?M[a1]:vector3<>(), M.size()?&M:0, magMomEquivalent)) //match position and spin
					die("Symmetries do not agree with atomic positions!\n");
			}
		}
}

void Symmetries::initAtomMaps()
{	const IonInfo& iInfo = e->iInfo;
	if(shouldPrintMatrices) logPrintf("\nMapping of atoms according to symmetries:\n");
	atomMap.resize(iInfo.species.size());
	
	for(unsigned sp = 0; sp < iInfo.species.size(); sp++)
	{	const SpeciesInfo& spInfo = *(iInfo.species[sp]);
		atomMap[sp].resize(spInfo.atpos.size());
		PeriodicLookup< vector3<> > plook(spInfo.atpos, e->gInfo.RTR);
		
		for(size_t a1=0; a1<spInfo.atpos.size(); a1++)
		{	if(shouldPrintMatrices) logPrintf("%s %3lu: ", spInfo.name.c_str(), a1);
			atomMap[sp][a1].resize(sym.size());
			
			for(unsigned iRot = 0; iRot<sym.size(); iRot++)
			{	size_t a2 = plook.find(sym[iRot] * spInfo.atpos[a1]);
				assert(a2 != string::npos);
				atomMap[sp][a1][iRot] = a2;
				if(not spInfo.constraints[a1].isEquivalent(spInfo.constraints[a2], e->gInfo.R*sym[iRot]*inv(e->gInfo.R)))
					die("Species %s atoms %lu and %lu are related by symmetry "
					"but have different move scale factors or inconsistent move constraints.\n\n",
						spInfo.name.c_str(), a1, a2);
				if(shouldPrintMatrices) logPrintf(" %3u", atomMap[sp][a1][iRot]);
			}
			if(shouldPrintMatrices) logPrintf("\n");
		}
	}
	logFlush();
}
