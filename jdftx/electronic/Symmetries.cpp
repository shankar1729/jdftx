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

#include <electronic/Symmetries.h>
#include <electronic/Everything.h>
#include <electronic/SpeciesInfo.h>
#include <electronic/ElecInfo.h>
#include <electronic/IonInfo.h>
#include <electronic/Basis.h>
#include <core/SphericalHarmonics.h>
#include <core/LatticeUtils.h>
#include <core/GridInfo.h>
#include <core/Thread.h>
#include <fluid/Euler.h>

static const int lMaxSpherical = 3;

Symmetries::Symmetries() : symSpherical(lMaxSpherical+1), symSpinAngle(lMaxSpherical+1), sup(vector3<int>(1,1,1)), isPertSup(false)
{	shouldPrintMatrices = false;
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
		default: //No symmetry (only operation is identity)
			sym.assign(1, SpaceGroupOp()); 
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
	for(const SpaceGroupOp& op: sym)
		if(op.rot==matrix3<int>(-1,-1,-1))
		{	invertList.resize(1); //inversion explicitly found in symmetry list, so remove from invertList
			break;
		}
	//Compile kpoint map:
	std::vector<unsigned long long>& kmap = ((Symmetries*)this)->kmap;
	kmap.assign(qnums.size(), ~0ULL); //list of source k-point and symmetry operation (ordered to prefer no inversion, earliest k-point and then earliest symmetry matrix)
	std::vector<int> isSymKmesh(sym.size(), true); //whether each symmetry matrix leaves the k-mesh invariant
	size_t iSrcStart, iSrcStop;
	TaskDivision(qnums.size(), mpiWorld).myRange(iSrcStart, iSrcStop);
	PeriodicLookup<QuantumNumber> plook(qnums, e->gInfo.GGT);
	for(size_t iSrc=iSrcStart; iSrc<iSrcStop; iSrc++)
		for(int invert: invertList)
			for(int iSym=0; iSym<int(sym.size()); iSym++)
			{	size_t iDest = plook.find(invert * qnums[iSrc].k * sym[iSym].rot);
				if(iDest != string::npos)
					kmap[iDest] = std::min(kmap[iDest], kmapPack(iSrc, invert, iSym));
				else
				{	if(invert>0) isSymKmesh[iSym] = false;
				}
			}
	//Sync map across processes
	mpiWorld->allReduce(kmap.data(), kmap.size(), MPIUtil::ReduceMin);
	mpiWorld->allReduce(isSymKmesh.data(), isSymKmesh.size(), MPIUtil::ReduceLAnd);
	//Print symmetry-incommensurate kmesh warning if necessary:
	size_t nSymKmesh = std::count(isSymKmesh.begin(), isSymKmesh.end(), true);
	if(nSymKmesh < sym.size()) //if even one of them is false
	{	logPrintf("\nWARNING: k-mesh symmetries are a subgroup of size %lu\n", nSymKmesh);
		if(shouldPrintMatrices)
		{	for(int iSym=0; iSym<int(sym.size()); iSym++)
				if(isSymKmesh[iSym])
				{	sym[iSym].rot.print(globalLog, " %2d ");
					sym[iSym].a.print(globalLog, " %lg ");
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
void Symmetries::symmetrize(ScalarField& x) const
{	if(sym.size()==1) return; // No symmetries, nothing to do
	complexScalarFieldTilde xTilde = J(Complex(x));
	symmetrize(xTilde);
	x = Real(I(xTilde));
}
void Symmetries::symmetrize(ScalarFieldTilde& x) const
{	if(sym.size()==1) return; // No symmetries, nothing to do
	complexScalarFieldTilde xComplex = Complex(x);
	symmetrize(xComplex);
	x = Real(xComplex);
}
void Symmetries::symmetrize(complexScalarFieldTilde& x) const
{	if(sym.size()==1) return; // No symmetries, nothing to do
	int nSymmClasses = symmIndex.nData() / sym.size(); //number of equivalence classes
	callPref(eblas_symmetrize)(nSymmClasses, sym.size(), symmIndex.dataPref(), symmMult.dataPref(), symmIndexPhase.dataPref(), x->dataPref());
}

//Symmetrize forces:
void Symmetries::symmetrize(IonicGradient& f) const
{	if(sym.size() <= 1) return;
	for(unsigned sp=0; sp<f.size(); sp++)
	{	std::vector<vector3<> > tempForces(f[sp].size());
		for(unsigned atom=0; atom<f[sp].size(); atom++)
			for(unsigned iRot=0; iRot<sym.size(); iRot++)
				tempForces[atom] += (~sym[iRot].rot) * f[sp][atomMap[sp][atom][iRot]];
		for(unsigned atom=0; atom<f[sp].size(); atom++)
			f[sp][atom] = tempForces[atom] / sym.size();
	}
}

//Symmetrize Ylm-basis matrices:
void Symmetries::symmetrizeSpherical(matrix& X, const SpeciesInfo* specie) const
{	//Find index of specie (so as to access atom map)
	unsigned sp = 0;
	for(sp=0; sp<e->iInfo.species.size(); sp++)
		if(e->iInfo.species[sp].get() == specie)
			break;
	int nAtoms = atomMap[sp].size();
	int spinorLength = e->eInfo.spinorLength();
	int l = (X.nRows()/(nAtoms*spinorLength)-1)/2; //matrix dimension = (2l+1)*nAtoms*spinorLength
	int orbCount = (2*l+1)*spinorLength;
	int nTot = orbCount*nAtoms;
	assert(X.nCols()==nTot);
	if(!l || sym.size()==1) return; //symmetries do nothing
	const std::vector<matrix>& sym_l = getSphericalMatrices(l, specie->isRelativistic());
	matrix result;
	for(unsigned iRot=0; iRot<sym_l.size(); iRot++)
	{	//Construct transformation matrix including atom maps:
		matrix m = zeroes(nTot, nTot);
		for(int atom=0; atom<nAtoms; atom++)
		{	int atomOut = atomMap[sp][atom][iRot];
			m.set(atomOut*orbCount,(atomOut+1)*orbCount, atom*orbCount,(atom+1)*orbCount, sym_l[iRot]);
		}
		//Apply
		result += m * X * dagger(m);
	}
	X = (1./sym_l.size()) * result;
}


const std::vector<SpaceGroupOp>& Symmetries::getMatrices() const
{	return sym;
}

const std::vector<matrix>& Symmetries::getSphericalMatrices(int l, bool relativistic) const
{	if(l>lMaxSpherical) die("l=%d > lMax=%d supported for density matrix symmetrization\n", l, lMaxSpherical);
	bool jSplit = relativistic && l>0;
	const std::vector< std::vector<matrix> >& cache = jSplit ? symSpinAngle : symSpherical;
	if(!cache[l].size()) //Not yet initialized, do so now:
	{	//Create a basis of unit vectors for which Ylm are linearly independent:
		int mCount = 2*l+1;
		int sCount = e->eInfo.spinorLength();
		int msCount = mCount * sCount;
		std::vector< vector3<> > nHat(mCount);
		nHat[0] = vector3<>(0,0,1);
		for(int m=1; m<=l; m++)
		{	double phi = 2./l; //chosen empirically to get small basis-matrix condition numbers for l<=3
			double theta = m*2./l;
			nHat[2*m-1] = vector3<>(sin(theta), 0, cos(theta));
			nHat[2*m-0] = vector3<>(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
		}
		//Optionally construct the spinor transformations:
		matrix Uspin;
		if(jSplit)
		{	Uspin = zeroes(msCount, msCount);
			Uspin.set(0,2*l,       0,msCount, dagger(SpeciesInfo::getYlmToSpinAngleMatrix(l, 2*l-1)));
			Uspin.set(2*l,msCount, 0,msCount, dagger(SpeciesInfo::getYlmToSpinAngleMatrix(l, 2*l+1)));
		}
		//Construct basis matrix at nHat:
		matrix bOrig = zeroes(msCount, msCount); complex* bOrigData = bOrig.data();
		for(int nIndex=0; nIndex<mCount; nIndex++)
			for(int m=-l; m<=l; m++)
				for(int s=0; s<sCount; s++)
					bOrigData[bOrig.index((l+m)*sCount+s,nIndex*sCount+s)] = Ylm(l, m, nHat[nIndex]);
		if(Uspin) bOrig = Uspin * bOrig;
		matrix bOrigInv = inv(bOrig);
		//For each rotation matrix, construct rotated basis, and deduce rotation matrix in l,m basis:
		std::vector<matrix>& out = (std::vector<matrix>&)cache[l];
		out.resize(sym.size());
		for(unsigned iRot=0; iRot<sym.size(); iRot++)
		{	matrix3<> rot = e->gInfo.R * sym[iRot].rot * inv(e->gInfo.R); //cartesian rotation matrix
			matrix bRot = zeroes(msCount, msCount); complex* bRotData = bRot.data();
			for(int nIndex=0; nIndex<mCount; nIndex++)
				for(int m=-l; m<=l; m++)
					for(int s=0; s<sCount; s++)
						bRotData[bRot.index((l+m)*sCount+s,nIndex*sCount+s)] = Ylm(l, m, rot * nHat[nIndex]);
			if(Uspin)
			{	matrix rotSpinor = getSpinorRotation(rot); //generate spinor rotation matrix
				bRot = Uspin * (tiledBlockMatrix(rotSpinor, mCount) * bRot);
			}
			out[iRot] = bRot * bOrigInv;
		}
	}
	return cache[l];
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

matrix Symmetries::getSpinorRotation(const matrix3<>& rot)
{	vector3<> euler = eulerFromMatrix(rot * (1./det(rot)));  //pure rotation (remove inversion which doesn't matter in the following)
	//R_z(alpha):
	matrix A = zeroes(2,2);
		complex cA = cis(0.5*euler[0]);
		A.set(0,0, cA);
		A.set(1,1, cA.conj());
	//R_y(beta):
	matrix B = zeroes(2,2);
		complex cB = cis(0.5*euler[1]);
		B.set(0,0,  cB.real()); B.set(0,1, cB.imag());
		B.set(1,0, -cB.imag()); B.set(1,1, cB.real());
	//R_z(gamma):
	matrix G = zeroes(2,2);
		complex cG = cis(0.5*euler[2]);
		G.set(0,0, cG);
		G.set(1,1, cG.conj());
	return A * B * G;
}

void Symmetries::calcSymmetries()
{
	//Find symmetries of bravais lattice
	matrix3<> Rreduced; matrix3<int> transmission;
	std::vector<matrix3<int>> symLattice = getSymmetries(e->gInfo.R, e->coulombParams.isTruncated(), &Rreduced, &transmission);
	logPrintf("\n");
	if(nrm2(Rreduced - e->gInfo.R) > symmThreshold * nrm2(Rreduced)) //i.e. if R != Rreduced
	{	logPrintf("Non-trivial transmission matrix:\n"); transmission.print(globalLog," %2d ");
		logPrintf("with reduced lattice vectors:\n"); Rreduced.print(globalLog," %12.6f ");
		logPrintf("\n");
	}
	logPrintf("Found %lu point-group symmetries of the bravais lattice\n", symLattice.size()); logFlush();

	//Find symmetries commensurate with atom positions:
	sym = findSpaceGroup(symLattice);
	logPrintf("Found %lu space-group symmetries with basis%s\n", sym.size(),
		sup==vector3<int>(1,1,1) ? "" : " (with translations restricted to unit cells)"); //clarify constraint in phonon case
	
	//Find symmetries commensurate with external electric field (if any):
	if(e->coulombParams.Efield.length_squared())
	{	std::vector<SpaceGroupOp> symNew;
		vector3<> RT_Efield_ramp, RT_Efield_wave;
		e->coulombParams.splitEfield(e->gInfo.R, RT_Efield_ramp, RT_Efield_wave);
		for(const SpaceGroupOp& op: sym)
		{	if((RT_Efield_ramp - RT_Efield_ramp * op.rot).length() > RT_Efield_ramp.length()*symmThreshold) continue; //does not leave E_ramp invariant
			if((RT_Efield_wave - RT_Efield_wave * op.rot).length() > RT_Efield_wave.length()*symmThreshold) continue; //does not leave E_wave invariant
			if(fabs(dot(op.a, RT_Efield_wave)) > RT_Efield_wave.length()*symmThreshold) continue; //no longer a symmetry because of translation along E_wave
			symNew.push_back(op); //leaves Efield invariant
		}
		sym = symNew;
		logPrintf("reduced to %lu space-group symmetries with electric field\n", sym.size());
	}
	
	//Make sure identity is the first symmetry
	sortSymmetries();

	//Print symmetry matrices
	if(shouldPrintMatrices)
	{	for(const SpaceGroupOp& op: sym)
		{	op.rot.print(globalLog, " %2d ");
			op.a.print(globalLog, " %lg ");
			logPrintf("\n");
		}
	}
	logFlush();
}

bool magMomEquivalent(const vector3<>& a, const vector3<>& b)
{	return (a-b).length_squared() < symmThresholdSq;
}

std::vector<SpaceGroupOp> Symmetries::findSpaceGroup(const std::vector< matrix3<int> >& symLattice) const
{	std::vector<SpaceGroupOp> spaceGroup;
	//Loop over lattice symmetries:
	for(const matrix3<int>& rot: symLattice)
	{	
		//Determine offsets after this rotation that map the structure onto itself 
		std::vector<vector3<>> aArr; //list of candidates for the offset
		bool firstAtom = true;
		for(auto sp: e->iInfo.species)
		{	const std::vector< vector3<> >* M = sp->initialMagneticMoments.size() ? &sp->initialMagneticMoments : 0;
			for(size_t a1=0; a1<sp->atpos.size(); a1++)
			{
				//Generate list of offsets that would work for current atom by itself:
				std::vector<vector3<>> aCur;
				PeriodicLookup<vector3<>> plook(aCur, (~e->gInfo.R) * e->gInfo.R);
				vector3<> pos1rot = rot*sp->atpos[a1]; //rotated version of a1 position
				vector3<> M1rot; if(M) M1rot = (e->eInfo.spinType==SpinVector ? rot*(*M)[a1] : (*M)[a1]); //original or rotated M[a1] depending on spin type
				for(size_t a2=0; a2<sp->atpos.size(); a2++)
					if( (!M) || magMomEquivalent(M1rot, (*M)[a2]) )
					{	vector3<> dpos = Diag(sup) * (sp->atpos[a2] - pos1rot); //note in unit cell coordinates (matters if this is a phonon supercell)
						for(int k=0; k<3; k++) dpos[k] -= floor(0.5+dpos[k]); //wrap offset to base cell
						if(plook.find(dpos) == string::npos) //keep offsets unique modulo unit cell (rather than supercell in the phonon case)
						{	plook.addPoint(aCur.size(), dpos);
							aCur.push_back(dpos);
						}
					}
				
				//Intersect current candidates with global list:
				if(firstAtom)
				{	aArr = aCur; //no previous list; intersection = current
					firstAtom = false;
					continue;
				}
				std::vector<vector3<>> aNext; //intersection results
				if(aCur.size())
				{	for(vector3<> a: aArr)
						if(plook.find(a) != string::npos)
							aNext.push_back(a);
				}
				aArr = aNext;
				
				if(!aArr.size()) break;
			}
			if(!aArr.size()) break;
		}
		
		//Special handling for system with no atoms:
		if(firstAtom)
		{	spaceGroup.push_back(SpaceGroupOp(rot, vector3<>())); //space group = point group
			continue;
		}
		
		//Refine offsets:
		for(vector3<>& a: aArr)
		{	a = inv(Diag(vector3<>(sup))) * a; //switch offset back to current cell coordinates (matters if this is a phonon supercell)
			vector3<> daSum; int nAtoms = 0.;
			for(auto sp: e->iInfo.species) //For each species
			{	PeriodicLookup< vector3<> > plook(sp->atpos, (~e->gInfo.R) * e->gInfo.R);
				const std::vector< vector3<> >* M = sp->initialMagneticMoments.size() ? &sp->initialMagneticMoments : 0;
				for(size_t a1=0; a1<sp->atpos.size(); a1++) //For each atom
				{	vector3<> pos1rot = rot*sp->atpos[a1] + a; //now including offset
					vector3<> M1rot; if(M) M1rot = (e->eInfo.spinType==SpinVector ? rot*(*M)[a1] : (*M)[a1]); //original or rotated M[a1] depending on spin type
					size_t a2 = plook.find(pos1rot, M1rot, M, magMomEquivalent); //match position and magentic moment
					assert(a2 != string::npos); //the above algorithm should guarantee this
					vector3<> da = sp->atpos[a2] - pos1rot;
					for(int k=0; k<3; k++) da[k] -= floor(0.5+da[k]);
					daSum += da; nAtoms++;
				}
			}
			a += daSum / nAtoms;
			spaceGroup.push_back(SpaceGroupOp(rot, a));
		}
	}
	return spaceGroup;
}

void Symmetries::initSymmIndex()
{	const GridInfo& gInfo = e->gInfo;
	if(sym.size()==1) return;

	std::vector<int> symmIndexVec, symmMultVec;
	std::vector<complex> symmIndexPhaseVec;
	symmIndexVec.reserve(gInfo.nr);
	symmMultVec.reserve(gInfo.nr / sym.size());
	symmIndexPhaseVec.reserve(gInfo.nr);
	std::vector<bool> done(gInfo.nr, false); //use full G-space for symmetrization
	//Loop over all points not already handled as an image of a previous one:
	{	const vector3<int>& S = gInfo.S;
		size_t iStart = 0, iStop = gInfo.nr;
		THREAD_fullGspaceLoop
		(	if(!done[i])
			{	std::set<int> orbit;
				//Loop over symmetry matrices:
				for(const SpaceGroupOp& op: sym)
				{	vector3<int> iG2 = iG * op.rot;
					complex phase = cis((-2*M_PI)*dot(iG,op.a));
					//project back into range:
					for(int k=0; k<3; k++)
					{	iG2[k] = positiveRemainder(iG2[k], S[k]);
						if(2*iG2[k]>S[k]) iG2[k]-=S[k];
					}
					int i2 = gInfo.fullGindex(iG2);
					symmIndexPhaseVec.push_back(phase);
					symmIndexVec.push_back(i2);
					done[i2] = true;
					orbit.insert(i2);
				}
				int multiplicity = sym.size()/orbit.size(); //number of times each point in orbit is covered
				if(multiplicity * orbit.size() != sym.size())
				{	die("\nSymmetry operations do not seem to form a group.\n"
						"This is most likely because the geometry has some border-line symmetries.\n"
						"Try either tightening or loosening the symmetry-threshold parameter.\n\n");
				}
				symmMultVec.push_back(multiplicity);
			}
		)
	}
	//Set the final pointers:
	int nSymmIndex = symmIndexVec.size();
	symmIndex.init(nSymmIndex);
	symmMult.init(symmMultVec.size());
	symmIndexPhase.init(nSymmIndex);
	memcpy(symmIndex.data(), &symmIndexVec[0], nSymmIndex*sizeof(int));
	memcpy(symmMult.data(), &symmMultVec[0], symmMultVec.size()*sizeof(int));
	memcpy(symmIndexPhase.data(), &symmIndexPhaseVec[0], nSymmIndex*sizeof(complex));
}

void Symmetries::sortSymmetries()
{	//Ensure first matrix is identity:
	SpaceGroupOp id;
	for(unsigned i=1; i<sym.size(); i++)
		if(sym[i].rot==id.rot && sym[i].a==id.a)
			std::swap(sym[0], sym[i]);
}


void Symmetries::checkFFTbox()
{	const vector3<int>& S = e->gInfo.S;
	for(unsigned iRot = 0; iRot<sym.size(); iRot++)
	{	//the mesh coordinate symmetry matrices are Diag(S) * m * Diag(inv(S))
		//and these must be integral for the mesh to be commensurate:
		matrix3<int> mMesh = Diag(S) * sym[iRot].rot;
		//Right-multiply by Diag(inv(S)) and ensure integer results:
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				if(mMesh(i,j) % S[j] == 0)
					mMesh(i,j) /= S[j];
				else
				{	logPrintf("FFT box not commensurate with symmetry matrix:\n");
					sym[iRot].rot.print(globalLog, " %2d ");
					sym[iRot].a.print(globalLog, " %lg ");
					die("FFT box not commensurate with symmetries.\n");
				}
	}
	
	//Check embedded truncation center:
	if(e->coulombParams.embed)
	{	const vector3<>& c = e->coulombParams.embedCenter;
		for(const SpaceGroupOp& op: sym)
		{	vector3<> cRot = op.rot*c + op.a;
			for(int dir=0; dir<3; dir++)
				if(!e->coulombParams.isTruncated()[dir])
					cRot[dir] = c[dir]; //don't need invariance in periodic directions
			if(circDistanceSquared(c, cRot) > symmThresholdSq)
			{	logPrintf("Coulomb truncation embedding center is not invariant under symmetry matrix:\n");
				op.rot.print(globalLog, " %2d ");
				op.a.print(globalLog, " %lg ");
				die("Coulomb truncation embedding center is not invariant under symmetries.\n");
			}
		}
	}
}


void Symmetries::checkSymmetries()
{	logPrintf("Checking manually specified symmetry matrices.\n");
	std::vector<SpaceGroupOp> symReduced; //reduced symmetries for a perturbed supercell
	for(const SpaceGroupOp& op: sym) //For each symmetry matrix
	{	bool isPertSym = true;
		for(auto sp: e->iInfo.species) //For each species
		{	PeriodicLookup< vector3<> > plook(sp->atpos, (~e->gInfo.R) * e->gInfo.R);
			const std::vector< vector3<> >* M = sp->initialMagneticMoments.size() ? &sp->initialMagneticMoments : 0;
			for(size_t a1=0; a1<sp->atpos.size(); a1++) //For each atom
			{	vector3<> M1rot; if(M) M1rot = (e->eInfo.spinType==SpinVector ? op.rot*(*M)[a1] : (*M)[a1]); //original or rotated M[a1] depending on spin type
				if(string::npos == plook.find(op.rot * sp->atpos[a1] + op.a, M1rot, M, magMomEquivalent)) //match position and spin
				{	if(isPertSup)
					{	isPertSym = false;
						break;
					}
					logPrintf("Ionic positions not invariant under symmetry matrix:\n");
					op.rot.print(globalLog, " %2d ");
					op.a.print(globalLog, " %lg ");
					die("Symmetries do not agree with atomic positions!\n");
				}
			}
		}
		if(isPertSup && isPertSym) symReduced.push_back(op);
	}
	if(isPertSup)
	{	logPrintf("Reduced %lu manually specified symmetries of unit cell to %lu symmetries of perturbed supercell.\n", sym.size(), symReduced.size());
		std::swap(sym, symReduced);
	}
}

void Symmetries::initAtomMaps()
{	const IonInfo& iInfo = e->iInfo;
	if(shouldPrintMatrices) logPrintf("\nMapping of atoms according to symmetries:\n");
	atomMap.resize(iInfo.species.size());
	double datposSqSum = 0.; int nAtomsTot = 0; //counters for atom symmetrization statistics
	
	for(unsigned sp = 0; sp < iInfo.species.size(); sp++)
	{	const SpeciesInfo& spInfo = *(iInfo.species[sp]);
		atomMap[sp].resize(spInfo.atpos.size());
		PeriodicLookup< vector3<> > plook(spInfo.atpos, (~e->gInfo.R) * e->gInfo.R);
		std::vector<vector3<> > datpos(spInfo.atpos.size()); //Displacements to exactly symmetrize atpos
		
		for(size_t a1=0; a1<spInfo.atpos.size(); a1++)
		{	if(shouldPrintMatrices) logPrintf("%s %3lu: ", spInfo.name.c_str(), a1);
			atomMap[sp][a1].resize(sym.size());
			
			for(unsigned iRot = 0; iRot<sym.size(); iRot++)
			{	vector3<> idealPos = sym[iRot].rot * spInfo.atpos[a1] + sym[iRot].a;
				size_t a2 = plook.find(idealPos);
				if(a2 == string::npos)
					die("Atom positions are marginally symmetric (errors comparable to detection threshold).\n"
						"Use command symmetry-threshold to either increase tolerance and include marginal\n"
						"symmetries, or reduce tolerance and exclude marginal symmetries, as appropriate.\n\n");
				atomMap[sp][a1][iRot] = a2;
				if(not spInfo.constraints[a1].isEquivalent(spInfo.constraints[a2], e->gInfo.R*sym[iRot].rot*inv(e->gInfo.R)))
					die("\nSpecies %s atoms %lu and %lu are related by symmetry "
					"but have different move scale factors or inconsistent move constraints.\n\n",
						spInfo.name.c_str(), a1, a2);
				if(shouldPrintMatrices) logPrintf(" %3u", atomMap[sp][a1][iRot]);
				
				//Add contributions to symmetrization displacements:
				vector3<> dat = idealPos - spInfo.atpos[a2];
				for(int j=0; j<3; j++) dat[j] -= floor(0.5+dat[j]); //wrap to [-0.5,0.5)
				datpos[a2] += (1./sym.size()) * dat;
			}
			if(shouldPrintMatrices) logPrintf("\n");
		}
		
		//Symmetrize atoms:
		for(size_t a=0; a<spInfo.atpos.size(); a++)
		{	((IonInfo&)e->iInfo).species[sp]->atpos[a] += datpos[a];
			datposSqSum += (e->gInfo.R * datpos[a]).length_squared();
			nAtomsTot++;
		}
	}
	
	//Print atom symmetrization statistics:
	logPrintf("Applied RMS atom displacement %lg bohrs to make symmetries exact.\n", sqrt(datposSqSum/nAtomsTot));
	logFlush();
}
