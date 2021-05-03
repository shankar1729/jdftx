/*-------------------------------------------------------------------
Copyright 2021 Ravishankar Sundararaman

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

#include <wannier/DefectSupercell.h>
#include <wannier/Wannier.h>
#include <commands/parser.h>
#include <core/LatticeUtils.h>
#include <core/WignerSeitz.h>
#include <core/ScalarFieldIO.h>
#include <libgen.h>

void DefectSupercell::initialize(const Wannier* wannier)
{	this->wannier = wannier;
	this->e = wannier->e;
	logPrintf("\n\n########### Initializing supercell for defect '%s' ###########\n\n", name.c_str());

	//Change directory:
	//--- first store present directory
	const int bufSize = 1024;
	char currentDir[bufSize];
	getcwd(currentDir, bufSize);
	//--- split inFile into directory and filename within:
	char dirBuf[2][bufSize];
	for(int k=0; k<2; k++)
		strcpy(dirBuf[k], inFile.c_str());
	const char* inFileDir = dirname(dirBuf[0]);
	const char* inFileName = basename(dirBuf[1]);
	//--- switch to inFileDir:
	logPrintf("Switching to directory '%s'\n", inFileDir);
	if(chdir(inFileDir) != 0)
		die("Failed to switch to defect input directory '%s'\n\n", inFileDir);
	
	//Read input file:
	logPrintf("Reading supercell input file '%s' ... ", inFileName);
	std::vector< std::pair<string,string> > input = readInputFile(inFileName);
	logPrintf("done."); logFlush();
	eSup = std::make_shared<WannierEverything>(); //WannierEverything needed here to handle default wannier commands (though unused)
	parse(input, *eSup);
	
	//Alter relevant settings and initialize:
	//--- Remove unnecessary initializations:
	eSup->eVars.skipWfnsInit = true; //skip because wavefunctions are set from unit cell calculation
	eSup->eVars.wfnsFilename.clear();
	eSup->eVars.eigsFilename.clear();
	eSup->eVars.fluidInitialStateFilename.clear();
	eSup->eInfo.initialFillingsFilename.clear();
	eSup->scfParams.historyFilename.clear();
	eSup->eVars.initLCAO = false; //only unit cell state needed
	//--- Make sure electronic potentials are read in:
	eSup->cntrl.fixed_H = true;
	eSup->eVars.nFilenamePattern.clear();
	eSup->eVars.VFilenamePattern = eSup->dump.format;
	//--- Initialize (only need grid, pseudopotentials and Vscloc):
	eSup->Everything::setup();
	
	//Switch back to original working directory:
	logPrintf("Switching back to directory '%s'\n", currentDir);
	if(chdir(currentDir) != 0)
		die("Failed to switch back to unit cell input directory '%s'\n\n", currentDir);
	
	//Determine / check supercells:
	logPrintf("\n---- Defect supercell ----\n");
	//--- determine input supercell
	matrix3<> supInTmp = inv(e->gInfo.R) * eSup->gInfo.R;
	supIn = round(vector3<>(supInTmp(0,0), supInTmp(1,1), supInTmp(2,2)));
	logPrintf("Input supercell: "); supIn.print(globalLog, " %d ");
	double err = nrm2(supInTmp - Diag(vector3<>(supIn)));
	if(err > symmThreshold)
	{	logPrintf("Linear combinations of unit cell in defect supercell:\n");
		supInTmp.print(globalLog, " %lf ");
		die("Defect supercell is not diagonal/integral within threshold %lg\n\n", symmThreshold);
	}
	//--- check output supercell
	logPrintf("Output supercell: "); supOut.print(globalLog, " %d ");
	for(int j=0; j<3; j++)
		if(!supOut[j] or (e->eInfo.kfold[j] % supOut[j]))
		{	die("kpoint folding %d is not a multiple of output supercell count %d for lattice direction %d.\n",
				e->eInfo.kfold[j], supOut[j], j);
		}
	
	//Map atom differences
	if(eSup->iInfo.coordsType == CoordsCartesian)
		xCenter = eSup->gInfo.invR * xCenter; //ensure lattice coordinates
	logPrintf("Center fractional coordinates: "); xCenter.print(globalLog, " %lf ");
	logPrintf("Atom changes from perfect to defective supercell:\n");
	logSuspend();
	wsSup = std::make_shared<WignerSeitz>(eSup->gInfo.R);
	logResume();
	atposRef.assign(eSup->iInfo.species.size(), std::vector<vector3<>>());
	U_rhoAtomRef.clear(); const matrix* U_rhoAtomUnit = e->eVars.U_rhoAtom.data();
	std::vector<bool> unitSpeciesFound(e->iInfo.species.size(), false);
	for(unsigned iSp=0; iSp<atposRef.size(); iSp++)
	{	SpeciesInfo& sp = *(eSup->iInfo.species[iSp]);
		//Map defect supercell atoms to Wigner-Seitz cell:
		for(vector3<>& x: sp.atpos)
			x = xCenter + wsSup->restrict(x - xCenter);
		//Find corresponding unit cell species and atoms:
		unsigned nAdded = sp.atpos.size(); //number of atoms of this species added relative to perfect supercell
		unsigned nMoved = 0; //number of atoms of this species moved from perfect supercell
		double rSqMoved = 0.; //total square distance of atom movement
		double rCutMove = 1.; //threshold displacement in bohrs above which atom is considered added/deleted rather than moved
		for(unsigned iSpUnit=0; iSpUnit<e->iInfo.species.size(); iSpUnit++)
			if(sp.name == e->iInfo.species[iSpUnit]->name)
			{	unitSpeciesFound[iSpUnit] = true;
				//Add reference atom positions for perfect supercell:
				const SpeciesInfo& spUnit = *(e->iInfo.species[iSpUnit]);
				matrix3<> invSupIn = inv(Diag(vector3<>(supIn)));
				vector3<int> iR;
				for(iR[0]=0; iR[0]<supIn[0]; iR[0]++)
				for(iR[1]=0; iR[1]<supIn[1]; iR[1]++)
				for(iR[2]=0; iR[2]<supIn[2]; iR[2]++)
				{	for(vector3<> pos: spUnit.atpos)
					{	vector3<> x = invSupIn * (pos + iR);
						//Find corresponding closest atom in defect supercell:
						double rSqSmallest = DBL_MAX; unsigned iClosest = -1;
						for(unsigned iAtom=0; iAtom<sp.atpos.size(); iAtom++)
						{	vector3<> dx = x - sp.atpos[iAtom];
							for(int dir=0; dir<3; dir++) dx[dir] -= floor(0.5 + dx[dir]); //minimum-image convention
							double rSq = eSup->gInfo.RTR.metric_length_squared(dx);
							if(rSq < rSqSmallest)
							{	rSqSmallest = rSq;
								iClosest = iAtom;
							}
						}
						//Count as a move if within threshold:
						if(rSqSmallest < std::pow(rCutMove,2))
						{	nAdded--;
							nMoved++;
							rSqMoved += rSqSmallest;
							const vector3<>& xRef = sp.atpos[iClosest];
							x = xRef + wsSup->restrict(x - xRef);
						}
						else x = xCenter + wsSup->restrict(x - xCenter);
						atposRef[iSp].push_back(x);
					}
				}
				//Repeat DFT+U matrices for unit cell repetitions, if needed:
				int nAtomsUnit = spUnit.atpos.size();
				int nMatricesPerAtom = spUnit.rhoAtom_nMatrices() / nAtomsUnit;
				for(int iMatrix=0; iMatrix<nMatricesPerAtom; iMatrix++)
				{	for(int iRep=0; iRep<supIn[0]*supIn[1]*supIn[2]; iRep++)
						U_rhoAtomRef.insert(U_rhoAtomRef.end(), U_rhoAtomUnit, U_rhoAtomUnit+nAtomsUnit);
					U_rhoAtomUnit += nAtomsUnit;
				}
				break;
			}
		unsigned nRemoved = atposRef[iSp].size() - nMoved;
		logPrintf("\t%s: added %u atoms, removed %u atoms, moved %u atoms",
			sp.name.c_str(), nAdded, nRemoved, nMoved);
		if(nMoved) logPrintf(" by %lf bohrs RMS", sqrt(rSqMoved/nMoved));
		logPrintf("\n");
		//Convert atoms to unit cell fratcional coordinates:
		matrix3<> DiagSupIn = Diag(vector3<>(supIn));
		for(vector3<>& x: sp.atpos) x = DiagSupIn * x;
		for(vector3<>& x: atposRef[iSp]) x = DiagSupIn * x;
	}
	for(unsigned iSpUnit=0; iSpUnit<e->iInfo.species.size(); iSpUnit++)
		if(not unitSpeciesFound[iSpUnit])
			logPrintf("WARNING: unit cell species %s not present in defect supercell.\n",
				e->iInfo.species[iSpUnit]->name.c_str());
		
	//Map Vscloc differences:
	vector3<int> Ssup = Diag(supIn) * e->gInfo.S; //expected supercell grid
	if(not (eSup->gInfo.S == Ssup))
		die("Supercell grid is not commensurate with unit cell.\n"
			"Specify fftbox explicitly in defect calculation.\n\n");
	//--- subtract unit cell potential:
	assert(e->eInfo.nDensities == eSup->eInfo.nDensities);
	bool potentialSubtraction = e->dump.potentialSubtraction;
	ScalarField dAtomic, dAtomicSup;
	#define PREPARE_dAtomic(target, e) \
		{	ScalarFieldTilde target##Tilde; \
			initZero(target##Tilde, e->gInfo); \
			for(auto sp: e->iInfo.species) \
				if(sp->atpos.size()) \
					sp->accumulateAtomicPotential(target##Tilde); \
			target = I(target##Tilde); \
		}
	if(potentialSubtraction)
	{	PREPARE_dAtomic(dAtomic, e) //Prepare atomic potential for unit cell
		PREPARE_dAtomic(dAtomicSup, eSup) //Prepare atomic potential for unit cell
	}
	int truncDir = e->coulombParams.geometry==CoulombParams::Slab ? e->coulombParams.iDir : -1;
	double Valign = 0., wAlign = 0.;
	for(int s=0; s<e->eInfo.nDensities; s++)
	{	const vector3<int>& S = e->gInfo.S;
		size_t iStart=0, iStop=e->gInfo.nr;
		vector3<> SsupInv(1./Ssup[0], 1./Ssup[1], 1./Ssup[2]);
		const double* dUnit = potentialSubtraction ? dAtomic->data() : NULL;
		const double* Vunit = e->eVars.Vscloc[s]->data();
		double* dSup = potentialSubtraction ? dAtomicSup->data() : NULL;
		double* Vsup = eSup->eVars.Vscloc[s]->data();
		double* nUnit = e->eVars.n[s]->data();
		THREAD_rLoop(
			vector3<int> ivSup;
			for(ivSup[0]=iv[0]; ivSup[0]<Ssup[0]; ivSup[0]+=S[0])
			for(ivSup[1]=iv[1]; ivSup[1]<Ssup[1]; ivSup[1]+=S[1])
			for(ivSup[2]=iv[2]; ivSup[2]<Ssup[2]; ivSup[2]+=S[2])
			{	size_t iSup = (ivSup[0]*Ssup[1] + ivSup[1])*Ssup[2] + ivSup[2];
				Vsup[iSup] -= Vunit[i];
				//Atomic corrections for alignment potential:
				if(potentialSubtraction and (s == 0)) //only subtract once (dAtomic is spin-independent)
					dSup[iSup] -= dUnit[i];
				double dVatom = potentialSubtraction ? (dSup[iSup] * e->gInfo.dV) : 0.; //Note: does not change Vscloc, only Valign
				//Contribution to alignment potential:
				vector3<> xSup(ivSup[0]*SsupInv[0], ivSup[1]*SsupInv[1], ivSup[2]*SsupInv[2]); //supercell lattice coordinates
				double bDist = wsSup->boundaryDistance(wsSup->restrict(xSup - xCenter), truncDir);
				double w = 0.5*erfc((bDist - alignWidth)/alignSmooth) * nUnit[i]; //weight for potential alignment
				Valign -= w * (Vsup[iSup] - dVatom); //note: alignment correction is unit cell - supercell
				wAlign += w;
			}
		)
	}
	Valign /= wAlign;
	ScalarField dV = potentialSubtraction ? (-dAtomicSup) : 0; //for output / debugging alignment only
	for(ScalarField& V: eSup->eVars.Vscloc)
	{	V += Valign; //apply alignment correction
		//Spin-averaged V difference (without dV factor):
		dV += V * (1./(eSup->eInfo.nDensities * e->gInfo.dV));
	}
	logPrintf("Long-range potential alignment correction: %lg Eh\n", Valign/e->gInfo.dV); //report without dV factor of Vscloc
	//--- write difference potential for manual verification:
	string fname = wannier->getFilename(Wannier::FilenameDump, "dV_"+name);
	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
	saveRawBinary(dV, fname.c_str());
	logPrintf("done; written on supercell grid: "); Ssup.print(globalLog, " %d ");
}

void DefectSupercell::project(const ColumnBundle& C, CachedProjections& proj)
{	for(int i=0; i<2; i++)
		proj.VdagC[i].resize(atposRef.size());
	const matrix* U_rhoAtomPtr[2] = { NULL, NULL };
	if(eSup->eInfo.hasU)
	{	for(int i=0; i<2; i++)
		{	proj.psiDagC[i].resize(atposRef.size());
			Urho[i].resize(atposRef.size());
		}
		U_rhoAtomPtr[0] = U_rhoAtomRef.data();
		U_rhoAtomPtr[1] = eSup->eVars.U_rhoAtom.data();
	}
	for(unsigned iSp=0; iSp<atposRef.size(); iSp++)
	{	SpeciesInfo& sp = *(eSup->iInfo.species[iSp]);
		for(int i=0; i<2; i++)
		{	std::swap(atposRef[iSp], sp.atpos);
			sp.sync_atpos();
			//--- First pass: atpos contains unit cell positions from atposRef
			//--- Second pass: atpos contains defect-supercell positions from eSup
			if(sp.atpos.size())
			{	proj.VdagC[i][iSp] = (*sp.getV(C)) ^ C;
				if(sp.rhoAtom_nMatrices())
				{	ColumnBundle psi;
					sp.rhoAtom_getV(C, U_rhoAtomPtr[i], psi, Urho[i][iSp]);
					proj.psiDagC[i][iSp] = psi ^ C;
					U_rhoAtomPtr[i] += sp.rhoAtom_nMatrices();
				}
			}
		}
	}
}

void DefectSupercell::bcast(CachedProjections& proj, int src) const
{	const int nBands = e->eInfo.nBands;
	for(int i=0; i<2; i++)
	{	proj.VdagC[i].resize(atposRef.size());
		if(eSup->eInfo.hasU)
			proj.psiDagC[i].resize(atposRef.size());
		for(unsigned iSp=0; iSp<atposRef.size(); iSp++)
		{	const SpeciesInfo& sp = *(eSup->iInfo.species[iSp]);
			const int nAtoms = (i ? sp.atpos.size() : atposRef[iSp].size());
			//Nonlocal PS projectors:
			const int nProj = sp.MnlAll.nRows() * nAtoms;
			matrix& m = proj.VdagC[i][iSp];
			if(nProj)
			{	if(not m) m.init(nProj, nBands);
				mpiWorld->bcastData(m, src);
			}
			//DFT+U atomic orbitals:
			if(eSup->eInfo.hasU)
			{	matrix& m = proj.psiDagC[i][iSp];
				int nOrbitals = m.nRows();
				mpiWorld->bcast(nOrbitals, src);
				if(nOrbitals)
				{	if(not m) m.init(nOrbitals, nBands);
					mpiWorld->bcastData(m, src);
				}
			}
		}
	}
}

matrix DefectSupercell::compute(const ColumnBundle& C1, const ColumnBundle& C2,
	const CachedProjections& proj1, const CachedProjections& proj2)
{
	//Local potential contributions:
	//--- collect effective dVscloc(q) on unit cell grid
	vector3<> q = C1.qnum->k - C2.qnum->k;
	vector3<> qSup2pi = (2*M_PI) * (Diag(vector3<>(supIn)) * q); //q in supercell recip. coords (and a factor of 2*pi)
	int nDensities = e->eInfo.nDensities;
	std::vector<complexScalarField> dVsclocq; nullToZero(dVsclocq, e->gInfo, nDensities);
	std::vector<complex*> dVq; for(complexScalarField& Vs: dVsclocq) dVq.push_back(Vs->data());
	std::vector<const double*> dVsup; for(const ScalarField& Vs: eSup->eVars.Vscloc) dVsup.push_back(Vs->data());
	{	const vector3<int>& S = e->gInfo.S;
		const vector3<int>& Ssup = eSup->gInfo.S;
		vector3<> SsupInv(1./Ssup[0], 1./Ssup[1], 1./Ssup[2]);
		size_t iStart=0, iStop=e->gInfo.nr;
		THREAD_rLoop(
			vector3<int> ivSup;
			for(ivSup[0]=iv[0]; ivSup[0]<Ssup[0]; ivSup[0]+=S[0])
			for(ivSup[1]=iv[1]; ivSup[1]<Ssup[1]; ivSup[1]+=S[1])
			for(ivSup[2]=iv[2]; ivSup[2]<Ssup[2]; ivSup[2]+=S[2])
			{	size_t iSup = (ivSup[0]*Ssup[1] + ivSup[1])*Ssup[2] + ivSup[2];
				vector3<> xSup(ivSup[0]*SsupInv[0], ivSup[1]*SsupInv[1], ivSup[2]*SsupInv[2]); //supercell lattice coordinates
				xSup = xCenter + wsSup->restrict(xSup - xCenter); //wrap to WS supercell 
				complex phase = cis(-dot(qSup2pi, xSup)); //Bloch phase
				for(int s=0; s<nDensities; s++)
					dVq[s][i] += phase * dVsup[s][iSup];
			}
		)
	}
	//--- compute corresponding matrix elements:
	matrix result = C1 ^ Idag_DiagV_I(C2, dVsclocq);
	
	//Nonlocal potential contributions:
	double projSign[2] = { -1., +1. }; //sign of the projector contribution (-1 for reference, +1 for defect)
	for(unsigned iSp=0; iSp<atposRef.size(); iSp++)
	{	SpeciesInfo& sp = *(eSup->iInfo.species[iSp]);
		for(int i=0; i<2; i++)
		{	if(not proj1.VdagC[i][iSp]) continue; //may not have any atoms for some eSup species when using unit cell atpos
			const int nAtoms = (i ? sp.atpos.size() : atposRef[iSp].size());
			axpy(projSign[i], dagger(proj1.VdagC[i][iSp]) * (tiledBlockMatrix(sp.MnlAll, nAtoms) * proj2.VdagC[i][iSp]), result);
			//DFT+U:
			if(eSup->eInfo.hasU and Urho[i][iSp])
				axpy(projSign[i], dagger(proj1.psiDagC[i][iSp]) * Urho[i][iSp] * proj2.psiDagC[i][iSp], result);
		}
	}
	
	return result;
}
