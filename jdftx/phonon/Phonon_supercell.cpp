/*-------------------------------------------------------------------
Copyright 2014 Ravishankar Sundararaman

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

#include <phonon/Phonon.h>
#include <commands/parser.h>
#include <commands/command.h>
#include <electronic/ColumnBundleTransform.h>

inline bool spinEqual(const QuantumNumber& qnum1, const QuantumNumber& qnum2) { return qnum1.spin == qnum2.spin; } //for k-point mapping (in spin polarized mode)

void Phonon::processPerturbation(const Perturbation& pert, string fnamePattern)
{
	//Start with eSupTemplate:
	eSup = std::make_shared<PhononEverything>(*this);
	eSup->cntrl.dragWavefunctions = false; //wavefunction-drag doesn't always play nice with setSupState (especially with relativity)
	logSuspend(); parse(input, *eSup); logResume();
	eSup->eInfo.kfold = eSupTemplate.eInfo.kfold;
	eSup->gInfo.S = eSupTemplate.gInfo.S;
	eSup->gInfo.R = eSupTemplate.gInfo.R;
	int nAtomsTot = 0;
	for(size_t sp=0; sp<e.iInfo.species.size(); sp++)
	{	const SpeciesInfo& spIn = *(eSupTemplate.iInfo.species[sp]);
		SpeciesInfo& spOut = *(eSup->iInfo.species[sp]);
		spOut.atpos = spIn.atpos;
		spOut.constraints = spIn.constraints;
		spOut.initialMagneticMoments.clear();
		for(SpeciesInfo::PlusU& Uparams: spOut.plusU)
			Uparams.Vext.resize(spOut.atpos.size());
		nAtomsTot += spOut.atpos.size();
	}
	if(eSup->symm.sym.size()) eSup->symm.sym = eSupTemplate.symm.sym; //use supercell version of manual symmetry matrices (if any)
	eSup->symm.isPertSup = true; //perturbed supercell; atom positions will reduce instead of check manual symmetries (if any)
	
	//Remove unit cell initial state settings (incompatible with supercell):
	eSup->eVars.skipWfnsInit = true; //skip because wavefunctions are set from unit cell calculation
	eSup->eVars.wfnsFilename.clear();
	eSup->eVars.eigsFilename.clear();
	eSup->eVars.fluidInitialStateFilename.clear();
	eSup->eInfo.initialFillingsFilename.clear();
	eSup->scfParams.historyFilename.clear();
	//ElecInfo:
	eSup->eInfo.nBands = nBandsOpt * prodSup;
	//ElecVars:
	eSup->eVars.initLCAO = false; //state will be initialized from unit cell anyway
	eSup->cntrl.convergeEmptyStates = false; //has no effect on any phonon results, so force-disable to save time
	
	//Instead read in appropriate supercell quantities:
	if(collectPerturbations)
	{	//Read perturbed Vscloc and fix H below (no minimize/SCF necessary):
		if(saveHsub) //only need Vscloc if outputting e-ph matrix elements
		{	eSup->eVars.VFilenamePattern = fnamePattern;
			eSup->cntrl.fixed_H = true;
		}
	}
	else
	{	//Read in supercell state if available:
		setAvailableFilenames(fnamePattern, *eSup);
		if(eSup->eVars.wfnsFilename.length())
			eSup->eVars.skipWfnsInit = false; //do read in wfns if available
	}
	eSup->dump.format = fnamePattern;
	
	//Apply perturbation and then setup (so that symmetries reflect perturbed state):
	vector3<> dxPert = inv(eSup->gInfo.R) * dr * pert.dir; //perturbation in lattice coordinates
	eSup->iInfo.species[pert.sp]->atpos[pert.at] += dxPert; //apply perturbation
	eSup->setup();
	if(dryRun)
	{	logPrintf("Dry run: supercell setup successful.\n");
		return;
	}
	
	//Map states:
	PeriodicLookup<QuantumNumber> plook(eSup->eInfo.qnums, eSup->gInfo.GGT);
	std::vector<int> nqPrev(eSup->eInfo.qnums.size(), 0);
	stateMap.clear();
	const Supercell& supercell = *(e.coulombParams.supercell);
	for(int iSpin=0; iSpin<nSpins; iSpin++)
		for(unsigned ik=0; ik<supercell.kmesh.size(); ik++)
		{	QuantumNumber qnum;
			qnum.k = supercell.kmesh[ik];
			qnum.spin = (nSpins==1 ? 0 : (iSpin==0 ? +1 : -1));
			vector3<> kSup = matrix3<>(Diag(sup)) * qnum.k; //qnum.k in supercell reciprocal lattice coords
			size_t qSup = plook.find(kSup, qnum, &(eSup->eInfo.qnums), spinEqual);
			if(qSup == string::npos) continue; //the corresponding supercell k-point must have been eliminated by symmetries
			//Add to stateMap:
			StateMapEntry sme;
			(Supercell::KmeshTransform&)sme = supercell.kmeshTransform[ik]; //copy base class properties
			sme.iReduced += iSpin*(e.eInfo.nStates/nSpins); //point to source k-point with appropriate spin
			sme.qSup = qSup;
			sme.nqPrev = nqPrev[qSup];
			sme.k = supercell.kmesh[ik];
			nqPrev[qSup]++;
			stateMap.push_back(sme);
		}
	for(int nq: nqPrev) assert(nq == prodSup); //each supercell k-point must be mapped to prod(sup) unit cell kpoints
	
	//Map corresponding basis objects:
	std::vector<SpaceGroupOp> sym = e.symm.getMatrices();
	for(int qSup=0; qSup<eSup->eInfo.nStates; qSup++)
		if(eSup->eInfo.isMine(qSup) || eSup->eInfo.qnums[qSup].k.length_squared()==0) //Make Gamma-point available on all processes
		{	ColumnBundleTransform::BasisWrapper basisSupWrapper(eSup->basis[qSup]);
			const vector3<>& kSup = eSup->eInfo.qnums[qSup].k;
			//Initialize index map for each unit cell k-point
			for(StateMapEntry& sme: stateMap) if(sme.qSup == qSup)
			{	const Basis& basis = e.basis[sme.iReduced];
				const vector3<>& k = e.eInfo.qnums[sme.iReduced].k;
				sme.transform = std::make_shared<ColumnBundleTransform>(k, basis, kSup, basisSupWrapper, nSpinor, sym[sme.iSym], sme.invert, Diag(sup));
			}
		}
	
	//Initialize state of supercell:
	std::vector<diagMatrix> Hsub0;
	if(!collectPerturbations || saveHsub)
		Hsub0 = setSupState();
	
	//Calculate energy and forces:
	IonicGradient dgrad_pert;
	if(collectPerturbations)
	{	dgrad_pert.init(eSup->iInfo);
		dgrad_pert.read(eSup->dump.getFilename("dforces").c_str());
		if(saveHsub) eSup->iInfo.augmentDensityGridGrad(eSup->eVars.Vscloc); //update Vscloc atom projections for ultrasoft psp's (needed by getPerturbedHsub)
	}
	else
	{	IonicGradient grad;
		IonicMinimizer(*eSup).compute(&grad, 0);
		dgrad_pert = (grad - grad0) * (1./dr);
		logPrintf("Energy change: %lg / unit cell\n", (relevantFreeEnergy(*eSup) - E0)/prodSup);
		logPrintf("RMS force: %lg\n", sqrt(dot(grad,grad)/(3*nAtomsTot)));
		//Output potentials and forces for future calculations:
		eSup->dump.insert(std::make_pair(DumpFreq_End, DumpVscloc));
		eSup->dump(DumpFreq_End, 0);
		dgrad_pert.write(eSup->dump.getFilename("dforces").c_str());
	}
	if(iPerturbation>=0) return; //remainder below not necessary except when doing a full calculation
	
	//Apply translational invariance correction:
	vector3<> fMean;
	for(const std::vector<vector3<>>& fArr: dgrad_pert)
		for(const vector3<>& f: fArr)
			fMean += f;
	fMean *= (1./nAtomsTot);
	for(std::vector<vector3<>>& fArr: dgrad_pert)
		for(vector3<>& f: fArr)
			f -= fMean;
	logPrintf("Applied translational invariance (net force) relative correction: %lg\n",
		sqrt(nAtomsTot*fMean.length_squared()/dot(dgrad_pert,dgrad_pert)));
	
	//Subspace hamiltonian change:
	std::vector<matrix> Hsub, dHsub_pert(nSpins);
	if(saveHsub)
	{	Hsub = getPerturbedHsub(pert, Hsub0);
		for(size_t s=0; s<Hsub.size(); s++)
			dHsub_pert[s] = (1./dr) * (Hsub[s] - Hsub0[s]);
	}
	
	//Accumulate results for all symmetric images of perturbation:
	const auto& atomMap = eSupTemplate.symm.getAtomMap();
	for(unsigned iSym=0; iSym<symSupCart.size(); iSym++)
	{
		//Figure out the mode that the rotated perturbation corresponds to:
		Mode mode;
		mode.sp = pert.sp; //rotations are not alchemists!
		mode.at = atomMap[pert.sp][pert.at][iSym];
		mode.dir = symSupCart[iSym] * pert.dir;
		
		//Reduce mode atom to fundamental unit cell if necessary:
		int nAtoms = e.iInfo.species[mode.sp]->atpos.size(); //per unit cell
		int unit = mode.at / nAtoms; //unit cell index of mapped atom
		mode.at -= nAtoms*unit; //mode.at is now in [0,nAtoms)
		vector3<int> cellOffset = -getCell(unit); //corresponding displacement in unit cell lattice coords
		
		//Find index of first mode that matches sp and at:
		unsigned iModeStart = 0;
		for(iModeStart=0; iModeStart<modes.size(); iModeStart++)
			if(mode.sp==modes[iModeStart].sp && mode.at==modes[iModeStart].at)
				break;
		assert(iModeStart+3 <= modes.size());
		
		//Accumulate dgrad contributions:
		for(unsigned sp2=0; sp2<eSup->iInfo.species.size(); sp2++)
		{	int nAtoms2 = e.iInfo.species[sp2]->atpos.size(); //per unit cell
			for(int at2=0; at2<nAtoms2*prodSup; at2++)
			{	int at2rot = atomMap[sp2][at2][iSym];
				int unit2rot = at2rot / nAtoms2;
				at2rot += nAtoms2 * (getUnit(getCell(unit2rot) + cellOffset) - unit2rot); //apply cellOffset
				vector3<> Frot = symSupCart[iSym] * dgrad_pert[sp2][at2]; //rotated force
				for(unsigned iMode2=iModeStart; iMode2<iModeStart+3; iMode2++)
					dgrad[iMode2][sp2][at2rot] += (pert.weight * dot(modes[iMode2].dir, mode.dir)) * Frot;
			}
		}
		
		//Accumulate Hsub contributions:
		if(saveHsub)
			for(int iSpin=0; iSpin<nSpins; iSpin++)
			{	//Fetch Hsub with rotations:
				matrix contrib = stateRot[iSpin][iSym].transform(dHsub_pert[iSpin]);
				//Collect k-vectors in order of the blocks of Hsub:
				int qSup = iSpin*(eSup->eInfo.nStates/nSpins); //Gamma point is always first in the list for each spin
				assert(eSup->eInfo.qnums[qSup].k.length_squared() == 0); //make sure that above is true
				std::vector< vector3<> > k; k.reserve(prodSup);
				for(const StateMapEntry& sme: stateMap) if(sme.qSup==qSup)
					k.push_back(sme.k);
				assert(int(k.size()) == prodSup);
				//Apply phase factors due to translations:
				int nBands = e.eInfo.nBands;
				for(unsigned ik1=0; ik1<k.size(); ik1++)
					for(unsigned ik2=0; ik2<k.size(); ik2++)
						contrib.set(ik1*nBands,(ik1+1)*nBands, ik2*nBands,(ik2+1)*nBands,
							contrib(ik1*nBands,(ik1+1)*nBands, ik2*nBands,(ik2+1)*nBands)
								* cis(-2*M_PI*dot(cellOffset, k[ik1]-k[ik2])));
				for(unsigned iMode2=iModeStart; iMode2<iModeStart+3; iMode2++)
					dHsub[iMode2][iSpin] += contrib * (pert.weight * dot(modes[iMode2].dir, mode.dir));
			}
	}
}

#define INITwfnsSup(C, nCols) \
	C.init(nCols, eSup->basis[qSup].nbasis * eSup->eInfo.spinorLength(), \
		&eSup->basis[qSup], &eSup->eInfo.qnums[qSup], isGpuEnabled());

std::vector<diagMatrix> Phonon::setSupState()
{	static StopWatch watch("phonon::setSupState"); watch.start();
	double scaleFac = 1./sqrt(prodSup); //to account for normalization
	
	//Calculate unperturbed subspace Hamiltonian (supercell Gamma point only, but all bands):
	std::vector<diagMatrix> Hsub0(nSpins);
	int nBands = e.eInfo.nBands;
	int nBandsSup = nBands * prodSup;
	for(int s=0; s<nSpins; s++)
	{	int qSup = s*(eSup->eInfo.nStates/nSpins); //Gamma point is always first in the list for each spin
		assert(eSup->eInfo.qnums[qSup].k.length_squared() == 0); //make sure that above is true
		Hsub0[s].resize(nBandsSup);
		if(eSup->eInfo.isMine(qSup))
		{	for(const StateMapEntry& sme: stateMap) if(sme.qSup==qSup)
			{	int nBandsPrev = nBands * sme.nqPrev;
				Hsub0[s].set(nBandsPrev,nBandsPrev+nBands, e.eVars.Hsub_eigs[sme.iReduced]);
			}
		}
		Hsub0[s].bcast(eSup->eInfo.whose(qSup));
	}
	if(collectPerturbations || eSup->eVars.wfnsFilename.length())
	{	//skip state initialization below if already read in, or if not needed since in collect mode:
		watch.stop();
		return Hsub0;
	}
	
	//Update supercell quantities:
	int nBandsOptSup = nBandsOpt * prodSup;
	for(int qSup=eSup->eInfo.qStart; qSup<eSup->eInfo.qStop; qSup++)
	{	//Prepare outputs:
		//--- Wavefuntions:
		ColumnBundle& Csup = eSup->eVars.C[qSup];
		if(Csup.nCols() != nBandsOptSup) INITwfnsSup(Csup, nBandsOptSup);
		Csup.zero(); //zero wavefunctions (since a scatter used below)
		//--- Fillings:
		diagMatrix& Fsup = eSup->eVars.F[qSup];
		Fsup.resize(nBandsOptSup);
		//--- Auxiliary Hamiltonian (if necessary):
		diagMatrix unused;
		diagMatrix& Haux_eigsSup = (e.eInfo.fillingsUpdate==ElecInfo::FillingsHsub) ? eSup->eVars.Haux_eigs[qSup] : unused;
		Haux_eigsSup.resize(nBandsOptSup);
		//Collect from unit cell to supercell:
		for(const StateMapEntry& sme: stateMap) if(sme.qSup==qSup)
		{	int nBandsPrev = nBandsOpt * sme.nqPrev;
			//Wavefunctions:
			const ColumnBundle& C = e.eVars.C[sme.iReduced];
			sme.transform->scatterAxpy(scaleFac, C.getSub(0,nBandsOpt), Csup,nBandsPrev,1);
			//Fillings:
			const diagMatrix& F = e.eVars.F[sme.iReduced];
			Fsup.set(nBandsPrev,nBandsPrev+nBandsOpt, F(0,nBandsOpt));
			//Auxiliary Hamiltonian (if necessary):
			if(e.eInfo.fillingsUpdate==ElecInfo::FillingsHsub)
			{	const diagMatrix& Haux_eigs = e.eVars.Haux_eigs[sme.iReduced];
				Haux_eigsSup.set(nBandsPrev,nBandsPrev+nBandsOpt, Haux_eigs(0,nBandsOpt));
			}
		}
		eSup->iInfo.project(Csup, eSup->eVars.VdagC[qSup]); //update wave function projections
	}
	
	//Update entropy contributions:
	if(eSup->eInfo.fillingsUpdate == ElecInfo::FillingsHsub)
	{	eSup->eInfo.updateFillingsEnergies(eSup->eVars.Haux_eigs, eSup->ener);
		eSup->eVars.HauxInitialized = true;
		eSup->eVars.Hsub_eigs = eSup->eVars.Haux_eigs;
	}
	
	watch.stop();
	return Hsub0;
}

std::vector<matrix> Phonon::getPerturbedHsub(const Perturbation& pert, const std::vector<diagMatrix>& Hsub0)
{	static StopWatch watch("phonon::getPerturbedHsub"); watch.start();
	double scaleFac = 1./sqrt(prodSup); //to account for normalization
	int nBands = e.eInfo.nBands;
	int nSpinor = e.eInfo.spinorLength();
	int nBandsSup = nBands * prodSup; //Note >= eSup->eInfo.nBands, depending on e.eInfo.nBands >= nBandsOpt
	int nqPrevStart, nqPrevStop; TaskDivision(prodSup, mpiWorld).myRange(nqPrevStart, nqPrevStop);
	//Get unperturbed projectors to account for ultrasoft augmentation in overlap (if any):
	SpeciesInfo& spPert = *(eSup->iInfo.species[pert.sp]); //species that has been perturbed
	std::vector<ColumnBundle> V0(nSpins); //unperturbed projectors of the perturbed atom
	int pStart = 0, pStop = 0; //range of projectors for perturbed atom
	std::vector<matrix> V0dagC(nSpins), VdagC(nSpins); //unperturbed and perturbed projections on perturbed atom
	if(spPert.QintAll)
	{	//Undo perturbation (needed to get unperturbed projectors below):
		vector3<> atpos0 = eSupTemplate.iInfo.species[pert.sp]->atpos[pert.at]; //unperturbed atom position
		std::swap(spPert.atpos[pert.at], atpos0);
		spPert.sync_atpos(); //Note: also clears cached projectors
		for(int s=0; s<nSpins; s++)
		{	int qSup = s*(eSup->eInfo.nStates/nSpins); //Gamma point is always first in the list for each spin
			pStart = spPert.QintAll.nRows() * pert.at; //perturbed atom projector starts here ...
			pStop = spPert.QintAll.nRows() * (pert.at+1); //... and ends here
			ColumnBundle Csup; INITwfnsSup(Csup, 1) //Dummy columnbundle for getV below
			V0[s] = spPert.getV(Csup)->getSub(pStart/nSpinor,pStop/nSpinor);
			V0dagC[s] = zeroes(pStop-pStart, nBandsSup);
			VdagC[s] = zeroes(pStop-pStart, nBandsSup);
		}
		//Restore perturbation:
		std::swap(spPert.atpos[pert.at], atpos0);
		spPert.sync_atpos();
	}
	std::vector<matrix> Hsub(nSpins);
	for(int s=0; s<nSpins; s++)
	{	int qSup = s*(eSup->eInfo.nStates/nSpins); //Gamma point is always first in the list for each spin
		assert(eSup->eInfo.qnums[qSup].k.length_squared() == 0); //make sure that above is true
		//Initialize outputs:
		Hsub[s] = zeroes(nBandsSup, nBandsSup);
		//Loop over supercell-commensurate unit cell k-points:
		for(const StateMapEntry& sme: stateMap)
			if(sme.qSup==qSup && sme.nqPrev>=nqPrevStart && sme.nqPrev<nqPrevStop) //MPI divide on nqPrev
			{	//Set supercell wavefunctions:
				const ColumnBundle& C = e.eVars.C[sme.iReduced];
				ColumnBundle& Csup = eSup->eVars.C[sme.qSup];
				if(Csup.nCols() != nBands) INITwfnsSup(Csup, nBands)
				Csup.zero();
				sme.transform->scatterAxpy(scaleFac, C, Csup,0,1);
				//Apply Hamiltonian:
				ColumnBundle HCsup; Energies ener;
				eSup->iInfo.project(Csup, eSup->eVars.VdagC[qSup]); //update wavefunction projections
				eSup->eVars.applyHamiltonian(qSup, eye(nBands), HCsup, ener, true);
				int start = nBands * sme.nqPrev;
				int stop = nBands * (sme.nqPrev+1);
				//Store projections for augmentation correction (if needed)
				if(spPert.QintAll)
				{	V0dagC[s].set(0,pStop-pStart, start,stop, V0[s] ^ Csup);
					VdagC[s].set(0,pStop-pStart, start,stop, eSup->eVars.VdagC[qSup][pert.sp](pStart,pStop, 0,Csup.nCols()));
				}
				//Second loop over supercell-commensurate unit cell k-points:
				for(const StateMapEntry& sme2: stateMap) if(sme2.qSup==qSup)
				{	const ColumnBundle& C2 = e.eVars.C[sme2.iReduced];
					ColumnBundle HC = C2.similar();
					HC.zero();
					sme2.transform->gatherAxpy(scaleFac, HCsup,0,1, HC);
					//Compute overlaps:
					int start2 = nBands * sme2.nqPrev;
					int stop2 = nBands * (sme2.nqPrev+1);
					matrix Hsub12 = HC^C2; //subspace Hamiltonian in the complex-conjugate convention of reduced C2
					if(sme2.invert<0) Hsub12 = conj(Hsub12); //switch conjugate convention back to that of supercell
					Hsub[s].set(start,stop, start2,stop2, Hsub12);
				}
			}
		Hsub[s].allReduce(MPIUtil::ReduceSum);
	}
	//Account for ultrasoft overlap augmentation (if any):
	if(spPert.QintAll)
	{	for(int s=0; s<nSpins; s++)
		{	VdagC[s].allReduce(MPIUtil::ReduceSum);
			V0dagC[s].allReduce(MPIUtil::ReduceSum);
			matrix dVdagC = VdagC[s] - V0dagC[s]; //change in projection of unperturbed states due to perturbation
			matrix contrib = dagger(dVdagC) * spPert.QintAll * VdagC[s] * Hsub0[s];
			Hsub[s] -= (contrib + dagger(contrib));
		}
	}
	watch.stop();
	return Hsub; 
}
