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

inline vector3<> getCoord(const QuantumNumber& qnum) { return qnum.k; } //for k-point mapping
inline bool spinEqual(const QuantumNumber& qnum1, const QuantumNumber& qnum2) { return qnum1.spin == qnum2.spin; } //for k-point mapping (in spin polarized mode)

void Phonon::processPerturbation(const Perturbation& pert)
{
	//Start with eSupTemplate:
	eSup = std::make_shared<Everything>();
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

	//Remove initial state settings (incompatible with supercell):
	eSup->eVars.wfnsFilename.clear();
	eSup->eVars.HauxFilename.clear();
	eSup->eVars.eigsFilename.clear();
	eSup->eVars.fluidInitialStateFilename.clear();
	eSup->eInfo.initialFillingsFilename.clear();
	eSup->scfParams.historyFilename.clear();
	//ElecInfo:
	eSup->eInfo.nBands = nBandsOpt * prodSup;
	//ElecVars:
	eSup->eVars.initLCAO = false; //state will be initialized from unit cell anyway
	
	//Apply perturbation and then setup (so that symmetries reflect perturbed state):
	double drSym = 0.;
	for(int iDir=0; iDir<3; iDir++)
		drSym = std::max(drSym, eSup->gInfo.R.column(iDir).length());
	drSym *= (10*symmThreshold); //ensure temporary perturbation is an order of magnitude larger than symmetry detection threshold
	vector3<> dxSym = inv(eSup->gInfo.R) * drSym * pert.dir;
	eSup->iInfo.species[pert.sp]->atpos[pert.at] += dxSym; //apply perturbation for detection of symmetry reduction
	eSup->setup();
	eSup->iInfo.species[pert.sp]->atpos[pert.at] -= dxSym; //restore unperturbed geometry
	eSup->iInfo.species[pert.sp]->sync_atpos();
	
	//Map states:
	PeriodicLookup<QuantumNumber> plook(eSup->eInfo.qnums, eSup->gInfo.GGT);
	std::vector<int> nqPrev(eSup->eInfo.qnums.size(), 0);
	stateMap.clear();
	const Supercell& supercell = *(e.coulombParams.supercell);
	for(int iSpin=0; iSpin<nSpins; iSpin++)
		for(unsigned ik=0; ik<supercell.kmesh.size(); ik++)
		{	QuantumNumber qnum;
			qnum.k = supercell.kmesh[ik];
			qnum.spin = (nSpins==1 ? 0 : (iSpin ? +1 : -1));
			vector3<> kSup = matrix3<>(Diag(sup)) * qnum.k; //qnum.k in supercell reciprocal lattice coords
			size_t qSup = plook.find(kSup, qnum, &(eSup->eInfo.qnums), spinEqual);
			if(qSup == string::npos) continue; //the corresponding supercell k-point must have been eliminated by symmetries
			vector3<> iGtmp = kSup - eSup->eInfo.qnums[qSup].k;
			//Add to stateMap:
			auto sme = std::make_shared<StateMapEntry>();
			(Supercell::KmeshTransform&)(*sme) = supercell.kmeshTransform[ik]; //copy base class properties
			sme->iReduced += iSpin*(e.eInfo.nStates/nSpins); //point to source k-point with appropriate spin
			sme->qSup = qSup;
			for(int j=0; j<3; j++)
				sme->iG[j] = round(iGtmp[j]);
			sme->nqPrev = nqPrev[qSup];
			nqPrev[qSup]++;
			stateMap.push_back(sme);
		}
	for(int nq: nqPrev) assert(nq == prodSup); //each supercell k-point must be mapped to prod(sup) unit cell kpoints
	
	//Map corresponding basis objects:
	std::vector< matrix3<int> > sym = e.symm.getMatrices();
	for(int qSup=eSup->eInfo.qStart; qSup<eSup->eInfo.qStop; qSup++)
	{	//Create lookup table for supercell indices:
		vector3<int> iGbox; //see Basis::setup(), this will bound the values of iG in the basis
		for(int i=0; i<3; i++)
			iGbox[i] = 1 + int(sqrt(2*eSup->cntrl.Ecut) * eSup->gInfo.R.column(i).length() / (2*M_PI));
		vector3<int> pitch;
		pitch[2] = 1;
		pitch[1] = pitch[2] * (2*iGbox[2]+1);
		pitch[0] = pitch[1] * (2*iGbox[1]+1);
		std::vector<int> supIndex(pitch[0] * (2*iGbox[0]+1), -1);
		const Basis& basisSup = eSup->basis[qSup];
		for(size_t n=0; n<basisSup.nbasis; n++)
			supIndex[dot(pitch,basisSup.iGarr[n]+iGbox)] = n;
		//Initialize index map for each unit cell k-point
		for(std::shared_ptr<StateMapEntry>& sme: stateMap) if(sme->qSup == qSup)
		{	const Basis& basis = e.basis[sme->iReduced];
			std::vector<int> indexVec(basis.nbasis);
			const matrix3<int> mRot = (~sym[sme->iSym]) * sme->invert; //effective transofrmation matrix
			for(size_t n=0; n<basis.nbasis; n++)
			{	vector3<int> iGtmp = basis.iGarr[n]; //original cell recip lattice coords
				iGtmp = mRot * iGtmp - sme->offset; //apply symmetry operation in unit cell 
				for(int j=0; j<3; j++) iGtmp[j] *= sup[j]; //to supercell reciprocal lattice coords
				iGtmp += sme->iG; //offset due to Bloch phase
				indexVec[n] = supIndex[dot(pitch,iGtmp+iGbox)]; //lookup from above table
			}
			assert(*std::min_element(indexVec.begin(), indexVec.end()) >= 0); //make sure all entries were found
			sme->setIndex(indexVec);
		}
	}
	
	//Initialize state of supercell:
	setSupState();
	IonicMinimizer imin(*eSup);
	IonicGradient grad0;
	imin.compute(&grad0); //compute initial forces and energy
	logPrintf("# Energy components:\n"); eSup->ener.print(); logPrintf("\n");
	double E0 = relevantFreeEnergy(*eSup);
	logPrintf("Supercell energy discrepancy: %lg / unit cell\n", E0/prodSup - relevantFreeEnergy(e));
	logPrintf("RMS force in initial configuration: %lg\n", sqrt(dot(grad0,grad0)/(3*nAtomsTot)));
	//subspace Hamiltonian of supercell Gamma point:
	std::vector<matrix> Hsub0;
	setSupState(&Hsub0);

	//Move to perturbed configuration:
	IonicGradient dir; dir.init(eSup->iInfo); //contains zeroes
	dir[pert.sp][pert.at] = pert.dir;
	imin.step(dir, dr);
	
	//Calculate energy and forces:
	IonicGradient grad, dgrad_pert;
	imin.compute(&grad);
	dgrad_pert = (grad - grad0) * (1./dr);
	logPrintf("Energy change: %lg / unit cell\n", (relevantFreeEnergy(*eSup) - E0)/prodSup);
	logPrintf("RMS force: %lg\n", sqrt(dot(grad,grad)/(3*nAtomsTot)));
	
	//Subspace hamiltonian change:
	std::vector<matrix> Hsub, dHsub_pert(nSpins);
	setSupState(&Hsub);
	for(size_t s=0; s<Hsub.size(); s++)
		dHsub_pert[s] = (1./dr) * (Hsub[s] - Hsub0[s]);
	
	//Restore atom position:
	bool dragWfns = false;
	std::swap(dragWfns, eSup->cntrl.dragWavefunctions); //disable wave function drag because state has already been restored to unperturbed version
	imin.step(dir, -dr);
	std::swap(dragWfns, eSup->cntrl.dragWavefunctions); //restore wave function drag flag
	
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
		
		//Accumulate dgrad constributions:
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
	}
}

void Phonon::setSupState(std::vector<matrix>* Hsub)
{
	int nBandsSup = e.eInfo.nBands * prodSup; //Note >= eSup->eInfo.nBands, depending on e.eInfo.nBands >= nBandsOpt
	
	//Zero wavefunctions and auxiliary Hamiltonia (since a scatter used below)
	for(int qSup=eSup->eInfo.qStart; qSup<eSup->eInfo.qStop; qSup++)
	{	ColumnBundle& Cq = eSup->eVars.C[qSup];
		eSup->eVars.Y[qSup].free(); //to save memory (will be regenerated below, after Hsub calculation)
		if(Cq.nCols() != nBandsSup)
			Cq.init(nBandsSup, Cq.colLength(), Cq.basis, Cq.qnum, isGpuEnabled());
		Cq.zero();
		if(e.eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux)
			eSup->eVars.B[qSup] = zeroes(nBandsSup,nBandsSup);
	}
	
	//Update supercell quantities:
	double scaleFac = 1./sqrt(prodSup); //to account for normalization
	for(const std::shared_ptr<StateMapEntry>& sme: stateMap) if(eSup->eInfo.isMine(sme->qSup))
	{	int nBandsPrev = e.eInfo.nBands * sme->nqPrev;
		//Wavefunctions:
		const ColumnBundle& C = e.eVars.C[sme->iReduced];
		ColumnBundle& Csup = eSup->eVars.C[sme->qSup];
		for(int b=0; b<e.eInfo.nBands; b++)
			for(int s=0; s<nSpinor; s++)
				callPref(eblas_scatter_zdaxpy)(sme->nIndices, scaleFac, sme->indexPref,
					C.dataPref() + C.index(b,s*C.basis->nbasis),
					Csup.dataPref() + Csup.index(b + nBandsPrev, s*Csup.basis->nbasis));
		//Fillings:
		const diagMatrix& F = e.eVars.F[sme->iReduced];
		diagMatrix& Fsup = eSup->eVars.F[sme->qSup];
		Fsup.resize(nBandsSup);
		Fsup.set(nBandsPrev,nBandsPrev+e.eInfo.nBands, F);
		//Auxiliary Hamiltonian (if necessary):
		if(e.eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux)
		{	const matrix& B = e.eVars.B[sme->iReduced];
			matrix& Bsup = eSup->eVars.B[sme->qSup];
			Bsup.set(nBandsPrev,nBandsPrev+e.eInfo.nBands, nBandsPrev,nBandsPrev+e.eInfo.nBands, B);
		}
	}
	
	//Compute gamma point Hamiltonian if requested:
	if(Hsub)
	{	Hsub->resize(nSpins);
		for(int s=0; s<nSpins; s++)
		{	int qSup = s*(eSup->eInfo.nStates/nSpins); //Gamma point is always first in the list for each spin
			assert(eSup->eInfo.qnums[qSup].k.length_squared() == 0); //make sure that above is true
			if(eSup->eInfo.isMine(qSup))
			{	ColumnBundle HC; Energies ener;
				eSup->iInfo.project(eSup->eVars.C[qSup], eSup->eVars.VdagC[qSup]); //update wavefunction projections
				eSup->eVars.applyHamiltonian(qSup, eye(nBandsSup), HC, ener, true);
				(*Hsub)[s] = eSup->eVars.Hsub[qSup] * prodSup; //account for scaling of wavefunctions above
			}
			else (*Hsub)[s].init(nBandsSup, nBandsSup);
			(*Hsub)[s].bcast(eSup->eInfo.whose(qSup));
		}
	}
	
	//Discard extra bands:
	for(int qSup=eSup->eInfo.qStart; qSup<eSup->eInfo.qStop; qSup++)
	{	ColumnBundle& Cq = eSup->eVars.C[qSup];
		ColumnBundle& Yq = eSup->eVars.Y[qSup];
		diagMatrix& Fq = eSup->eVars.F[qSup];
		matrix& Bq = eSup->eVars.B[qSup];
		Yq = Cq.similar(eSup->eInfo.nBands);
		diagMatrix Ftmp(eSup->eInfo.nBands);
		matrix Btmp = zeroes(eSup->eInfo.nBands, eSup->eInfo.nBands);
		for(int nqPrev=0; nqPrev<prodSup; nqPrev++)
		{	int offsIn = nqPrev * e.eInfo.nBands;
			int offsOut = nqPrev * nBandsOpt;
			callPref(eblas_copy)(Yq.data()+Yq.index(offsOut,0), Cq.data()+Cq.index(offsIn,0), nBandsOpt*Yq.colLength());
			Ftmp.set(offsOut,offsOut+nBandsOpt, Fq(offsIn,offsIn+nBandsOpt));
			if(e.eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux)
				Btmp.set(offsOut,offsOut+nBandsOpt, offsOut,offsOut+nBandsOpt, Bq(offsIn,offsIn+nBandsOpt, offsIn,offsIn+nBandsOpt));
		}
		Cq = Yq;
		std::swap(Fq, Ftmp);
		if(e.eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux)
			std::swap(Bq, Btmp);
		eSup->iInfo.project(Cq, eSup->eVars.VdagC[qSup]); //update wave function projections
	}
	
	//Update entropy contributions:
	if(eSup->eInfo.fillingsUpdate != ElecInfo::ConstantFillings)
		eSup->eInfo.updateFillingsEnergies(eSup->eVars.F, eSup->ener);
	
	if(e.eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux)
		eSup->eVars.HauxInitialized = true;
}

//------------- class Phonon::StateMapEntry --------------

Phonon::StateMapEntry::StateMapEntry() : nIndices(0), indexPref(0), index(0)
#ifdef GPU_ENABLED
, indexGpu(0)
#endif
{

}

Phonon::StateMapEntry::~StateMapEntry()
{
	if(index) delete[] index;
	#ifdef GPU_ENABLED
	if(indexGpu) cudaFree(indexGpu);
	#endif
}

void Phonon::StateMapEntry::setIndex(const std::vector<int>& indexVec)
{	nIndices = indexVec.size();
	//CPU version
	if(index) delete[] index;
	index = new int[nIndices];
	eblas_copy(index, indexVec.data(), nIndices);
	indexPref = index;
	//GPU version
	#ifdef GPU_ENABLED
	if(indexGpu) cudaFree(indexGpu);
	cudaMalloc(&indexGpu, sizeof(int)*nIndices); gpuErrorCheck();
	cudaMemcpy(indexGpu, index, sizeof(int)*nIndices, cudaMemcpyHostToDevice); gpuErrorCheck();
	indexPref = indexGpu;
	#endif
}
