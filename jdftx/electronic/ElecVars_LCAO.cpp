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

#include <electronic/Everything.h>
#include <electronic/ElecVars.h>
#include <electronic/ColumnBundle.h>
#include <electronic/ElecMinimizer.h>
#include <core/matrix.h>
#include <core/ScalarFieldIO.h>

void printSymmetryError(const matrix& m, const char* name)
{	double num = nrm2(m - dagger(m));
	double den = nrm2(m);
	logPrintf("%s symmetry error = %le (rel), %le (abs)\n", name, num/den, num);
}

struct LCAOminimizer : Minimizable<ElecGradient> //Uses only the Haux entries of ElecGradient
{	int nBands;
	ElecVars& eVars;
	const Everything& e;
	const ElecInfo& eInfo;
	const ExCorr* exCorr;
	std::vector<matrix> HniSub;
	std::vector<matrix> rotPrev; //Accumulated rotations of the wavefunctions
	
	LCAOminimizer(ElecVars& eVars, const Everything& e)
	: eVars(eVars), e(e), eInfo(e.eInfo), HniSub(eInfo.nStates), rotPrev(eInfo.nStates)
	{
	}
	
	void step(const ElecGradient& dir, double alpha)
	{	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{	assert(dir.Haux[q]);
			//Move aux along dir after transforming dir to match rotations:
			matrix Haux = eVars.Haux_eigs[q], Haux_evecs;
			axpy(alpha, dagger(rotPrev[q])*dir.Haux[q]*rotPrev[q], Haux);
			//Adjust rotations to make Haux diagonal again:
			Haux.diagonalize(Haux_evecs, eVars.Haux_eigs[q]);
			rotPrev[q] = rotPrev[q] * Haux_evecs;
			eVars.C[q] = eVars.C[q] * Haux_evecs;
			for(unsigned sp=0; sp<e.iInfo.species.size(); sp++)
				if(eVars.VdagC[q][sp]) eVars.VdagC[q][sp] = eVars.VdagC[q][sp] * Haux_evecs;
		}
	}
	
	double compute(ElecGradient* grad, ElecGradient* Kgrad)
	{	Energies ener = e.ener;
		if(grad) grad->init(e);
		if(Kgrad) Kgrad->init(e);
		
		//Update fillings (Aux algorithm, fixed N only):
		double Bz, mu = eInfo.findMu(eVars.Haux_eigs, eInfo.nElectrons, Bz);
		double dmuNum[2] = {0.,0.}, dmuDen[2]={0.,0.};
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
			eVars.F[q] = eInfo.smear(eInfo.muEff(mu,Bz,q), eVars.Haux_eigs[q]);
		eInfo.updateFillingsEnergies(eVars.Haux_eigs, ener);
		
		//Update density and density-matices if needed:
		eVars.n = eVars.calcDensity();
		if(eInfo.hasU) e.iInfo.rhoAtom_calc(eVars.F, eVars.C, eVars.rhoAtom);
		
		//Update local potential:
		eVars.EdensityAndVscloc(ener, exCorr);
		if(grad) e.iInfo.augmentDensityGridGrad(eVars.Vscloc);
		
		//Wavefunction dependent parts:
		ener.E["NI"] = 0.;
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{	const QuantumNumber& qnum = eInfo.qnums[q];
			
			//KE and Nonlocal pseudopotential from precomputed subspace matrix:
			matrix HniRot = dagger(rotPrev[q]) * HniSub[q] * rotPrev[q];
			ener.E["NI"] += qnum.weight * trace(eVars.F[q] * HniRot).real();
		
			//Gradient and subspace Hamiltonian:
			if(grad)
			{	ColumnBundle HCq = Idag_DiagV_I(eVars.C[q], eVars.Vscloc); //Accumulate Idag Diag(Vscloc) I C
				if(eInfo.hasU) e.iInfo.rhoAtom_grad(eVars.C[q], eVars.U_rhoAtom, HCq); //Contribution via atomic density matrices (DFT+U)
				std::vector<matrix> HVdagCq(e.iInfo.species.size());
				e.iInfo.augmentDensitySphericalGrad(qnum, eVars.VdagC[q], HVdagCq); //Contribution via pseudopotential density augmentation
				e.iInfo.projectGrad(HVdagCq, eVars.C[q], HCq);
				eVars.Hsub[q] = HniRot + (eVars.C[q]^HCq);
				eVars.Hsub[q].diagonalize(eVars.Hsub_evecs[q], eVars.Hsub_eigs[q]);
				//N/M constraint contributions to gradient:
				diagMatrix fprime = eInfo.smearPrime(eInfo.muEff(mu,Bz,q), eVars.Haux_eigs[q]);
				double w = eInfo.qnums[q].weight;
				int sIndex = eInfo.qnums[q].index();
				dmuNum[sIndex] += w * trace(fprime * (diag(eVars.Hsub[q])-eVars.Haux_eigs[q]));
				dmuDen[sIndex] += w * trace(fprime);
			}
		}
		mpiWorld->allReduce(ener.E["NI"], MPIUtil::ReduceSum);
		
		//Final gradient propagation to auxiliary Hamiltonian:
		if(grad) 
		{	mpiWorld->allReduce(dmuNum, 2, MPIUtil::ReduceSum);
			mpiWorld->allReduce(dmuDen, 2, MPIUtil::ReduceSum);
			double dmuContrib, dBzContrib;
			if((eInfo.spinType==SpinZ) and std::isnan(eInfo.Bz))
			{	//Fixed N and M (effectively independent constraints on Nup and Ndn)
				double dmuContribUp = dmuNum[0]/dmuDen[0];
				double dmuContribDn = dmuNum[1]/dmuDen[1];
				dmuContrib = 0.5*(dmuContribUp + dmuContribDn);
				dBzContrib = 0.5*(dmuContribUp - dmuContribDn);
			}
			else
			{	//Fixed N only
				dmuContrib = (dmuNum[0]+dmuNum[1])/(dmuDen[0]+dmuDen[1]);
				dBzContrib = 0.;
			}
			for(int q=eInfo.qStart; q<eInfo.qStop; q++)
			{	const QuantumNumber& qnum = eInfo.qnums[q];
				matrix gradF = eVars.Hsub[q]-eVars.Haux_eigs[q] - eye(nBands)*eInfo.muEff(dmuContrib,dBzContrib,q); //gradient w.r.t fillings
				grad->Haux[q] = qnum.weight * eInfo.smearGrad(eInfo.muEff(mu,Bz,q), eVars.Haux_eigs[q], gradF);
				if(Kgrad) Kgrad->Haux[q] = -gradF; //Drop the fermiPrime factors and state weights in preconditioned gradient
				//Transform gradients back to original rotation (which CG remains in):
				grad->Haux[q] = rotPrev[q] * grad->Haux[q] * dagger(rotPrev[q]);
				if(Kgrad) Kgrad->Haux[q] = rotPrev[q] * Kgrad->Haux[q] * dagger(rotPrev[q]);
			}
		}
		if(!std::isnan(e.eInfo.mu)) ener.E["minusMuN"] = -e.eInfo.mu*e.eInfo.nElectrons; //Fix energy printout for fixed-mu calculation
		return ener.F();
	}

	double sync(double x) const { mpiWorld->bcast(x); return x; } //!< All processes minimize together; make sure scalars are in sync to round-off error
	
	bool report(int iter)
	{	eInfo.smearReport();
		return false;
	}
};


int ElecVars::LCAO()
{	const ElecInfo& eInfo = e->eInfo;
	const IonInfo& iInfo = e->iInfo;
	
	//Count total atomic orbitals:
	int nAtomic = iInfo.nAtomicOrbitals();
	if(nAtomic)
	{	logPrintf("linear combination of atomic orbitals\n");
		if(nAtomic < eInfo.nBands)
			logPrintf("Note: number of bands (%d) exceeds available atomic orbitals (%d)\n",
				eInfo.nBands, nAtomic);
	}
	else
	{	logPrintf("(No orbitals for LCAO)  ");
		return 0;
	}
	
	LCAOminimizer lcao(*this, *e);
	
	//Check exchange-correlation functional, and replace with PBE if not strictly (semi-)local
	ExCorr exCorrPBE; //defaults to gga-PBE
	if(e->exCorr.exxFactor() || e->exCorr.needsKEdensity() || (!e->exCorr.hasEnergy()) || e->exCorr.orbitalDep) //Hybrid, meta-GGA, potential-only or orbital-dependent respectively
	{	logPrintf("Initializing semi-local functional for LCAO:\n");
		exCorrPBE.setup(*e);
		lcao.exCorr = &exCorrPBE;
	}
	else lcao.exCorr = &e->exCorr;
	
	//Temporarily disable the fluid (which is yet to be initialized)
	FluidType fluidTypeTemp = FluidNone;
	std::swap(fluidParams.fluidType, fluidTypeTemp);
	
	//Backup original fillings:
	std::vector<diagMatrix> Forig = F;
	
	//Get orthonormal atomic orbitals and non-interacting part of subspace Hamiltonian:
	lcao.nBands = std::max(nAtomic+1, std::max(eInfo.nBands, int(ceil(1+eInfo.nElectrons/eInfo.qWeightSum))));
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
	{	C[q] = iInfo.getAtomicOrbitals(q, false, lcao.nBands-nAtomic);
		if(nAtomic<lcao.nBands) C[q].randomize(nAtomic, lcao.nBands); //Randomize extra columns if any
		orthonormalize(q);
		//Non-interacting Hamiltonian:
		ColumnBundle HniCq = -0.5*L(C[q]);
		std::vector<matrix> HVdagCq(iInfo.species.size());
		iInfo.EnlAndGrad(eInfo.qnums[q], eye(lcao.nBands), VdagC[q], HVdagCq); //non-local pseudopotentials
		iInfo.projectGrad(HVdagCq, C[q], HniCq);
		lcao.HniSub[q] = C[q]^HniCq;
		lcao.rotPrev[q] = eye(lcao.nBands);
		F[q].resize(lcao.nBands, 0.);
	}
	
	//Get electron density obtained by adding those of the atoms:
	if(!e->cntrl.fixed_H)
	{	ScalarFieldTildeArray nTilde(n.size());
		for(auto sp: e->iInfo.species)
			if(sp->atpos.size())  // Check for unused species
				sp->accumulateAtomicDensity(nTilde);
		nullToZero(nTilde, e->gInfo);
		for(unsigned s=0; s<n.size(); s++)
		{	n[s] = I(nTilde[s]);
			nTilde[s] = 0; //free memory
			e->symm.symmetrize(n[s]);
		}
	}
	
	//Select a multi-pass method to use eVars.F, if it is non-default
	int nPasses = 1;
	if(!e->cntrl.fixed_H)
	{	if(eInfo.initialFillingsFilename.length()) nPasses = 2; //custom fillings
		if(eInfo.Qinitial || eInfo.Minitial) nPasses = 2; //net charges modified
	}
	
	//Compute local-potential at the atomic reference density (for pass 1) and resulting density based on pass 1 (for pass 2):
	for(int pass=0; pass<nPasses; pass++)
	{	if(!(e->cntrl.fixed_H && e->eVars.VFilenamePattern.length())) //compute Vscloc unless it has been read in
		{	Energies ener;
			EdensityAndVscloc(ener, lcao.exCorr);
		}
		iInfo.augmentDensityInit();
		iInfo.augmentDensityGridGrad(Vscloc); //Update Vscloc projections on ultrasoft pseudopotentials
		
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{
			//Calculate subspace Hamiltonian:
			ColumnBundle HCq = Idag_DiagV_I(C[q], Vscloc); //local self-consistent potential
			std::vector<matrix> HVdagCq(iInfo.species.size());
			iInfo.augmentDensitySphericalGrad(eInfo.qnums[q], VdagC[q], HVdagCq); //ultrasoft augmentation
			iInfo.projectGrad(HVdagCq, C[q], HCq);
			Hsub[q] = dagger(lcao.rotPrev[q]) * lcao.HniSub[q] * lcao.rotPrev[q] + (C[q]^HCq);
			
			//Switch to eigenvectors of Hsub:
			Hsub[q].diagonalize(Hsub_evecs[q], Hsub_eigs[q]);
			C[q] = C[q] * Hsub_evecs[q];
			for(unsigned sp=0; sp<iInfo.species.size(); sp++)
				if(VdagC[q][sp]) VdagC[q][sp] = VdagC[q][sp] * Hsub_evecs[q]; 
			lcao.rotPrev[q] = lcao.rotPrev[q] * Hsub_evecs[q];
			Hsub[q] = Hsub_eigs[q];
			Hsub_evecs[q] = eye(lcao.nBands);
		}
		
		if(pass+1<nPasses) n = calcDensity(); //only needed here for a subsequent pass
	}
	
	//Subspace minimize:
	Haux_eigs = Hsub_eigs;
	if(!e->cntrl.fixed_H) //current Hsub suffices for fixed_H
	{	MinimizeParams mp;
		mp.nDim = eInfo.nStates * lcao.nBands*lcao.nBands;
		mp.fpLog = globalLog;
		mp.linePrefix = "LCAOMinimize: ";
		mp.energyLabel = relevantFreeEnergyName(*e);
		mp.energyFormat = "%+.16f";
		mp.energyDiffThreshold = lcaoTol;
		mp.nIterations = (lcaoIter>=0) ? lcaoIter : ( eInfo.fillingsUpdate==ElecInfo::FillingsHsub ? 30 : 3 );
		lcao.minimize(mp);
	}
	
	//Cut wavefunctions and subspace Hamiltonia back down to size:
	if(eInfo.nBands<lcao.nBands)
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{	Hsub[q] = Hsub[q](0,eInfo.nBands, 0,eInfo.nBands);
			Hsub[q].diagonalize(Hsub_evecs[q], Hsub_eigs[q]);
			C[q] = C[q].getSub(0,eInfo.nBands);
			Haux_eigs[q].resize(eInfo.nBands);
		}
	
	//Transition fillings :
	if(eInfo.fillingsUpdate==ElecInfo::FillingsHsub)
	{	//Hsub fillings: use same eigenvalues, but recalculate to account for reduced band count:
		double Bz, mu = eInfo.findMu(Haux_eigs, eInfo.nElectrons, Bz);
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
			F[q] = eInfo.smear(eInfo.muEff(mu,Bz,q), Haux_eigs[q]);
	}
	else
	{	//Constant fillings: remove Haux_eigs and restore F to original:
		Haux_eigs.clear();
		F = Forig;
	}
	HauxInitialized = true;
	
	std::swap(fluidParams.fluidType, fluidTypeTemp); //Restore the fluid type
	return eInfo.nBands;
}
