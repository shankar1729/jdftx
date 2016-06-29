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
#include <electronic/matrix.h>
#include <electronic/operators.h>
#include <electronic/ElecMinimizer.h>
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
	std::vector<matrix> Haux; //Auxiliary hamiltonian (minimizer state)
	std::vector<matrix> rotPrev; //Accumulated rotations of the wavefunctions
	ElecGradient Kgrad;
	
	LCAOminimizer(ElecVars& eVars, const Everything& e)
	: eVars(eVars), e(e), eInfo(e.eInfo), HniSub(eInfo.nStates), Haux(eInfo.nStates), rotPrev(eInfo.nStates)
	{	Kgrad.init(e);
	}
	
	void step(const ElecGradient& dir, double alpha)
	{	assert(dir.Haux.size() == Haux.size());
		for(unsigned q=0; q<dir.Haux.size(); q++)
			if(dir.Haux[q]) axpy(alpha, dir.Haux[q], Haux[q]);
	}
	
	double compute(ElecGradient* grad)
	{	Energies ener = e.ener;
		if(grad) grad->init(e);
		
		//Simplified version of ElecVars::orthonormalize()
		std::vector<matrix> B_evecs(eInfo.nStates);
		std::vector<diagMatrix> B_eigs(eInfo.nStates);
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{	Haux[q].diagonalize(B_evecs[q], B_eigs[q]);
			matrix rot = B_evecs[q];
			if(rotPrev[q])
				rot = dagger(rotPrev[q]) * rot; //account for previous rotation
			eVars.C[q] = eVars.C[q] * rot;
			for(unsigned sp=0; sp<e.iInfo.species.size(); sp++)
				if(eVars.VdagC[q][sp]) eVars.VdagC[q][sp] = eVars.VdagC[q][sp] * rot; 
			rotPrev[q] = B_evecs[q]; //total rotation thus far
		}
		
		//Update fillings (Aux algorithm, fixed N only):
		double Bz, mu = eInfo.findMu(B_eigs, eInfo.nElectrons, Bz);
		double dmuNum[2] = {0.,0.}, dmuDen[2]={0.,0.};
		std::vector<diagMatrix> F(eInfo.nStates);
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
			F[q] = eInfo.fermi(eInfo.muEff(mu,Bz,q), B_eigs[q]);
		eInfo.updateFillingsEnergies(F, ener);
		
		//Update density:
		for(ScalarField& ns: eVars.n) ns=0;
		e.iInfo.augmentDensityInit();
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{	eVars.n += eInfo.qnums[q].weight * diagouterI(F[q], eVars.C[q], eVars.n.size(), &e.gInfo);
			e.iInfo.augmentDensitySpherical(eInfo.qnums[q], F[q], eVars.VdagC[q]); //pseudopotential contribution
		}
		e.iInfo.augmentDensityGrid(eVars.n);
		for(ScalarField& ns: eVars.n)
		{	nullToZero(ns, e.gInfo);
			e.symm.symmetrize(ns);
			ns->allReduce(MPIUtil::ReduceSum);
		}
		if(eInfo.hasU) e.iInfo.rhoAtom_calc(F, eVars.C, eVars.rhoAtom);
		
		//Update local potential:
		eVars.EdensityAndVscloc(ener, exCorr);
		if(grad) e.iInfo.augmentDensityGridGrad(eVars.Vscloc);
		
		//Wavefunction dependent parts:
		std::vector<ColumnBundle> HC(eInfo.nStates);
		std::vector< std::vector<matrix> > HVdagC(eInfo.nStates, std::vector<matrix>(e.iInfo.species.size()));

		ener.E["NI"] = 0.;
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{	const QuantumNumber& qnum = eInfo.qnums[q];
			
			//KE and Nonlocal pseudopotential from precomputed subspace matrix:
			matrix HniRot = dagger(B_evecs[q]) * HniSub[q] * B_evecs[q];
			ener.E["NI"] += qnum.weight * trace(F[q] * HniRot).real();
		
			//Gradient and subspace Hamiltonian:
			if(grad)
			{	HC[q] += Idag_DiagV_I(eVars.C[q], eVars.Vscloc); //Accumulate Idag Diag(Vscloc) I C
				if(eInfo.hasU) e.iInfo.rhoAtom_grad(eVars.C[q], eVars.U_rhoAtom, HC[q]); //Contribution via atomic density matrices (DFT+U)
				e.iInfo.augmentDensitySphericalGrad(qnum, F[q], eVars.VdagC[q], HVdagC[q]); //Contribution via pseudopotential density augmentation
				e.iInfo.projectGrad(HVdagC[q], eVars.C[q], HC[q]);
				eVars.Hsub[q] = HniRot + (eVars.C[q]^HC[q]);
				eVars.Hsub[q].diagonalize(eVars.Hsub_evecs[q], eVars.Hsub_eigs[q]);
				//N/M constraint contributions to gradient:
				diagMatrix fprime = eInfo.fermiPrime(eInfo.muEff(mu,Bz,q), B_eigs[q]);
				double w = eInfo.qnums[q].weight;
				int sIndex = eInfo.qnums[q].index();
				dmuNum[sIndex] += w * trace(fprime * (diag(eVars.Hsub[q])-B_eigs[q]));
				dmuDen[sIndex] += w * trace(fprime);
			}
		}
		mpiUtil->allReduce(ener.E["NI"], MPIUtil::ReduceSum);
		
		//Final gradient propagation to auxiliary Hamiltonian:
		if(grad) 
		{	mpiUtil->allReduce(dmuNum, 2, MPIUtil::ReduceSum);
			mpiUtil->allReduce(dmuDen, 2, MPIUtil::ReduceSum);
			double dmuContrib, dBzContrib;
			if(eInfo.Mconstrain)
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
				matrix gradF = eVars.Hsub[q]-B_eigs[q] - eye(nBands)*eInfo.muEff(dmuContrib,dBzContrib,q); //gradient w.r.t fillings
				grad->Haux[q] = qnum.weight * dagger_symmetrize(B_evecs[q] * eInfo.fermiGrad(eInfo.muEff(mu,Bz,q), B_eigs[q], gradF) * dagger(B_evecs[q]));
				//Drop the fermiPrime factors and state weights in preconditioned gradient:
				Kgrad.Haux[q] = dagger_symmetrize(B_evecs[q] * (-gradF) * dagger(B_evecs[q]));
			}
		}
		return ener.F();
	}
	
	ElecGradient precondition(const ElecGradient& grad)
	{	return Kgrad;
	}

	double sync(double x) const { mpiUtil->bcast(x); return x; } //!< All processes minimize together; make sure scalars are in sync to round-off error
	
	bool report(int iter)
	{	eInfo.printFermi();
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
		
		//Set initial auxiliary hamiltonian to the subspace Hamiltonian:
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{	ColumnBundle HCq = Idag_DiagV_I(C[q], Vscloc); //local self-consistent potential
			std::vector<matrix> HVdagCq(iInfo.species.size());
			iInfo.augmentDensitySphericalGrad(eInfo.qnums[q], eye(lcao.nBands), VdagC[q], HVdagCq); //ultrasoft augmentation
			iInfo.projectGrad(HVdagCq, C[q], HCq);
			lcao.Haux[q] = lcao.HniSub[q] + (C[q]^HCq);
		}
		
		if(nPasses==2 && pass==0) //update the density for next pass
		{	for(ScalarField& ns: n) ns=0;
			iInfo.augmentDensityInit();
			for(int q=eInfo.qStart; q<eInfo.qStop; q++)
			{	matrix Bq_evecs; diagMatrix Bq_eigs; lcao.Haux[q].diagonalize(Bq_evecs, Bq_eigs);
				C[q] = C[q] * Bq_evecs;
				for(unsigned sp=0; sp<iInfo.species.size(); sp++)
					if(VdagC[q][sp]) VdagC[q][sp] = VdagC[q][sp] * Bq_evecs; 
				lcao.rotPrev[q] = Bq_evecs;
				diagMatrix Fq = F[q]; Fq.resize(lcao.nBands, 0.); //length-corrected eVars.F
				n += eInfo.qnums[q].weight * diagouterI(Fq, C[q], n.size(), &e->gInfo);
				iInfo.augmentDensitySpherical(eInfo.qnums[q], Fq, VdagC[q]); //pseudopotential contribution
			}
			iInfo.augmentDensityGrid(n);
			for(ScalarField& ns: n)
			{	nullToZero(ns, e->gInfo);
				e->symm.symmetrize(ns);
				ns->allReduce(MPIUtil::ReduceSum);
			}
		}
	}
	
	//Subspace minimize:
	MinimizeParams mp;
	mp.nDim = eInfo.nStates * lcao.nBands*lcao.nBands;
	mp.fpLog = globalLog;
	mp.linePrefix = "LCAOMinimize: ";
	mp.energyLabel = "F";
	mp.energyFormat = "%+.16f";
	mp.energyDiffThreshold = lcaoTol;
	mp.nIterations = (lcaoIter>=0) ? lcaoIter : ( eInfo.fillingsUpdate==ElecInfo::FillingsHsub ? 30 : 3 );
	if(e->cntrl.fixed_H) { Hsub = lcao.Haux; } //bypass subspace iteration
	else lcao.minimize(mp);
	
	//Set wavefunctions to eigenvectors:
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
	{	Hsub[q].diagonalize(Hsub_evecs[q], Hsub_eigs[q]);
		if(eInfo.nBands<lcao.nBands)
		{	Hsub_evecs[q] = Hsub_evecs[q](0,lcao.nBands, 0,eInfo.nBands); //drop extra eigenvectors
			Hsub_eigs[q] = Hsub_eigs[q](0,eInfo.nBands); //drop extra eigenvalues
		}
		C[q] = C[q] * Hsub_evecs[q];
		Hsub_evecs[q] = eye(eInfo.nBands);
		Hsub[q] = Hsub_eigs[q];
		if(eInfo.fillingsUpdate==ElecInfo::FillingsHsub)
			Haux_eigs[q] = Hsub_eigs[q];
	}
	if(eInfo.fillingsUpdate==ElecInfo::FillingsHsub)
	{	double Bz, mu = eInfo.findMu(Hsub_eigs, eInfo.nElectrons, Bz);
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
			F[q] = eInfo.fermi(eInfo.muEff(mu,Bz,q), Hsub_eigs[q]);
	}
	HauxInitialized = true;
	
	std::swap(fluidParams.fluidType, fluidTypeTemp); //Restore the fluid type
	return eInfo.nBands;
}
