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
#include <core/DataIO.h>

struct LCAOminimizer : Minimizable<ElecGradient> //Uses only the B entries of ElecGradient
{	int nBands;
	ElecVars& eVars;
	const Everything& e;
	const ExCorr* exCorr;
	std::vector<matrix> HniSub;
	std::vector<matrix> B; //Auxiliary hamiltonian (minimizer state)
	std::vector< std::vector<matrix> > VdagY;
	ElecGradient Kgrad;
	
	LCAOminimizer(ElecVars& eVars, const Everything& e)
	: eVars(eVars), e(e), HniSub(e.eInfo.nStates), B(e.eInfo.nStates),
		VdagY(e.eInfo.nStates, std::vector<matrix>(e.iInfo.species.size()))
	{	Kgrad.init(e);
	}
	
	void step(const ElecGradient& dir, double alpha)
	{	assert(dir.B.size() == B.size());
		for(unsigned q=0; q<dir.Y.size(); q++)
			if(dir.B[q]) axpy(alpha, dir.B[q], B[q]);
	}
	
	double compute(ElecGradient* grad)
	{	Energies ener = e.ener;
		if(grad) grad->init(e);
		
		//Simplified version of ElecVars::orthonormalize()
		std::vector<matrix> B_evecs(e.eInfo.nStates);
		std::vector<diagMatrix> B_eigs(e.eInfo.nStates);
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		{	B[q].diagonalize(B_evecs[q], B_eigs[q]);
			eVars.C[q] = eVars.Y[q] * B_evecs[q];
			for(unsigned sp=0; sp<e.iInfo.species.size(); sp++)
				if(VdagY[q][sp]) eVars.VdagC[q][sp] = VdagY[q][sp] * B_evecs[q]; 
		}
		
		//Update fillings (Aux algorithm, fixed N only):
		double mu = e.eInfo.findMu(B_eigs, e.eInfo.nElectrons), dmuNumDen[2] = { 0., 0. };
		std::vector<diagMatrix> F(e.eInfo.nStates);
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
			F[q] = e.eInfo.fermi(mu, B_eigs[q]);
		e.eInfo.updateFillingsEnergies(F, ener);
		
		//Update density:
		for(DataRptr& ns: eVars.n) ns=0;
		e.iInfo.augmentDensityInit();
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		{	eVars.n[e.eInfo.qnums[q].index()] += e.eInfo.qnums[q].weight * diagouterI(F[q], eVars.C[q], &e.gInfo);
			e.iInfo.augmentDensitySpherical(e.eInfo.qnums[q], F[q], eVars.VdagC[q]); //pseudopotential contribution
		}
		e.iInfo.augmentDensityGrid(eVars.n);
		for(DataRptr& ns: eVars.n)
		{	nullToZero(ns, e.gInfo);
			e.symm.symmetrize(ns);
			ns->allReduce(MPIUtil::ReduceSum);
		}
		
		//Update local potential:
		eVars.EdensityAndVscloc(ener, exCorr);
		if(grad) e.iInfo.augmentDensityGridGrad(eVars.Vscloc);
		
		//Wavefunction dependent parts:
		std::vector<ColumnBundle> HC(e.eInfo.nStates);
		std::vector< std::vector<matrix> > HVdagC(e.eInfo.nStates, std::vector<matrix>(e.iInfo.species.size()));

		ener.E["U"] = e.iInfo.computeU(F, eVars.C, grad ? &HC : 0);
		ener.E["NI"] = 0.;
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		{	const QuantumNumber& qnum = e.eInfo.qnums[q];
			
			//KE and Nonlocal pseudopotential from precomputed subspace matrix:
			matrix HniRot = dagger(B_evecs[q]) * HniSub[q] * B_evecs[q];
			ener.E["NI"] += qnum.weight * trace(F[q] * HniRot).real();
		
			//Gradient and subspace Hamiltonian:
			if(grad)
			{	HC[q] += Idag_DiagV_I(eVars.C[q], eVars.Vscloc[qnum.index()]); //Accumulate Idag Diag(Vscloc) I C
				e.iInfo.augmentDensitySphericalGrad(qnum, F[q], eVars.VdagC[q], HVdagC[q]); //Contribution via pseudopotential density augmentation
				e.iInfo.projectGrad(HVdagC[q], eVars.C[q], HC[q]);
				eVars.Hsub[q] = HniRot + (eVars.C[q]^HC[q]);
				//Nconstraint contributions to gradient:
				diagMatrix fprime = e.eInfo.fermiPrime(mu, B_eigs[q]);
				dmuNumDen[0] += qnum.weight * trace(fprime * (diag(eVars.Hsub[q])-B_eigs[q]));
				dmuNumDen[1] += qnum.weight * trace(fprime);
			}
		}
		mpiUtil->allReduce(ener.E["NI"], MPIUtil::ReduceSum);
		
		//Final gradient propagation to auxiliary Hamiltonian:
		if(grad) 
		{	mpiUtil->allReduce(dmuNumDen, 2, MPIUtil::ReduceSum);
			matrix dmuContrib = eye(nBands) * (dmuNumDen[0]/dmuNumDen[1]); //contribution due to Nconstraint via the mu gradient 
			for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
			{	const QuantumNumber& qnum = e.eInfo.qnums[q];
				matrix gradF = eVars.Hsub[q]-B_eigs[q]-dmuContrib; //gradient w.r.t fillings
				grad->B[q] = qnum.weight * dagger_symmetrize(B_evecs[q] * e.eInfo.fermiGrad(mu, B_eigs[q], gradF) * dagger(B_evecs[q]));
				//Drop the fermiPrime factors and state weights in preconditioned gradient:
				Kgrad.B[q] = dagger_symmetrize(B_evecs[q] * (-gradF) * dagger(B_evecs[q]));
			}
		}
		return ener.F();
	}
	
	ElecGradient precondition(const ElecGradient& grad)
	{	return Kgrad;
	}

	double sync(double x) const { mpiUtil->bcast(x); return x; } //!< All processes minimize together; make sure scalars are in sync to round-off error
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
	lcao.nBands = std::max(nAtomic+1, std::max(eInfo.nBands, int(ceil(1+eInfo.nElectrons/2))));
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
	{	Y[q] = iInfo.getAtomicOrbitals(q, lcao.nBands-nAtomic);
		if(nAtomic<lcao.nBands) Y[q].randomize(nAtomic, lcao.nBands); //Randomize extra columns if any
		matrix YtoC = invsqrt(Y[q]^O(Y[q], &lcao.VdagY[q]));
		Y[q] = Y[q] * YtoC; //orthonormalize
		iInfo.project(Y[q], lcao.VdagY[q], &YtoC);
		//Non-interacting Hamiltonian:
		ColumnBundle HniYq = -0.5*L(Y[q]);
		std::vector<matrix> HVdagYq(iInfo.species.size());
		iInfo.EnlAndGrad(eInfo.qnums[q], eye(lcao.nBands), lcao.VdagY[q], HVdagYq); //non-local pseudopotentials
		iInfo.projectGrad(HVdagYq, Y[q], HniYq);
		lcao.HniSub[q] = Y[q]^HniYq;
	}
	
	//Get electron density obtained by adding those of the atoms:
	if(!e->cntrl.fixed_H)
	{	DataGptrCollection nTilde(n.size());
		for(auto sp: e->iInfo.species)
			if(sp->atpos.size())  // Check for unused species
				sp->accumulateAtomicDensity(nTilde);
		for(unsigned s=0; s<n.size(); s++)
		{	n[s] = I(nTilde[s]);
			nTilde[s] = 0; //free memory
			e->symm.symmetrize(n[s]);
		}
	}
	
	//Select a multi-pass method to use eVars.F, if it is non-default
	int nPasses = 1;
	if(!e->cntrl.fixed_H && (eInfo.customFillings.size() || eInfo.initialFillingsFilename.length())) nPasses = 2; //custom fillings
	for(double q: eInfo.qNet) if(q) nPasses = 2; //net charges modified
	
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
		{	ColumnBundle HYq = Idag_DiagV_I(Y[q], Vscloc[eInfo.qnums[q].index()]); //local self-consistent potential
			std::vector<matrix> HVdagYq(iInfo.species.size());
			iInfo.augmentDensitySphericalGrad(eInfo.qnums[q], eye(lcao.nBands), lcao.VdagY[q], HVdagYq); //ultrasoft augmentation
			iInfo.projectGrad(HVdagYq, Y[q], HYq);
			lcao.B[q] = lcao.HniSub[q] + (Y[q]^HYq);
		}
		
		if(nPasses==2 && pass==0) //update the density for next pass
		{	for(DataRptr& ns: n) ns=0;
			iInfo.augmentDensityInit();
			for(int q=eInfo.qStart; q<eInfo.qStop; q++)
			{	matrix Bq_evecs; diagMatrix Bq_eigs; lcao.B[q].diagonalize(Bq_evecs, Bq_eigs);
				C[q] = Y[q] * Bq_evecs;
				for(unsigned sp=0; sp<iInfo.species.size(); sp++)
					if(lcao.VdagY[q][sp]) VdagC[q][sp] = lcao.VdagY[q][sp] * Bq_evecs; 
				diagMatrix Fq = F[q]; Fq.resize(lcao.nBands, 0.); //length-corrected eVars.F
				n[eInfo.qnums[q].index()] += eInfo.qnums[q].weight * diagouterI(Fq, C[q], &e->gInfo);
				iInfo.augmentDensitySpherical(eInfo.qnums[q], Fq, VdagC[q]); //pseudopotential contribution
			}
			iInfo.augmentDensityGrid(n);
			for(DataRptr& ns: n)
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
	mp.energyDiffThreshold = lcaoTol;
	mp.dirUpdateScheme = MinimizeParams::LBFGS;
	mp.history = 10;
	mp.alphaTincreaseFactor = 10;
	mp.nIterations = (lcaoIter>=0) ? lcaoIter : ( eInfo.subspaceRotation ? 30 : 3 );
	if(e->cntrl.fixed_H) { C = Y; Hsub = lcao.B; } //bypass subspace iteration
	else lcao.minimize(mp);
	
	//Set wavefunctions to eigenvectors:
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
	{	Hsub[q].diagonalize(Hsub_evecs[q], Hsub_eigs[q]);
		if(eInfo.nBands<lcao.nBands)
		{	Hsub_evecs[q] = Hsub_evecs[q](0,lcao.nBands, 0,eInfo.nBands); //drop extra eigenvectors
			Hsub_eigs[q] = Hsub_eigs[q](0,eInfo.nBands); //drop extra eigenvalues
		}
		Y[q] = C[q] * Hsub_evecs[q]; C[q].free();
		if(eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux)
		{	matrix Bq_evecs; diagMatrix Bq_eigs;
			lcao.B[q].diagonalize(Bq_evecs, Bq_eigs);
			B[q] = dagger(Hsub_evecs[q]) * Bq_eigs * Hsub_evecs[q];
			HauxInitialized = true;
		}
		Hsub_evecs[q] = eye(eInfo.nBands);
		Hsub[q] = Hsub_eigs[q];
	}
	if(eInfo.fillingsUpdate!=ElecInfo::ConstantFillings)
	{	double mu = eInfo.findMu(Hsub_eigs, eInfo.nElectrons);
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
			F[q] = eInfo.fermi(mu, Hsub_eigs[q]);
	}
	
	std::swap(fluidParams.fluidType, fluidTypeTemp); //Restore the fluid type
	return eInfo.nBands;
}
