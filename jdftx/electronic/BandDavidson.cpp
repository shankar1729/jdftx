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

#include <electronic/BandDavidson.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>

BandDavidson::BandDavidson(Everything& e, int q): e(e), eVars(e.eVars), eInfo(e.eInfo), q(q)
{	assert(e.cntrl.fixed_H); // Check whether the electron Hamiltonian is fixed
}

void BandDavidson::minimize()
{	//Use the same working set as the CG minimizer:
	ColumnBundle& C = eVars.C[q];
	std::vector<matrix>& VdagC = eVars.VdagC[q];
	matrix& Hsub = eVars.Hsub[q];
	matrix& Hsub_evecs = eVars.Hsub_evecs[q];
	diagMatrix& Hsub_eigs = eVars.Hsub_eigs[q];
	const QuantumNumber& qnum = eInfo.qnums[q];
	int nBandsOut = eInfo.nBands; //number of final output bands desired
	int nBandsMax = ceil(e.cntrl.davidsonBandRatio * nBandsOut);
	if(2*nBandsMax >= int(C.basis->nbasis))
		die_alone("Cannot use Davidson eigenvalue algorithm when 2 x nBands x davidsonBandRatio > nBasis.\n"
			"Reduce nBands or davidsonBandRatio, increase nBasis (Ecut) or use elec-eigen-algo CG.\n\n");
	
	//Initial subspace eigenvalue problem:
	ColumnBundle HC;
	diagMatrix I = eye(nBandsOut);
	Energies ener; //not really used here
	eVars.applyHamiltonian(q, I, HC, ener, true);
	Hsub = C^HC;
	Hsub.diagonalize(Hsub_evecs, Hsub_eigs);
	//--- switch C to subspace eigenbasis:
	C = C * Hsub_evecs;
	HC = HC * Hsub_evecs;
	e.iInfo.project(C, VdagC, &Hsub_evecs);
	double Eband = qnum.weight * trace(Hsub_eigs);
	logPrintf("BandDavidson: Iter: %3d  Eband: %+.15lf\n", 0, Eband); fflush(globalLog);
	
	const MinimizeParams& mp = e.elecMinParams;
	int iter=1;
	for(; iter<=mp.nIterations; iter++)
	{	int nBands = C.nCols();
		//Compute subspace expansion:
		diagMatrix KEref = (-0.5) * diagDot(C, L(C)); //Update reference KE for preconditioning:
		ColumnBundle Cexp = HC; Cexp -= O(C) * Hsub_eigs; //Calculate residual of current eigenvector guesses
		precond_inv_kinetic_band(Cexp, KEref); //Davidson approximate inverse (using KE as the diagonal)
		//Drop converged eigenpairs and approximately normalize subspace expansion (for avoiding roundoff issues only):
		diagMatrix CexpNorm = diagDot(Cexp, Cexp);
		double CexpNormCut = std::max(mp.energyDiffThreshold/nBands, 1e-15*Cexp.colLength());
		{	//Drop columns whose norm falls below above cutoff
			complex* CexpData = Cexp.dataPref();
			int bOut = 0;
			for(int b=0; b<nBands; b++)
			{	if(CexpNorm[b]<CexpNormCut) continue;
				CexpNorm[bOut] = 1/sqrt(CexpNorm[b]);
				if(bOut<b) callPref(eblas_copy)(CexpData+Cexp.index(bOut,0), CexpData+Cexp.index(b,0), Cexp.colLength());
				bOut++;
			}
			if(!bOut) //This is unlikely, but just in case (to avoid zero column matrices below)
			{	logPrintf("BandDavidson: Converged (dEband<%le)\n", mp.energyDiffThreshold);
				break;
			}
			if(bOut<nBands)
			{	Cexp = Cexp.getSub(0,bOut);
				CexpNorm = CexpNorm(0,bOut);
			}
		}
		Cexp = Cexp * CexpNorm;
		int nBandsNew = Cexp.nCols();
		int nBandsBig = nBands + nBandsNew;
		//Expansion subspace overlaps:
		matrix bigOsub(nBandsBig, nBandsBig);
		std::vector<matrix> VdagCexp;
		{	ColumnBundle OCexp = O(Cexp, &VdagCexp);
			matrix rotExisting = eye(Cexp.nCols());
			e.iInfo.project(Cexp, VdagCexp, &rotExisting);
			matrix CdagOCexp = C ^ OCexp;
			bigOsub.set(0,nBands, 0,nBands, eye(nBands)); //since C's are already orthonormal
			bigOsub.set(nBands,nBandsBig, nBands,nBandsBig, Cexp ^ OCexp);
			bigOsub.set(0,nBands, nBands,nBandsBig, CdagOCexp);
			bigOsub.set(nBands,nBandsBig, 0,nBands, dagger(CdagOCexp));
		}
		//Expansion subspace Hamiltonian:
		matrix bigHsub(nBandsBig, nBandsBig);
		ColumnBundle HCexp;
		{	matrix HsubExp; diagMatrix HsubExp_eigs;
			#define SWAP_C_Cexp \
				std::swap(C, Cexp); \
				std::swap(VdagC, VdagCexp); \
				std::swap(Hsub, HsubExp); \
				std::swap(Hsub_eigs, HsubExp_eigs);
			SWAP_C_Cexp //Temporarily swap C and Cexp
			eVars.applyHamiltonian(q, eye(nBandsNew), HCexp, ener, true); //Hamiltonian always operates on C, where we put Cexp 
			SWAP_C_Cexp  //Restore C and Cexp to correct places
			matrix CdagHCexp = C  ^ HCexp;
			bigHsub.set(0,nBands, 0,nBands, Hsub_eigs);
			bigHsub.set(nBands,nBandsBig, nBands,nBandsBig, HsubExp);
			bigHsub.set(0,nBands, nBands,nBandsBig, CdagHCexp);
			bigHsub.set(nBands,nBandsBig, 0,nBands, dagger(CdagHCexp));
		}
		//Solve expanded subspace generalized eigenvalue problem:
		matrix bigU = invsqrt(bigOsub);
		bigHsub = dagger_symmetrize(dagger(bigU) * bigHsub * bigU); //switch to the symmetrically-orthonormalized basis
		matrix bigHsub_evecs; diagMatrix bigHsub_eigs;
		bigHsub.diagonalize(bigHsub_evecs, bigHsub_eigs);
		matrix rot = bigU * bigHsub_evecs; //rotation from [C,Cexp] to the expanded subspace eigenbasis
		int nBandsNext = std::min(nBandsMax, nBandsBig); //number of bands to retain for next iteration
		matrix Crot = rot(0,nBands, 0,nBandsNext); //contribution of C to lowest nBandsNext eigenvectors
		matrix CexpRot = rot(nBands,nBandsBig, 0,nBandsNext); //contribution of Cexp to lowest nBandsNext eigenvectors
		//Update C to optimum nBands subspace from [C,Cexp]
		C = C*Crot + Cexp*CexpRot;
		HC = HC*Crot + HCexp*CexpRot;
		Hsub_eigs = bigHsub_eigs(0,nBandsNext);
		for(size_t sp=0; sp<VdagC.size(); sp++) if(VdagC[sp])
			VdagC[sp] = VdagC[sp]*Crot + VdagCexp[sp]*CexpRot;
		//Print and test convergence
		double EbandPrev = Eband;
		Eband = qnum.weight * trace(Hsub_eigs(0,nBandsOut));
		double dEband = Eband - EbandPrev;
		logPrintf("BandDavidson: Iter: %3d  Eband: %+.15lf  dEband: %le  t[s]: %9.2lf\n", iter, Eband, dEband, clock_sec()); fflush(globalLog);
		if(dEband<0 and fabs(dEband)<mp.energyDiffThreshold)
		{	logPrintf("BandDavidson: Converged (dEband<%le)\n", mp.energyDiffThreshold);
			break;
		}
	}
	if(iter>mp.nIterations)
		logPrintf("BandDavidson: None of the convergence criteria satisfied after %d iterations.\n", mp.nIterations);
	fflush(globalLog);
	
	//Update final quantities:
	if(C.nCols() != nBandsOut)
	{	//reduce outputs to size:
		C = C.getSub(0, nBandsOut);
		for(size_t sp=0; sp<VdagC.size(); sp++) if(VdagC[sp])
			VdagC[sp] = VdagC[sp](0,VdagC[sp].nRows(), 0,nBandsOut);
		Hsub_eigs = Hsub_eigs(0,nBandsOut);
	}
	Hsub = Hsub_eigs;
	Hsub_evecs = I;
}
