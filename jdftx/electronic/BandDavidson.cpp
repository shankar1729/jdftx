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
#include <electronic/operators.h>

BandDavidson::BandDavidson(Everything& e, int qActive): qActive(qActive), e(e), eVars(e.eVars), eInfo(e.eInfo)
{	assert(e.cntrl.fixed_H); // Check whether the electron Hamiltonian is fixed
}

void BandDavidson::minimize()
{	//Use the same working set as the CG minimizer:
	ColumnBundle& C = eVars.C[qActive]; //Use as expansion subset
	ColumnBundle Y = C.similar();       //Use to store current guesses to eigenvectors
	std::vector<matrix>& VdagC = eVars.VdagC[qActive];
	const QuantumNumber& qnum = eInfo.qnums[qActive];
	int nBandsOut = eInfo.nBands; //number of final output bands desired
	int nBandsMax = ceil(e.cntrl.davidsonBandRatio * nBandsOut);
	
	//Initial subspace eigenvalue problem:
	ColumnBundle OY = O(Y, &VdagC), HY;
	//--- compute H * Y
	C = Y; //Note Hamiltonian is always applied to C (as this is what makes sense in the CG context)
	e.iInfo.project(C, VdagC);
	diagMatrix I = eye(nBandsOut);
	Energies ener; //not really used here
	eVars.applyHamiltonian(qActive, I, HY, ener, true);
	//--- solve subspace eigenvalue problem
	matrix U = invsqrt(Y^OY); //symmetric orthonormalization
	matrix Hsub = dagger(U) * (Y^HY) * U;
	matrix Hsub_evecs; diagMatrix Hsub_eigs;
	Hsub.diagonalize(Hsub_evecs, Hsub_eigs);
	//--- switch Y to subspace eigenbasis:
	matrix initialRot = U * Hsub_evecs;
	Y = Y * initialRot;
	OY = OY * initialRot;
	HY = HY * initialRot;
	double Eband = qnum.weight * trace(Hsub_eigs);
	logPrintf("BandDavidson: Iter: %3d  Eband: %+.15lf\n", 0, Eband); fflush(globalLog);
	
	const MinimizeParams& mp = e.elecMinParams;
	int iter=1;
	for(; iter<=mp.nIterations; iter++)
	{	int nBands = Y.nCols();
		//Compute subspace expansion:
		diagMatrix KEref = (-0.5) * diagDot(Y, L(Y)); //Update reference KE for preconditioning:
		C = HY; C -= OY * Hsub_eigs; //Set C to residual of current eigenvector guesses
		precond_inv_kinetic_band(C, KEref); //Davidson approximate inverse (using KE as the diagonal)
		//Drop converged eigenpairs and approximately normalize subspace expansion (for avoiding roundoff issues only):
		diagMatrix Cnorm = diagDot(C, C);
		double CnormCut = std::max(mp.energyDiffThreshold/nBands, 1e-15*C.colLength());
		{	//Drop columns whose norm falls below above cutoff
			complex* Cdata = C.dataPref();
			int bOut = 0;
			for(int b=0; b<nBands; b++)
			{	if(Cnorm[b]<CnormCut) continue;
				Cnorm[bOut] = 1/sqrt(Cnorm[b]);
				if(bOut<b) callPref(eblas_copy)(Cdata+C.index(bOut,0), Cdata+C.index(b,0), C.colLength());
				bOut++;
			}
			if(!bOut) //This is unlikely, but just in case (to avoid zero column matrices below)
			{	logPrintf("BandDavidson: Converged (dEband<%le)\n", mp.energyDiffThreshold);
				break;
			}
			if(bOut<nBands)
			{	C = C.getSub(0,bOut);
				Cnorm = Cnorm(0,bOut);
			}
		}
		C = C * Cnorm;
		int nBandsNew = C.nCols();
		int nBandsBig = nBands + nBandsNew;
		//Compute overlap and Hamiltonian on expansion:
		ColumnBundle OC = O(C, &VdagC), HC;
		e.iInfo.project(C, VdagC);
		eVars.applyHamiltonian(qActive, eye(nBandsNew), HC, ener, true);
		//Setup matrices for expanded subspace generalized eigenvalue problem:
		matrix bigOsub(nBandsBig, nBandsBig);
		{	bigOsub.set(0,nBands, 0,nBands, eye(nBands)); //since Y's are already orthonormal
			bigOsub.set(nBands,nBandsBig, nBands,nBandsBig, C^OC);
			matrix YdagOC = Y^OC;
			bigOsub.set(0,nBands, nBands,nBandsBig, YdagOC);
			bigOsub.set(nBands,nBandsBig, 0,nBands, dagger(YdagOC));
		}
		matrix bigHsub(nBandsBig, nBandsBig);
		{	bigHsub.set(0,nBands, 0,nBands, Hsub_eigs);
			bigHsub.set(nBands,nBandsBig, nBands,nBandsBig, C^HC);
			matrix YdagHC = Y^HC;
			bigHsub.set(0,nBands, nBands,nBandsBig, YdagHC);
			bigHsub.set(nBands,nBandsBig, 0,nBands, dagger(YdagHC));
		}
		//Solve expanded subspace generalized eigenvalue problem:
		matrix bigU = invsqrt(bigOsub);
		bigHsub = dagger_symmetrize(dagger(bigU) * bigHsub * bigU); //switch to the symmetrically-orthonormalized basis
		matrix bigHsub_evecs; diagMatrix bigHsub_eigs;
		bigHsub.diagonalize(bigHsub_evecs, bigHsub_eigs);
		matrix rot = bigU * bigHsub_evecs; //rotation from [Y,C] to the expanded subspace eigenbasis
		int nBandsNext = std::min(nBandsMax, nBandsBig); //number of bands to retain for next iteration
		matrix Yrot = rot(0,nBands, 0,nBandsNext); //contribution of Y to lowest nBands eigenvectors
		matrix Crot = rot(nBands,nBandsBig, 0,nBandsNext); //contribution of C to lowest nBands eigenvectors
		//Update Y to optimum nBands subspace from [Y,C]
		Y = Y*Yrot + C*Crot;
		OY = OY*Yrot + OC*Crot;
		HY = HY*Yrot + HC*Crot;
		Hsub_eigs = bigHsub_eigs(0,nBandsNext);
		//Print and test convergence
		double EbandPrev = Eband;
		Eband = qnum.weight * trace(Hsub_eigs(0,nBandsOut));
		double dEband = Eband - EbandPrev;
		logPrintf("BandDavidson: Iter: %3d  Eband: %+.15lf  dEband: %le\n", iter, Eband, dEband); fflush(globalLog);
		if(dEband<0 and fabs(dEband)<mp.energyDiffThreshold)
		{	logPrintf("BandDavidson: Converged (dEband<%le)\n", mp.energyDiffThreshold);
			break;
		}
	}
	if(iter>mp.nIterations)
		logPrintf("BandDavidson: None of the convergence criteria satisfied after %d iterations.\n", mp.nIterations);
	fflush(globalLog);
	
	//Update final quantities:
	C = (Y.nCols()==nBandsOut ? Y : Y.getSub(0, nBandsOut)); //already orthonormal and in eigenbasis
	e.iInfo.project(C, VdagC);
	eVars.Hsub_eigs[qActive] = Hsub_eigs(0,nBandsOut);
	eVars.Hsub[qActive] = eVars.Hsub_eigs[qActive];
	eVars.Hsub_evecs[qActive] = I;
}
