/*-------------------------------------------------------------------
Copyright 2022 Ravishankar Sundararaman

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

#include <electronic/Dump_internal.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/BandMinimizer.h>
#include <electronic/BandDavidson.h>
#include <core/LatticeUtils.h>


DumpCprime::DumpCprime(double dk, double degeneracyThreshold, double vThreshold, bool realSpaceTruncated)
: dk(dk), degeneracyThreshold(degeneracyThreshold), vThreshold(vThreshold), realSpaceTruncated(realSpaceTruncated)
{
}

void DumpCprime::dump(Everything& e) const
{
	logPrintf("\n---- dC/dk calculation with dk = %lg a0^-1 and degeneracyThreshold = %le Eh ----\n", dk, degeneracyThreshold);
	bool needL = e.dump.count(std::make_pair(DumpFreq_End, DumpL));
	bool needQ = e.dump.count(std::make_pair(DumpFreq_End, DumpQ));
	bool needBerry = e.dump.count(std::make_pair(DumpFreq_End, DumpBerry));
	int nBands = e.eInfo.nBands;
	
	//Compute dC/dk, Berry Curvature Omega_k = i <dC/dk| X |dC/dk>, r*[r,H] moments and hence L or Q moments for each reduced k:
	std::vector<matrix> L(e.eInfo.nStates), Q(e.eInfo.nStates), BerryCurvature(e.eInfo.nStates);
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{	//Compute dC/dk and r*[r,H] moments using dC/dk:
		ColumnBundle Cprime_i[3]; matrix3<matrix> rrH;
		for(int iDir=0; iDir<3; iDir++)
		{	matrix CprimeOC;
			Cprime_i[iDir] = getCprime(e, q, iDir, CprimeOC); //Compute idC/dk (effectively r*C)
			if (needL || needQ){
				matrix CprimeHC = CprimeOC * e.eVars.Hsub_eigs[q]; //since C are eigenvectors
				rrH.set_row(iDir, e.iInfo.rHcommutator(Cprime_i[iDir], e.eVars.C[q], CprimeHC)); //effectively (r*C) ^ p . C
			}
		}
		if (needBerry){
			BerryCurvature[q] = zeroes(nBands, nBands * 3);
			const ColumnBundle Y[3] = { Cprime_i[0], Cprime_i[1], Cprime_i[2] };
			for (int kDir = 0; kDir < 3; kDir++){
				int iDir = (kDir + 1) % 3;
				int jDir = (kDir + 2) % 3;
				// BerryCurvature = i <dC/dk| X |dC/dk> = i <idC/dk| X |idC/dk>
				matrix BerryCurvatureqk = complex(0, 1) * ((Y[iDir] ^ Y[jDir]) - (Y[jDir] ^ Y[iDir]));
				BerryCurvature[q].set(0, nBands, kDir*nBands, (kDir + 1)*e.eInfo.nBands, BerryCurvatureqk);
			}
		}
		//Extract L if needed:
		if(needL)
		{	L[q] = zeroes(nBands, nBands*3);
			for(int kDir=0; kDir<3; kDir++)
			{	int iDir = (kDir + 1) % 3;
				int jDir = (kDir + 2) % 3;
				matrix Lqk = complex(0, -1) * (rrH(iDir, jDir) - rrH(jDir, iDir));
				L[q].set(0, nBands, kDir*nBands, (kDir+1)*e.eInfo.nBands, dagger_symmetrize(Lqk));
			}
		}
		//Extract Q if needed:
		if(needQ)
		{	Q[q] = zeroes(nBands, nBands*5);
			//xy, yz and zx components:
			for(int iDir=0; iDir<3; iDir++)
			{	int jDir = (iDir + 1) % 3;
				matrix Qqij = complex(0, -1) * (rrH(iDir, jDir) + rrH(jDir, iDir)); //Includes + (i <-> j)
				Q[q].set(0, nBands, iDir*nBands, (iDir+1)*e.eInfo.nBands, dagger_symmetrize(Qqij));
			}
			//xx - r^2/3 and yy - r^2/3 components:
			for(int iDir=0; iDir<2; iDir++)
			{	matrix Qqiir = complex(0, -2) * (rrH(iDir, iDir) - (1./3)*trace(rrH));  //factor of 2 from the + (i <-> j)
				Q[q].set(0, nBands, (iDir+3)*nBands, (iDir+4)*e.eInfo.nBands, dagger_symmetrize(Qqiir));
			}
		}
	}
	logPrintf("\n---- dC/dk complete ----\n\n"); logFlush();

	//Output L if requested:
	if (needBerry)
	{	string fname = e.dump.getFilename("BerryCurvature");
		logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
		e.eInfo.write(BerryCurvature, fname.c_str(), nBands, nBands * 3);
		logPrintf("done.\n");
	}

	if(needL)
	{	string fname = e.dump.getFilename("L");
		logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
		e.eInfo.write(L, fname.c_str(), nBands, nBands * 3);
		logPrintf("done.\n");
	}

	//Output Q if requested:
	if(needQ)
	{	string fname = e.dump.getFilename("Q");
		logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
		e.eInfo.write(Q, fname.c_str(), nBands, nBands*5);
		logPrintf("done.\n");
	}
}


ColumnBundle DumpCprime::getCprime(Everything& e, int q, int iDir, matrix& CprimeOC) const
{	//Check if direction is non-periodic:
	vector3<> ei; ei[iDir] = -1.; //unit vector along iDir (sign to compensate V = -E.r for Efield potential)
	std::swap(ei, e.coulombParams.Efield); //use Efield mechanism in Coulomb to check/create ramp potential
	vector3<> RT_ei_ramp, RT_ei_wave;
	e.coulombParams.splitEfield(e.gInfo.R, RT_ei_ramp, RT_ei_wave);
	if(realSpaceTruncated and (RT_ei_wave.length_squared() < symmThresholdSq))
	{	//No wave => all ramp => non-periodic: use ramp to compute ri * C:
		logPrintf("\n---- Multiplying %c for quantum number: ", "xyz"[iDir]);
		e.eInfo.kpointPrint(globalLog, q, true); logPrintf(" ----\n");
		ScalarFieldArray ri(1);
		ri[0] = (1./e.gInfo.nr) * e.coulomb->getEfieldPotential(); //scale factor converts Idag to J
		std::swap(ei, e.coulombParams.Efield); //restore Efield
		ColumnBundle Cprime = Idag_DiagV_I(e.eVars.C[q], ri);
		CprimeOC = Cprime ^ O(e.eVars.C[q]);
		return Cprime;
	}
	else std::swap(ei, e.coulombParams.Efield); //restore Efield
	//Handle periodic case via dC/dk:
	vector3<> dkVec = inv(e.gInfo.GT).column(iDir) * dk; //fractional k perturbation
	matrix Vi = complex(0,-1) * e.iInfo.rHcommutator(e.eVars.C[q], iDir, e.eVars.Hsub_eigs[q]); //velocity matrix
	matrix CpOC, CmOC;
	ColumnBundle Cp = getCpert(e, q, dkVec, dk * Vi, CpOC);
	ColumnBundle Cm = getCpert(e, q, -dkVec, -dk * Vi, CmOC);
	complex prefac(0., 0.5/dk); //prefactor for i d/dk in central difference formula
	CprimeOC = prefac.conj() * (CpOC - CmOC);
	return prefac * (Cp - Cm);
}


ColumnBundle DumpCprime::getCpert(Everything& e, int q, vector3<> dkVec, const matrix& dkDotV, matrix& CpertOC) const
{
	//Perform a fixed-H calculation at this q with perturbed k:
	ColumnBundle Cq = e.eVars.C[q]; //backup orginal wavefunction at this q
	std::vector<matrix> VdagCq = e.eVars.VdagC[q]; //corresponding projections
	QuantumNumber qnum = e.eInfo.qnums[q];
	qnum.k += dkVec; //perturbed quantum number
	std::swap(qnum, e.eInfo.qnums[q]); //change quantum number in e, but keep same basis
	e.iInfo.project(e.eVars.C[q], e.eVars.VdagC[q]);
	//--- Simplified from bandMinimize() for a single q and with no EXX loop
	{	bool fixed_H = true; std::swap(fixed_H, e.cntrl.fixed_H); //remember fixed_H flag and temporarily set it to true
		matrix Hsub_evecs_q; std::swap(Hsub_evecs_q, e.eVars.Hsub_evecs[q]);
		diagMatrix Hsub_eigs_q; std::swap(Hsub_eigs_q, e.eVars.Hsub_eigs[q]);
		logPrintf("\n---- Minimization of quantum number: "); e.eInfo.kpointPrint(globalLog, q, true); logPrintf(" ----\n");
		switch(e.cntrl.elecEigenAlgo)
		{	case ElecEigenCG: 
			{	BandMinimizer(e, q).minimize(e.elecMinParams);
				//Set eigenvectors:
				e.eVars.C[q] = e.eVars.C[q] * e.eVars.Hsub_evecs[q];
				for(matrix& VdagCq_sp: e.eVars.VdagC[q])
					if(VdagCq_sp) VdagCq_sp = VdagCq_sp * e.eVars.Hsub_evecs[q];
				e.eVars.Hsub[q] = e.eVars.Hsub_eigs[q]; //now diagonal
				e.eVars.Hsub_evecs[q] = eye(e.eInfo.nBands);
				break;
			}
			case ElecEigenDavidson: { BandDavidson(e, q).minimize(false); break; }
		}
		std::swap(fixed_H, e.cntrl.fixed_H); //restore fixed_H flag
		std::swap(Hsub_evecs_q, e.eVars.Hsub_evecs[q]);
		std::swap(Hsub_eigs_q, e.eVars.Hsub_eigs[q]);
	}
	std::swap(qnum, e.eInfo.qnums[q]); //restore quantum number
	std::swap(Cq, e.eVars.C[q]); //restore wavefunction (Cq is now perturbed version)
	std::swap(VdagCq, e.eVars.VdagC[q]); //restore corresponding projections
	
	//Align phases / unitary rotations:
	CpertOC = Cq ^ O(e.eVars.C[q]);
	matrix U = fixUnitary(CpertOC, e.eVars.Hsub_eigs[q], dkDotV);
	CpertOC = dagger(U) * CpertOC;
	return Cq * U;
}


//Fix phase / unitary rotations within degenerate subspaces of E,
//given overlap matrix CpertOC between the perturbed and original
//Bloch functions, and the perturbation Hamiltonian dH.
//Output is a unitary matrix, which will be block diagonal in degenerate subspaces.
matrix DumpCprime::fixUnitary(const matrix& CpertOC, const diagMatrix& E, const matrix& dH) const
{	int N = E.nRows();
	assert(CpertOC.nRows() == N);
	assert(dH.nRows() == N);
	matrix U = zeroes(N, N);
	double EpertPrev = NAN; //last (largest) perturbed energy of previous degenerate subspace
	for(int bStart = 0; bStart < N;)
	{	int bStop = bStart;
		while(bStop < N and (E[bStop] < E[bStart] + degeneracyThreshold))
			bStop++;
		int g = bStop - bStart; //degenerate subspace size
		
		//Resolve degeneracy (partly) by perturbation:
		matrix dHsub = dH(bStart, bStop, bStart, bStop);
		diagMatrix dE; matrix dHevecs;
		dHsub.diagonalize(dHevecs, dE);
		matrix Osub = CpertOC(bStart, bStop, bStart, bStop) * dHevecs; //Cpert^O(C) with degeneracies resolved by dH
		
		//Best align wavefunctions within remaining degenerate sub-subspaces:
		matrix Usub = zeroes(g, g);
		for(int aStart=0; aStart < g;)
		{	int aStop = aStart;
			while(aStop < g and (dE[aStop] < dE[aStart] + vThreshold*dk))
				aStop++;
			matrix OsubSub = Osub(aStart, aStop, aStart, aStop);
			OsubSub = OsubSub * invsqrt(dagger(OsubSub) * OsubSub); //make exactly unitary
			Usub.set(aStart, aStop, aStart, aStop, OsubSub);
			aStart = aStop;
		}
		Usub = Usub * dagger(dHevecs); //rotate back into original linear combinations of C
		U.set(bStart, bStop, bStart, bStop, Usub);
		
		//Check for overlapping subspaces:
		if(bStart)
		{	double EpertMin = E[bStart] + dE[0]; //lowest energy of current subspace
			if(EpertMin < EpertPrev + degeneracyThreshold)
				logPrintf("WARNING: degenerate subspaces overlap upon perturbation near E = %lg\n", EpertMin);
		}
		EpertPrev = E[bStop-1] + dE[g-1];
		
		bStart = bStop;
	}
	return U;
}
