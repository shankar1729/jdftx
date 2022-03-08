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


DumpCprime::DumpCprime(double dk, double degeneracyThreshold)
: dk(dk), degeneracyThreshold(degeneracyThreshold)
{
}

void DumpCprime::dump(Everything& e) const
{
	logPrintf("\n---- dC/dk calculation with dk = %lg a0^-1 and degeneracyThreshold = %le Eh ----\n", dk, degeneracyThreshold);
	bool needL = e.dump.count(std::make_pair(DumpFreq_End, DumpL));
	bool needQ = e.dump.count(std::make_pair(DumpFreq_End, DumpQ));
	int nBands = e.eInfo.nBands;
	
	//Compute dC/dk, r*[r,H] moments and hence L or Q moments for each reduced k:
	std::vector<matrix> L(e.eInfo.nStates), Q(e.eInfo.nStates);
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{	//Compute r*[r,H] moments using dC/dk:
		matrix3<matrix> rrH;
		for(int iDir=0; iDir<3; iDir++)
		{	matrix CprimeOC;
			ColumnBundle Cprime_i = getCprime(e, q, iDir, CprimeOC); //Compute idC/dk (effectively r*C)
			matrix CprimeHC = CprimeOC * e.eVars.Hsub_eigs[q]; //since C are eigenvectors
			rrH.set_row(iDir, e.iInfo.rHcommutator(Cprime_i, e.eVars.C[q], CprimeHC)); //effectively (r*C) ^ p . C
		}
		//Extract L if needed:
		if(needL)
		{	L[q] = zeroes(nBands, nBands*3);
			for(int kDir=0; kDir<3; kDir++)
			{	int iDir = (kDir + 1) % 3;
				int jDir = (kDir + 2) % 3;
				L[q].set(0, nBands, kDir*nBands, (kDir+1)*e.eInfo.nBands,
					complex(0, -1) * (rrH(iDir, jDir) - rrH(jDir, iDir)));
			}
		}
		//Extract Q if needed:
		if(needQ)
		{	Q[q] = zeroes(nBands, nBands*5);
			//xy, yz and zx components:
			for(int iDir=0; iDir<3; iDir++)
			{	int jDir = (iDir + 1) % 3;
				Q[q].set(0, nBands, iDir*nBands, (iDir+1)*e.eInfo.nBands,
					complex(0, -1) * (rrH(iDir, jDir) - dagger(rrH(jDir, iDir)))); //ri pj + pi rj
			}
			//xx - r^2/3 and yy - r^2/3 components:
			matrix traceTerm = complex(0., -1./3) * trace(rrH); //sum(ri pi)/3
			traceTerm += dagger(traceTerm); //+= sum(pi ri)/3
			for(int iDir=0; iDir<2; iDir++)
				Q[q].set(0, nBands, (iDir+3)*nBands, (iDir+4)*e.eInfo.nBands,
					complex(0, -1) * (rrH(iDir, iDir) - dagger(rrH(iDir, iDir))) //ri pi + pi ri
					- traceTerm);
		}
	}
	logPrintf("\n---- dC/dk complete ----\n\n"); logFlush();

	//Output L if requested:
	if(needL)
	{	string fname = e.dump.getFilename("L");
		logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
		e.eInfo.write(L, fname.c_str(), nBands, nBands*3);
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
{	vector3<> dkVec = inv(e.gInfo.GT).column(iDir) * dk; //fractional k perturbation
	matrix CpOC, CmOC;
	ColumnBundle Cp = getCpert(e, q, dkVec, CpOC);
	ColumnBundle Cm = getCpert(e, q, -dkVec, CmOC);
	complex prefac(0., 0.5/dk); //prefactor for i d/dk in central difference formula
	CprimeOC = prefac.conj() * (CpOC - CmOC);
	return prefac * (Cp - Cm);
}


ColumnBundle DumpCprime::getCpert(Everything& e, int q, vector3<> dkVec, matrix& CpertOC) const
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
	matrix Osub = Cq ^ O(e.eVars.C[q]);
	matrix U = fixUnitary(Osub, e.eVars.Hsub_eigs[q]);
	CpertOC = dagger(U) * Osub;
	return Cq * U;
}


//Fix phase / unitary rotations within degenerate subspaces of E, given overlap matrix O.
//Output is a unitary matrix, which will be block diagonal in degenrate subspaces.
matrix DumpCprime::fixUnitary(const matrix& O, const diagMatrix& E) const
{	int N = E.nRows();
	assert(O.nRows() == N);
	matrix U = zeroes(N, N);
	for(int bStart = 0; bStart < N-1;)
	{	int bStop = bStart;
		while(bStop < N and (E[bStop] < E[bStart] + degeneracyThreshold))
			bStop++;
		if(bStop - bStart > 1) //degeneracy
		{	matrix Osub = O(bStart, bStop, bStart, bStop);
			U.set(bStart, bStop, bStart, bStop, Osub * invsqrt(dagger(Osub) * Osub));
		}
		else //separated band
		{	complex Obb = O(bStart, bStart);
			U.set(bStart, bStart,  Obb / Obb.abs());
		}
		bStart = bStop;
	}
	return U;
}
