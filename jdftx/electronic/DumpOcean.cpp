/*-------------------------------------------------------------------
Copyright 2015 Kathleen Schwarz

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

#include <electronic/Dump.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/matrix.h>
#include <electronic/operators.h>
#include <electronic/ElecMinimizer.h>
#include <electronic/ColumnBundleTransform.h>
#include <core/Thread.h>
#include <core/Operators.h>
#include <core/LatticeUtils.h>
#include <config.h>
#include <map>

void Dump::dumpOcean()
{
	const ElecInfo &eInfo = e->eInfo;
	const ElecVars &eVars = e->eVars;

	if(e->eInfo.isNoncollinear())
		die("OCEAN output is not supported for noncollinear spin modes.\n");
	
	int nSpins = (eInfo.spinType==SpinNone ? 1 : 2);
	int nkPoints = eInfo.nStates / nSpins;
	double sqrtvol=sqrt(e->gInfo.detR);

	logPrintf("\nDumping Kohn-Sham orbitals for the OCEAN code:\n");
	for(int ik=0; ik<nkPoints; ik++)
	{
		if(!mpiUtil->isHead()) //Non-head processes only need to ship wavefunctions to head
		{
			for(int s=0; s<nSpins; s++)
			{	int q = ik + nkPoints*s; //net quantum number
				if(eInfo.isMine(q))
					eVars.C[q].send(0, s); //use spin as a tag
			}
		}
		else //All the writing happens from head
		{
			//Open a file for each k-point:			
			std::string fname = "kpoint" + std::to_string(ik);
			logPrintf("\tDumping '%s' ... ", fname.c_str()); logFlush();
			FILE *fp = fopen(fname.c_str(), "wb");
			if(!fp) die_alone("Error opening %s for writing.\n", fname.c_str());
			
			//Initialize k-point transformation from JDFTx to OCEAN convention:
			const Basis& basis = e->basis[ik];
			const QuantumNumber& qnum = eInfo.qnums[ik];
			//--- quantum number
			QuantumNumber qnumOut = qnum;
			logPrintf("at k = [ ");
			for(int j=0; j<3; j++)
			{	qnumOut.k[j] -= floor(qnumOut.k[j]); //wrap to [0,1)
				if(fabs(qnumOut.k[j] - 1.) < symmThreshold)
					qnumOut.k[j] -= 1.; //make sure components ~ 0.9999 wrap to 0
				logPrintf("%+.6f ", qnumOut.k[j]);
			}
			logPrintf("] "); logFlush();
			//--- basis
			Basis basisOut;
			{	const GridInfo& gInfoBasis = e->gInfoWfns ? *(e->gInfoWfns) : e->gInfo;
				const vector3<> kbasis = (e->cntrl.basisKdep==BasisKpointDep) ? qnumOut.k : vector3<>();
				logSuspend();
				basisOut.setup(gInfoBasis, e->iInfo, e->cntrl.Ecut, kbasis);
				logResume();
			}
			//--- transformation indices
			ColumnBundleTransform transform(qnum.k, basis, qnumOut.k, basisOut, 1, matrix3<int>(1,1,1), +1);
			
			//Write header:
			int nbasis = basisOut.nbasis; //cast nbasis to int
			fwrite(&nbasis, sizeof(int), 1, fp); //Number of G-vectors
			fwrite(basisOut.iGarr, sizeof(int), 3*nbasis, fp); //List of G-vectors
			
			//Collect relevant wavefunctions:
			std::vector<ColumnBundle> Ck(nSpins);
			for(int s=0; s<nSpins; s++)
			{	int q = ik + nkPoints*s; //net quantum number
				//Get wavefunction at JDFTx k-point:
				ColumnBundle CkTemp;
				const ColumnBundle* CkIn;
				if(eInfo.isMine(q))
					CkIn = &eVars.C[q];
				else
				{	CkIn = &CkTemp;
					CkTemp.init(eInfo.nBands, nbasis, &basis, &qnum);
					CkTemp.recv(eInfo.whose(q), s);
				}
				//Convert to OCEAN k-point (and add normalization factor):
				Ck[s].init(eInfo.nBands, nbasis, &basisOut, &qnumOut);
				Ck[s].zero();
				transform.scatterAxpy(sqrtvol, *CkIn, Ck[s], 0, 1);
			}
			
			//Write wavefunctions:
			std::vector<complex> phases(2);
			phases[0] = complex(+1, 0); //get real part below
			phases[1] = complex(0, -1); //get imag part below
			for(complex phase: phases) //loop over real/imag part
				for(const ColumnBundle& Cq: Ck) //loop over spin
					(phase * Cq).write_real(fp); //outer loop over bands, inner loop over G-vectors
			
			fclose(fp);
			logPrintf("done.\n"); logFlush();
		}
	}
}
