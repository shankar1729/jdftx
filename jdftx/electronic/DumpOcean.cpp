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
		//Output filename for each k-point:
		std::string fname = "kpoint" + std::to_string(ik);
		logPrintf("\tDumping '%s' ... ", fname.c_str()); logFlush();
		
		//Initialize k-point transformation from JDFTx to OCEAN convention:
		const Basis& basis = e->basis[ik];
		const QuantumNumber& qnum = eInfo.qnums[ik];
		vector3<> kOut = qnum.k;
		logPrintf("at k = [ ");
		for(int j=0; j<3; j++)
		{	//Wrap to [0,1)
			kOut[j] -= floor(kOut[j]);
			if(fabs(kOut[j] - 1.) < symmThreshold)
				kOut[j] -= 1.; //make sure components ~ 0.9999 wrap to 0
			logPrintf("%+.6f ", kOut[j]);
		}
		logPrintf("] "); logFlush();
		vector3<int> kOffset = round(qnum.k - kOut);
		
		if(eInfo.isMine(ik)) //All the writing happens from the process that owns the first spin channel
		{
			//Open the file:
			FILE *fp = fopen(fname.c_str(), "wb");
			if(!fp) die_alone("Error opening %s for writing.\n", fname.c_str());
			
			//Write header:
			int nbasis = basis.nbasis; //cast nbasis to int
			fwrite(&nbasis, sizeof(int), 1, fp); //Number of G-vectors
			for(int n=0; n<nbasis; n++)
			{	vector3<int> iGout = basis.iGarr[n] + kOffset;
				fwrite(&iGout, sizeof(int), 3, fp); //List of G-vectors
			}
			
			//Collect relevant wavefunctions:
			std::vector<ColumnBundle> CkTemp(nSpins);
			std::vector<const ColumnBundle*> Ck(nSpins);
			for(int s=0; s<nSpins; s++)
			{	int q = ik + nkPoints*s; //net quantum number
				if(eInfo.isMine(q))
					Ck[s] = &eVars.C[q];
				else
				{	Ck[s] = &CkTemp[s];
					CkTemp[s].init(eInfo.nBands, nbasis, &basis, &qnum);
					CkTemp[s].recv(eInfo.whose(q), s);
				}
			}
			
			//Write wavefunctions:
			std::vector<complex> phases(2);
			phases[0] = complex(+sqrtvol, 0); //normalize and get real part below
			phases[1] = complex(0, -sqrtvol); //normalize and get imag part below
			for(complex phase: phases) //loop over real/imag part
				for(const ColumnBundle* Cq: Ck) //loop over spin
					(phase * (*Cq)).write_real(fp); //outer loop over bands, inner loop over G-vectors
			
			fclose(fp);
			
		}
		else //Send wavefunctions to the owner of the first spin channel at this k
		{
			for(int s=0; s<nSpins; s++)
			{	int q = ik + nkPoints*s; //net quantum number
				if(eInfo.isMine(q))
					eVars.C[q].send(eInfo.whose(ik), s); //use spin as a tag
			}
		}
		logPrintf("done.\n"); logFlush();
	}
}
