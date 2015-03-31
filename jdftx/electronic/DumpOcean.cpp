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
			std::string fname = "kpoint" + std::to_string (ik);
			logPrintf("\tDumping '%s' ... ", fname.c_str()); logFlush();
			FILE *fp = fopen(fname.c_str(), "wb");
			if(!fp) die_alone("Error opening %s for writing.\n", fname.c_str());
			
			//Write header:
			int nbasis = e->basis[ik].nbasis; //cast nbasis to int
			fwrite(&nbasis, sizeof(int), 1, fp); //Number of G-vectors
			fwrite(e->basis[ik].iGarr, sizeof(int), 3*nbasis, fp); //List of G-vectors
			
			//Collect relevant wavefunctions:
			std::vector<ColumnBundle> CkTemp(nSpins);
			std::vector<const ColumnBundle*> Ck(nSpins);
			for(int s=0; s<nSpins; s++)
			{	int q = ik + nkPoints*s; //net quantum number
				if(eInfo.isMine(q))
					Ck[s] = &eVars.C[q];
				else
				{	Ck[s] = &CkTemp[s];
					CkTemp[s].init(eInfo.nBands, e->basis[q].nbasis, &e->basis[q], &eInfo.qnums[q]);
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
			logPrintf("done.\n"); logFlush();
		}
	}
	
}
