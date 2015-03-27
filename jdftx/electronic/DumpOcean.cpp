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
	const GridInfo& gInfo = e->gInfo;

	if(e->eInfo.isNoncollinear())
		die("OCEAN output may not be supported for noncollinear spin modes.\n");
	
	std::string fname; ofstream ofs;

	int nSpins = (eInfo.spinType==SpinNone ? 1 : 2);
	int nkPoints = eInfo.nStates / nSpins;

	for(int ik=0; ik<nkPoints; ik++)
	  {
	    
	    vector3<> k(eInfo.qnums[ik].k * gInfo.G); //cartesian k-vector
	    
	    //make a string for the output file with ik in it:
	    std::string kfile = "kpoint" + std::to_string (ik);
	    
	    
	    if(!mpiUtil->isHead()) fname = "/dev/null";
	    else 
	      {  fname=kfile;
		
		//open the file just to delete old contents
		FILE *fp2 = fopen(fname.c_str(),"wb");
		if(!fp2) die("Error opening %s for writing.\n", fname.c_str());
		fclose(fp2);
	      }

	    FILE *fp = fopen(fname.c_str(),"ab");
	    if(!fp) die("Error opening %s for writing.\n", fname.c_str());
	    
	    
	    
	    for(int s=0; s<nSpins; s++)
	      {
		int q = ik + nkPoints*s; //net quantum number
		
		//Get relevant wavefunctions 
		ColumnBundle CqTemp; 
		const ColumnBundle* Cq=0; 
		if(mpiUtil->isHead())
		  {	if(eInfo.isMine(q))
		      {	Cq = &eVars.C[q];
		      }
		    else
		      {	Cq = &CqTemp;
			CqTemp.init(eInfo.nBands, e->basis[q].nbasis, &e->basis[q], &eInfo.qnums[q]);
			CqTemp.recv(eInfo.whose(q));
		      }
		  }
		else
		  {	if(eInfo.isMine(q))
		      {	eVars.C[q].send(0);
		      }
		    continue; //only head performs computation below
		  }
		const Basis& basis = e->basis[q];
		fwrite(&basis.nbasis,sizeof(int),1, fp);
		
		for(size_t i=0; i<basis.nbasis; i++)
		  fwrite(&basis.iGarr[i],sizeof(int),3, fp);
		
		Cq->write_real(fp);
		((*Cq) * complex(0,1)).write_real(fp);
		fclose(fp);
	      }
	  }
	logPrintf("done.\n"); logFlush();
}
