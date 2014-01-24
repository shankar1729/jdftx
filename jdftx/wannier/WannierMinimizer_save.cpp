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

#include <wannier/WannierMinimizer.h>
#include <electronic/operators.h>
#include <core/LatticeUtils.h>

void WannierMinimizer::saveMLWF()
{	for(int iSpin=0; iSpin<nSpins; iSpin++)
		saveMLWF(iSpin);
}

void WannierMinimizer::saveMLWF(int iSpin)
{	
	//Check for initial state:
	string fname = wannier.getFilename(true, "mlwfU", &iSpin);
	bool readInitialMats = false;
	FILE* fp = fopen(fname.c_str(), "r");
	if(fp)
	{	logPrintf("Reading initial matrices from %s (ignoring trial projections).\n", fname.c_str());
		for(size_t ik=0; ik<kMesh.size(); ik++)
		{	kMesh[ik].U0.init(nCenters, nCenters);
			kMesh[ik].U0.read(fp);
		}
		fclose(fp);
		readInitialMats = true;
	}
	
	//Compute the overlap matrices and initial rotations for current group of centers:
	for(int jProcess=0; jProcess<mpiUtil->nProcesses(); jProcess++)
	{	//Send/recv wavefunctions to other processes:
		Cother.assign(e.eInfo.nStates, ColumnBundle());
		if(jProcess == mpiUtil->iProcess()) //send
		{	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
				((ColumnBundle&)e.eVars.C[q]).bcast(jProcess);
		}
		else //recv
		{	for(int q=e.eInfo.qStartOther(jProcess); q<e.eInfo.qStopOther(jProcess); q++)
			{	Cother[q].init(nBands, e.basis[q].nbasis, &e.basis[q], &e.eInfo.qnums[q]);
				Cother[q].bcast(jProcess);
			}
		}
		
		for(size_t ik=0; ik<kMesh.size(); ik++) if(isMine_q(ik,iSpin))
		{	ColumnBundle Ci = getWfns(kMesh[ik].point, iSpin); //Bloch functions at ik
			ColumnBundle OCi = O(Ci);
			//Overlap with neighbours:
			for(EdgeFD& edge: kMesh[ik].edge)
				if(whose_q(edge.ik,iSpin)==jProcess)
				{	//Pick up result from reverse edge if it has already been computed:
					bool foundReverse = false;
					if(jProcess==mpiUtil->iProcess() && edge.ik<ik)
					{	auto neighbourEdges = kMesh[edge.ik].edge;
						for(const EdgeFD& reverseEdge: neighbourEdges)
							if(reverseEdge.ik==ik)
							{	edge.M0 = dagger(reverseEdge.M0);
								foundReverse = true;
								break;
							}
					}
					//Compute overlap if reverse edge not yet computed:
					if(!foundReverse)
						edge.M0 = OCi ^ getWfns(edge.point, iSpin);
				}
			if(!jProcess) //Do only once (will get here multiple times for local wfns)
			{	//Initial rotation:
				if(!readInitialMats)
				{	matrix CdagOg = OCi ^ trialWfns(kMesh[ik].point);
					kMesh[ik].U0 = CdagOg * invsqrt(dagger(CdagOg) * CdagOg);
				}
			}
		}
	}
	
	//Broadcast overlaps and initial rotations:
	for(size_t ik=0; ik<kMesh.size(); ik++)
	{	for(EdgeFD& edge: kMesh[ik].edge)
		{	if(!isMine_q(ik,iSpin)) edge.M0 = zeroes(nCenters, nCenters);
			edge.M0.bcast(whose_q(ik,iSpin));
			if(!isMine(ik)) edge.M0 = matrix(); //not needed any more on this process
		}
		if(!readInitialMats)
		{	if(!isMine_q(ik,iSpin)) kMesh[ik].U0 = zeroes(nCenters, nCenters);
			kMesh[ik].U0.bcast(whose_q(ik,iSpin));
		}
		kMesh[ik].B = zeroes(nCenters, nCenters);
	}
	
	//Apply initial rotations to the overlap matrices:
	for(size_t ik=ikStart; ik<ikStop; ik++)
		for(EdgeFD& edge: kMesh[ik].edge)
			edge.M0 = dagger(kMesh[ik].U0) * edge.M0 * kMesh[edge.ik].U0;
	
	//Minimize:
	minimize(wannier.minParams);
	
	//List the centers:
	logPrintf("Centers in %s coords:\n", e.iInfo.coordsType==CoordsCartesian ? "cartesian" : "lattice");
	for(int n=0; n<nCenters; n++)
	{	vector3<> rCoords = e.iInfo.coordsType==CoordsCartesian
			? rExpect[n] : e.gInfo.invR * rExpect[n]; //r in coordinate system of choice
		logPrintf("\t[ %lg %lg %lg ] spread: %lg bohrs\n", rCoords[0], rCoords[1], rCoords[2], sqrt(rSqExpect[n] - rExpect[n].length_squared()));
	}
	logFlush();
	
	//Save the matrices:
	fname = wannier.getFilename(false, "mlwfU", &iSpin);
	logPrintf("Dumping '%s' ... ", fname.c_str());
	if(mpiUtil->isHead())
	{	FILE* fp = fopen(fname.c_str(), "w");
		for(const auto& kMeshEntry: kMesh)
			(kMeshEntry.U0 * kMeshEntry.V).write(fp);
		fclose(fp);
	}
	logPrintf("done.\n"); logFlush();

	//Compute supercell wavefunctions:
	logPrintf("Computing supercell wavefunctions ... "); logFlush();
	ColumnBundle Csuper(nCenters, basisSuper.nbasis, &basisSuper, &qnumSuper, isGpuEnabled());
	Csuper.zero();
	for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
		Csuper += getWfns(kMesh[i].point, iSpin, true) * (kMesh[i].U0 * kMesh[i].V * wk);
	Csuper.allReduce(MPIUtil::ReduceSum);
	Csuper = translate(Csuper, vector3<>(.5,.5,.5)); //center in supercell
	logPrintf("done.\n"); logFlush();
	
	//Save supercell wavefunctions:
	for(int n=0; n<nCenters; n++)
	{	//Generate filename
		ostringstream varName;
		varName << n << ".mlwf";
		string fname = wannier.getFilename(false, varName.str(), &iSpin);
		logPrintf("Dumping '%s':\n", fname.c_str());
		//Convert to real space and remove phase:
		complexDataRptr psi = I(Csuper.getColumn(n));
		if(qnumSuper.k.length_squared() > symmThresholdSq)
			multiplyBlochPhase(psi, qnumSuper.k);
		complex* psiData = psi->data();
		double meanPhase, sigmaPhase, rmsImagErr;
		removePhase(gInfoSuper.nr, psiData, meanPhase, sigmaPhase, rmsImagErr);
		logPrintf("\tPhase = %lf +/- %lf\n", meanPhase, sigmaPhase); logFlush();
		logPrintf("\tRMS imaginary part = %le (after phase removal)\n", rmsImagErr);
		logFlush();
		//Write real part of supercell wavefunction to file:
		if(mpiUtil->isHead())
		{	FILE* fp = fopen(fname.c_str(), "wb");
			if(!fp) die("Failed to open file '%s' for binary write.\n", fname.c_str());
			for(int i=0; i<gInfoSuper.nr; i++)
				fwrite(psiData+i, sizeof(double), 1, fp);
			fclose(fp);
		}
	}

	return; //TODO: Add Hsub read-in so that the following works
	//Save Hamiltonian in Wannier basis:
	const vector3<int>& nSuper = wannier.supercell;
	vector3<int> offsetSuper(-nSuper[0]/2, -nSuper[1]/2, -nSuper[2]/2); //supercell offset
	std::vector<matrix> Hwannier(nSuper[0]*nSuper[1]*nSuper[2]);
	for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
	{	//Fetch Hamiltonian for subset of bands in center:
		matrix Hsub(nCenters, nCenters);
		{	complex* HsubData = Hsub.data();
			const matrix& HsubFull = e.eVars.Hsub[kMesh[i].point.q + iSpin*qCount];
			const complex* HsubFullData = HsubFull.data();
			for(int c1=0; c1<nCenters; c1++)
				for(int c2=0; c2<nCenters; c2++)
					HsubData[Hsub.index(c1,c2)] = HsubFullData[HsubFull.index(wannier.bStart+c1,wannier.bStart+c2)];
		}
		//Apply MLWF-optimizd rotation:
		matrix U = kMesh[i].U0 * kMesh[i].V;
		Hsub = dagger(U) * Hsub * U;
		//Accumulate with each requested Bloch phase
		vector3<int> iSuper;
		std::vector<matrix>::iterator HwannierIter = Hwannier.begin();
		for(iSuper[0]=0; iSuper[0]<nSuper[0]; iSuper[0]++)
		for(iSuper[1]=0; iSuper[1]<nSuper[1]; iSuper[1]++)
		for(iSuper[2]=0; iSuper[2]<nSuper[2]; iSuper[2]++)
			*(HwannierIter++) += (wk * cis(2*M_PI*dot(kMesh[i].point.k, offsetSuper+iSuper))) * Hsub;
	}
	for(matrix& H: Hwannier) H.allReduce(MPIUtil::ReduceSum);
	//-- save to file
	fname = wannier.getFilename(false, "mlwfH", &iSpin);
	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
	if(mpiUtil->isHead())
	{	FILE* fp = fopen(fname.c_str(), "wb");
		if(!fp) die("Failed to open file '%s' for binary write.\n", fname.c_str());
		for(matrix& H: Hwannier) H.write_real(fp);
		fclose(fp);
	}
	logPrintf("done.\n"); logFlush();
}
