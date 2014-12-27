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
#include <core/DataIO.h>

void WannierMinimizer::saveMLWF()
{	for(int iSpin=0; iSpin<nSpins; iSpin++)
		saveMLWF(iSpin);
}

void WannierMinimizer::saveMLWF(int iSpin)
{	
	//Load initial rotations if necessary:
	if(wannier.loadRotations)
	{	//Read U
		string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfU", &iSpin);
		FILE* fp = fopen(fname.c_str(), "r");
		if(!fp) die("Could not open '%s' for reading.\n", fname.c_str());
		logPrintf("Reading initial rotations from '%s' ... ", fname.c_str());
		for(auto& ke: kMesh)
		{	ke.U.init(nBands, nCenters);
			ke.U.read(fp);
		}
		fclose(fp);
		logPrintf("done.\n"); logFlush();
		//Read U2:
		fname = wannier.getFilename(Wannier::FilenameDump, "mlwfU2", &iSpin);
		fp = fopen(fname.c_str(), "r");
		if(!fp) die("Could not open '%s' for reading.\n", fname.c_str());
		logPrintf("Reading initial outer rotations from '%s' ... ", fname.c_str());
		for(auto& ke: kMesh)
		{	ke.U2.init(nCenters, nCenters);
			ke.U2.read(fp);
		}
		fclose(fp);
		logPrintf("done.\n"); logFlush();
		logPrintf("NOTE: ignoring trial orbitals since we are resuming a previous WannierMinimize.\n");
	}
	
	//Compute the initial rotations for current group of centers:
	for(size_t ik=0; ik<kMesh.size(); ik++) if(isMine_q(ik,iSpin))
	{	KmeshEntry& ke = kMesh[ik];
		
		//Band ranges:
		int bStart=0, bStop=0, bFixedStart=0, bFixedStop=0;
		if(wannier.outerWindow)
		{	const std::vector<double>& eigs = e.eVars.Hsub_eigs[ke.point.iReduced + iSpin*qCount];
			bStart = 0;
			while(bStart<nBands && eigs[bStart]<wannier.eOuterMin)
				bStart++;
			bStop = bStart;
			while(bStop<nBands && eigs[bStop]<=wannier.eOuterMax)
				bStop++;
			if(bStop-bStart < nCenters)
				die("Number of bands within outer window = %d less than nCenters = %d at k = [ %lg %lg %lg ]\n",
					bStop-bStart, nCenters, ke.point.k[0], ke.point.k[1], ke.point.k[2]);
			//Optionally range for inner window:
			if(wannier.innerWindow)
			{	bFixedStart = bStart;
				while(bFixedStart<bStop && eigs[bFixedStart]<wannier.eInnerMin)
					bFixedStart++;
				bFixedStop = bFixedStart;
				while(bFixedStop<bStop && eigs[bFixedStop]<=wannier.eInnerMax)
					bFixedStop++;
				if(bFixedStop-bFixedStart > nCenters)
					die("Number of bands within inner window = %d exceeds nCenters = %d at k = [ %lg %lg %lg ]\n",
						bFixedStop-bFixedStart, nCenters, ke.point.k[0], ke.point.k[1], ke.point.k[2]);
			}
			else bFixedStart = bFixedStop = bStart; //fixed interval is empty
		}
		else //fixed bands
		{	bFixedStart = bStart = wannier.bStart;
			bFixedStop  = bStop  = wannier.bStart + nCenters;
		}
		
		//Initial rotation of bands to get to Wannier subspace:
		ke.nFixed = bFixedStop - bFixedStart;
		int nFree = nCenters - ke.nFixed;
		ke.nIn = (nFree>0) ? (bStop-bStart) : nCenters; //number of bands contributing to Wannier subspace	
		if(wannier.loadRotations)
		{	//Factorize U (nBands x nCenters) into U1 (nBands x nIn) and U2 (nCenters x nCenters):
			//--- check unitarity:
			const double tol = 1e-6 * nCenters;
			if(nrm2(dagger(ke.U) * ke.U - eye(nCenters)) > tol) die("Initial matrices U are not unitary.\n");
			if(nrm2(dagger(ke.U2) * ke.U2 - eye(nCenters)) > tol) die("Initial matrices U2 are not unitary.\n");
			//--- compute and check U1:
			ke.U1 = zeroes(nBands, ke.nIn);
			ke.U1.set(0,nBands, 0,nCenters, ke.U * dagger(ke.U2));
			if( (bStart>0 && nrm2(ke.U1(0,bStart, 0,nCenters))>tol) 
				|| (bStop>nBands && nrm2(ke.U1(bStop,nBands, 0,nCenters))>tol) )
				die("Initial matrices are incompatible with current outer window / band selection.\n");
			if( ke.nFixed>0
				&& ( (nrm2(ke.U1(bFixedStart,bFixedStop, 0,ke.nFixed) - eye(ke.nFixed))>tol)
				|| (bStart<bFixedStart && nrm2(ke.U1(bStart,bFixedStop, 0,ke.nFixed))>tol)
				|| (bFixedStop<bStop && nrm2(ke.U1(bFixedStop,bStop, 0,ke.nFixed))>tol)
				|| (nFree>0 && nrm2(ke.U1(bFixedStart,bFixedStop, ke.nFixed,nCenters))>tol) ) )
				die("Initial matrices are incompatible with current inner window.\n");
			//--- Compute extra linearly indep columns of U1 (if necessary):
			if(ke.nIn > nCenters)
			{	matrix U, Vdag; diagMatrix S;
				ke.U1.svd(U, S, Vdag);
				ke.U1 = U(0,nBands, 0,ke.nIn) * Vdag;
			}
		}
		else
		{	//Determine from trial orbitals:
			matrix CdagG = getWfns(ke.point, iSpin) ^ trialWfns(ke.point);
			ke.U1 = zeroes(nBands, ke.nIn);
			//--- Pick up fixed bands directly
			if(ke.nFixed > 0)
				ke.U1.set(bFixedStart,bFixedStop, 0,ke.nFixed, eye(ke.nFixed));
			//--- Pick up best linear combination of remaining bands (if any)
			if(nFree > 0)
			{	//Create overlap matrix with contribution from fixed bands projected out:
				matrix CdagGFree;
				if(ke.nFixed > 0)
				{	//SVD the fixed band contribution to the trial space:
					matrix U, Vdag; diagMatrix S;
					CdagG(bFixedStart,bFixedStop, 0,nCenters).svd(U, S, Vdag);
					//Project out the fixed bands (use only the zero singular values)
					CdagGFree = CdagG * dagger(Vdag(ke.nFixed,nCenters, 0,nCenters));
				}
				else CdagGFree = CdagG;
				//Truncate to non-zero rows:
				int nLo = bFixedStart-bStart;
				int nHi = bStop-bFixedStop;
				int nOuter = nLo+nHi;
				matrix CdagGFreeNZ = zeroes(nOuter, nOuter);
				if(nLo>0) CdagGFreeNZ.set(0,nLo, 0,nFree, CdagGFree(bStart,bFixedStart, 0,nFree));
				if(nHi>0) CdagGFreeNZ.set(nLo,nOuter, 0,nFree, CdagGFree(bFixedStop,bStop, 0,nFree));
				//SVD to get best linear combinations first:
				matrix U, Vdag; diagMatrix S;
				CdagGFreeNZ.svd(U, S, Vdag);
				//Convert left space from non-zero to all bands:
				matrix Upad = zeroes(nBands, nOuter);
				if(nLo>0) Upad.set(bStart,bFixedStart, 0,nOuter, U(0,nLo, 0,nOuter));
				if(nHi>0) Upad.set(bFixedStop,bStop, 0,nOuter, U(nLo,nOuter, 0,nOuter));
				//Include this combination in U1:
				ke.U1.set(0,nBands, ke.nFixed,ke.nIn, Upad * Vdag);
			}
			
			//Optimal initial rotation within Wannier subspace:
			matrix WdagG = dagger(ke.U1(0,nBands, 0,nCenters)) * CdagG;
			ke.U2 = WdagG * invsqrt(dagger(WdagG) * WdagG);
		}
	}
	suspendOperatorThreading();
	
	//Broadcast initial rotations:
	for(size_t ik=0; ik<kMesh.size(); ik++)
	{	KmeshEntry& ke = kMesh[ik];
		mpiUtil->bcast(ke.nIn, whose_q(ik,iSpin));
		mpiUtil->bcast(ke.nFixed, whose_q(ik,iSpin));
		if(!isMine_q(ik,iSpin))
		{	ke.U1 = zeroes(nBands, ke.nIn);
			ke.U2 = zeroes(nCenters, nCenters);
		}
		ke.U1.bcast(whose_q(ik,iSpin));
		ke.U2.bcast(whose_q(ik,iSpin));
		if(isMine(ik))
			ke.B = zeroes(nCenters, ke.nIn);
		else //No longer need sub-matrices on this process
		{	ke.U1 = matrix();
			ke.U2 = matrix();
			ke.B = matrix();
		}
		ke.U = zeroes(nBands, nCenters);
	}
	
	//Sub-class specific initialization:
	initialize(iSpin);
	
	//Minimize:
	double Omega = minimize(wannier.minParams);
	double OmegaI = getOmegaI();
	logPrintf("\nOptimum spread:\n\tOmega:  %.15le\n\tOmegaI: %.15le\n", Omega, OmegaI);
	
	//List the centers:
	logPrintf("\nCenters in %s coords:\n", e.iInfo.coordsType==CoordsCartesian ? "cartesian" : "lattice");
	for(int n=0; n<nCenters; n++)
	{	vector3<> rCoords = e.iInfo.coordsType==CoordsCartesian
			? rExpect[n] : e.gInfo.invR * rExpect[n]; //r in coordinate system of choice
		logPrintf("\t[ %lg %lg %lg ] spread: %lg bohr^2\n", rCoords[0], rCoords[1], rCoords[2], rSqExpect[n] - rExpect[n].length_squared());
	}
	logFlush();
	
	//Save the matrices:
	if(mpiUtil->isHead())
	{	//Write U:
		string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfU", &iSpin);
		logPrintf("Dumping '%s' ... ", fname.c_str());
		FILE* fp = fopen(fname.c_str(), "w");
		for(const auto& ke: kMesh) ke.U.write(fp);
		fclose(fp);
		logPrintf("done.\n"); logFlush();
		//Write U2:
		fname = wannier.getFilename(Wannier::FilenameDump, "mlwfU2", &iSpin);
		logPrintf("Dumping '%s' ... ", fname.c_str());
		fp = fopen(fname.c_str(), "w");
		for(unsigned ik=0; ik<kMesh.size(); ik++)
		{	const KmeshEntry& ke = kMesh[ik];
			matrix U2out;
			if(isMine(ik)) U2out = ke.U2 * ke.V2;
			else
			{	U2out = zeroes(nCenters,nCenters);
				U2out.recv(whose(ik)); //recieve from another process (see below)
			}
			U2out.write(fp);
		}
		fclose(fp);
		logPrintf("done.\n"); logFlush();
	}
	else
	{	for(unsigned ik=ikStart; ik<ikStop; ik++)
		{	const KmeshEntry& ke = kMesh[ik];
			(ke.U2 * ke.V2).send(0); //send to head for output (see above)
		}
	}
	
	bool realPartOnly = !e.eInfo.isNoncollinear(); //cannot save only real part in noncollinear calculations
	
	if(wannier.saveWfns || wannier.saveWfnsRealSpace)
	{	resumeOperatorThreading();
		//--- Compute supercell wavefunctions:
		logPrintf("Computing supercell wavefunctions ... "); logFlush();
		ColumnBundle Csuper(nCenters, basisSuper.nbasis*nSpinor, &basisSuper, &qnumSuper, isGpuEnabled());
		Csuper.zero();
		for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
		{	const KmeshEntry& ki = kMesh[i];
			axpyWfns(ki.point.weight, ki.U, ki.point, iSpin, Csuper);
		}
		Csuper.allReduce(MPIUtil::ReduceSum);
		Csuper = translate(Csuper, vector3<>(.5,.5,.5)); //center in supercell
		logPrintf("done.\n"); logFlush();
		
		//--- Save supercell wavefunctions in reciprocal space:
		if(mpiUtil->isHead() && wannier.saveWfns)
		{	string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfC", &iSpin);
			logPrintf("Dumping '%s'... ", fname.c_str()); logFlush();
			Csuper.write(fname.c_str());
			logPrintf("done.\n"); logFlush();
			//Header:
			fname = fname + ".header";
			logPrintf("Dumping '%s'... ", fname.c_str()); logFlush();
			FILE* fp = fopen(fname.c_str(), "w");
			fprintf(fp, "%d %lu #nColumns, columnLength\n", nCenters, basisSuper.nbasis);
			for(int i=0; i<3; i++)
				for(int j=0; j<3; j++)
					fprintf(fp, "%.15g ", gInfoSuper.GT(i,j));
			fprintf(fp, "#GT row-major (G col-major), iGarr follows:\n");
			const vector3<int>* iGarr = basisSuper.iGarr;
			for(size_t j=0; j<basisSuper.nbasis; j++)
				fprintf(fp, "%d %d %d\n", iGarr[j][0], iGarr[j][1], iGarr[j][2]);
			fclose(fp);
			logPrintf("done.\n"); logFlush();
		}
		
		//--- Save supercell wavefunctions in real space:
		if(mpiUtil->isHead() && wannier.saveWfnsRealSpace) for(int n=0; n<nCenters; n++) for(int s=0; s<nSpinor; s++)
		{	//Generate filename
			ostringstream varName;
			varName << (nSpinor*n+s) << ".mlwf";
			string fname = wannier.getFilename(Wannier::FilenameDump, varName.str(), &iSpin);
			logPrintf("Dumping '%s':\n", fname.c_str());
			//Convert to real space and optionally remove phase:
			complexDataRptr psi = I(Csuper.getColumn(n,s));
			if(qnumSuper.k.length_squared() > symmThresholdSq)
				multiplyBlochPhase(psi, qnumSuper.k);
			if(realPartOnly)
			{	complex* psiData = psi->data();
				double meanPhase, sigmaPhase, rmsImagErr;
				removePhase(gInfoSuper.nr, psiData, meanPhase, sigmaPhase, rmsImagErr);
				logPrintf("\tPhase = %lf +/- %lf\n", meanPhase, sigmaPhase); logFlush();
				logPrintf("\tRMS imaginary part = %le (after phase removal)\n", rmsImagErr);
				logFlush();
				//Write real part of supercell wavefunction to file:
				FILE* fp = fopen(fname.c_str(), "wb");
				if(!fp) die("Failed to open file '%s' for binary write.\n", fname.c_str());
				for(int i=0; i<gInfoSuper.nr; i++)
					fwrite(psiData+i, sizeof(double), 1, fp);
				fclose(fp);
			}
			else saveRawBinary(psi, fname.c_str());
		}
		suspendOperatorThreading();
	}
	
	//Save Hamiltonian in Wannier basis:
	std::vector<matrix> Hwannier(iCellMap.size(), zeroes(nCenters,nCenters));
	for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
	{	//Fetch Hamiltonian for subset of bands in center:
		matrix Hsub = e.eVars.Hsub_eigs[kMesh[i].point.iReduced + iSpin*qCount];
		//Apply MLWF-optimized rotation:
		Hsub = dagger(kMesh[i].U) * Hsub * kMesh[i].U;
		//Accumulate with each requested Bloch phase
		std::vector<matrix>::iterator HwannierIter = Hwannier.begin();
		for(auto cell: iCellMap)
			*(HwannierIter++) += (cell.second * kMesh[i].point.weight * cis(2*M_PI*dot(kMesh[i].point.k, cell.first))) * Hsub;
	}
	for(matrix& H: Hwannier) H.allReduce(MPIUtil::ReduceSum);
	//-- save to file
	if(mpiUtil->isHead())
	{	string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfH", &iSpin);
		logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
		FILE* fp = fopen(fname.c_str(), "wb");
		if(!fp) die("Failed to open file '%s' for binary write.\n", fname.c_str());
		double normTot=0., normIm=0.;
		for(matrix& H: Hwannier)
		{	if(realPartOnly)
			{	H.write_real(fp);
				//Collect imaginary part:
				normTot += pow(nrm2(H), 2);
				eblas_dscal(H.nData(), 0., ((double*)H.data()), 2); //zero out real parts
				normIm += pow(nrm2(H), 2);
			}
			else H.write(fp);
		}
		fclose(fp);
		if(realPartOnly) logPrintf("done. Relative discarded imaginary part: %le\n", sqrt(normIm/normTot)); else logPrintf("done.\n"); logFlush();
	}
	resumeOperatorThreading();
	
	//Save momenta in Wannier basis:
	if(wannier.saveMomenta)
	{	//--- compute momentum matrix elements of Bloch states:
		std::vector< std::vector<matrix> > pBloch(3, std::vector<matrix>(e.eInfo.nStates));
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
			if(e.eInfo.qnums[q].spin==iSpin)
				for(int iDir=0; iDir<3; iDir++)
					pBloch[iDir][q] = e.gInfo.detR * (e.eVars.C[q] ^ D(e.eVars.C[q], iDir)); //note factor of iota dropped to make it real (and anti-symmetric)
		//--- convert to Wannier basis:
		std::vector< std::vector<matrix> > pWannier(3, std::vector<matrix>(iCellMap.size(), zeroes(nCenters,nCenters)));
		for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
		{	matrix pSub[3]; vector3<complex*> pSubData;
			for(int iDir=0; iDir<3; iDir++)
			{	matrix pBlochCur = pBloch[iDir][kMesh[i].point.iReduced + iSpin*qCount];
				if(kMesh[i].point.invert<0) //apply complex conjugate:
					callPref(eblas_dscal)(pBlochCur.nData(), -1., ((double*)pBlochCur.dataPref())+1, 2);
				pSub[iDir] = dagger(kMesh[i].U) * pBlochCur * kMesh[i].U;
				pSubData[iDir] = pSub[iDir].data();
			}
			//Apply spatial transformation:
			matrix3<complex> rot(inv(e.gInfo.R * sym[kMesh[i].point.iSym] * e.gInfo.invR)); //cartesian symmetry matrix
			for(size_t j=0; j<pSub[0].nData(); j++)
				storeVector(rot * loadVector(pSubData,j), pSubData,j);
			//Accumulate with each requested Bloch phase
			for(int iDir=0; iDir<3; iDir++)
			{	std::vector<matrix>::iterator pWannierIter = pWannier[iDir].begin();
				for(auto cell: iCellMap)
					*(pWannierIter++) += (cell.second * kMesh[i].point.weight * cis(2*M_PI*dot(kMesh[i].point.k, cell.first))) * pSub[iDir];
			}
		}
		for(int iDir=0; iDir<3; iDir++) for(matrix& p: pWannier[iDir]) p.allReduce(MPIUtil::ReduceSum);
		//--- save to file
		if(mpiUtil->isHead())
			for(int iDir=0; iDir<3; iDir++)
			{	string fname = wannier.getFilename(Wannier::FilenameDump, string("mlwfP")+"xyz"[iDir], &iSpin);
				logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
				FILE* fp = fopen(fname.c_str(), "wb");
				if(!fp) die("Failed to open file '%s' for binary write.\n", fname.c_str());
				double normTot=0., normIm=0.;
				for(matrix& p: pWannier[iDir])
				{	if(realPartOnly)
					{	p.write_real(fp);
						//Collect imaginary part
						normTot += pow(nrm2(p), 2);
						eblas_dscal(p.nData(), 0., ((double*)p.data()), 2); //zero out real parts
						normIm += pow(nrm2(p), 2);
					}
					else p.write(fp);
				}
				fclose(fp);
				if(realPartOnly) logPrintf("done. Relative discarded imaginary part: %le\n", sqrt(normIm/normTot)); else logPrintf("done.\n"); logFlush();
			}
	}
	
	//Electron-phonon matrix elements:
	if(wannier.phononSup.length_squared())
	{	//--- generate list of commensurate k-points in order present in the unit cell calculation
		int prodPhononSup = wannier.phononSup[0] * wannier.phononSup[1] * wannier.phononSup[2];
		std::vector<int> ikArr; ikArr.reserve(prodPhononSup);
		for(unsigned ik=0; ik<kMesh.size(); ik++)
		{	vector3<> kSup = kMesh[ik].point.k * Diag(wannier.phononSup);
			double roundErr; round(kSup, &roundErr);
			if(roundErr < symmThreshold) //integral => commensurate with supercell
				ikArr.push_back(ik);
		}
		assert(int(ikArr.size()) == prodPhononSup);
		//---- generate pairs of commensurate k-points along with pointer to Wannier rotation
		struct KpointPair { vector3<> k1, k2; int ik1, ik2; };
		std::vector<KpointPair> kpointPairs; //pairs of k-points in the same order as matrices in phononHsub
		for(int ik1: ikArr)
		for(int ik2: ikArr)
		{	KpointPair kPair;
			kPair.k1 = kMesh[ik1].point.k;
			kPair.k2 = kMesh[ik2].point.k;
			kPair.ik1 = ik1;
			kPair.ik2 = ik2;
			kpointPairs.push_back(kPair);
		}
		int iPairStart = (kpointPairs.size() * mpiUtil->iProcess()) / mpiUtil->nProcesses();
		int iPairStop = (kpointPairs.size() * (mpiUtil->iProcess()+1)) / mpiUtil->nProcesses();
		//--- convert phononHsub from Bloch to wannier for each nuclear displacement mode
		string fnameIn = wannier.getFilename(Wannier::FilenameInit, "phononHsub", &iSpin);
		string fnameOut = wannier.getFilename(Wannier::FilenameDump, "mlwfHePh", &iSpin);
		logPrintf("Wannierizing '%s' to '%s' ... ", fnameIn.c_str(), fnameOut.c_str()); logFlush();
		//--- --- open input file
		MPIUtil::File fpIn;
		size_t matSizeIn = nBands*nBands * sizeof(complex);
		size_t modeStrideIn = kpointPairs.size() * matSizeIn; //input data size per phonon mode
		mpiUtil->fopenRead(fpIn, fnameIn.c_str(), nPhononModes * modeStrideIn);
		//--- --- open output file
		FILE* fpOut = 0;
		if(mpiUtil->isHead()) fpOut = fopen(fnameOut.c_str(), "w");
		bool isOutOpen = fpOut; mpiUtil->bcast(isOutOpen);
		if(!isOutOpen) die("Could not open '%s' for writing.\n", fnameOut.c_str());
		//--- --- loop over modes
		double normTot=0., normIm=0.;
		for(int iMode=0; iMode<nPhononModes; iMode++)
		{	//Read phononHsub and apply Wannier rotations:
			std::vector<matrix> phononHsub(kpointPairs.size());
			mpiUtil->fseek(fpIn, iMode*modeStrideIn + iPairStart*matSizeIn, SEEK_SET);
			double kPairWeight = 1./(prodPhononSup*prodPhononSup);
			for(int iPair=iPairStart; iPair<iPairStop; iPair++)
			{	matrix Hsub(nBands, nBands);
				mpiUtil->fread(Hsub.data(), sizeof(complex), nBands*nBands, fpIn);
				const KpointPair& pair = kpointPairs[iPair];
				phononHsub[iPair] = kPairWeight * (dagger(kMesh[pair.ik1].U) * Hsub * kMesh[pair.ik2].U); //save with Wannier rotations and k-integration weights
			}
			//Transform to real space (on phononCellMap squared)
			for(const auto& entry1: phononCellMap)
			for(const auto& entry2: phononCellMap)
			{	const vector3<int>& iR1 = entry1.first;
				const vector3<int>& iR2 = entry2.first;
				matrix HePh = zeroes(nCenters, nCenters);
				for(int iPair=iPairStart; iPair<iPairStop; iPair++)
				{	const KpointPair& pair = kpointPairs[iPair];
					HePh += cis(2*M_PI*(dot(pair.k2,iR2) - dot(pair.k1,iR1))) * phononHsub[iPair];
				}
				HePh *= entry1.second * entry2.second; //nclude cell weights due to boundary symmetrization
				HePh.allReduce(MPIUtil::ReduceSum); //collect contributions over all pairs
				if(mpiUtil->isHead())
				{	if(realPartOnly)
					{	HePh.write_real(fpOut);
						//Collect imaginary part
						normTot += pow(nrm2(HePh), 2);
						eblas_dscal(HePh.nData(), 0., ((double*)HePh.data()), 2); //zero out real parts
						normIm += pow(nrm2(HePh), 2);
					}
					else HePh.write(fpOut);
				}
			}
		}
		if(mpiUtil->isHead()) fclose(fpOut);
		mpiUtil->fclose(fpIn);
		if(mpiUtil->isHead() && realPartOnly)
			logPrintf("done. Relative discarded imaginary part: %le\n", sqrt(normIm/normTot));
		else logPrintf("done.\n");
		logFlush();
	}
}
