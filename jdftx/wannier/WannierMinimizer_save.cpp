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
		logPrintf("Reading initial rotations from '%s' ... ", fname.c_str()); logFlush();
		FILE* fp = fopen(fname.c_str(), "r");
		if(!fp) die("could not open '%s' for reading.\n", fname.c_str());
		for(auto& ke: kMesh)
		{	ke.U.init(nBands, nCenters);
			ke.U.read(fp);
		}
		fclose(fp);
		logPrintf("done.\n"); logFlush();
		logPrintf("NOTE: ignoring trial orbitals since we are resuming a previous WannierMinimize.\n");
	}
	
	//Load frozen rotations if necessary:
	std::vector<matrix> Ufrozen(kMesh.size());
	if(nFrozen)
	{	//Read U
		logPrintf("Reading frozen rotations from '%s' ... ", wannier.frozenUfilename.c_str()); logFlush();
		FILE* fp = fopen(wannier.frozenUfilename.c_str(), "r");
		if(!fp) die("could not open '%s' for reading.\n", wannier.frozenUfilename.c_str());
		for(size_t ik=0; ik<kMesh.size(); ik++)
		{	Ufrozen[ik].init(nBands, nFrozen);
			Ufrozen[ik].read(fp);
		}
		fclose(fp);
		logPrintf("done.\n"); logFlush();
	}
	
	//Compute / check the initial rotations:
	ostringstream ossErr;
	for(size_t ik=0; ik<kMesh.size(); ik++) if(isMine_q(ik,iSpin))
	{	KmeshEntry& ke = kMesh[ik];
		string kString;
		{	ostringstream kOss;
			kOss << " at k = [ " << ke.point.k[0] << ' ' << ke.point.k[1] << ' ' << ke.point.k[2] << " ]\n";
			kString = kOss.str();
		}
		//Band ranges:
		int bStart=0, bStop=0, bFixedStart=0, bFixedStop=0;
		const std::vector<double>& eigs = e.eVars.Hsub_eigs[ke.point.iReduced + iSpin*qCount];
		if(wannier.outerWindow)
		{	bStart = 0;
			while(bStart<nBands && eigs[bStart]<wannier.eOuterMin)
				bStart++;
			bStop = bStart;
			while(bStop<nBands && eigs[bStop]<=wannier.eOuterMax)
				bStop++;
			if(bStop-bStart < nCenters)
			{	ossErr << "Number of bands within outer window = " << bStop-bStart
					<< " less than nCenters = " << nCenters << kString;
				break;
			}
			//Optionally range for inner window:
			if(wannier.innerWindow)
			{	bFixedStart = bStart;
				while(bFixedStart<bStop && eigs[bFixedStart]<wannier.eInnerMin)
					bFixedStart++;
				bFixedStop = bFixedStart;
				while(bFixedStop<bStop && eigs[bFixedStop]<=wannier.eInnerMax)
					bFixedStop++;
				if(bFixedStop-bFixedStart > nCenters)
				{	ossErr << "Number of bands within inner window = " << bStop-bStart
						<< " exceeds nCenters = " << nCenters << kString;
					break;
				}
			}
			else bFixedStart = bFixedStop = bStart; //fixed interval is empty
		}
		else //fixed bands
		{	bFixedStart = bStart = wannier.bStart;
			bFixedStop  = bStop  = wannier.bStart + nCenters;
		}
		
		//Check compatibility of frozen window with band range:
		if(nFrozen)
		{	const double tol = 1e-6 * nFrozen;
			if( (bFixedStart>0 && nrm2(Ufrozen[ik](0,bFixedStart, 0,nFrozen)) > tol)
			 || (bFixedStop<nBands && nrm2(Ufrozen[ik](bFixedStop,nBands, 0,nFrozen)) > tol) )
			{	ossErr << "Frozen rotations are incompatible with current outer window" << kString;
				break;
			}
		}
		
		//Initial rotation of bands to get to Wannier subspace:
		ke.nFixed = bFixedStop - bFixedStart;
		int nFree = nCenters - ke.nFixed;
		ke.nIn = (nFree>0) ? (bStop-bStart) : nCenters; //number of bands contributing to Wannier subspace
		const double tol = 1e-6 * nCenters; //tolerance for checking unitarity etc.
		//--- Create fixed part of U1 (if any):
		ke.U1 = zeroes(nBands, ke.nIn);
		if(ke.nFixed > 0)
		{	matrix U1fixed = eye(ke.nFixed);
			if(nFrozen)
			{	if(ke.nFixed < nFrozen)
				{	ossErr << "Number of bands within " << (wannier.innerWindow ? "inner window" : "fixed band range")
						<< " = " << bStop-bStart << " less than nFrozen = " << nFrozen << kString;
					break;
				}
				//Rotate the fixed subspace so that the frozen centers come first:
				matrix UfrozenNZ = Ufrozen[ik](bFixedStart,bFixedStop, 0,nFrozen); //drop the zero rows (should be zero, else SVD check below will fail)
				U1fixed.set(0,ke.nFixed, 0,nFrozen, UfrozenNZ);
				//fill in the null space of UfrozenNZ using an SVD
				matrix U, Vdag; diagMatrix S;
				UfrozenNZ.svd(U, S, Vdag);
				if(nrm2(S(0,nFrozen)-eye(nFrozen))>tol)
				{	ossErr << "Frozen rotations are incompatible with current outer window" << kString;
					break;
				}
				if(nFrozen<ke.nFixed)
					U1fixed.set(0,ke.nFixed, nFrozen,ke.nFixed, U(0,ke.nFixed, nFrozen,ke.nFixed));
			}
			ke.U1.set(bFixedStart,bFixedStop, 0,ke.nFixed, U1fixed);
		}
		if(wannier.loadRotations)
		{	//Factorize U (nBands x nCenters) into U1 (nBands x nIn) and U2 (nCenters x nCenters):
			//--- check unitarity:
			if(nrm2(dagger(ke.U) * ke.U - eye(nCenters)) > tol)
			{	ossErr << "Initial rotations U are not unitary" << kString;
				break;
			}
			//--- check compatibility with frozen rotations:
			if(nFrozen && nrm2(ke.U(0,nBands, 0,nFrozen) - Ufrozen[ik])>tol)
			{	ossErr << "Initial rotations U are incompatible with frozen rotations" << kString;
				break;
			}
			//--- determine the free columns of U1:
			if(nFree > 0)
			{	matrix U1free = ke.U;
				//project out the fixed columns (if any):
				if(ke.nFixed > 0)
				{	U1free -= ke.U1*(dagger(ke.U1) * ke.U);
					//SVD to find the remaining non-null columns:
					matrix U, Vdag; diagMatrix S;
					U1free.svd(U, S, Vdag);
					if(nrm2(S(0,nFree)-eye(nFree))>tol || nrm2(S(nFree,nCenters))>tol)
					{	ossErr << "Initial rotations are incompatible with current inner window" << kString;
						break;
					}
					U1free = U(0,nBands, 0,nFree); //discard null subspace
				}
				ke.U1.set(0,nBands, ke.nFixed,nCenters, U1free);
			}
			//--- check U1:
			if(nrm2((dagger(ke.U1) * ke.U1)(0,nCenters, 0,nCenters) - eye(nCenters)) > tol)
			{	ossErr << "Initial rotations U1 are not unitary" << kString;
				break;
			}
			if( (bStart>0 && nrm2(ke.U1(0,bStart, 0,nCenters))>tol) 
				|| (bStop<nBands && nrm2(ke.U1(bStop,nBands, 0,nCenters))>tol) )
			{	ossErr << "Initial rotations are incompatible with current outer window / band selection" << kString;
				break;
			}
			//--- calculate and check U2:
			ke.U2 = eye(nCenters);
			ke.U2.set(nFrozen,nCenters, nFrozen,nCenters,
				dagger(ke.U1(0,nBands, nFrozen,nCenters)) * ke.U(0,nBands, nFrozen,nCenters));
			if(nrm2(dagger(ke.U2) * ke.U2 - eye(nCenters)) > tol)
			{	ossErr << "Initial rotations U2 are not unitary" << kString;
				break;
			}
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
			int nNew = nCenters - nFrozen; //number of new centers
			//--- Pick up best linear combination of remaining bands (if any)
			if(nFree > 0)
			{	//Create overlap matrix with contribution from fixed bands projected out:
				matrix CdagGFree;
				if(ke.nFixed > 0)
				{	//SVD the fixed band contribution to the trial space:
					matrix U, Vdag; diagMatrix S;
					CdagG(bFixedStart,bFixedStop, 0,nNew).svd(U, S, Vdag);
					//Project out the fixed bands (use only the zero singular values)
					CdagGFree = CdagG * dagger(Vdag(nNew-nFree,nNew, 0,nNew));
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
			matrix WdagG = dagger(ke.U1(0,nBands, nFrozen,nCenters)) * CdagG;
			ke.U2 = eye(nCenters);
			ke.U2.set(nFrozen,nCenters, nFrozen,nCenters, WdagG * invsqrt(dagger(WdagG) * WdagG));
		}
		//Make initial rotations exactly unitary:
		ke.U1 = fixUnitary(ke.U1);
		ke.U2 = fixUnitary(ke.U2);
	}
	mpiUtil->checkErrors(ossErr);
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
	}
	
	//Output range of centers that span each band
	{	string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfBandContrib", &iSpin);
		logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
		//Find range of centers that span each band:
		std::vector<int> nMin(nBands, nCenters), nMax(nBands, 0);
		for(const KmeshEntry& ke: kMesh)
			for(int n=0; n<nCenters; n++)
				for(int b=0; b<nBands; b++)
					if(ke.U(b,n).norm() > 1e-12)
					{	nMin[b] = std::min(nMin[b], n);
						nMax[b] = std::max(nMax[b], n);
					}
		//Find energy range that is exact for each range of centers:
		std::map< std::pair<int,int>, std::pair<double,double> > nRangeToErange;
		nRangeToErange[std::make_pair(-1,-1)] = std::make_pair(NAN,NAN);
		for(int b=0; b<nBands; b++)
			if(nMin[b] <= nMax[b])
			{	std::pair<int,int> nRange(nMin[b], nMax[b]);
				if(nRangeToErange.find(nRange) == nRangeToErange.end())
				{	double eMin =-INFINITY, eMax = +INFINITY;
					for(size_t ik=0; ik<kMesh.size(); ik++) if(isMine_q(ik,iSpin))
					{	KmeshEntry& ke = kMesh[ik];
						const std::vector<double>& eigs = e.eVars.Hsub_eigs[ke.point.iReduced + iSpin*qCount];
						matrix Utrunc = ke.U(0,nBands, nMin[b],nMax[b]+1);
						diagMatrix overlap = diag(Utrunc * dagger(Utrunc));
						double eMin_k = +INFINITY, eMax_k = -INFINITY;
						for(int iBand=0; iBand<nBands; iBand++)
							if(fabs(overlap[iBand] - 1.) < 1e-6)
							{	eMin_k = std::min(eMin_k, iBand ? eigs[iBand-1] : -INFINITY);
								eMax_k = std::max(eMax_k, (iBand<nBands) ? eigs[iBand+1] : +INFINITY);
							}
						eMin = std::max(eMin, eMin_k);
						eMax = std::min(eMax, eMax_k);
					}
					mpiUtil->allReduce(eMin, MPIUtil::ReduceMax);
					mpiUtil->allReduce(eMax, MPIUtil::ReduceMin);
					if(eMin >= eMax) { eMin = eMax = NAN; }
					nRangeToErange[nRange] = std::make_pair(eMin, eMax);
				}
			}
			else { nMin[b] = nMax[b] = -1; } //unused band
		if(mpiUtil->isHead())
		{	FILE* fp = fopen(fname.c_str(), "w");
			for(int b=0; b<nBands; b++)
			{	const std::pair<double,double> eRange = nRangeToErange[std::make_pair(nMin[b],nMax[b])];
				fprintf(fp, "%d %d   %+10.5lf %+10.5lf\n", nMin[b], nMax[b], eRange.first, eRange.second);
			}
			fclose(fp);
		}
		logPrintf("done.\n"); logFlush();
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
	int nqMine = 0;
	for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
		nqMine++;
	nqMine = std::max(nqMine,1); //avoid zero size matrices below
	matrix phase = zeroes(nqMine, iCellMap.size());
	{	matrix HwannierTilde = zeroes(nCenters*nCenters, nqMine);
		int iqMine = 0;
		for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
		{	//Apply MLWF-optimized unitary rotations to Hamiltonian at current k:
			matrix Hsub = dagger(kMesh[i].U) * e.eVars.Hsub_eigs[kMesh[i].point.iReduced + iSpin*qCount] * kMesh[i].U;
			callPref(eblas_copy)(HwannierTilde.dataPref()+HwannierTilde.index(0,iqMine), Hsub.dataPref(), Hsub.nData());
			//Calculate required phases:
			int iCell = 0;
			for(auto cell: iCellMap)
				phase.set(iqMine, iCell++, cell.second * kMesh[i].point.weight * cis(2*M_PI*dot(kMesh[i].point.k, cell.first)));
			iqMine++;
		}
		//Fourier transform to Wannier space and save
		matrix Hwannier = HwannierTilde * phase;
		Hwannier.allReduce(MPIUtil::ReduceSum);
		dumpMatrix(Hwannier, "mlwfH", realPartOnly, iSpin);
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
		matrix pWannierTilde = zeroes(nCenters*nCenters*3, nqMine);
		int iqMine = 0;
		for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
		{	matrix pSub(nCenters*nCenters, 3);
			for(int iDir=0; iDir<3; iDir++)
			{	matrix pSubDir = pBloch[iDir][kMesh[i].point.iReduced + iSpin*qCount];
				if(kMesh[i].point.invert<0) //apply complex conjugate:
					callPref(eblas_dscal)(pSubDir.nData(), -1., ((double*)pSubDir.dataPref())+1, 2);
				pSubDir = dagger(kMesh[i].U) * pSubDir * kMesh[i].U; //apply MLWF-optimized rotations
				callPref(eblas_copy)(pSub.dataPref()+pSub.index(0,iDir), pSubDir.dataPref(), pSubDir.nData());
			}
			//Store with spatial transformation:
			matrix rot(inv(e.gInfo.R * sym[kMesh[i].point.iSym] * e.gInfo.invR)); //cartesian symmetry matrix
			pSub = pSub * dagger(rot);
			callPref(eblas_copy)(pWannierTilde.dataPref()+pWannierTilde.index(0,iqMine), pSub.dataPref(), pSub.nData());
			iqMine++;
		}
		//Fourier transform to Wannier space and save
		matrix pWannier = pWannierTilde * phase;
		pWannier.allReduce(MPIUtil::ReduceSum);
		dumpMatrix(pWannier, "mlwfP", realPartOnly, iSpin);
	}
	
	//Momentum-squared matrix element in Wannier basis (for CEDA):
	if(wannier.saveMomenta && wannier.ceda)
	{	int dirPairs[6][2] = { {0,0}, {1,1}, {2,2}, {1,2}, {2,0}, {0,1} }; //pairs of directions in stored order for P^2 matrix elements
		//--- compute momentum squared matrix elements of Bloch states:
		std::vector< std::vector<matrix> > pSqBloch(6, std::vector<matrix>(e.eInfo.nStates));
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
			if(e.eInfo.qnums[q].spin==iSpin)
			{	for(int iDirPair=0; iDirPair<6; iDirPair++)
					pSqBloch[iDirPair][q] = e.gInfo.detR * (e.eVars.C[q] ^ DD(e.eVars.C[q], dirPairs[iDirPair][0], dirPairs[iDirPair][1])); //factor of iota dropped for consistency with above
			}
		//--- convert to Wannier basis:
		matrix pSqWannierTilde = zeroes(nCenters*nCenters*6, nqMine);
		int iqMine = 0;
		for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
		{	matrix pSqSub(nCenters*nCenters, 6);
			for(int iDirPair=0; iDirPair<6; iDirPair++)
			{	matrix pSqSubDir = pSqBloch[iDirPair][kMesh[i].point.iReduced + iSpin*qCount];
				if(kMesh[i].point.invert<0) //apply complex conjugate:
					callPref(eblas_dscal)(pSqSubDir.nData(), -1., ((double*)pSqSubDir.dataPref())+1, 2);
				pSqSubDir = dagger(kMesh[i].U) * pSqSubDir * kMesh[i].U; //apply MLWF-optimized rotations
				callPref(eblas_copy)(pSqSub.dataPref()+pSqSub.index(0,iDirPair), pSqSubDir.dataPref(), pSqSubDir.nData());
			}
			//Initialize rotation:
			matrix3<> rot = inv(e.gInfo.R * sym[kMesh[i].point.iSym] * e.gInfo.invR); //cartesian symmetry matrix
			matrix rotSq(6,6); //corresponding transformation of symmetric rank-2 tensor
			for(int iDirPair=0; iDirPair<6; iDirPair++)
			{	//Construct unrotated Cartesian tensor for this direction pair
				matrix3<> e;
				e(dirPairs[iDirPair][0], dirPairs[iDirPair][1]) = 1.;
				e = 0.5*(e + ~e); //symmetrize
				//Rotate it:
				e = rot * e * (~rot);
				//Extract direction pairs from rotated tensor:
				for(int jDirPair=0; jDirPair<6; jDirPair++)
					rotSq.set(jDirPair, iDirPair, e(dirPairs[jDirPair][0], dirPairs[jDirPair][1]));
			}
			//Store with spatial transformation:
			pSqSub = pSqSub * dagger(rotSq);
			callPref(eblas_copy)(pSqWannierTilde.dataPref()+pSqWannierTilde.index(0,iqMine), pSqSub.dataPref(), pSqSub.nData());
			iqMine++;
		}
		//Fourier transform to Wannier space and save
		matrix pSqWannier = pSqWannierTilde * phase;
		pSqWannier.allReduce(MPIUtil::ReduceSum);
		dumpMatrix(pSqWannier, "mlwfPsq", realPartOnly, iSpin);
	}
	
	//Electron-electron linewidths:
	{	string fname = wannier.getFilename(Wannier::FilenameInit, "ImSigma_ee", &iSpin);
		if(fileSize(fname.c_str()) >= 0)
		{	//Read Bloch version:
			logPrintf("Reading '%s' ... ", fname.c_str()); logFlush();
			std::vector<diagMatrix> ImSigma_ee;
			e.eInfo.read(ImSigma_ee, fname.c_str());
			logPrintf("done.\n"); logFlush();
			//Fill in linewidths for states out of range (to ensure smoothness in Wannier interpolation)
			//--- determine weights in fit:
			std::vector<diagMatrix> fitWeight(e.eInfo.nStates);
			for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
			{	fitWeight[q].reserve(nBands);
				for(double IS: ImSigma_ee[q]) //first pass: whether ImSigma is valid or not
					fitWeight[q].push_back(std::isnan(IS) ? 0. : 1.);
				for(int b=0; b<nBands; b++) if(fitWeight[q][b])
				{	double minEdist = 100.; //keep finite so that no error below if all ImSigma's already defined
					for(int b2=0; b2<nBands; b2++)
						if(!fitWeight[q][b2])
							minEdist = std::min(minEdist, fabs(e.eVars.Hsub_eigs[q][b]-e.eVars.Hsub_eigs[q][b2]));
					fitWeight[q][b] = 1./hypot(minEdist, 0.1); //cap the weights at some reasonable amount
				}
			}
			//--- construct matrices for polynomial fit:
			int order = 3; //quadratic
			int nRows = e.eInfo.nStates * nBands;
			matrix Lhs = zeroes(nRows, order);
			matrix rhs = zeroes(nRows, 1);
			diagMatrix w(nRows, 0.);
			int row = e.eInfo.qStart * nBands;
			for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
				for(int b=0; b<nBands; b++)
				{	double E = e.eVars.Hsub_eigs[q][b], Epow = 1.;
					for(int n=0; n<order; n++)
					{	Lhs.set(row,n, Epow);
						Epow *= E;
					}
					rhs.set(row,0, fitWeight[q][b] ? ImSigma_ee[q][b] : 0.); //Handle NANs correctly
					w[row] = fitWeight[q][b];
					row++;
				}
			Lhs.allReduce(MPIUtil::ReduceSum);
			rhs.allReduce(MPIUtil::ReduceSum);
			w.allReduce(MPIUtil::ReduceSum);
			//--- weighted least squares polynomial fit
			matrix rhsFit = Lhs * (inv(dagger(Lhs)*w*Lhs) * (dagger(Lhs)*w*rhs));
			//--- fill in missing values
			row = e.eInfo.qStart * nBands;
			for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
				for(int b=0; b<nBands; b++)
				{	if(std::isnan(ImSigma_ee[q][b]))
						ImSigma_ee[q][b] =  rhsFit(row,0).real();
					row++;
				}
			//Wannierize:
			matrix ImSigma_eeWannierTilde = zeroes(nCenters*nCenters, nqMine);
			int iqMine = 0;
			for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
			{	matrix ImSigma_eeSub = dagger(kMesh[i].U) * ImSigma_ee[kMesh[i].point.iReduced + iSpin*qCount] * kMesh[i].U;
				callPref(eblas_copy)(ImSigma_eeWannierTilde.dataPref()+ImSigma_eeWannierTilde.index(0,iqMine), ImSigma_eeSub.dataPref(), ImSigma_eeSub.nData());
				iqMine++;
			}
			matrix ImSigma_eeWannier = ImSigma_eeWannierTilde * phase;
			ImSigma_eeWannier.allReduce(MPIUtil::ReduceSum);
			dumpMatrix(ImSigma_eeWannier, "mlwfImSigma_ee", realPartOnly, iSpin);
		}
	}
	
	//Electron-phonon linewidths:
	{	string fname = wannier.getFilename(Wannier::FilenameInit, "ImSigma_ePh", &iSpin);
		off_t fsize = fileSize(fname.c_str());
		if(fsize >= 0)
		{	//Read from file:
			int nBandsIn = fsize / (sizeof(double) * e.eInfo.nStates);
			if(int(nBandsIn * e.eInfo.nStates * sizeof(double)) != fsize)
				die("Length of file '%s' = %ld is not a multiple of nStates = %d doubles.\n", fname.c_str(), fsize, e.eInfo.nStates);
			logPrintf("Reading '%s' ... ", fname.c_str()); logFlush();
			std::vector<diagMatrix> ImSigma_ePh;
			e.eInfo.read(ImSigma_ePh, fname.c_str(), nBandsIn);
			logPrintf("done.\n"); logFlush();
			//Fill in extra bands if necessary:
			if(nBandsIn < nBands)
				for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
					for(int b=nBandsIn; b<nBands; b++)
						ImSigma_ePh[q].push_back(ImSigma_ePh[q].back());
			//Convert to log for better interpolation:
			for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
				for(double& IS: ImSigma_ePh[q])
					IS = log(IS);
			//Wannierize:
			matrix ImSigma_ePhWannierTilde = zeroes(nCenters*nCenters, nqMine);
			int iqMine = 0;
			for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
			{	matrix ImSigma_ePhSub = dagger(kMesh[i].U) * ImSigma_ePh[kMesh[i].point.iReduced + iSpin*qCount] * kMesh[i].U;
				callPref(eblas_copy)(ImSigma_ePhWannierTilde.dataPref()+ImSigma_ePhWannierTilde.index(0,iqMine), ImSigma_ePhSub.dataPref(), ImSigma_ePhSub.nData());
				iqMine++;
			}
			matrix ImSigma_ePhWannier = ImSigma_ePhWannierTilde * phase;
			ImSigma_ePhWannier.allReduce(MPIUtil::ReduceSum);
			dumpMatrix(ImSigma_ePhWannier, "mlwfImSigma_ePh", realPartOnly, iSpin);
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
		//--- generate pairs of commensurate k-points along with pointer to Wannier rotation
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
		int iPairStart, iPairStop;
		TaskDivision(kpointPairs.size(), mpiUtil).myRange(iPairStart, iPairStop);
		//--- calculate Fourier transform phase (with integration weights)
		int nPairsMine = std::max(1, iPairStop-iPairStart); //avoid zero size matrices below
		matrix phase = zeroes(nPairsMine, std::pow(phononCellMap.size(),2));
		double kPairWeight = 1./(prodPhononSup*prodPhononSup);
		for(int iPair=iPairStart; iPair<iPairStop; iPair++)
		{	const KpointPair& pair = kpointPairs[iPair];
			int iCellPair = 0;
			for(const auto& entry1: phononCellMap)
			for(const auto& entry2: phononCellMap)
			{	const vector3<int>& iR1 = entry1.first;
				const vector3<int>& iR2 = entry2.first;
				double weight = kPairWeight * entry1.second * entry2.second;
				phase.set(iPair-iPairStart, iCellPair++, weight * cis(2*M_PI*(dot(pair.k2,iR2) - dot(pair.k1,iR1))) );
			}
		}
		//--- read Bloch electron - nuclear displacement matrix elements
		string fnameIn = wannier.getFilename(Wannier::FilenameInit, "phononHsub", &iSpin);
		logPrintf("Reading '%s' ... ", fnameIn.c_str()); logFlush();
		MPIUtil::File fpIn;
		size_t matSizeIn = nBands*nBands * sizeof(complex);
		size_t modeStrideIn = kpointPairs.size() * matSizeIn; //input data size per phonon mode
		mpiUtil->fopenRead(fpIn, fnameIn.c_str(), nPhononModes * modeStrideIn);
		matrix HePhTilde = zeroes(nCenters*nCenters*nPhononModes, nPairsMine);
		for(int iMode=0; iMode<nPhononModes; iMode++)
		{	//Read phononHsub and store with Wannier rotations:
			mpiUtil->fseek(fpIn, iMode*modeStrideIn + iPairStart*matSizeIn, SEEK_SET);
			for(int iPair=iPairStart; iPair<iPairStop; iPair++)
			{	const KpointPair& pair = kpointPairs[iPair];
				matrix phononHsub(nBands, nBands);
				mpiUtil->fread(phononHsub.data(), sizeof(complex), phononHsub.nData(), fpIn); //read from file
				phononHsub = (dagger(kMesh[pair.ik1].U) * phononHsub * kMesh[pair.ik2].U); //apply Wannier rotations
				callPref(eblas_copy)(HePhTilde.dataPref() + HePhTilde.index(0,iPair-iPairStart) + nCenters*nCenters*iMode,
					phononHsub.dataPref(), phononHsub.nData());
			}
		}
		logPrintf("done.\n"); logFlush();
		//--- convert phononHsub from Bloch to wannier for each nuclear displacement mode
		matrix HePh = HePhTilde * phase;
		HePhTilde = 0; //free memory
		HePh.allReduce(MPIUtil::ReduceSum);
		//--- save output:
		dumpMatrix(HePh, "mlwfHePh", realPartOnly, iSpin);
	}
	
	//Process subclass-specific outputs, if any:
	saveExtra(iSpin);
}


void WannierMinimizer::dumpMatrix(const matrix& H, string varName, bool realPartOnly, int iSpin) const
{	if(mpiUtil->isHead())
		H.dump(wannier.getFilename(Wannier::FilenameDump, varName, &iSpin).c_str(), realPartOnly);
}
