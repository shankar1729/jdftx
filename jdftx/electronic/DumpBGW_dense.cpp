/*-------------------------------------------------------------------
Copyright 2019 Ravishankar Sundararaman

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

#ifdef HDF5_ENABLED //BGW output requires HDF5

#include <electronic/DumpBGW_internal.h>
#include <electronic/ColumnBundle.h>
#include <electronic/ExactExchange.h>

#ifndef SCALAPACK_ENABLED //Dense solve requires ScaLAPACK

void BGW::denseWriteWfn(hid_t gidWfns)
{
	die("\n\nDense solve for extra bands for BGW (nBandsDense > 0) requires ScaLAPACK.\n\n");
}

#else //SCALAPACK_ENABLED

extern "C"
{
	void blacs_pinfo_(int* mypnum, int* nprocs);
	void blacs_get_(const int* icontxt, const int* what, int* val);
	void blacs_gridinit_(const int* icontxt, const char* layout, const int* nprow, const int* npcol);
	void blacs_gridinfo_(const int* icontxt, int* nprow, int* npcol, int* myprow, int* mypcol);
	void blacs_gridexit_(const int* icontxt);
	void blacs_exit_(const int* cont);
	
	void descinit_(int* desc, const int* m, const int* n, const int* mb, const int* nb,
		const int* irsrc, const int* icsrc, const int* ictxt, const int* lld, int* info);
	int numroc_(const int* n, const int* nb, const int* iproc, const int* srcproc, const int* nprocs);
	
	void pzheevx_(const char *jobz, const char *range, const char *uplo, const int* n, complex* a, const int* ia, const int* ja, int* desca,
		double* vl, double* vu, int* il, int* iu, double* abstol, int* m, int* nz, double* w, double* orfac,
		complex* z, const int* iz, const int* jz, int* descz, complex* work, int* lwork, double* rwork, int* lrwork,
		int* iwork, int* liwork, int* ifail, int* iclustr, double* gap, int* info);
	
	void pzgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
		const complex* alpha, const complex* a, const int* ia, const int* ja, const int* desca,
		const complex* b, const int* ib, const int* jb, const int* descb, const complex* beta,
		complex* c, const int* ic, const int* jc, const int* descc);
	
	void pzgemr2d_(const int* m, const int* n,
		const complex* a, const int* ia, const int* ja, const int* desca,
		complex* b, const int* ib, const int* jb, const int* descb,
		const int* ictxt);
}


//Return list of indices in a given dimension (row or column) that belong to me in block-cyclic distribution
std::vector<int> distributedIndices(int nTotal, int blockSize, int iProcDim, int nProcsDim)
{	int zero = 0;
	int nMine = numroc_(&nTotal, &blockSize, &iProcDim, &zero, &nProcsDim);
	std::vector<int> myIndices; myIndices.reserve(nMine);
	int blockStride = blockSize * nProcsDim;
	int nBlocksMineMax = (nTotal + blockStride - 1) / blockStride;
	for(int iBlock=0; iBlock<nBlocksMineMax; iBlock++)
	{	int iStart = iProcDim*blockSize + iBlock*blockStride;
		int iStop = std::min(iStart+blockSize, nTotal);
		for(int i=iStart; i<iStop; i++)
			myIndices.push_back(i);
	}
	assert(int(myIndices.size()) == nMine);
	return myIndices;
}


//Create a vector by indexing an attay i.e. return v[index] in octave/numpy notation
template<typename T> std::vector<T> indexVector(const T* v, const std::vector<int>& index)
{	std::vector<T> result;
	result.reserve(index.size());
	for(int i: index)
		result.push_back(v[i]);
	return result;
}


void transformV(int q, matrix& V, const matrix& evecs, std::vector<matrix>& Vsub,
	const ElecInfo& eInfo, int desc[9], int nProcsRow, int nProcsCol,
	int nRows, int nEigs, int nRowsMine, int nColsMine, int blockSize,
	int iProcRow, int iProcCol)
{
	static StopWatch watch("scalapackTransformV"); watch.start();
	matrix Vevecs = zeroes(nRowsMine, nColsMine); //local part of V * evecs
	//--- parallel matrix multiplies:
	complex alpha(1.,0.), beta(0.,0.); int one = 1;
	pzgemm_("N", "N", &nRows, &nRows, &nRows, &alpha,
		V.data(), &one, &one, desc,
		evecs.data(), &one, &one, desc, &beta,
		Vevecs.data(), &one, &one, desc); //VEvecs = V * evecs
	pzgemm_("C", "N", &nRows, &nRows, &nRows, &alpha,
		evecs.data(), &one, &one, desc,
		Vevecs.data(), &one, &one, desc, &beta,
		V.data(), &one, &one, desc); //VNew = evecx' * VEvecs
	Vevecs = 0; //cleanup memory
	
	//Collect required portion of V on a single process (the one that owns state q):
	if(eInfo.isMine(q))
	{	Vsub[q] = zeroes(nEigs, nEigs);
		for(int jProcRow=0; jProcRow<nProcsRow; jProcRow++)
			for(int jProcCol=0; jProcCol<nProcsCol; jProcCol++)
			{	int jProcess = jProcRow * nProcsCol + jProcCol;
				//Determine indices within nEigs for this block:
				std::vector<int> jEigRowsMine = distributedIndices(nEigs, blockSize, jProcRow, nProcsRow);
				std::vector<int> jEigColsMine = distributedIndices(nEigs, blockSize, jProcCol, nProcsCol);
				
				if((!jEigRowsMine.size()) or (!jEigColsMine.size()))
					continue; //nothing from this block
				//Get current block: locally or from another process
				matrix Vcur;
				if(jProcess == mpiWorld->iProcess())
					Vcur = V(0,jEigRowsMine.size(), 0,jEigColsMine.size()); //local
				else
				{	Vcur.init(jEigRowsMine.size(), jEigColsMine.size());
					mpiWorld->recvData(Vcur, jProcess, 0);
				}
				//Distribute to full matrix:
				const complex* inData = Vcur.data();
				for(int c: jEigColsMine)
					for(int r: jEigRowsMine)
						Vsub[q].set(r, c, *(inData++));
			}
	}
	else
	{	std::vector<int> iEigRowsMine = distributedIndices(nEigs, blockSize, iProcRow, nProcsRow);
		std::vector<int> iEigColsMine = distributedIndices(nEigs, blockSize, iProcCol, nProcsCol);
		if(iEigRowsMine.size() and iEigColsMine.size())
			mpiWorld->sendData(matrix(V(0,iEigRowsMine.size(), 0,iEigColsMine.size())), eInfo.whose(q), 0);
	}
	watch.stop();
}


//Solve wavefunctions using ScaLAPACK and write to hdf5 file:
void BGW::denseWriteWfn(hid_t gidWfns)
{	static StopWatch watchDiag("scalapackDiagonalize"), watchSetup("scalapackMatrixSetup"), watchIO("scalapackWriteHDF5");
	logPrintf("\n");
	nBands = bgwp.nBandsDense;
	nBandsV = bgwp.nBandsV ? bgwp.nBandsV : nBands;
	if(nBandsV > nBands) die("nBandsV = %d must be less than nBandsDense = %d\n", nBandsV, nBands);
	if(nSpinor > 1) die("\nDense diagonalization not yet implemented for spin-orbit / vector-spin modes.\n");
	for(const auto& sp: e.iInfo.species)
		if(sp->isUltrasoft())
			 die("\nDense diagonalization not supported for ultrasoft pseudopotentials.\n");
	//Calculate squarest possible process grid:
	int nProcesses = mpiWorld->nProcesses();
	int nProcsRow = int(round(sqrt(nProcesses)));
	while(nProcesses % nProcsRow) nProcsRow--;
	int nProcsCol = nProcesses / nProcsRow;

	//Initialize BLACS process grid:
	int blacsContext, blacsContextCol, iProcRow, iProcCol;
	{	int unused=-1, what=0;
		blacs_get_(&unused, &what, &blacsContext);
		blacs_gridinit_(&blacsContext, "Row-major", &nProcsRow, &nProcsCol);
		blacs_gridinfo_(&blacsContext, &nProcsRow, &nProcsCol, &iProcRow, &iProcCol);
		assert(mpiWorld->iProcess() == iProcRow * nProcsCol + iProcCol); //this mapping is assumed below, so check
		//Initialize trivial context for column-split output:
		int one = 1;
		blacs_get_(&unused, &what, &blacsContextCol);
		blacs_gridinit_(&blacsContextCol, "Row-major", &one, &nProcesses);
	}
	logPrintf("\tInitialized %d x %d process BLACS grid.\n", nProcsRow, nProcsCol);
	
	//Create dataset (must happen on all processes together):
	hsize_t nGtot = nBasisPrev.back() + nBasis.back();
	hsize_t dims[4] = { hsize_t(nBands), hsize_t(nSpins*nSpinor), nGtot, 2 };
	hid_t sid = H5Screate_simple(4, dims, NULL);
	hid_t did = H5Dcreate(gidWfns, "coeffs", H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t plid = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plid, H5FD_MPIO_INDEPENDENT);
	H5Sclose(sid);
	H5Fflush(gidWfns, H5F_SCOPE_GLOBAL);
	
	//Loop over states:
	for(int iSpin=0; iSpin<nSpins; iSpin++)
	for(int ik=0; ik<nReducedKpts; ik++)
	{	int q = ik + iSpin*nReducedKpts;
		const Basis& basis = e.basis[q];
		const QuantumNumber& qnum = eInfo.qnums[q];
		logPrintf("\tDiagonalizing state ");
		eInfo.kpointPrint(globalLog, q, true);
		
		//Initialize matrix distribution:
		const int nRows = basis.nbasis * nSpinor; //Hamiltonian dimension
		logPrintf(" with dimension %d at t[s]: %.2lf\n", nRows, clock_sec()); logFlush();
		const int blockSize = bgwp.blockSize; //block dimensions
		const int nEigs = bgwp.nBandsDense; //number of eigenvalues/eigenvectors requested
		if(nRows < blockSize * (std::max(nProcsRow, nProcsCol) - 1))
			die("\tNo data on some processes: reduce blockSize or # processes.\n");
		watchSetup.start();
		std::vector<int> iRowsMine = distributedIndices(nRows, blockSize, iProcRow, nProcsRow); //indices of rows on current process
		std::vector<int> iColsMine = distributedIndices(nRows, blockSize, iProcCol, nProcsCol); //indices of cols on current process
		int nRowsMine = iRowsMine.size();
		int nColsMine = iColsMine.size();
		int descH[9];
		{	int zero=0, info;
			descinit_(descH, &nRows, &nRows, &blockSize, &blockSize, &zero, &zero, &blacsContext, &nRowsMine, &info); assert(info==0);
		}
		
		//Create basis objects restricted to row and column ranges of current process:
		Basis basisRow, basisCol;
		basisRow.setup(*basis.gInfo, *basis.iInfo, indexVector(basis.index.data(), iRowsMine));
		basisCol.setup(*basis.gInfo, *basis.iInfo, indexVector(basis.index.data(), iColsMine));
		ColumnBundle Crow(1, nRowsMine, &basisRow, &qnum); //dummy ColumnBundle for SpeciesInfo::getV() etc.
		ColumnBundle Ccol(1, nColsMine, &basisCol, &qnum); //dummy ColumnBundle for SpeciesInfo::getV() etc.
		
		//Initialize Hamiltonian matrix:
		matrix H = zeroes(nRowsMine, nColsMine);
		matrix Vxc, Vxx;
		if(bgwp.saveVxc) Vxc = zeroes(nRowsMine, nColsMine);
		//--- Kinetic, potential and kinetic-potential contributions:
		{	complex* Hdata = H.data();
			complex* VxcData = (bgwp.saveVxc ? Vxc.data() : NULL);
			const vector3<int>* iGarr = basis.iGarr.data();
			//Prepare potential in reciprocal space:
			complexScalarFieldTilde Vtilde = Complex(J(eVars.Vscloc[qnum.index()]*(1./e.gInfo.dV)));
			const complex* Vdata = Vtilde->data();
			//Prepare potential in reciprocal space:
			complexScalarFieldTilde VxcTilde = Complex(J(eVars.Vxc[qnum.index()]));
			const complex* VxcTildeData = VxcTilde->data();
			//Prepare kinetic potential in reciprocal space (if needed):
			complexScalarFieldTilde VtauTilde; const complex* VtauData = 0;
			if(e.exCorr.needsKEdensity() && eVars.Vtau[qnum.index()])
			{	VtauTilde = Complex(J(eVars.Vtau[qnum.index()]));
				VtauData = VtauTilde->data();
			}
			for(int jCur=0; jCur<nColsMine; jCur++)
			{	int j = iColsMine[jCur]; //global column index
				for(int iCur=0; iCur<nRowsMine; iCur++)
				{	int i = iRowsMine[iCur]; //global row index
					//Kinetic (diagonal) contributions:
					if(i == j)
						*Hdata += 0.5*gInfo.GGT.metric_length_squared(qnum.k + iGarr[i]);
					//Potential contributions:
					size_t diffIndex = e.gInfo.fullGindex(iGarr[i] - iGarr[j]); //wrapped difference index into potential arrays
					*Hdata += Vdata[diffIndex];
					if(VxcData) *VxcData += VxcTildeData[diffIndex]; //Hartree etc. not included here
					if(VtauData)
					{	complex VtauCur = VtauData[diffIndex] * (0.5*dot(iGarr[i]+qnum.k, gInfo.GGT * (iGarr[j]+qnum.k)));
						*Hdata += VtauCur;
						if(VxcData) *VxcData += VtauCur;
					}
					//Increment pointers:
					Hdata++;
					if(VxcData) VxcData++;
				}
			}
		}
		//--- Nonlocal pseudopotential contributions:
		for(const auto& sp: e.iInfo.species)
		{	int nAtoms = sp->atpos.size();
			int nProj = sp->MnlAll.nRows(); //number of projectors per atom
			if(!nAtoms || !nProj) continue; //unused species or purely local psp
			const std::shared_ptr<ColumnBundle> Vrow = sp->getV(Crow);
			const std::shared_ptr<ColumnBundle> Vcol = sp->getV(Ccol);
			matrix Vrow_a(Vrow->colLength(), nProj);
			matrix Vcol_a(Vcol->colLength(), nProj);
			for(int a=0; a<nAtoms; a++)
			{	callPref(eblas_copy)(Vrow_a.dataPref(), Vrow->dataPref() + a*Vrow_a.nData(), Vrow_a.nData());
				callPref(eblas_copy)(Vcol_a.dataPref(), Vcol->dataPref() + a*Vcol_a.nData(), Vcol_a.nData());
				H += (1./e.gInfo.detR) * Vrow_a * sp->MnlAll * dagger(Vcol_a);
			}
			sp->sync_atpos(); //free cached projectors
		}
		//--- DFT+U contributions:
		const matrix* U_rhoPtr = e.eVars.U_rhoAtom.data();
		for(const auto& sp: e.iInfo.species)
		{	int nAtoms = sp->atpos.size();
			if(!nAtoms || !sp->plusU.size()) continue; //unused species or purely local psp
			for(unsigned iU=0; iU<sp->plusU.size(); iU++)
			{	const SpeciesInfo::PlusU& Uparams = sp->plusU[iU];
				int orbCount = (2*Uparams.l+1) * nSpinor;
				ColumnBundle OpsiRow(Crow.similar(orbCount*nAtoms)); sp->setAtomicOrbitals(OpsiRow, true, Uparams.n, Uparams.l);
				ColumnBundle OpsiCol(Ccol.similar(orbCount*nAtoms)); sp->setAtomicOrbitals(OpsiCol, true, Uparams.n, Uparams.l);
				matrix OpsiRow_a(OpsiRow.colLength(), orbCount);
				matrix OpsiCol_a(OpsiCol.colLength(), orbCount);
				for(int s=0; s<qnum.index(); s++) U_rhoPtr += nAtoms;
				for(int a=0; a<nAtoms; a++)
				{	callPref(eblas_copy)(OpsiRow_a.dataPref(), OpsiRow.dataPref() + a*OpsiRow_a.nData(), OpsiRow_a.nData());
					callPref(eblas_copy)(OpsiCol_a.dataPref(), OpsiCol.dataPref() + a*OpsiCol_a.nData(), OpsiCol_a.nData());
					matrix Ucontrib = (1./(e.gInfo.detR*eInfo.spinWeight)) * OpsiRow_a * (*(U_rhoPtr++)) * dagger(OpsiCol_a);
					H += Ucontrib;
					if(Vxc) Vxc += Ucontrib; //DFT+U is logically an ex-corr extension
				}
				for(int s=qnum.index()+1; s<nSpins; s++) U_rhoPtr += nAtoms;
			}
		}
		//--- EXX contributions:
		if(e.exCorr.exxFactor())
		{	matrix HXX = zeroes(nRowsMine, nColsMine);
			e.exx->addHamiltonian(e.exCorr.exxFactor(), e.exCorr.exxRange(), q, HXX, iRowsMine, iColsMine);
			H += HXX;
			if(Vxc) Vxc += HXX;
			//Save for Vxx if needed / applicable:
			if(bgwp.saveVxx and (not e.exCorr.exxRange())) //must be bare exchange
				Vxx = (1./e.exCorr.exxFactor()) * HXX; //remove hybrid scale factor
		}
		if(bgwp.saveVxx and (not Vxx))
		{	Vxx = zeroes(nRowsMine, nColsMine);
			e.exx->addHamiltonian(1., 0., q, Vxx, iRowsMine, iColsMine);
		}
		matrix evecs(nRowsMine, nColsMine);
		diagMatrix eigs(nRows);
		watchSetup.stop();
		
		//Call the scalapack diagonalizer:
		watchDiag.start();
		int eigStart=1, eigStop=nEigs; //find lowest nEigs eigs
		double vlUnused=0., vuUnused=0.; //eigenvalue ranges (unused)
		double absTol = -1.; //use default accuracy
		double orFac = 5e-7; //default orthogonalization threshold
		int nEigsFound, nEvecsFound; //number of calculated eigenvalues and eigenvectors (output)
		int lwork = -1, lrwork = -1, liwork = -1, one = 1;
		ManagedArray<complex> work; 
		ManagedArray<double> rwork, clusterGaps; 
		ManagedArray<int> iwork, iFail, iClusters;
		work.init(1); rwork.init(1); iwork.init(1); //These will be sized after workspace query
		iFail.init(nRows);
		iClusters.init(2*nProcesses);
		clusterGaps.init(nProcesses);
		int info = 0;
		for(int pass=0; pass<2; pass++) //first pass is workspace query, next pass is actual calculation
		{	pzheevx_("V", "I", "U", &nRows, H.data(), &one, &one, descH,
				&vlUnused, &vuUnused, &eigStart, &eigStop, &absTol, &nEigsFound, &nEvecsFound, eigs.data(), &orFac,
				evecs.data(), &one, &one, descH, work.data(), &lwork, rwork.data(), &lrwork,
				iwork.data(), &liwork, iFail.data(), iClusters.data(), clusterGaps.data(), &info);
			if(info < 0)
			{	int errCode = -info;
				if(errCode < 100) die("\tError in argument# %d to pzheevx.\n", errCode)
				else die("\tError in entry %d of argument# %d to pzheevx.\n", errCode%100, errCode/100)
			}
			if(info > 0)
			{	ostringstream err;
				if(info & 0x01) err << "\tSome eigenvectors failed to converge in pzheevx.\n";
				if(info & 0x02) err << "\tSome eigenvectors could not be orthogonalized in pzheevx.\n";
				if(info & 0x04) err << "\tInsufficeint space to compute eigenvectors in pzheevx.\n";
				if(info & 0x08) err << "\tFailed to compute tridiagonal-matrix eigenvalues in pzheevx.\n";
				die("%s", err.str().c_str());
			}
			if(pass) break; //done
			//After first-pass, use results of work-space query to allocate:
			lwork = int(work.data()[0].real()) + bgwp.clusterSize*nRows; work.init(lwork);
			lrwork = int(rwork.data()[0]) + bgwp.clusterSize*nRows; rwork.init(lrwork);
			liwork = int(iwork.data()[0]); iwork.init(liwork);
		}
		watchDiag.stop();
		
		//Select needed eigenvalues:
		eigs.resize(nEigs); //rest zero
		if(eInfo.isMine(q))
		{	E[q] = eigs;
			F[q].resize(nBands, 0.); //update to nBandsDense (padded with zeroes)
		}
		
		//Write the wavefunctions:
		watchIO.start();
		//--- Switch to column-split form:
		int colBlockSize = ceildiv(nEigs, nProcesses); //such that 1D block cyclic is just column-split
		int descIn[9], descOut[9];
		{	int zero=0, info;
			//Input descriptor: same as descH, but only use first nEigs columns (of evecs)
			descinit_(descIn, &nRows, &nEigs, &blockSize, &blockSize, &zero, &zero, &blacsContext, &nRowsMine, &info);
			assert(info==0);
			//Output descriptor: full rows, and contiguous column blocks:
			descinit_(descOut, &nRows, &nEigs, &nRows, &colBlockSize, &zero, &zero, &blacsContextCol, &nRows, &info);
			assert(info==0);
		}
		int iFullColStart = std::min(mpiWorld->iProcess() * colBlockSize, nEigs);
		int iFullColStop = std::min((mpiWorld->iProcess()+1) * colBlockSize, nEigs);
		int nFullColsMine = iFullColStop - iFullColStart;
		matrix buf(nRows, nFullColsMine);
		pzgemr2d_(&nRows, &nEigs,
			evecs.data(), &one, &one, descIn,
			buf.data(), &one, &one, descOut,
			&blacsContext);
		//Select destination in hdf5 file and wwrite from buffer:
		hsize_t offset[4] = { hsize_t(iFullColStart), hsize_t(iSpin), hsize_t(nBasisPrev[ik]), 0 };
		hsize_t count[4] = { hsize_t(nFullColsMine), 1, hsize_t(nRows), 2 };
		sid = H5Dget_space(did);
		H5Sselect_hyperslab(sid, H5S_SELECT_SET, offset, NULL, count, NULL);
		hid_t sidMem = H5Screate_simple(4, count, NULL);
		H5Dwrite(did, H5T_NATIVE_DOUBLE, sidMem, sid, plid, buf.data());
		H5Sclose(sidMem);
		H5Fflush(gidWfns, H5F_SCOPE_GLOBAL);
		buf = 0; //cleanup
		watchIO.stop();
		
		//Truncate eigenvectors down to those needed for V output:
		
		//Transform Vxc to eigenbasis:
		if(bgwp.saveVxc)
			transformV(q, Vxc, evecs, VxcSub, eInfo, descH, nProcsRow, nProcsCol,
				nRows, nBandsV, nRowsMine, nColsMine, blockSize, iProcRow, iProcCol);
		
		if(bgwp.saveVxx)
			//Transform Vxx to eigenbasis:
			transformV(q, Vxx, evecs, VxxSub, eInfo, descH, nProcsRow, nProcsCol,
				nRows, nBandsV, nRowsMine, nColsMine, blockSize, iProcRow, iProcCol);
	}
	H5Pclose(plid);
	H5Dclose(did);
	
	//Update fillings if necessary:
	if(eInfo.fillingsUpdate == ElecInfo::FillingsHsub)
	{	double Bz, mu = eInfo.findMu(E, eInfo.nElectrons, Bz);
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
			F[q] = eInfo.smear(eInfo.muEff(mu,Bz,q), E[q]);
		logPrintf("\t"); eInfo.smearReport();
	}
	logPrintf("\t");
}

#endif //SCALAPACK_ENABLED
#endif //HDF5_ENABLED
