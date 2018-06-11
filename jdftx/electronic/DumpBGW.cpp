/*-------------------------------------------------------------------
Copyright 2017 Ravishankar Sundararaman

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
#include <electronic/Dump_internal.h>

#ifndef HDF5_ENABLED
void Dump::dumpBGW()
{	die("BerkeleyGW output requires HDF5 support.\n");
}
#else

#include <core/H5io.h>

//! Helper class for DumpBGW
class BGW
{
	const Everything& e;
	const BGWparams& bgwp;
	const GridInfo& gInfo;
	const ElecInfo& eInfo;
	const ElecVars& eVars;
	const int nSpins; //!< number of spin channels
	const int nSpinor; //!< number of spinor component
	const int nReducedKpts; //!< number of reduced k-points
	std::vector<int> nBasis; //!< number of basis functions per reduced k-point
	std::vector<int> nBasisPrev; //!< cumulative number of basis functions before reduced k-point
	int nBasisMax; //!< maximum number of basis functions for any reduced k-point
	std::vector<vector3<>> k; //!< k-points in BGW convention [0,1)
	std::vector<vector3<int>> kOffset; //!< k offsets to switch from JDFTx to BGW convention
	std::vector<double> wk; //!< k-point weights
	
	int nBands; //!< eInfo.nBands, or overridden by nBandsDense
	std::vector<diagMatrix> E, F; //!< eigenvalues and fillings
	std::vector<matrix> VxcSub; //!< exchange-correlation matrix elements
	
public:
	BGW(const Everything& e, const BGWparams& bgwp);
	void writeWfn() const; //!< Write wavefunction file
	void denseWriteWfn(hid_t gidWfns); //!< Solve wavefunctions using ScaLAPACK and write to hdf5 file
	void writeVxc() const; //!< Write exchange-correlation matrix elements
	void writeChiFluid() const; //!< Write fluid polarizability
private:
	hid_t openHDF5(string fname) const; //!< Open HDF5 file for collective access
	void writeHeaderMF(hid_t fid) const; //!< Write common HDF5 header specifying the mean-field claculation for BGW outputs
};


void Dump::dumpBGW()
{	if(!bgwParams) bgwParams = std::make_shared<BGWparams>(); //default parameters
	BGW bgw(*e, *bgwParams);
	bgw.writeWfn();
	bgw.writeVxc();
	if(e->eVars.fluidSolver)
		bgw.writeChiFluid();
}


//------------ Implementation of class DumpBGW -----------------

BGW::BGW(const Everything& e, const BGWparams& bgwp)
: e(e), bgwp(bgwp), gInfo(e.gInfo), eInfo(e.eInfo), eVars(e.eVars),
nSpins(eInfo.nSpins()), nSpinor(eInfo.spinorLength()), nReducedKpts(eInfo.nStates/nSpins)
{
	//nBasis arrays:
	nBasis.assign(nReducedKpts, 0);
	nBasisPrev.assign(nReducedKpts, 0);
	for(int q=0; q<nReducedKpts; q++)
	{	nBasis[q] = e.basis[q].nbasis;
		if(q+1<nReducedKpts) nBasisPrev[q+1] = nBasisPrev[q] + nBasis[q];
	}
	nBasisMax = *std::max_element(nBasis.begin(), nBasis.end());
	
	//k-points:
	k.resize(nReducedKpts);
	kOffset.resize(nReducedKpts);
	wk.resize(nReducedKpts);
	for(int q=0; q<nReducedKpts; q++)
	{	k[q] = eInfo.qnums[q].k;
		//Switch k to BGW convention of [0,1)
		for(int iDir=0; iDir<3; iDir++)
		{	kOffset[q][iDir] = int(floor(k[q][iDir]));
			k[q][iDir] -= kOffset[q][iDir];
		}
		wk[q] = eInfo.qnums[q].weight / eInfo.spinWeight; //Set sum(wk) = 1 in all cases
	}
	
	//DFT nBands; corresponding E and F:
	nBands = eInfo.nBands;
	E = eVars.Hsub_eigs;
	F = eVars.F;
	VxcSub.resize(eInfo.nStates);
}


//Write wavefunction file for BGW
void BGW::writeWfn() const
{	//Open file:
	hid_t fid = openHDF5(e.dump.getFilename("bgw.wfn.h5"));
	
	//Wavefunction group:
	hid_t gidWfns = h5createGroup(fid, "wfns");
	//--- G-vectors:
	std::vector<vector3<int>> iGarr;
	for(int q=0; q<nReducedKpts; q++)
		for(const vector3<int>& iG: e.basis[q].iGarr)
			iGarr.push_back(iG + kOffset[q]); //shift iG to account for change in k convention from JDFTx to BGW
	hsize_t dimsGwfns[2] = { iGarr.size(), 3 };
	h5writeVector(gidWfns, "gvecs", &iGarr[0][0], dimsGwfns, 2);
	//--- Coefficients:
	if(bgwp.nBandsDense)
	{	//Write results of a ScaLAPACK solve:
		((BGW*)this)->denseWriteWfn(gidWfns);
	}
	else //Default output of bands from usual totalE / bandstructure calculation
	{	//Create dataset (must happen on all processes together):
		hsize_t dims[4] = { hsize_t(nBands), hsize_t(nSpins*nSpinor), iGarr.size(), 2 };
		hid_t sid = H5Screate_simple(4, dims, NULL);
		hid_t did = H5Dcreate(gidWfns, "coeffs", H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		hid_t plid = H5Pcreate(H5P_DATASET_XFER);
		H5Sclose(sid);
		//Loop over k, bands and spin/spinors:
		hsize_t offset[4] = { 0, 0, 0, 0 };
		hsize_t count[4] = { 1, 1, 1, 2 };
		std::vector<complex> buffer(*std::max_element(nBasis.begin(), nBasis.end()));
		double volScaleFac = sqrt(gInfo.detR);
		for(int iSpin=0; iSpin<nSpins; iSpin++)
		for(int iSpinor=0; iSpinor<nSpinor; iSpinor++)
		{	offset[1] = iSpin*nSpinor + iSpinor;
			for(int ik=0; ik<nReducedKpts; ik++)
			{	int q=iSpin*nReducedKpts+ik;
				if(!eInfo.isMine(q)) continue;
				count[2] = nBasis[ik];
				offset[2] = nBasisPrev[ik];
				hid_t sidMem = H5Screate_simple(4, count, NULL);
				for(int b=0; b<nBands; b++)
				{	offset[0] = b;
					sid = H5Dget_space(did);
					H5Sselect_hyperslab(sid, H5S_SELECT_SET, offset, NULL, count, NULL);
					//Copy to buffer and scale:
					eblas_copy(buffer.data(), eVars.C[q].data()+eVars.C[q].index(b, iSpinor*nBasis[ik]), nBasis[ik]);
					eblas_zdscal(nBasis[ik], volScaleFac, buffer.data(), 1);
					//Write buffer to HDF5:
					H5Dwrite(did, H5T_NATIVE_DOUBLE, sidMem, sid, plid, buffer.data());
				}
				H5Sclose(sidMem);
			}
		}
		H5Pclose(plid);
		H5Dclose(did);
	}
	H5Gclose(gidWfns);
	
	//Write common header at end (so as to use updated eigenvalues and fillings from scalapack solve, if any)
	writeHeaderMF(fid);
	
	//Close file:
	H5Fclose(fid);
	logPrintf("Done.\n"); logFlush();
}


//Write exchange-correlation matrix elements for BGW
void BGW::writeVxc() const
{
	//Open file:
	string fname = e.dump.getFilename("bgw.vxc.dat");
	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
	
	//Calculate X-C matrix elements:
	std::vector<matrix>& VxcSub = ((BGW*)this)->VxcSub;
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
	{	if(VxcSub[q]) continue; //already calculated (dense version)
		ColumnBundle HCq = gInfo.dV * Idag_DiagV_I(eVars.C[q], eVars.Vxc);
		if(e.exCorr.needsKEdensity() && eVars.Vtau[eInfo.qnums[q].index()]) //metaGGA KE potential
		{	for(int iDir=0; iDir<3; iDir++)
				HCq -= (0.5*gInfo.dV) * D(Idag_DiagV_I(D(eVars.C[q],iDir), eVars.Vtau), iDir);
		}
		if(e.eInfo.hasU) //Contribution via atomic density matrix projections (DFT+U)
			e.iInfo.rhoAtom_grad(eVars.C[q], eVars.U_rhoAtom, HCq);
		VxcSub[q] = eVars.C[q] ^ HCq;
	}
	
	//Output from head
	if(mpiWorld->isHead())
	{	FILE* fp = fopen(fname.c_str(), "w");
		if(!fp) die_alone("failed to open for writing.\n");
		for(int ik=0; ik<nReducedKpts; ik++)
		{	fprintf(fp, "%.9f %.9f %.9f %4d %4d\n", k[ik][0], k[ik][1], k[ik][2],
				nBands*nSpins, nBands*nBands*nSpins);
			for(int iSpin=0; iSpin<nSpins; iSpin++)
			{	int q=iSpin*nReducedKpts+ik;
				if(!eInfo.isMine(q))
				{	VxcSub[q] = zeroes(nBands, nBands);
					VxcSub[q].recv(eInfo.whose(q), q);
				}
				//Diagonal elements:
				for(int b=0; b<nBands; b++)
				{	complex V_eV = VxcSub[q](b,b) / eV; //convert to eV
					fprintf(fp, "%4d %4d %+14.9f %+14.9f\n", iSpin+1, b+1, V_eV.real(), V_eV.imag());
				}
				//Off-diagonal elements:
				for(int b2=0; b2<nBands; b2++)
				{	for(int b1=0; b1<nBands; b1++)
					{	complex V_eV = VxcSub[q](b1,b2) / eV; //convert to eV
						fprintf(fp, "%4d %4d %4d %+14.9f %+14.9f\n", iSpin+1, b1+1, b2+1, V_eV.real(), V_eV.imag());
					}
				}
				if(!eInfo.isMine(q)) VxcSub[q] = 0; //cleanup temporary memory
			}
		}
		fclose(fp);
	}
	else //send ones tored on other processes to head
	{	for(int q=0; q<eInfo.nStates; q++)
			if(eInfo.isMine(q))
				VxcSub[q].send(0, q);
	}
	logPrintf("Done.\n");
}


//Write fluid polarizability
void BGW::writeChiFluid() const
{	assert(eVars.fluidSolver); //should only be called in solvated cases
	
	//Open file and write common header:
	hid_t fid = openHDF5(e.dump.getFilename("bgw.chiFluid.h5"));
	writeHeaderMF(fid);
	
	//---------- Epsilon header ----------
	hid_t gidHeader = h5createGroup(fid, "eps_header");
	h5writeScalar(gidHeader, "versionnumber", 1);
	h5writeScalar(gidHeader, "flavor", 2);  //=> complex wavefunctions
	
	//----- Params group:
	hid_t gidParams = h5createGroup(gidHeader, "params");
	h5writeScalar(gidParams, "matrix_type", 2); //=> chi
	h5writeScalar(gidParams, "has_advanced", 0); //=> retarded response only
	h5writeScalar(gidParams, "nmatrix", nSpins); //=> polarizability per spin
	h5writeScalar(gidParams, "matrix_flavor", 2); //=> complex matrices
	h5writeScalar(gidParams, "icutv", 0); //TODO: what is this?
	h5writeScalar(gidParams, "ecuts", e.cntrl.Ecut); //TODO: what Ecut should I use here?
	h5writeScalar(gidParams, "nband", 0); //No bands used for fluid polarizability
	h5writeScalar(gidParams, "efermi", 0.); //Not meaningful for fluid polarizability
	
	//---- q-points group:
	die("Not yet implemented!\n"); //TODO
	
	//Close file:
	H5Fclose(fid);
	logPrintf("Done.\n"); logFlush();
}


//Open HDF5 file for collective access:
hid_t BGW::openHDF5(string fname) const
{	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
	//Create MPI file access across all processes:
	hid_t plid = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plid, MPI_COMM_WORLD, MPI_INFO_NULL);
	//Open file with MPI access:
	hid_t fid = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plid);
	if(fid<0) die("Could not open/create output HDF5 file '%s'\n", fname.c_str());
	H5Pclose(plid);
	return fid;
}


//Write common HDF5 header specifying the mean-field claculation for BGW outputs
void BGW::writeHeaderMF(hid_t fid) const
{
	hid_t gidHeader = h5createGroup(fid, "mf_header");
	h5writeScalar(gidHeader, "versionnumber", 1);
	h5writeScalar(gidHeader, "flavor", 2);  //=> complex wavefunctions
	
	//---------- kpoint group ----------
	hid_t gidKpts = h5createGroup(gidHeader, "kpoints");
	vector3<> kShift = Diag(eInfo.qnums[0].k) * eInfo.kFoldingCount(); //k-shift used before folding
	h5writeScalar(gidKpts, "nspin", nSpins);
	h5writeScalar(gidKpts, "nspinor", nSpinor);
	h5writeScalar(gidKpts, "nrk", nReducedKpts);
	h5writeScalar(gidKpts, "mnband", nBands);
	h5writeScalar(gidKpts, "ngkmax", nBasisMax);
	h5writeScalar(gidKpts, "ecutwfc", e.cntrl.Ecut/Ryd);
	h5writeVector(gidKpts, "kgrid", &eInfo.kFoldingCount()[0], 3);
	h5writeVector(gidKpts, "shift", &kShift[0], 3);
	h5writeVector(gidKpts, "ngk", nBasis);
	
	//--- occupied band ranges:
	const double Fcut = 0.5;
	std::vector<int> ifmin(eInfo.nStates, 1), ifmax(eInfo.nStates, 1);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
	{	for(int b=0; b<nBands; b++)
		{	if(F[q][b] < Fcut) break;
			ifmax[q] = b+1;
		}
	}
	mpiWorld->allReduce(ifmax.data(), ifmax.size(), MPIUtil::ReduceMax);
	hsize_t dimsKspin[2] = { hsize_t(nSpins), hsize_t(nReducedKpts) };
	h5writeVector(gidKpts, "ifmin", ifmin.data(), dimsKspin, 2);
	h5writeVector(gidKpts, "ifmax", ifmax.data(), dimsKspin, 2);
	
	//--- kpoints and weights:
	hsize_t dimsK[2] = { hsize_t(nReducedKpts), 3 };
	h5writeVector(gidKpts, "rk", &k[0][0], dimsK, 2);
	h5writeVector(gidKpts, "w", wk);
	
	//--- eigenvalues and occupations:
	std::vector<double> Eall, Fall;
	Eall.reserve(eInfo.nStates*nBands);
	Fall.reserve(eInfo.nStates*nBands);
	for(int q=0; q<eInfo.nStates; q++)
	{	diagMatrix Ecur(nBands), Fcur(nBands);
		if(eInfo.isMine(q)) { Ecur = E[q]*(1./Ryd); Fcur = F[q]; }
		Ecur.bcast(eInfo.whose(q)); Eall.insert(Eall.end(), Ecur.begin(), Ecur.end());
		Fcur.bcast(eInfo.whose(q)); Fall.insert(Fall.end(), Fcur.begin(), Fcur.end());
	}
	hsize_t dimsKspinBands[3] = { hsize_t(nSpins), hsize_t(nReducedKpts), hsize_t(nBands) };
	h5writeVector(gidKpts, "el", Eall.data(), dimsKspinBands, 3);
	h5writeVector(gidKpts, "occ", Fall.data(), dimsKspinBands, 3);
	Eall.clear(); Fall.clear();
	H5Gclose(gidKpts);
	
	//---------- G-space group ----------
	hid_t gidGspace = h5createGroup(gidHeader, "gspace");
	std::vector<vector3<int>> iGarr(gInfo.nr);
	{	const vector3<int>& S = gInfo.S;
		size_t iStart = 0, iStop = gInfo.nr;
		THREAD_fullGspaceLoop(
			iGarr[i] = iG;
			for(int iDir=0; iDir<3; iDir++)
				if(2*iGarr[i][iDir]==gInfo.S[iDir])
					iGarr[i][iDir]-=gInfo.S[iDir]; //[-S/2,S/2) in BGW (rather than (-S/2,S/2] in JDFTx)
		)
	}
	hsize_t dimsG[2] = { hsize_t(gInfo.nr), 3 };
	h5writeScalar(gidGspace, "ng", gInfo.nr);
	h5writeScalar(gidGspace, "ecutrho", std::max(e.cntrl.EcutRho, 4*e.cntrl.Ecut)/Ryd);
	h5writeVector(gidGspace, "FFTgrid", &gInfo.S[0], 3);
	h5writeVector(gidGspace, "components", &iGarr[0][0], dimsG, 2);
	H5Gclose(gidGspace);
	
	//---------- Symmetries group ----------
	hid_t gidSymm = h5createGroup(gidHeader, "symmetry");
	const std::vector<SpaceGroupOp> ops = e.symm.getMatrices();
	std::vector<matrix3<int> > rots(ops.size());
	std::vector<vector3<> > trans(ops.size());
	for(size_t iSym=0; iSym<ops.size(); iSym++)
	{	matrix3<int> rotInv = det(ops[iSym].rot) * adjugate(ops[iSym].rot); //since |det(rot)| = 1
		rots[iSym] = rotInv; //BGW uses inverse convention
		trans[iSym] = (2*M_PI)*ops[iSym].a; //BGW used 2*pi times translation
	}
	hsize_t dimsRot[3] = { ops.size(), 3, 3 };
	hsize_t dimsTrans[2] = { ops.size(), 3 };
	h5writeScalar(gidSymm, "ntran", int(ops.size()));
	h5writeScalar(gidSymm, "cell_symmetry", 0);
	h5writeVector(gidSymm, "mtrx", &rots[0](0,0), dimsRot, 3);
	h5writeVector(gidSymm, "tnp", &trans[0][0], dimsTrans, 2);
	H5Gclose(gidSymm);
	
	//---------- Crystal group ----------
	hid_t gidCrystal = h5createGroup(gidHeader, "crystal");
	hsize_t dims33[2] = { 3, 3 };
	h5writeScalar(gidCrystal, "celvol", gInfo.detR);
	h5writeScalar(gidCrystal, "recvol", fabs(det(gInfo.G)));
	h5writeScalar(gidCrystal, "alat", 1);
	h5writeScalar(gidCrystal, "blat", 1);
	h5writeVector(gidCrystal, "avec", &gInfo.RT(0,0), dims33, 2); //BGW lattice vectors in rows
	h5writeVector(gidCrystal, "bvec", &gInfo.G(0,0), dims33, 2); //BGW recip lattice vectors in rows
	h5writeVector(gidCrystal, "adot", &gInfo.RTR(0,0), dims33, 2);
	h5writeVector(gidCrystal, "bdot", &gInfo.GGT(0,0), dims33, 2);
	//--- collect atoms:
	std::vector<vector3<>> apos; std::vector<int> atyp;
	for(const auto& sp: e.iInfo.species)
	{	apos.insert(apos.end(), sp->atpos.begin(), sp->atpos.end());
		atyp.insert(atyp.end(), sp->atpos.size(), sp->atomicNumber);
	}
	hsize_t dimsApos[2] = { apos.size(), 3 };
	h5writeScalar(gidCrystal, "nat", int(apos.size()));
	h5writeVector(gidCrystal, "atyp", atyp);
	h5writeVector(gidCrystal, "apos", &apos[0][0], dimsApos, 2);
	H5Gclose(gidCrystal);

	H5Gclose(gidHeader);
}


//--------------- ScaLAPACK solve of bands -------------
#ifdef SCALAPACK_ENABLED
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
	
	void zlacpy_(const char* uplo, const int* m, const int* n, const complex* a, const int* lda, complex* b, const int* ldb);
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
#endif

//Solve wavefunctions using ScaLAPACK and write to hdf5 file:
void BGW::denseWriteWfn(hid_t gidWfns)
{	static StopWatch watchDiag("scalapackDiagonalize"), watchSetup("scalapackMatrixSetup"), watchVxc("scalapackVxcTransform");
	logPrintf("\n");
#ifdef SCALAPACK_ENABLED
	nBands = bgwp.nBandsDense;
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
	int blacsContext, iProcRow, iProcCol;
	{	int unused=-1, what=0;
		blacs_get_(&unused, &what, &blacsContext);
		blacs_gridinit_(&blacsContext, "Row-major", &nProcsRow, &nProcsCol);
		blacs_gridinfo_(&blacsContext, &nProcsRow, &nProcsCol, &iProcRow, &iProcCol);
		assert(mpiWorld->iProcess() == iProcRow * nProcsCol + iProcCol); //this mapping is assumed below, so check
	}
	logPrintf("\tInitialized %d x %d process BLACS grid.\n", nProcsRow, nProcsCol);
	
	//Create dataset (must happen on all processes together):
	hsize_t nGtot = nBasisPrev.back() + nBasis.back();
	hsize_t dims[4] = { hsize_t(nBands), hsize_t(nSpins*nSpinor), nGtot, 2 };
	hid_t sid = H5Screate_simple(4, dims, NULL);
	hid_t did = H5Dcreate(gidWfns, "coeffs", H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t plid = H5Pcreate(H5P_DATASET_XFER);
	H5Sclose(sid);
	
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
		logPrintf(" with dimension %d\n", nRows); logFlush();
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
		matrix Vxc = zeroes(nRowsMine, nColsMine);
		//--- Kinetic, potential and kinetic-potential contributions:
		{	complex* Hdata = H.data();
			complex* VxcData = Vxc.data();
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
					*VxcData += VxcTildeData[diffIndex]; //Hartree etc. not included here
					if(VtauData)
					{	complex VtauCur = VtauData[diffIndex] * (0.5*dot(iGarr[i]+qnum.k, gInfo.GGT * (iGarr[j]+qnum.k)));
						*Hdata += VtauCur;
						*VxcData += VtauCur;
					}
					//Increment pointers:
					Hdata++;
					VxcData++;
				}
			}
		}
		//--- Nonlocal pseudopotential contributions:
		for(const auto& sp: e.iInfo.species)
		{	int nAtoms = sp->atpos.size();
			int nProj = sp->MnlAll.nRows(); //number of projectors per atom
			if(!nAtoms || !nProj) continue; //unused species or purely local psp
			const ColumnBundle& Vrow = *(sp->getV(Crow));
			const ColumnBundle& Vcol = *(sp->getV(Ccol));
			matrix Vrow_a(Vrow.colLength(), nProj);
			matrix Vcol_a(Vcol.colLength(), nProj);
			for(int a=0; a<nAtoms; a++)
			{	callPref(eblas_copy)(Vrow_a.dataPref(), Vrow.dataPref() + a*Vrow_a.nData(), Vrow_a.nData());
				callPref(eblas_copy)(Vcol_a.dataPref(), Vcol.dataPref() + a*Vcol_a.nData(), Vcol_a.nData());
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
					Vxc += Ucontrib; //DFT+U is logically an ex-corr extension
				}
				for(int s=qnum.index()+1; s<nSpins; s++) U_rhoPtr += nAtoms;
			}
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
			{	string err;
				if(info & 0x01) err += "\tSome eigenvectors failed to converge in pzheevx.\n";
				if(info & 0x02) err += "\tSome eigenvectors could not be orthogonalized in pzheevx.\n";
				if(info & 0x04) err += "\tInsufficeint space to compute eigenvectors in pzheevx.\n";
				if(info & 0x08) err += "\tFailed to compute tridiagonal-matrix eigenvalues in pzheevx.\n";
				die("%s", err.c_str());
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
		std::vector<int> iEigColsMine = iColsMine;
		for(int jMine=0; jMine<nColsMine; jMine++)
			if(iEigColsMine[jMine] >= nEigs)
			{	iEigColsMine.resize(jMine); //remaining columns irrelevant
				break;
			}
		int nEigColsMine = iEigColsMine.size();
		if(eInfo.isMine(q))
		{	E[q] = eigs;
			F[q].resize(nBands, 0.); //update to nBandsDense (padded with zeroes)
		}
		
		//Write the wavefunctions:
		hsize_t offset[4] = { 0, hsize_t(iSpin), 0, 0 };
		hsize_t count[4] = { 1, 1, 1, 2 };
		matrix buf(blockSize, blockSize);
		int nBlockRows = (nRowsMine+blockSize-1)/blockSize;
		int nBlockCols = (nEigColsMine+blockSize-1)/blockSize;
		for(int jBlock=0; jBlock<nBlockCols; jBlock++)
		{	int jStart = jBlock*blockSize; //start index in local matrix
			int jCount = std::min(blockSize, nEigColsMine-jStart); //number
			int jOffset = iEigColsMine[jStart]; //start index in global matrix
			for(int iBlock=0; iBlock<nBlockRows; iBlock++)
			{	int iStart = iBlock*blockSize; //start index in local matrix
				int iCount = std::min(blockSize, nRowsMine-iStart); //number
				int iOffset = iRowsMine[iStart]; //start index in global matrix
				//Copy block from local matrix to contiguous buffer:
				zlacpy_("A", &iCount, &jCount, evecs.data()+nRowsMine*jStart+iStart, &nRowsMine, buf.data(), &iCount);
				//Select destination in hdf5 file and wwrite from buffer:
				count[0] = jCount; offset[0] = jOffset; //band index
				count[2] = iCount; offset[2] = iOffset+nBasisPrev[ik]; //G-vector index
				sid = H5Dget_space(did);
				H5Sselect_hyperslab(sid, H5S_SELECT_SET, offset, NULL, count, NULL);
				hid_t sidMem = H5Screate_simple(4, count, NULL);
				H5Dwrite(did, H5T_NATIVE_DOUBLE, sidMem, sid, plid, buf.data());
				H5Sclose(sidMem);
			}
		}
		
		//Transform Vxc to eigenbasis:
		watchVxc.start();
		matrix VxcEvecs = zeroes(nRowsMine, nColsMine); //local part of Vxc * evecs
		//--- parallel matrix multiplies:
		complex alpha(1.,0.), beta(0.,0.);
		pzgemm_("N", "N", &nRows, &nRows, &nRows, &alpha,
			Vxc.data(), &one, &one, descH,
			evecs.data(), &one, &one, descH, &beta,
			VxcEvecs.data(), &one, &one, descH); //VxcEvecs = Vxc * evecs
		pzgemm_("C", "N", &nRows, &nRows, &nRows, &alpha,
			evecs.data(), &one, &one, descH,
			VxcEvecs.data(), &one, &one, descH, &beta,
			Vxc.data(), &one, &one, descH); //VxcNew = evecx' * VxcEvecs
		VxcEvecs = 0; //cleanup memory
		
		//Collect required portion of Vxc on a single process (the one that owns state q):
		if(eInfo.isMine(q))
		{	VxcSub[q] = zeroes(nEigs, nEigs);
			for(int jProcRow=0; jProcRow<nProcsRow; jProcRow++)
				for(int jProcCol=0; jProcCol<nProcsCol; jProcCol++)
				{	int jProcess = jProcRow * nProcsCol + jProcCol;
					//Determine indices within nEigs for this block:
					std::vector<int> jEigRowsMine = distributedIndices(nEigs, blockSize, jProcRow, nProcsRow);
					std::vector<int> jEigColsMine = distributedIndices(nEigs, blockSize, jProcCol, nProcsCol);
					
					if((!jEigRowsMine.size()) or (!jEigColsMine.size()))
						continue; //nothing from this block
					//Get current block: locally or from another process
					matrix VxcCur;
					if(jProcess == mpiWorld->iProcess())
						VxcCur = Vxc(0,jEigRowsMine.size(), 0,jEigColsMine.size()); //local
					else
					{	VxcCur.init(jEigRowsMine.size(), jEigColsMine.size());
						VxcCur.recv(jProcess);
					}
					//Distribute to full matrix:
					const complex* inData = VxcCur.data();
					for(int c: jEigColsMine)
						for(int r: jEigRowsMine)
							VxcSub[q].set(r, c, *(inData++));
				}
		}
		else
		{	std::vector<int> iEigRowsMine = distributedIndices(nEigs, blockSize, iProcRow, nProcsRow);
			if(iEigRowsMine.size() and iEigColsMine.size())
				Vxc(0,iEigRowsMine.size(), 0,iEigColsMine.size()).send(eInfo.whose(q));
		}
		watchVxc.stop();
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
#else
	die("\nDense solve for extra bands for BGW (nBandsDense > 0) requires ScaLAPACK.\n\n");
#endif
}

#endif
