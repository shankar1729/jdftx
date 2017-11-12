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
	
public:
	BGW(const Everything& e);
	void writeWfn() const; //!< Write wavefunction file
	void writeVxc() const; //!< Write exchange-correlation matrix elements
	void writeChiFluid() const; //!< Write fluid polarizability
private:
	hid_t openHDF5(string fname) const; //!< Open HDF5 file for collective access
	void writeHeaderMF(hid_t fid) const; //!< Write common HDF5 header specifying the mean-field claculation for BGW outputs
};


void Dump::dumpBGW()
{	BGW bgw(*e);
	bgw.writeWfn();
	bgw.writeVxc();
	if(e->eVars.fluidSolver)
		bgw.writeChiFluid();
}


//------------ Implementation of class DumpBGW -----------------

BGW::BGW(const Everything& e)
: e(e), gInfo(e.gInfo), eInfo(e.eInfo), eVars(e.eVars),
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
}


//Write wavefunction file for BGW
void BGW::writeWfn() const
{	//Open file and write common header:
	hid_t fid = openHDF5(e.dump.getFilename("bgw.wfn.h5"));
	writeHeaderMF(fid);
	
	//Wavefucntion group:
	hid_t gidWfns = h5createGroup(fid, "wfns");
	//--- G-vectors:
	std::vector<vector3<int>> iGarr;
	for(int q=0; q<nReducedKpts; q++)
		for(const vector3<int>& iG: e.basis[q].iGarr)
			iGarr.push_back(iG + kOffset[q]); //shift iG to account for change in k convention from JDFTx to BGW
	hsize_t dimsGwfns[2] = { iGarr.size(), 3 };
	h5writeVector(gidWfns, "gvecs", &iGarr[0][0], dimsGwfns, 2);
	//--- Coefficients:
	{	//Create dataset (must happen on all processes together):
		hsize_t dims[4] = { hsize_t(eInfo.nBands), hsize_t(nSpins*nSpinor), iGarr.size(), 2 };
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
				for(int b=0; b<eInfo.nBands; b++)
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
		
	//Close file:
	H5Fclose(fid);
	logPrintf("done\n"); logFlush();
}


//Write exchange-correlation matrix elements for BGW
void BGW::writeVxc() const
{
	//Open file:
	string fname = e.dump.getFilename("bgw.vxc.dat");
	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
	
	//Calculate X-C matrix elements:
	std::vector<matrix> VxcSub(eInfo.nStates);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
	{	ColumnBundle HCq = gInfo.dV * Idag_DiagV_I(eVars.C[q], eVars.Vxc);
		if(e.exCorr.needsKEdensity() && eVars.Vtau[eInfo.qnums[q].index()]) //metaGGA KE potential
		{	for(int iDir=0; iDir<3; iDir++)
				HCq -= (0.5*gInfo.dV) * D(Idag_DiagV_I(D(eVars.C[q],iDir), eVars.Vtau), iDir);
		}
		if(e.eInfo.hasU) //Contribution via atomic density matrix projections (DFT+U)
			e.iInfo.rhoAtom_grad(eVars.C[q], eVars.U_rhoAtom, HCq);
		VxcSub[q] = eVars.C[q] ^ HCq;
	}
	
	//Make available on all processes
	for(int q=0; q<eInfo.nStates; q++)
	{	if(!eInfo.isMine(q)) VxcSub[q] = zeroes(eInfo.nBands, eInfo.nBands);
		VxcSub[q].bcast(eInfo.whose(q));
	}
	
	//Output from head
	if(mpiUtil->isHead())
	{	FILE* fp = fopen(fname.c_str(), "w");
		if(!fp) die_alone("failed to open for writing.\n");
		for(int ik=0; ik<nReducedKpts; ik++)
		{	fprintf(fp, "%.9f %.9f %.9f %4d %4d\n", k[ik][0], k[ik][1], k[ik][2], eInfo.nBands*nSpins, 0);
			for(int iSpin=0; iSpin<nSpins; iSpin++)
			{	int q=iSpin*nReducedKpts+ik;
				for(int b=0; b<eInfo.nBands; b++)
				{	complex V_eV = VxcSub[q](b,b) / eV; //convert to eV
					fprintf(fp, "%4d %4d %+14.9f %+14.9f\n", iSpin+1, b+1, V_eV.real(), V_eV.imag());
				}
			}
		}
		fclose(fp);
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
	logPrintf("done\n"); logFlush();
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
	h5writeScalar(gidKpts, "mnband", eInfo.nBands);
	h5writeScalar(gidKpts, "ngkmax", nBasisMax);
	h5writeScalar(gidKpts, "ecutwfc", e.cntrl.Ecut/Ryd);
	h5writeVector(gidKpts, "kgrid", &eInfo.kFoldingCount()[0], 3);
	h5writeVector(gidKpts, "shift", &kShift[0], 3);
	h5writeVector(gidKpts, "ngk", nBasis);
	
	//--- occupied band ranges:
	const double Fcut = 1e-6;
	std::vector<int> ifmin(eInfo.nStates, 1), ifmax(eInfo.nStates, 1);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
	{	for(int b=0; b<eInfo.nBands; b++)
		{	if(eVars.F[q][b] < Fcut) break;
			ifmax[q] = b+1;
		}
	}
	mpiUtil->allReduce(ifmax.data(), ifmax.size(), MPIUtil::ReduceMax);
	hsize_t dimsKspin[2] = { hsize_t(nSpins), hsize_t(nReducedKpts) };
	h5writeVector(gidKpts, "ifmin", ifmin.data(), dimsKspin, 2);
	h5writeVector(gidKpts, "ifmax", ifmax.data(), dimsKspin, 2);
	
	//--- kpoints and weights:
	hsize_t dimsK[2] = { hsize_t(nReducedKpts), 3 };
	h5writeVector(gidKpts, "rk", &k[0][0], dimsK, 2);
	h5writeVector(gidKpts, "w", wk);
	
	//--- eigenvalues and occupations:
	std::vector<double> Eall, Fall;
	Eall.reserve(eInfo.nStates*eInfo.nBands);
	Fall.reserve(eInfo.nStates*eInfo.nBands);
	for(int q=0; q<eInfo.nStates; q++)
	{	diagMatrix Ecur(eInfo.nBands), Fcur(eInfo.nBands);
		if(eInfo.isMine(q)) { Ecur = eVars.Hsub_eigs[q]*(1./Ryd); Fcur = eVars.F[q]; }
		Ecur.bcast(eInfo.whose(q)); Eall.insert(Eall.end(), Ecur.begin(), Ecur.end());
		Fcur.bcast(eInfo.whose(q)); Fall.insert(Fall.end(), Fcur.begin(), Fcur.end());
	}
	hsize_t dimsKspinBands[3] = { hsize_t(nSpins), hsize_t(nReducedKpts), hsize_t(eInfo.nBands) };
	h5writeVector(gidKpts, "el", Eall.data(), dimsKspinBands, 3);
	h5writeVector(gidKpts, "occ", Fall.data(), dimsKspinBands, 3);
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
	h5writeVector(gidCrystal, "avec", &gInfo.R(0,0), dims33, 2);
	h5writeVector(gidCrystal, "bvec", &gInfo.GT(0,0), dims33, 2);
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

#endif
