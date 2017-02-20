/*-------------------------------------------------------------------
Copyright 2016 Ravishankar Sundararaman

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
#include "Everything.h"

#ifndef HDF5_ENABLED
void Dump::dumpBGW()
{	assert(!"BerkeleyGW output requires HDF5 support.");
}
#else

#include <core/H5io.h>

void Dump::dumpBGW()
{	//Get filename:
	string fname = getFilename("bgw.wfn.h5");
	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
	
	//Prepare data:
	const ElecInfo& eInfo = e->eInfo;
	const ElecVars& eVars = e->eVars;
	int nSpins = eInfo.nSpins();
	int nSpinor = eInfo.spinorLength();
	int nReducedKpts = eInfo.nStates/nSpins;
	vector3<> kShift = Diag(eInfo.qnums[0].k) * eInfo.kFoldingCount(); //k-shift used before folding
	//--- nBasis array:
	std::vector<int> nBasis(nReducedKpts);
	for(int q=0; q<nReducedKpts; q++)
		nBasis[q] = e->basis[q].nbasis;
	int nBasisMax = *std::max_element(nBasis.begin(), nBasis.end());
	
	//Open file:
	hid_t plid = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plid, MPI_COMM_WORLD, MPI_INFO_NULL);
	hid_t fid = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plid);
	if(fid<0) die("Could not open/create output HDF5 file '%s'\n", fname.c_str());
	H5Pclose(plid);
	
	//======== Header ========
	hid_t gidHeader = h5createGroup(fid, "mf_header");
	h5writeScalar(gidHeader, "versionnumber", 1);
	h5writeScalar(gidHeader, "flavor", 2);  //=> complex wavefunctions
	//----- kpoint related -----
	hid_t gidKpts = h5createGroup(gidHeader, "kpoints");
	h5writeScalar(gidKpts, "nspin", nSpins);
	h5writeScalar(gidKpts, "nspinor", nSpinor);
	h5writeScalar(gidKpts, "nrk", nReducedKpts);
	h5writeScalar(gidKpts, "mnband", eInfo.nBands);
	h5writeScalar(gidKpts, "ngkmax", nBasisMax);
	h5writeScalar(gidKpts, "ecutwfc", e->cntrl.Ecut);
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
	std::vector<vector3<>> k(nReducedKpts);
	std::vector<double> wk(nReducedKpts);
	for(int q=0; q<nReducedKpts; q++)
	{	k[q] = eInfo.qnums[q].k;
		wk[q] = eInfo.qnums[q].weight;
	}
	hsize_t dimsK[2] = { hsize_t(nReducedKpts), 3 };
	h5writeVector(gidKpts, "rk", &k[0][0], dimsK, 2);
	h5writeVector(gidKpts, "w", wk);
	//--- eigenvalues and occupations:
	std::vector<double> Eall, Fall;
	Eall.reserve(eInfo.nStates*eInfo.nBands);
	Fall.reserve(eInfo.nStates*eInfo.nBands);
	for(int q=0; q<eInfo.nStates; q++)
	{	diagMatrix Ecur(eInfo.nBands), Fcur(eInfo.nBands);
		if(eInfo.isMine(q)) { Ecur = eVars.Hsub_eigs[q]; Fcur = eVars.F[q]; }
		Ecur.bcast(eInfo.whose(q)); Eall.insert(Eall.end(), Ecur.begin(), Ecur.end());
		Fcur.bcast(eInfo.whose(q)); Fall.insert(Fall.end(), Fcur.begin(), Fcur.end());
	}
	hsize_t dimsKspinBands[3] = { hsize_t(nSpins), hsize_t(nReducedKpts), hsize_t(eInfo.nBands) };
	h5writeVector(gidKpts, "el", Eall.data(), dimsKspinBands, 3);
	h5writeVector(gidKpts, "occ", Fall.data(), dimsKspinBands, 3);
	H5Gclose(gidKpts);
	//----- G-space related -----
	hid_t gidGspace = h5createGroup(gidHeader, "gspace");
	const GridInfo& gInfo = e->gInfo;
	std::vector<vector3<int>> iGarr(gInfo.nr);
	{	const vector3<int>& S = gInfo.S;
		size_t iStart = 0, iStop = gInfo.nr;
		THREAD_fullGspaceLoop( iGarr[i] = iG; )
	}
	hsize_t dimsG[2] = { hsize_t(gInfo.nr), 3 };
	h5writeScalar(gidGspace, "ng", gInfo.nr);
	h5writeScalar(gidGspace, "ecutrho", std::max(e->cntrl.EcutRho, 4*e->cntrl.Ecut));
	h5writeVector(gidGspace, "FFTgrid", &gInfo.S[0], 3);
	h5writeVector(gidGspace, "components", &iGarr[0][0], dimsG, 2);
	H5Gclose(gidGspace);
	//----- symmetries related -----
	hid_t gidSymm = h5createGroup(gidHeader, "symmetry");
	const std::vector<SpaceGroupOp> ops = e->symm.getMatrices();
	std::vector<matrix3<int> > rots(ops.size());
	std::vector<vector3<> > trans(ops.size());
	for(size_t iSym=0; iSym<ops.size(); iSym++)
	{	rots[iSym] = ops[iSym].rot;
		trans[iSym] = ops[iSym].a;
	}
	hsize_t dimsRot[3] = { ops.size(), 3, 3 };
	hsize_t dimsTrans[2] = { ops.size(), 3 };
	h5writeScalar(gidSymm, "ntran", int(ops.size()));
	h5writeScalar(gidSymm, "cell_symmetry", 0);
	h5writeVector(gidSymm, "mtrx", &rots[0](0,0), dimsRot, 3);
	h5writeVector(gidSymm, "tnp", &trans[0][0], dimsTrans, 2);
	H5Gclose(gidSymm);
	//----- crystal related -----
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
	for(const auto& sp: e->iInfo.species)
	{	apos.insert(apos.end(), sp->atpos.begin(), sp->atpos.end());
		atyp.insert(atyp.end(), sp->atpos.size(), sp->atomicNumber);
	}
	hsize_t dimsApos[2] = { apos.size(), 3 };
	h5writeScalar(gidCrystal, "nat", int(apos.size()));
	h5writeVector(gidCrystal, "atyp", atyp);
	h5writeVector(gidCrystal, "apos", &apos[0][0], dimsApos, 2);
	H5Gclose(gidCrystal);
	H5Gclose(gidHeader);
	
	//Close file:
	H5Fclose(fid);
	logPrintf("done\n"); logFlush();
}

#endif