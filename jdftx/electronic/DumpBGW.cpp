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
#include <electronic/DumpBGW_internal.h>
#include <electronic/ColumnBundle.h>
#include <electronic/ExactExchange.h>
#include <core/CoulombIsolated.h>

#ifndef HDF5_ENABLED
void Dump::dumpBGW()
{	die("BerkeleyGW output requires HDF5 support.\n");
}
#else

void Dump::dumpBGW()
{	if(!bgwParams) bgwParams = std::make_shared<BGWparams>(); //default parameters
	BGW bgw(*e, *bgwParams);
	bgw.writeWfn();
	if(bgwParams->saveVxc) bgw.writeVxc();
	if(bgwParams->saveVxx) bgw.writeVxx();
	if(bgwParams->rpaExx)
	{	double EXX_RPA = (*e->exx)(1., 0., e->eVars.F, e->eVars.C, NULL, NULL, true, &e->eVars.Hsub_eigs);
		logPrintf("\n EXX(RPA) = %25.16lf\n\n", EXX_RPA);
		logFlush();
	}
	if(e->eVars.fluidSolver and bgwParams->EcutChiFluid)
	{	bgw.writeChiFluid(true); //write q0
		bgw.writeChiFluid(false); //write q != q0
	}
	if(bgwParams->Ecut_rALDA)
	{	bgw.write_rALDA(true); //write q0
		bgw.write_rALDA(false); //write q != q0
	}
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
	
	//G-vectors for charge-density mesh:
	iGrho.resize(gInfo.nr);
	{	const vector3<int>& S = gInfo.S;
		size_t iStart = 0, iStop = gInfo.nr;
		THREAD_fullGspaceLoop(
			iGrho[i] = iG;
			for(int iDir=0; iDir<3; iDir++)
				if(2*iGrho[i][iDir]==gInfo.S[iDir])
					iGrho[i][iDir]-=gInfo.S[iDir]; //[-S/2,S/2) in BGW (rather than (-S/2,S/2] in JDFTx)
		)
	}
	
	
	//DFT nBands; corresponding E and F:
	nBands = eInfo.nBands;
	E = eVars.Hsub_eigs;
	F = eVars.F;
	VxcSub.resize(eInfo.nStates);
	VxxSub.resize(eInfo.nStates);
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
{	logPrintf("Computing Vxc ... "); logFlush();
	std::vector<matrix>& VxcSub = ((BGW*)this)->VxcSub;
	if(e.exCorr.exxFactor()) e.exx->prepareHamiltonian(e.exCorr.exxRange(), e.eVars.F, e.eVars.C);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
	{	if(VxcSub[q]) continue; //already calculated (dense version)
		//Calculate X-C matrix elements:
		ColumnBundle HCq = gInfo.dV * Idag_DiagV_I(eVars.C[q], eVars.Vxc);
		if(e.exCorr.needsKEdensity() && eVars.Vtau[eInfo.qnums[q].index()]) //metaGGA KE potential
		{	for(int iDir=0; iDir<3; iDir++)
				HCq -= (0.5*gInfo.dV) * D(Idag_DiagV_I(D(eVars.C[q],iDir), eVars.Vtau), iDir);
		}
		if(e.eInfo.hasU) //Contribution via atomic density matrix projections (DFT+U)
			e.iInfo.rhoAtom_grad(eVars.C[q], eVars.U_rhoAtom, HCq);
		if(e.exCorr.exxFactor()) //Exact-exchange contributions:
			e.exx->applyHamiltonian(e.exCorr.exxFactor(), e.exCorr.exxRange(), q, eVars.F[q], eVars.C[q], HCq);
		VxcSub[q] = eVars.C[q] ^ HCq;
	}
	logPrintf("Done.\n");
	writeV(VxcSub, e.dump.getFilename("bgw.vxc.dat"));
}


//Write exact exchange matrix elements for BGW
void BGW::writeVxx() const
{	if(((not bgwp.saveVxc) or e.exCorr.exxRange() or (not e.exCorr.exxFactor())) //writeVxc() did not prepare bare exchange
			and (not bgwp.nBandsDense)) //and neither did the dense diagonalize
		e.exx->prepareHamiltonian(0., e.eVars.F, e.eVars.C);
	logPrintf("Computing Vxx ... "); logFlush();
	std::vector<matrix>& VxxSub = ((BGW*)this)->VxxSub;
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
	{	if(VxxSub[q]) continue; //already calculated (dense version)
		//Calculate Vxx matrix elements:
		ColumnBundle HCq = eVars.C[q].similar(); HCq.zero();
		e.exx->applyHamiltonian(1., 0., q, eVars.F[q], eVars.C[q], HCq);
		VxxSub[q] = eVars.C[q] ^ HCq;
	}
	logPrintf("Done.\n");
	writeV(VxxSub, e.dump.getFilename("bgw.vxx.dat"));
	
	//Report EXX energy computed from matrix elements:
	double EXX = 0;
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		EXX += 0.5 * eInfo.qnums[q].weight * trace(diag(VxxSub[q]) * eVars.F[q]);
	mpiWorld->allReduce(EXX, MPIUtil::ReduceSum);
	logPrintf("\n      EXX = %25.16lf\n\n", EXX);
	logFlush();
}


void BGW::writeV(std::vector<matrix>& Vsub, string fname) const
{	//Output from head
	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
	if(mpiWorld->isHead())
	{	FILE* fp = fopen(fname.c_str(), "w");
		if(!fp) die_alone("failed to open for writing.\n");
		for(int ik=0; ik<nReducedKpts; ik++)
		{	fprintf(fp, "%.9f %.9f %.9f %4d %4d\n", k[ik][0], k[ik][1], k[ik][2],
				nBands*nSpins, bgwp.offDiagV ? nBands*nBands*nSpins : 0);
			for(int iSpin=0; iSpin<nSpins; iSpin++)
			{	int q=iSpin*nReducedKpts+ik;
				if(!eInfo.isMine(q))
				{	Vsub[q] = zeroes(nBands, nBands);
					mpiWorld->recvData(Vsub[q], eInfo.whose(q), q);
				}
				//Diagonal elements:
				for(int b=0; b<nBands; b++)
				{	complex V_eV = Vsub[q](b,b) / eV; //convert to eV
					fprintf(fp, "%4d %4d %+14.9f %+14.9f\n", iSpin+1, b+1, V_eV.real(), V_eV.imag());
				}
				//Off-diagonal elements:
				if(bgwp.offDiagV)
					for(int b2=0; b2<nBands; b2++)
					{	for(int b1=0; b1<nBands; b1++)
						{	complex V_eV = Vsub[q](b1,b2) / eV; //convert to eV
							fprintf(fp, "%4d %4d %4d %+14.9f %+14.9f\n", iSpin+1, b1+1, b2+1, V_eV.real(), V_eV.imag());
						}
					}
				if(!eInfo.isMine(q)) Vsub[q] = 0; //cleanup temporary memory
			}
		}
		fclose(fp);
	}
	else //send ones stored on other processes to head
	{	for(int q=0; q<eInfo.nStates; q++)
			if(eInfo.isMine(q))
				mpiWorld->sendData(Vsub[q], 0, q);
	}
	logPrintf("Done.\n");
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
	mpiWorld->allReduceData(ifmax, MPIUtil::ReduceMax);
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
		mpiWorld->bcastData(Ecur, eInfo.whose(q)); Eall.insert(Eall.end(), Ecur.begin(), Ecur.end());
		mpiWorld->bcastData(Fcur, eInfo.whose(q)); Fall.insert(Fall.end(), Fcur.begin(), Fcur.end());
	}
	hsize_t dimsKspinBands[3] = { hsize_t(nSpins), hsize_t(nReducedKpts), hsize_t(nBands) };
	h5writeVector(gidKpts, "el", Eall.data(), dimsKspinBands, 3);
	h5writeVector(gidKpts, "occ", Fall.data(), dimsKspinBands, 3);
	Eall.clear(); Fall.clear();
	H5Gclose(gidKpts);
	
	//---------- G-space group ----------
	hid_t gidGspace = h5createGroup(gidHeader, "gspace");
	hsize_t dimsG[2] = { hsize_t(gInfo.nr), 3 };
	h5writeScalar(gidGspace, "ng", gInfo.nr);
	h5writeScalar(gidGspace, "ecutrho", std::max(e.cntrl.EcutRho, 4*e.cntrl.Ecut)/Ryd);
	h5writeVector(gidGspace, "FFTgrid", &gInfo.S[0], 3);
	h5writeVector(gidGspace, "components", &iGrho[0][0], dimsG, 2);
	H5Gclose(gidGspace);
	
	//---------- Symmetries group ----------
	hid_t gidSymm = h5createGroup(gidHeader, "symmetry");
	size_t nSymmMax = 48;
	const std::vector<SpaceGroupOp> ops = e.symm.getMatrices();
	std::vector<matrix3<int> > rots(nSymmMax);
	std::vector<vector3<> > trans(nSymmMax);
	for(size_t iSym=0; iSym<std::min(ops.size(),nSymmMax); iSym++)
	{	matrix3<int> rotInv = det(ops[iSym].rot) * adjugate(ops[iSym].rot); //since |det(rot)| = 1
		rots[iSym] = rotInv; //BGW uses inverse convention
		trans[iSym] = (2*M_PI)*ops[iSym].a; //BGW used 2*pi times translation
	}
	hsize_t dimsRot[3] = { nSymmMax, 3, 3 };
	hsize_t dimsTrans[2] = { nSymmMax, 3 };
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
	h5writeScalar(gidCrystal, "alat", 1.);
	h5writeScalar(gidCrystal, "blat", 1.);
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


void BGW::write_rALDA(bool write_q0) const
{
	if((not write_q0) and k.size()==1) return; //only q0 present, so nothing to do
	if(nSpinor != 1) die("rALDA not supported for spinorial calculations");

	//Open file and write common header:
	hid_t fid = openHDF5(e.dump.getFilename(write_q0 ? "bgw.fxc0_rALDA.h5" : "bgw.fxc_rALDA.h5"));
	std::vector<vector3<>> q;
	std::vector<complex> freq;
	std::vector<std::vector<vector3<int>>> iGarr;
	std::vector<int> nBasis;
	int nBasisMax = 0;
	writeHeaderEps(fid, write_q0, "rALDA", q, freq, iGarr, nBasis, nBasisMax);
	
	//Prepare q-cutoff and XC contribution at each point on the grid:
	const double nCut = 1e-12; //regularize n below this value
	ScalarField nTot = eVars.get_nTot(); //spin-summed
	ScalarFieldArray xcData(2);
	nullToZero(xcData, gInfo);
	const double* n = nTot->data();
	double* qSqCut = xcData[0]->data();
	double* fxLDA = xcData[1]->data();
	for(int i=0; i<gInfo.nr; i++)
	{	double n_i = std::max(n[i], nCut); //regularized density
		double kFsq = std::pow((3*M_PI*M_PI)*n_i, 2./3); //Fermi wave-vector squared
		qSqCut[i] = 4*kFsq;
		fxLDA[i] = -M_PI/kFsq;
	}
	
	//Prepare truncation:
	//NOTE: only non-embedded spherical truncation supported as of now
	double Rc = 0.0;
	if(e.coulombParams.geometry != CoulombParams::Periodic)
	{	if((e.coulombParams.geometry == CoulombParams::Spherical) and (not e.coulombParams.embed))
		{	Rc = std::static_pointer_cast<CoulombSpherical>(e.coulomb)->Rc;
			logPrintf("rALDA Coulomb kernel will be spherical-truncated with Rc = %lg\n", Rc);
		}
		else logPrintf(
			"\nWARNING: only non-embedded spherical kernel supported for rALDA.\n"
			"Falling back to periodic coulomb kernel for rALDA output.\n\n");
	}
	else logPrintf("rALDA Coulomb kernel will be periodic.\n");
	
	//------------ Calculate and write matrices --------------
	//--- Process grid:
	int nProcesses = mpiWorld->nProcesses();
	int nProcsRow = int(round(sqrt(nProcesses)));
	while(nProcesses % nProcsRow) nProcsRow--;
	int nProcsCol = nProcesses / nProcsRow;
	int iProcRow = mpiWorld->iProcess() / nProcsCol;
	int iProcCol = mpiWorld->iProcess() % nProcsCol;
	logPrintf("\tInitialized %d x %d process grid.\n", nProcsRow, nProcsCol);
	//--- Determine matrix division:
	int rowStart = (nBasisMax * iProcRow) / nProcsRow;
	int rowStop = (nBasisMax * (iProcRow+1)) / nProcsRow;
	int nRowsMine = rowStop - rowStart;
	int colStart = (nBasisMax * iProcCol) / nProcsCol;
	int colStop = (nBasisMax * (iProcCol+1)) / nProcsCol;
	int nColsMine = colStop - colStart;
	if(std::max(nProcsRow, nProcsCol) > nBasisMax)
		die("\tNo data on some processes: reduce # processes.\n");
	//--- create dataset:
	hid_t gidMats = h5createGroup(fid, "mats");
	hsize_t dims[6] = { hsize_t(q.size()), hsize_t(nSpins), hsize_t(freq.size()),
		hsize_t(nBasisMax), hsize_t(nBasisMax), 2 };
	hid_t sid = H5Screate_simple(6, dims, NULL);
	hid_t did = H5Dcreate(gidMats, "matrix", H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t plid = H5Pcreate(H5P_DATASET_XFER);
	H5Sclose(sid);
	//--- loop over q:
	logPrintf("\tState:"); logFlush();
	for(int iq=0; iq<int(q.size()); iq++)
	{	std::vector<matrix> buf(nSpins, zeroes(nRowsMine, nColsMine));
		for(int iColMine=0; iColMine<nColsMine; iColMine++)
		{	int iCol = colStart + iColMine;
			if(iCol >= nBasis[iq]) continue;
			for(int iRowMine=0; iRowMine<nRowsMine; iRowMine++)
			{	int iRow = rowStart + iRowMine;
				if(iRow >= nBasis[iq]) continue;
				//Current wave vectors:
				vector3<int> iGdiff = iGarr[iq][iCol] - iGarr[iq][iRow]; //(G'-G) in recip coords for G in row, G' in col
				double qSqSym = sqrt(
					gInfo.GGT.metric_length_squared(q[iq] + iGarr[iq][iCol]) *
					gInfo.GGT.metric_length_squared(q[iq] + iGarr[iq][iRow]));
				//Collect contributions to XC and short-ranged coulomb-screening parts separately:
				complex fxc, coulombS;
				//--- loop over grid points
				const vector3<int>& S = gInfo.S;
				vector3<> iGdiffByS = inv(Diag(vector3<>(S))) * iGdiff; //elementwise iGdiff / S
				size_t iStart=0, iStop=gInfo.nr;
				THREAD_rLoop(
					complex phase = cis(2*M_PI*dot(iGdiffByS, iv)); //exp(-i(G-G').r) for G in row, G' in col
					if(qSqSym < qSqCut[i])
						fxc += phase * fxLDA[i];
					else
					{	double truncation = (Rc ? (1.0 - cos(Rc * sqrt(qSqSym))) : 1.0);
						coulombS -= phase * (truncation * 4*M_PI/qSqSym);  //difference between screened and bare Coulomb parts
					}
				)
				//Set spin diagonal and off-diagonal components:
				double prefac = 1./gInfo.nr; //scale factor to convert above sum to unit cell average
				for(int iSpin=0; iSpin<nSpins; iSpin++)
					buf[iSpin].set(iRowMine, iColMine, prefac*(
						coulombS //coulomb in both diagonal and off-diagonal
						+ (iSpin==0 ? nSpins*fxc : 0))); //fxc only in spin-diagonal part
			}
		}
		for(int iSpin=0; iSpin<nSpins; iSpin++)
		{	//Write results:
			hsize_t offset[6] = { hsize_t(iq), hsize_t(iSpin), hsize_t(0), hsize_t(colStart), hsize_t(rowStart), 0 };
			hsize_t count[6] = { 1, 1, 1, hsize_t(nColsMine), hsize_t(nRowsMine), 2 };
			sid = H5Dget_space(did);
			H5Sselect_hyperslab(sid, H5S_SELECT_SET, offset, NULL, count, NULL);
			hid_t sidMem = H5Screate_simple(6, count, NULL);
			H5Dwrite(did, H5T_NATIVE_DOUBLE, sidMem, sid, plid, buf[iSpin].data());
			H5Sclose(sidMem);
		}
		logPrintf(" %d", iq+1); logFlush();
	}
	logPrintf("\n"); logFlush();
	//--- close dataset:
	H5Pclose(plid);
	H5Dclose(did);
	H5Gclose(gidMats);
	
	//Close file:
	H5Fclose(fid);
	logPrintf("\tDone.\n"); logFlush();
}

#endif //HDF5_ENABLED
