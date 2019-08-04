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
#include <fluid/FluidSolver.h>
#include <core/SphericalHarmonics.h>
#include <core/LatticeUtils.h>


//Sort comparator for iG sorting
struct IndexedComparator
{	double* keyA; int* keyB; //sort first by keyA (KE), and then by keyB (grid index) to break ties
	IndexedComparator(double* keyA, int* keyB) : keyA(keyA), keyB(keyB) {}
	bool operator()(const int& i1, const int& i2) const
	{	double A1=keyA[i1], A2=keyA[i2];
		if(fabs(A1-A2) > symmThresholdSq) //tolerance for comparing double precision
			return A1 < A2; //unequal at specified tolerance
		else //check second key since first key is equal at specified tolerance
			return keyB[i1] < keyB[i2];
	}
};


template<typename T> void applyIndex(std::vector<T>& arr, const std::vector<int>& index)
{	assert(arr.size() == index.size());
	std::vector<T> arrCopy(arr);
	for(size_t i=0; i<index.size(); i++)
		arr[i] = arrCopy[index[i]];
}


inline void initYlmAndWeights(const  vector3<>& q, const std::vector<vector3<int>> iGarr, int start, int stop,
	const matrix3<>& GT, const std::vector<FluidSolver::SusceptibilityTerm>& susceptibility,
	std::vector<diagMatrix>& w, std::vector<std::vector<diagMatrix>>& YlmArr)
{
	int nTerms = susceptibility.size();
	int nMine = stop - start;
	w.assign(nTerms, diagMatrix(nMine, 1.)); //initialize to local
	YlmArr.resize(nTerms);
	//Calculate G-vector magnitudes and directions:
	diagMatrix Gmag; Gmag.reserve(nMine);
	std::vector<vector3<>> Gvec; Gvec.reserve(nMine);
	for(int iG=start; iG<stop; iG++)
	{	Gvec.push_back(GT * (iGarr[iG] + q));
		Gmag.push_back(Gvec.back().length());
	}
	for(int iTerm=0; iTerm<nTerms; iTerm++)
	{	const FluidSolver::SusceptibilityTerm& term = susceptibility[iTerm];
		//Override weight function if needed:
		if(term.w)
		{	for(int iMine=0; iMine<nMine; iMine++)
				w[iTerm][iMine] = (*term.w)(Gmag[iMine]);
		}
		//Calculate spherical harmonics
		YlmArr[iTerm].assign(2*term.l+1, diagMatrix(nMine));
		for(int m=-term.l; m<=term.l; m++)
			for(int iMine=0; iMine<nMine; iMine++)
				YlmArr[iTerm][m+term.l][iMine] = Ylm(term.l, m, Gvec[iMine]);
	}
}


//Translate CoulombParams::Geometry to BerkeleyGW flag:
inline int get_icutv(CoulombParams::Geometry geom)
{	switch(geom)
	{	case CoulombParams::Periodic: return 0; //no truncation
		case CoulombParams::Spherical: return 1; //spherical truncation - 0D systems
		case CoulombParams::Isolated: return 5; //cell box truncation - 0D systems
		case CoulombParams::Cylindrical: return 4; //cell wire truncation - 1D systems
		case CoulombParams::Wire: return 4; //cell  wire truncation - 1D systems
		case CoulombParams::Slab: return 6; //cell slab truncation - 2D systems
		default: return 0; //never encountered (only to prevent warning)
	}
}


//Write fluid polarizability
void BGW::writeChiFluid(bool write_q0) const
{	assert(eVars.fluidSolver); //should only be called in solvated cases
	if((not write_q0) and k.size()==1) return; //only q0 present, so nothing to do
	
	//Open file and write common header:
	hid_t fid = openHDF5(e.dump.getFilename(write_q0 ? "bgw.chi0Fluid.h5" : "bgw.chiFluid.h5"));
	writeHeaderMF(fid);
	logPrintf("\n");
	
	//---------- Epsilon header ----------
	hid_t gidHeader = h5createGroup(fid, "eps_header");
	h5writeScalar(gidHeader, "versionnumber", 1);
	h5writeScalar(gidHeader, "flavor", 2);  //=> complex wavefunctions
	
	//--- Params group:
	hid_t gidParams = h5createGroup(gidHeader, "params");
	h5writeScalar(gidParams, "matrix_type", 2); //=> chi
	h5writeScalar(gidParams, "has_advanced", 0); //=> retarded response only
	h5writeScalar(gidParams, "nmatrix", nSpins); //=> polarizability per spin
	h5writeScalar(gidParams, "matrix_flavor", 2); //=> complex matrices
	h5writeScalar(gidParams, "icutv", get_icutv(e.coulombParams.geometry)); //truncation geometry
	h5writeScalar(gidParams, "ecuts", bgwp.EcutChiFluid/Ryd); //output cutoff for polarizability
	h5writeScalar(gidParams, "nband", 0); //No bands used for fluid polarizability
	h5writeScalar(gidParams, "efermi", 0.); //Not meaningful for fluid polarizability
	H5Gclose(gidParams);
	
	//--- q-points group:
	hid_t gidQpoints = h5createGroup(gidHeader, "qpoints");
	std::vector<vector3<>> q;
	if(write_q0)
		q.assign(1, bgwp.q0); //only q0
	else
	{	//same as k-mesh except Gamma-point removed (since separate q0):
		bool foundGamma = false;
		for(const vector3<>& qCur: k)
			if(qCur.length_squared() < symmThresholdSq)
				foundGamma = true;
			else
				q.push_back(qCur);
		if(not foundGamma)
			die("Fluid polarizability output for BGW only supported for Gamma-centered mesh.\n");
	}
	int nq = q.size();
	h5writeScalar(gidQpoints, "nq", nq);
	hsize_t dimsQ[2] = { hsize_t(nq), 3 };
	h5writeVector(gidQpoints, "qpts", &q[0][0], dimsQ, 2);
	h5writeVector(gidQpoints, "qgrid", &eInfo.kFoldingCount()[0], 3);
	h5writeVector(gidQpoints, "qpt_done", std::vector<int>(nq, 1));
	H5Gclose(gidQpoints);
	
	//--- freq group:
	//------ Create frequency grid in eV:
	std::vector<complex> freq;
	for(double freqRe=0.; freqRe<bgwp.freqReMax_eV+symmThreshold; freqRe+=bgwp.freqReStep_eV)
		freq.push_back(complex(freqRe, bgwp.freqBroaden_eV));
	if(bgwp.freqPlasma)
	{	//GW imaginary freq. grid
		double freqPlasma_eV = bgwp.freqPlasma/eV;
		for(int i=0; i<bgwp.freqNimag; i++)
			freq.push_back(complex(0., i*freqPlasma_eV/(bgwp.freqNimag-i)));
	}
	else
	{	//RPA imaginary freq. grid
		for(int i=0; i<bgwp.freqNimag; i++)
			freq.push_back(complex(0., (1./eV)*tan(0.5*M_PI*(1.-(i+1)*1./bgwp.freqNimag))));
	}
	//------ Output frequency grid in eV:
	hid_t gidFreq = h5createGroup(gidHeader, "freqs");
	h5writeScalar(gidFreq, "freq_dep", 2); //=> full-frequency
	h5writeScalar(gidFreq, "nfreq", int(freq.size())); //number of frequencies
	h5writeScalar(gidFreq, "nfreq_imag", bgwp.freqNimag); //number of imaginary frequencies
	hsize_t dimsFreq[2] = { hsize_t(freq.size()), 2 };
	h5writeVector(gidFreq, "freqs", &(freq[0].real()), dimsFreq, 2);
	H5Gclose(gidFreq);
	//------ Convert to Hartrees for internal storage:
	for(complex& f: freq)
		f *= eV;
	logPrintf("\tInitialized frequency grid of size %d.\n", int(freq.size()));
	
	//--- gspace group:
	//------ initialize list of G vectors for each q:
	int nG = iGrho.size();
	std::vector<double> KE(nq * nG); //in Ryd, in rho order
	std::vector<int> indEpsToRho(nq * nG);
	std::vector<int> indRhoToEps(nq * nG);
	std::vector<std::vector<vector3<int>>> iGarr(nq);
	std::vector<int> nBasis(nq);
	int nBasisMax = 0;
	{	//Grid index used for tie-breaking
		vector3<int> iGstride(e.gInfo.S[1]*e.gInfo.S[2], e.gInfo.S[2], 1);
		std::vector<int> gridInd(nG); 
		for(int g=0; g<nG; g++)
			gridInd[g] = dot(iGrho[g], iGstride);
		for(int iq=0; iq<nq; iq++)
		{	double* KEcur = &KE[nG*iq]; //KE in Ryd with q
			int* sortIndex = &indEpsToRho[nG*iq];
			int* sortIndexInv = &indRhoToEps[nG*iq];
			//Calculate KE with q and sort by it (eps sorting):
			vector3<> qCur = write_q0 ? vector3<>() : q[iq]; //don't account for small q0 offset in sorting
			for(int g=0; g<nG; g++)
			{	KEcur[g] = e.gInfo.GGT.metric_length_squared(iGrho[g] + qCur); //no 0.5 since in Ryd
				sortIndex[g] = g;
			}
			std::sort(sortIndex, sortIndex+nG, IndexedComparator(KEcur, gridInd.data()));
			//Collect basis functions within G-sphere in KE order:
			double EcutRy = 2.*bgwp.EcutChiFluid;
			for(int g=0; (g<nG && KEcur[sortIndex[g]]<=EcutRy); g++)
				iGarr[iq].push_back(iGrho[sortIndex[g]]);
			nBasis[iq] = iGarr[iq].size();
			nBasisMax = std::max(nBasisMax, nBasis[iq]);
			//Invert sort index and convert to 1-based:
			for(int g=0; g<nG; g++)
				sortIndexInv[sortIndex[g]++] = g+1; //switches both from 0 to 1-based indexing
		}
	}
	//------ output g-space header to file:
	hid_t gidGspace = h5createGroup(gidHeader, "gspace");
	h5writeVector(gidGspace, "nmtx", nBasis); //matrix dimensions for each q
	h5writeScalar(gidGspace, "nmtx_max", nBasisMax); //maximum value of above over all q
	hsize_t dimsInd[2] = { hsize_t(nq), hsize_t(nG) };
	h5writeVector(gidGspace, "ekin", KE.data(), dimsInd, 2);
	h5writeVector(gidGspace, "gind_eps2rho", indEpsToRho.data(), dimsInd, 2);
	h5writeVector(gidGspace, "gind_rho2eps", indRhoToEps.data(), dimsInd, 2);
	H5Gclose(gidGspace);
	KE.clear();
	indEpsToRho.clear();
	indRhoToEps.clear();
	H5Gclose(gidHeader); //end eps header
	
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
	//--- Get susceptibility description from fluidSolver:
	std::vector<FluidSolver::SusceptibilityTerm> susceptibility;
	ScalarFieldTildeArray sTilde;
	e.eVars.fluidSolver->getSusceptibility(freq, susceptibility, sTilde, bgwp.elecOnly);
	std::vector<const complex*> sTildeData; for(ScalarFieldTilde& sT: sTilde) sTildeData.push_back(sT->data());
	//--- create dataset:
	hid_t gidMats = h5createGroup(fid, "mats");
	hsize_t dims[6] = { hsize_t(q.size()), hsize_t(nSpins*nSpinor),
		hsize_t(freq.size()), hsize_t(nBasisMax), hsize_t(nBasisMax), 2 };
	hid_t sid = H5Screate_simple(6, dims, NULL);
	hid_t did = H5Dcreate(gidMats, "matrix", H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t plid = H5Pcreate(H5P_DATASET_XFER);
	H5Sclose(sid);
	//--- loop over q and frequencies:
	matrix buf(nRowsMine, nColsMine);
	for(int iq=0; iq<int(q.size()); iq++)
	{	logPrintf("\tState %d, freq:", iq+1); logFlush();
		//Initialize frequency-independent terms by row and column:
		int nTerms = susceptibility.size();
		std::vector<diagMatrix> wRow, wCol; //weight functions w(G)
		std::vector<std::vector<diagMatrix>> YlmRow, YlmCol; //Ylm(Ghat) * G^l
		initYlmAndWeights(q[iq], iGarr[iq], rowStart, rowStop, gInfo.GT, susceptibility, wRow, YlmRow);
		initYlmAndWeights(q[iq], iGarr[iq], colStart, colStop, gInfo.GT, susceptibility, wCol, YlmCol);
		//Loop over frequencies:
		for(int iFreq=0; iFreq<int(freq.size()); iFreq++)
		{	buf.zero();
			for(int iColMine=0; iColMine<nColsMine; iColMine++)
			{	int iCol = colStart + iColMine;
				if(iCol >= nBasis[iq]) continue;
				for(int iRowMine=0; iRowMine<nRowsMine; iRowMine++)
				{	int iRow = rowStart + iRowMine;
					if(iRow >= nBasis[iq]) continue;
					//Determine index into shape arrays:
					vector3<int> iGdiff = iGarr[iq][iCol] - iGarr[iq][iRow];
					for(int dir=0; dir<3; dir++)
						iGdiff[dir] = positiveRemainder(iGdiff[dir], gInfo.S[dir]); //wrap positive
					bool conj = false;
					if(iGdiff[2] > gInfo.S[2]/2) //account for half G space
					{	conj = true; //pick up from -G with complex conjugate
						for(int dir=0; dir<3; dir++)
							iGdiff[dir] = iGdiff[dir] ? gInfo.S[dir]-iGdiff[dir] : 0; //G -> -G
					}
					int diffIndex = gInfo.halfGindex(iGdiff);
					//Collect contributions to chi:
					complex result;
					for(int iTerm=0; iTerm<nTerms; iTerm++)
					{	const FluidSolver::SusceptibilityTerm& term = susceptibility[iTerm];
						//Spherical harmonics (div.grad for l=1) factor:
						double Ydot = 0.;
						int mCount = 2*term.l+1;
						for(int lm=0; lm<mCount; lm++)
							Ydot += YlmRow[iTerm][lm][iRowMine] * YlmCol[iTerm][lm][iColMine];
						complex sTildeCur = sTildeData[term.iSite][diffIndex];
						if(conj) sTildeCur = sTildeCur.conj();
						result -= term.prefactor[iFreq] * sTildeCur
							* (wRow[iTerm][iRowMine] * wCol[iTerm][iColMine] * Ydot * (4*M_PI)/mCount);
					}
					buf.set(iRowMine, iColMine, result);
				}
			}
			//Write results:
			for(int iSpin=0; iSpin<nSpins*nSpinor; iSpin++)
			{	hsize_t offset[6] = { hsize_t(iq), hsize_t(iSpin), hsize_t(iFreq), hsize_t(colStart), hsize_t(rowStart), 0 };
				hsize_t count[6] = { 1, 1, 1, hsize_t(nColsMine), hsize_t(nRowsMine), 2 };
				sid = H5Dget_space(did);
				H5Sselect_hyperslab(sid, H5S_SELECT_SET, offset, NULL, count, NULL);
				hid_t sidMem = H5Screate_simple(6, count, NULL);
				H5Dwrite(did, H5T_NATIVE_DOUBLE, sidMem, sid, plid, buf.data());
				H5Sclose(sidMem);
			}
			logPrintf(" %d", iFreq+1); logFlush();
		}
		logPrintf("\n"); logFlush();
	}
	//--- close dataset:
	H5Pclose(plid);
	H5Dclose(did);
	H5Gclose(gidMats);
	
	//Close file:
	H5Fclose(fid);
	logPrintf("\tDone.\n"); logFlush();
}

#endif //HDF5_ENABLED