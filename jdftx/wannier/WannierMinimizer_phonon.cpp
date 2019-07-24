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

#include <wannier/WannierMinimizer.h>
#include <core/WignerSeitz.h>

//Output electron-phonon related quantities in Wannier basis:
void WannierMinimizer::saveMLWF_phonon(int iSpin)
{
	//Generate list of commensurate k-points in order present in the unit cell calculation
	prodPhononSup = wannier.phononSup[0] * wannier.phononSup[1] * wannier.phononSup[2];
	std::vector<int> ikArr; ikArr.reserve(prodPhononSup);
	for(unsigned ik=0; ik<kMesh.size(); ik++)
	{	vector3<> kSup = kMesh[ik].point.k * Diag(wannier.phononSup);
		double roundErr; round(kSup, &roundErr);
		if(roundErr < symmThreshold) //integral => commensurate with supercell
			ikArr.push_back(ik);
	}
	assert(int(ikArr.size()) == prodPhononSup);
	
	//Generate and write electron-phonon cell map:
	xAtoms.clear();
	for(const auto& sp: e.iInfo.species)
		xAtoms.insert(xAtoms.end(), sp->atpos.begin(), sp->atpos.end());
	assert(3*int(xAtoms.size()) == nPhononModes);
	ePhCellMap = getCellMap(
		e.gInfo.R, e.gInfo.R * Diag(wannier.phononSup),
		e.coulombParams.isTruncated(), xAtoms, xExpect, wannier.rSmooth,
		wannier.getFilename(Wannier::FilenameDump, "mlwfCellMapPh", &iSpin));
	//--- Corresponding cell map:
	if(mpiWorld->isHead())
	{	string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfCellWeightsPh", &iSpin);
		logPrintf("Dumping '%s'... ", fname.c_str()); logFlush();
		FILE* fp = fopen(fname.c_str(), "w");
		if(!fp) die_alone("could not open file for writing.\n");
		for(const auto& entry: ePhCellMap)
			entry.second.write_real(fp);
		fclose(fp);
		logPrintf("done.\n"); logFlush();
	}
	//--- Re-organize by unique cells in supercell:
	struct UniqueCell
	{	vector3<int> iR; //unique cell in supercell
		std::vector<std::pair<vector3<int>,matrix>> cells; //equivalent ones with weights
	};
	std::vector<UniqueCell> uniqueCells(prodPhononSup);
	vector3<int> stride(wannier.phononSup[1]*wannier.phononSup[2], wannier.phononSup[2], 1);
	for(const auto& entry: ePhCellMap)
	{	//Compute unique index:
		vector3<int> iR = entry.first;
		for(int iDir=0; iDir<3; iDir++)
			iR[iDir] = positiveRemainder(iR[iDir], wannier.phononSup[iDir]);
		int uniqIndex = dot(iR, stride);
		//Set/collect unique cell properties:
		uniqueCells[uniqIndex].iR = iR;
		uniqueCells[uniqIndex].cells.push_back(entry);
	}
	
	//Generate pairs of commensurate k-points along with pointer to Wannier rotation
	kpointPairs.clear();
	for(int ik1: ikArr)
	for(int ik2: ikArr)
	{	KpointPair kPair;
		kPair.k1 = kMesh[ik1].point.k;
		kPair.k2 = kMesh[ik2].point.k;
		kPair.ik1 = ik1;
		kPair.ik2 = ik2;
		kpointPairs.push_back(kPair);
	}
	TaskDivision(kpointPairs.size(), mpiWorld).myRange(iPairStart, iPairStop);
	
	//Read Bloch electron - nuclear displacement matrix elements
	string fnameIn = wannier.getFilename(Wannier::FilenameInit, "phononHsub", &iSpin);
	logPrintf("Reading '%s' ... ", fnameIn.c_str()); logFlush();
	MPIUtil::File fpIn;
	size_t matSizeIn = nBands*nBands * sizeof(complex);
	size_t modeStrideIn = kpointPairs.size() * matSizeIn; //input data size per phonon mode
	mpiWorld->fopenRead(fpIn, fnameIn.c_str(), nPhononModes * modeStrideIn);
	std::vector<std::vector<matrix>> phononHsub(nPhononModes, std::vector<matrix>(kpointPairs.size()));
	for(int iMode=0; iMode<nPhononModes; iMode++)
	{	//Read phononHsub and store with Wannier rotations:
		mpiWorld->fseek(fpIn, iMode*modeStrideIn + iPairStart*matSizeIn, SEEK_SET);
		for(int iPair=iPairStart; iPair<iPairStop; iPair++)
		{	matrix& phononHsubCur = phononHsub[iMode][iPair];
			phononHsubCur.init(nBands, nBands);
			mpiWorld->freadData(phononHsubCur, fpIn); //read from file
		}
	}
	
	//Apply translational invariance correction:
	std::vector<diagMatrix> Hsub_eigs = e.eVars.Hsub_eigs; //make available on all processes
	for(int q=0; q<e.eInfo.nStates; q++)
	{	Hsub_eigs[q].resize(nBands);
		mpiWorld->bcastData(Hsub_eigs[q], e.eInfo.whose(q));
	}
	for(int ik: ikArr)
	{	if(not isMine_q(ik, iSpin))
			DblochMesh[ik].init(nBands*nBands, 3);
		mpiWorld->bcastData(DblochMesh[ik], whose_q(ik, iSpin));
	}
	double nrmTot = 0., nrmCorr = 0.;
	for(int iPair=iPairStart; iPair<iPairStop; iPair++)
	{	const KpointPair& pair = kpointPairs[iPair];
		if(pair.ik1 == pair.ik2) //only Gamma-point phonons
		{	const diagMatrix& E = Hsub_eigs[kMesh[pair.ik1].point.iReduced + iSpin*qCount];
			const matrix& U = kMesh[pair.ik1].U;
			int nAtoms = nPhononModes/3;
			assert(nAtoms*3 == nPhononModes);
			for(int iDir=0; iDir<3; iDir++)
			{	//Calculate matrix element due to uniform translation of all atoms:
				matrix phononHsubMean;
				for(int iAtom=0; iAtom<nAtoms; iAtom++)
				{	int iMode = 3*iAtom + iDir;
					phononHsubMean += (1./(nAtoms*invsqrtM[iMode])) * phononHsub[iMode][iPair];
				}
				nrmTot += std::pow(nrm2(dagger(U) * phononHsubMean * U), 2);
				//Subtract the expected matrix elements from the sum rule connecting it to momentum matrix elements:
				complex* Hdata = phononHsubMean.data();
				const complex* Ddata = DblochMesh[pair.ik1].data() + DblochMesh[pair.ik1].index(0,iDir);
				double normFac = 1./(nAtoms*prodPhononSup);
				for(int b2=0; b2<nBands; b2++)
					for(int b1=0; b1<nBands; b1++)
						*(Hdata++) -= *(Ddata++) * (E[b1]-E[b2]) * normFac;
				nrmCorr += std::pow(nrm2(dagger(U) * phononHsubMean * U), 2);
				//Apply correction:
				for(int iAtom=0; iAtom<nAtoms; iAtom++)
				{	int iMode = 3*iAtom + iDir;
					phononHsub[iMode][iPair] -= invsqrtM[iMode] * phononHsubMean;
				}
			}
		}
	}
	DblochMesh.clear();
	mpiWorld->allReduce(nrmCorr, MPIUtil::ReduceSum);
	mpiWorld->allReduce(nrmTot, MPIUtil::ReduceSum);
	logPrintf("done. Translation invariance correction: %le\n", sqrt(nrmCorr/nrmTot)); logFlush();
	
	//Read quantities required for polar subtraction:
	std::vector<vector3<>> Zeff; 
	std::shared_ptr<LongRangeSum> lrs;
	if (wannier.polar)
	{	string fnameZeff = wannier.getFilename(Wannier::FilenameInit, "Zeff");
		logPrintf("\n"); logFlush();
		Zeff = readArrayVec3(fnameZeff);
		string fnameEps = wannier.getFilename(Wannier::FilenameInit, "epsInf"); 
		std::vector<vector3<>> eps = readArrayVec3(fnameEps);
		matrix3<> epsInf; epsInf.set_rows(eps[0], eps[1], eps[2]);
		lrs = std::make_shared<LongRangeSum>(e.gInfo.R, epsInf);
	}

	//Apply Wannier rotations
	logPrintf("Applying Wannier rotations ... "); logFlush();
	int nPairsMine = std::max(1, iPairStop-iPairStart); //avoid zero size matrices below
	matrix HePhTilde = zeroes(nCenters*nCenters*nPhononModes, nPairsMine);
	for(int iMode=0; iMode<nPhononModes; iMode++) //in Cartesian atom displacement basis
		for(int iPair=iPairStart; iPair<iPairStop; iPair++)
		{	const KpointPair& pair = kpointPairs[iPair];
			matrix& phononHsubCur = phononHsub[iMode][iPair]; //in Bloch basis
			phononHsubCur = (dagger(kMesh[pair.ik1].U) * phononHsubCur * kMesh[pair.ik2].U); //apply Wannier rotations
			//Subtract polar part in wannier-rotated version
			if(wannier.polar)
			{	vector3<> q = kMesh[pair.ik1].point.k - kMesh[pair.ik2].point.k;
				complex gLij =  complex(0,1)
					* ((4*M_PI) * invsqrtM[iMode] / (e.gInfo.detR * prodPhononSup))
					*  (*lrs)(q, Zeff[iMode], xAtoms[iMode/3]);
				for(int b=0; b<nCenters; b++)   
					phononHsubCur.data()[phononHsubCur.index(b,b)] -= gLij; //diagonal correction
			}
			callPref(eblas_copy)(HePhTilde.dataPref() + HePhTilde.index(0,iPair-iPairStart) + nCenters*nCenters*iMode,
				phononHsubCur.dataPref(), phononHsubCur.nData());
		}
	logPrintf("done.\n"); logFlush();
	
	//Wannierize and output one cell fixed at a time to minimize memory usage:
	string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfHePh", &iSpin);
	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
	FILE* fp = 0;
	if(mpiWorld->isHead())
	{	fp = fopen(fname.c_str(), "wb");
		if(!fp) die_alone("Error opening %s for writing.\n", fname.c_str());
	}
	matrix phase = zeroes(HePhTilde.nCols(), prodPhononSup);
	double kPairWeight = 1./prodPhononSup;
	double nrm2totSq = 0., nrm2imSq = 0.;
	std::map<vector3<int>, matrix> Hsum;
	for(const UniqueCell& cell1: uniqueCells)
	{	//calculate Fourier transform phase (with integration weights):
		for(int iPair=iPairStart; iPair<iPairStop; iPair++)
		{	const KpointPair& pair = kpointPairs[iPair];
			int iCell2 = 0;
			for(const UniqueCell& cell2: uniqueCells)
				phase.set(iPair-iPairStart, iCell2++, 
					kPairWeight * cis(2*M_PI*(dot(pair.k1,cell1.iR) - dot(pair.k2,cell2.iR))) );
		}
		//convert phononHsub from Bloch to wannier for each nuclear displacement mode:
		matrix H = HePhTilde * phase;
		mpiWorld->allReduceData(H, MPIUtil::ReduceSum);
		//write results for unique cells to file:
		if(realPartOnly)
		{	if(mpiWorld->isHead()) H.write_real(fp);
			nrm2totSq += std::pow(nrm2(H), 2); 
			nrm2imSq += std::pow(callPref(eblas_dnrm2)(H.nData(), ((double*)H.dataPref())+1, 2), 2); //look only at imaginary parts with a stride of 2
		}
		else { if(mpiWorld->isHead()) H.write(fp); }
		//collect sum rule, accounting for cell map and weights:
		const complex* Hdata = H.dataPref();
		int nDataPerCell = nPhononModes*nCenters*nCenters;
		for(const UniqueCell& cell2: uniqueCells)
		{	for(const auto& entry1: cell1.cells)
			{	for(const auto& entry2: cell2.cells)
				{	//copy H for current cells:
					std::vector<complex> HcopyVec(Hdata, Hdata+nDataPerCell);
					complex* Hcopy = HcopyVec.data();
					//prepare for sum rule collection:
					vector3<int> iRdiff = entry2.first - entry1.first;
					matrix& HsumCur = Hsum[iRdiff];
					if(not HsumCur) HsumCur = zeroes(nCenters*nCenters, 3);
					//Loop over atoms and directions:
					int iMode = 0;
					for(size_t iAtom=0; iAtom<xAtoms.size(); iAtom++)
					{	matrix w1i = entry1.second(iAtom,iAtom+1, 0,xExpect.size());
						matrix w2i = entry2.second(iAtom,iAtom+1, 0,xExpect.size());
						matrix w = transpose(w1i) * w2i;
						complex* HsumData = HsumCur.dataPref(); //summed over atoms, but not over directions
						for(int iVector=0; iVector<3; iVector++)
						{	callPref(eblas_zmul)(w.nData(), w.dataPref(), 1, Hcopy, 1); //apply weights
							callPref(eblas_zaxpy)(w.nData(), 1./invsqrtM[iMode], Hcopy, 1, HsumData, 1); //collect sum rule
							Hcopy += w.nData();
							HsumData += w.nData();
							iMode++;
						}
					}
				}
			}
			Hdata += nDataPerCell;
		}
	}
	if(mpiWorld->isHead()) fclose(fp);
	if(realPartOnly)
		logPrintf("done. Relative discarded imaginary part: %le\n", sqrt(nrm2imSq / nrm2totSq));
	else
		logPrintf("done.\n");
	
	//Write sum rule matrices and cell map:
	if(mpiWorld->isHead())
	{
		//Matrices:
		string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfHePhSum", &iSpin);
		logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
		FILE* fp = fopen(fname.c_str(), "wb");
		if(!fp) die_alone("Error opening %s for writing.\n", fname.c_str());
		double nrm2totSq = 0., nrm2imSq = 0.;
		for(const auto& entry: Hsum)
		{	const matrix& M = entry.second;
			if(realPartOnly)
			{	M.write_real(fp);
				nrm2totSq += std::pow(nrm2(M), 2); 
				nrm2imSq += std::pow(callPref(eblas_dnrm2)(M.nData(), ((double*)M.dataPref())+1, 2), 2); //look only at imaginary parts with a stride of 2
			}
			else M.write(fp);
		}
		if(realPartOnly)
			logPrintf("done. Relative discarded imaginary part: %le\n", sqrt(nrm2imSq / nrm2totSq));
		else
			logPrintf("done.\n");
		
		//Cell map:
		fname = wannier.getFilename(Wannier::FilenameDump, "mlwfCellMapPhSum", &iSpin);
		writeCellMap(Hsum, e.gInfo.R, fname);
	}
}
