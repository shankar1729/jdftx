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
	//Generate phonon cell map (will match that used by phonon):
	xAtoms.clear();
	for(const auto& sp: e.iInfo.species)
		xAtoms.insert(xAtoms.end(), sp->atpos.begin(), sp->atpos.end());
	assert(3*int(xAtoms.size()) == nPhononModes);
	phononCellMap = getCellMap(
		e.gInfo.R, e.gInfo.R * Diag(wannier.phononSup),
		e.coulombParams.isTruncated(), xAtoms, xAtoms, wannier.rSmooth); //phonon force-matrix cell map
	
	//Read phonon force matrix:
	std::map<vector3<int>, matrix> phononOmegaSq;
	if(mpiWorld->isHead())
	{	string fname = wannier.getFilename(Wannier::FilenameInit, "phononOmegaSq");
		logPrintf("Reading '%s' ... ", fname.c_str()); logFlush();
		FILE* fp = fopen(fname.c_str(), "r");
		for(const auto iter: phononCellMap)
		{	matrix omegaSqCur(nPhononModes, nPhononModes);
			omegaSqCur.read_real(fp);
			phononOmegaSq[iter.first] = omegaSqCur;
		}
		fclose(fp);
		logPrintf("done.\n"); logFlush();
	}
	
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
	
	//Wrap atoms to WS cell if necessary:
	std::vector<vector3<int>> dxAtoms;
	if(wannier.wrapWS)
	{	logPrintf("\nWrapping phonon basis to Wigner-Seitz cell:\n");
		//Calculate offsets and new atom positions:
		std::vector<vector3<>> xAtomsNew;
		for(const vector3<>& xAt: xAtoms)
		{	vector3<> x = xAt - e.coulombParams.embedCenter; //lattice coordinates w.r.t Coulomb center
			vector3<> xWS = ws->restrict(x);
			vector3<int> dx = round(xWS-x);
			logPrintf("\t[ %3d %3d %3d ]\n", dx[0], dx[1], dx[2]);
			dxAtoms.push_back(dx);
			xAtomsNew.push_back(xAt + dx);
		}
		//Construct offset force matrix:
		std::map<vector3<int>, matrix> phononOmegaSqNew;
		for(const auto& iter: phononOmegaSq)
			for(size_t iAtom=0; iAtom<xAtoms.size(); iAtom++)
			for(size_t jAtom=0; jAtom<xAtoms.size(); jAtom++)
			{	vector3<int> iRnew = iter.first - (dxAtoms[jAtom] - dxAtoms[iAtom]); //iR + x_j - x_i stays invariant
				matrix& oSqCur = phononOmegaSqNew[iRnew];
				if(!oSqCur) oSqCur = zeroes(nPhononModes, nPhononModes);
				oSqCur.set(3*iAtom, 3*iAtom+3, 3*jAtom, 3*jAtom+3,
					iter.second(3*iAtom, 3*iAtom+3, 3*jAtom, 3*jAtom+3) );
			}
		//Construct new cell map:
		std::map<vector3<int>,matrix> phononCellMapNew = getCellMap(
			e.gInfo.R, e.gInfo.R * Diag(wannier.phononSup),
			e.coulombParams.isTruncated(), xAtomsNew, xAtomsNew, wannier.rSmooth);
		//Reduce force matrix to it:
		double nrmSqKept = 0., nrmSqDropped = 0.;
		for(auto iter=phononOmegaSqNew.begin(); iter!=phononOmegaSqNew.end();)
		{	double nrmSqCur = std::pow(nrm2(iter->second), 2);
			if(phononCellMapNew.find(iter->first) == phononCellMapNew.end())
			{	//Not in new cell map; remove
				nrmSqDropped += nrmSqCur;
				iter = phononOmegaSqNew.erase(iter);
			}
			else
			{	//In new cell map; keep
				nrmSqKept += nrmSqCur;
				iter++;
			}
		}
		logPrintf("\tRelative error in phonon force matrix remapping: %le\n\n",
			sqrt(nrmSqDropped / (nrmSqKept + nrmSqDropped)));
		//Replace originals:
		xAtoms = xAtomsNew;
		phononOmegaSq = phononOmegaSqNew;
		phononCellMap = phononCellMapNew;
	}
	
	//Generate electron-phonon cell map:
	ePhCellMap = getCellMap(
		e.gInfo.R, e.gInfo.R * Diag(wannier.phononSup),
		e.coulombParams.isTruncated(), xAtoms, xExpect, wannier.rSmooth); //e-ph elements cell map
	//--- Add ePhCellMap cells missing in phononCellMap:
	for(const auto iter: ePhCellMap)
		if(phononCellMap.find(iter.first) == phononCellMap.end())
		{	phononCellMap[iter.first] = zeroes(xAtoms.size(), xAtoms.size());
			phononOmegaSq[iter.first] = zeroes(nPhononModes, nPhononModes);
		}
	//--- Add phononCellMap cells missing in ePhCellMap:
	for(const auto iter: phononCellMap)
		if(ePhCellMap.find(iter.first) == ePhCellMap.end())
			ePhCellMap[iter.first] = zeroes(xAtoms.size(), xExpect.size());
	
	//Output force matrix on unified phonon cellMap:
	if(mpiWorld->isHead())
	{	string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfOmegaSqPh", &iSpin);
		logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
		FILE* fp = fopen(fname.c_str(), "w");
		for(const auto iter: phononOmegaSq)
			iter.second.write_real(fp);
		fclose(fp);
		logPrintf("done.\n"); logFlush();
	}
	
	//Output unified phonon cellMap:
	if(mpiWorld->isHead())
	{	string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfCellMapPh", &iSpin);
		logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
		FILE* fp = fopen(fname.c_str(), "w");
		fprintf(fp, "#i0 i1 i2  x y z  (integer lattice combinations, and cartesian offsets)\n");
		for(const auto& entry: phononCellMap)
		{	const vector3<int>& i = entry.first;
			vector3<> r = e.gInfo.R * i;
			fprintf(fp, "%+2d %+2d %+2d  %+11.6lf %+11.6lf %+11.6lf\n", i[0], i[1], i[2], r[0], r[1], r[2]);
		}
		fclose(fp);
		logPrintf("done.\n"); logFlush();
	}
	
	//Output phonon cellMapSq:
	if(mpiWorld->isHead())
	{	string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfCellMapSqPh", &iSpin);
		logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
		FILE* fp = fopen(fname.c_str(), "w");
		fprintf(fp, "#i0 i1 i2  i0' i1' i2'   (integer lattice combinations for pairs of sites)\n");
		for(const auto& entry1: phononCellMap)
		for(const auto& entry2: phononCellMap)
		{	const vector3<int>& iR1 = entry1.first;
			const vector3<int>& iR2 = entry2.first;
			fprintf(fp, "%+2d %+2d %+2d  %+2d %+2d %+2d\n", iR1[0], iR1[1], iR1[2], iR2[0], iR2[1], iR2[2]);
		}
		fclose(fp);
		logPrintf("done.\n"); logFlush();
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
			//Translate for phonon basis wrapping (if any):
			if(wannier.wrapWS)
				phononHsubCur *= cis(-2*M_PI*dot(dxAtoms[iMode/3], kpointPairs[iPair].k1 - kpointPairs[iPair].k2));
		}
	}
	
	//Apply translational invariance correction:
	std::vector<diagMatrix> Hsub_eigs = e.eVars.Hsub_eigs; //make available on all processes
	for(int q=0; q<e.eInfo.nStates; q++)
	{	Hsub_eigs[q].resize(nBands);
		mpiWorld->bcastData(Hsub_eigs[q], e.eInfo.whose(q));
	}
	double nrmTot = 0., nrmCorr = 0.;
	for(int iPair=iPairStart; iPair<iPairStop; iPair++)
	{	const KpointPair& pair = kpointPairs[iPair];
		if(pair.ik1 == pair.ik2) //only Gamma-point phonons
		{	const diagMatrix& E = Hsub_eigs[kMesh[pair.ik1].point.iReduced + iSpin*qCount];
			int nAtoms = nPhononModes/3;
			assert(nAtoms*3 == nPhononModes);
			for(int iDir=0; iDir<3; iDir++)
			{	//Calculate matrix element due to uniform translation of all atoms:
				matrix phononHsubMean;
				for(int iAtom=0; iAtom<nAtoms; iAtom++)
				{	int iMode = 3*iAtom + iDir;
					phononHsubMean += (1./(nAtoms*invsqrtM[iMode])) * phononHsub[iMode][iPair];
					nrmTot += std::pow(nrm2(phononHsub[iMode][iPair])/invsqrtM[iMode], 2);
				}
				//Restrict correction to degenerate subspaces:
				for(int b1=0; b1<nBands; b1++)
					for(int b2=0; b2<nBands; b2++)
						if(fabs(E[b1]-E[b2]) > 1e-4)
							phononHsubMean.set(b1,b2, 0.);
				nrmCorr += nAtoms*std::pow(nrm2(phononHsubMean), 2);
				//Apply correction:
				for(int iAtom=0; iAtom<nAtoms; iAtom++)
				{	int iMode = 3*iAtom + iDir;
					phononHsub[iMode][iPair] -= invsqrtM[iMode] * phononHsubMean;
				}
			}
		}
	}
	logPrintf("done. Translation invariance correction: %le\n", sqrt(nrmCorr/nrmTot)); logFlush();
	
	//Apply Wannier rotations
	logPrintf("Applying Wannier rotations ... "); logFlush();
	int nPairsMine = std::max(1, iPairStop-iPairStart); //avoid zero size matrices below
	matrix HePhTilde = zeroes(nCenters*nCenters*nPhononModes, nPairsMine);
	for(int iMode=0; iMode<nPhononModes; iMode++)
		for(int iPair=iPairStart; iPair<iPairStop; iPair++)
		{	const KpointPair& pair = kpointPairs[iPair];
			matrix& phononHsubCur = phononHsub[iMode][iPair];
			phononHsubCur = (dagger(kMesh[pair.ik1].U) * phononHsubCur * kMesh[pair.ik2].U); //apply Wannier rotations
			callPref(eblas_copy)(HePhTilde.dataPref() + HePhTilde.index(0,iPair-iPairStart) + nCenters*nCenters*iMode,
				phononHsubCur.dataPref(), phononHsubCur.nData());
		}
	logPrintf("done.\n"); logFlush();
	
	//Save Wanneirized:
	dumpWannierizedPh(HePhTilde, 1, "mlwfHePh", realPartOnly, iSpin);
}


void WannierMinimizer::dumpWannierizedPh(const matrix& Htilde, int nMatrices, string varName, bool realPartOnly, int iSpin) const
{
	//Wannierize and output one cell fixed at a time to minimize memory usage:
	string fname = wannier.getFilename(Wannier::FilenameDump, varName, &iSpin);
	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
	FILE* fp = 0;
	if(mpiWorld->isHead())
	{	fp = fopen(fname.c_str(), "wb");
		if(!fp) die_alone("Error opening %s for writing.\n", fname.c_str());
	}
	matrix phase = zeroes(Htilde.nCols(), phononCellMap.size());
	double kPairWeight = 1./prodPhononSup;
	double nrm2totSq = 0., nrm2imSq = 0.;
	for(const auto& entry1: ePhCellMap)
	{	//calculate Fourier transform phase (with integration weights):
		for(int iPair=iPairStart; iPair<iPairStop; iPair++)
		{	const KpointPair& pair = kpointPairs[iPair];
			int iCell2 = 0;
			for(const auto& entry2: ePhCellMap)
			{	const vector3<int>& iR1 = entry1.first;
				const vector3<int>& iR2 = entry2.first;
				phase.set(iPair-iPairStart, iCell2++, kPairWeight * cis(2*M_PI*(dot(pair.k1,iR1) - dot(pair.k2,iR2))) );
			}
		}
		//convert phononHsub from Bloch to wannier for each nuclear displacement mode:
		matrix HePh = Htilde * phase;
		mpiWorld->allReduceData(HePh, MPIUtil::ReduceSum);
		//apply cell weights:
		complex* HePhData = HePh.dataPref();
		for(const auto& entry2: ePhCellMap)
			for(size_t iAtom=0; iAtom<xAtoms.size(); iAtom++)
			{	matrix w1i = entry1.second(iAtom,iAtom+1, 0,xExpect.size());
				matrix w2i = entry2.second(iAtom,iAtom+1, 0,xExpect.size());
				matrix w = transpose(w1i) * w2i;
				for(int iVector=0; iVector<3; iVector++)
				{	callPref(eblas_zmul)(w.nData(), w.dataPref(), 1, HePhData, 1);
					HePhData += w.nData();
				}
			}
		if(realPartOnly)
		{	if(mpiWorld->isHead()) HePh.write_real(fp);
			nrm2totSq += std::pow(nrm2(HePh), 2); 
			nrm2imSq += std::pow(callPref(eblas_dnrm2)(HePh.nData(), ((double*)HePh.dataPref())+1, 2), 2); //look only at imaginary parts with a stride of 2
		}
		else { if(mpiWorld->isHead()) HePh.write(fp); }
	}
	if(mpiWorld->isHead()) fclose(fp);
	if(realPartOnly)
		logPrintf("done. Relative discarded imaginary part: %le\n", sqrt(nrm2imSq / nrm2totSq));
	else
		logPrintf("done.\n");
}
