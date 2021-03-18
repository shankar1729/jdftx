/*-------------------------------------------------------------------
Copyright 2021 Ravishankar Sundararaman

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

//Shared helper for defect and phonon matrix element Wannierization routines:
void WannierMinimizer::initializeCommensurate(const vector3<int>& sup, 
	const std::vector<vector3<>>& xCenters, string suffix, int iSpin, int& prodSup,
	std::vector<int>& ikArr, std::map<vector3<int>,matrix>& cellMap, std::vector<UniqueCell>& uniqueCells)
{
	//Generate list of commensurate k-points in order present in the unit cell calculation
	prodSup = sup[0] * sup[1] * sup[2];
	ikArr.clear();
	for(unsigned ik=0; ik<kMesh.size(); ik++)
	{	vector3<> kSup = kMesh[ik].point.k * Diag(sup);
		double roundErr; round(kSup, &roundErr);
		if(roundErr < symmThreshold) //integral => commensurate with supercell
			ikArr.push_back(ik);
	}
	assert(int(ikArr.size()) == prodSup);
	
	//Generate and write cell map:
	cellMap = getCellMap(
		e.gInfo.R, e.gInfo.R * Diag(sup),
		e.coulombParams.isTruncated(), xCenters, xExpect, wannier.rSmooth,
		wannier.getFilename(Wannier::FilenameDump, "mlwfCellMap"+suffix, &iSpin));
	//--- Corresponding cell weights:
	if(mpiWorld->isHead())
	{	string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfCellWeights"+suffix, &iSpin);
		logPrintf("Dumping '%s'... ", fname.c_str()); logFlush();
		FILE* fp = fopen(fname.c_str(), "w");
		if(!fp) die_alone("could not open file for writing.\n");
		for(const auto& entry: cellMap)
			entry.second.write_real(fp);
		fclose(fp);
		logPrintf("done.\n"); logFlush();
	}
	//--- Re-organize by unique cells in supercell:
	uniqueCells.assign(prodSup, UniqueCell());
	vector3<int> stride(sup[1]*sup[2], sup[2], 1);
	for(const auto& entry: cellMap)
	{	//Compute unique index:
		vector3<int> iR = entry.first;
		for(int iDir=0; iDir<3; iDir++)
			iR[iDir] = positiveRemainder(iR[iDir], sup[iDir]);
		int uniqIndex = dot(iR, stride);
		//Set/collect unique cell properties:
		uniqueCells[uniqIndex].iR = iR;
		uniqueCells[uniqIndex].cells.push_back(entry);
	}
}


//Output defect-related quantities in Wannier basis:
void WannierMinimizer::saveMLWF_defect(int iSpin, DefectSupercell& ds)
{
	//Get commensurate k-points, cell map and unqiue cells for defect:
	std::vector<vector3<>> xCenter(1, Diag(vector3<>(ds.supIn)) * ds.xCenter); //in unit cell coords
	int prodSup;
	std::vector<int> ikArr;
	std::map<vector3<int>,matrix> defectCellMap;
	std::vector<UniqueCell> uniqueCells;
	initializeCommensurate(ds.supOut, xCenter, "D_"+ds.name, iSpin, prodSup, ikArr, defectCellMap, uniqueCells);
	
	//Pairs of commensurate k-points in order they will be stored:
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
	int ikArrStart, ikArrStop; //MPI division for working on commensurate k-points
	TaskDivision(prodSup, mpiWorld).myRange(ikArrStart, ikArrStop);
	int iPairStart=ikArrStart*prodSup, iPairStop=ikArrStop*prodSup; //Divide by ik1 and handle all ik2 for each
	
	//Compute defect matrix elements in reciprocal space:
	int nPairsMine = std::max(1, iPairStop-iPairStart); //avoid zero size matrices below
	matrix HDtilde = zeroes(nCenters*nCenters, nPairsMine);
	//TODO
	
	//Wannierize and output one cell fixed at a time to minimize memory usage:
	string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfHD_"+ds.name, &iSpin);
	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
	FILE* fp = 0;
	if(mpiWorld->isHead())
	{	fp = fopen(fname.c_str(), "wb");
		if(!fp) die_alone("Error opening %s for writing.\n", fname.c_str());
	}
	matrix phase = zeroes(HDtilde.nCols(), prodSup);
	double kPairWeight = 1./prodSup;
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
		matrix H = HDtilde * phase;
		mpiWorld->allReduceData(H, MPIUtil::ReduceSum);
		//write results for unique cells to file:
		if(realPartOnly)
		{	if(mpiWorld->isHead()) H.write_real(fp);
			nrm2totSq += std::pow(nrm2(H), 2); 
			nrm2imSq += std::pow(callPref(eblas_dnrm2)(H.nData(), ((double*)H.dataPref())+1, 2), 2); //look only at imaginary parts with a stride of 2
		}
		else { if(mpiWorld->isHead()) H.write(fp); }
	}
	if(mpiWorld->isHead()) fclose(fp);
	if(realPartOnly)
		logPrintf("done. Relative discarded imaginary part: %le\n", sqrt(nrm2imSq / nrm2totSq));
	else
		logPrintf("done.\n");
}
