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
	TaskDivision ikArrDiv(prodSup, mpiWorld);
	ikArrDiv.myRange(ikArrStart, ikArrStop);
	int iPairStart=ikArrStart*prodSup, iPairStop=ikArrStop*prodSup; //Divide by ik1 and handle all ik2 for each
	
	//Compute defect matrix elements in reciprocal space:
	int nPairsMine = std::max(1, iPairStop-iPairStart); //avoid zero size matrices below
	matrix HDtilde = zeroes(nCenters*nCenters, nPairsMine);
	//--- get wavefunctions for each commensurate k (split by ikArr as above):
	std::vector<ColumnBundle> C(prodSup);
	std::vector<MPIUtil::Request> requests;
	for(int ikIndex=0; ikIndex<prodSup; ikIndex++)
	{	int ik = ikArr[ikIndex];
		ColumnBundle& Ck = C[ikIndex];
		if(isMine_q(ik, iSpin))
		{	Ck = getWfns(kMesh[ik].point, iSpin);
			int dest = ikArrDiv.whose(ikIndex);
			if(dest != mpiWorld->iProcess())
			{	requests.push_back(MPIUtil::Request());
				mpiWorld->sendData(Ck, dest, ikIndex, &requests.back());
			}
		}
		else if(ikArrDiv.isMine(ikIndex))
		{	int src = whose_q(ik, iSpin);
			Ck.init(nBands, basis.nbasis*nSpinor, &basis, &kMesh[ik].point, isGpuEnabled());
			requests.push_back(MPIUtil::Request());
			mpiWorld->recvData(Ck, src, ikIndex, &requests.back());
		}
	}
	mpiWorld->waitAll(requests);
	std::vector<DefectSupercell::CachedProjections> proj(prodSup);
	for(int ikIndex=0; ikIndex<prodSup; ikIndex++)
		if(ikArrDiv.isMine(ikIndex))
			ds.project(C[ikIndex], proj[ikIndex]); //cache projections
		else
			C[ikIndex] = 0; //clean up un-needed
	//--- loop over all ikIndex2
	logPrintf("Computing matrix elements for defect '%s' ...  ", ds.name.c_str()); logFlush(); 
	int nPairsInterval = std::max(1, int(round(nPairsMine/20.))); //interval for reporting progress
	int nPairsDone = 0;
	for(int ikIndex2=0; ikIndex2<prodSup; ikIndex2++)
	{	int ik2 = ikArr[ikIndex2];
		//Make C2 available on all processes:
		ColumnBundle C2; DefectSupercell::CachedProjections proj2;
		if(ikArrDiv.isMine(ikIndex2))
		{	C2 = C[ikIndex2];
			proj2 = proj[ikIndex2];
		}
		else C2.init(nBands, basis.nbasis*nSpinor, &basis, &kMesh[ik2].point, isGpuEnabled());
		mpiWorld->bcastData(C2, ikArrDiv.whose(ikIndex2));
		ds.bcast(proj2, ikArrDiv.whose(ikIndex2));
		//Compute matrix element with all local C1:
		for(int ikIndex1=ikArrStart; ikIndex1<ikArrStop; ikIndex1++)
		{	int ik1 = ikArr[ikIndex1];
			//Compute matrix elements in eigenbasis and apply Wannier rotations:
			const ColumnBundle& C1 = C[ikIndex1];
			const DefectSupercell::CachedProjections& proj1 = proj[ikIndex1];
			matrix HDcur = dagger(kMesh[ik1].U) * ds.compute(C1, C2, proj1, proj2) * kMesh[ik2].U;
			//Store at appropriate location in global array:
			int iPairMine = (ikIndex1-ikArrStart)*prodSup + ikIndex2;
			callPref(eblas_copy)(HDtilde.dataPref() + HDtilde.index(0,iPairMine), HDcur.dataPref(), HDcur.nData());
			//Print progress:
			nPairsDone++;
			if(nPairsDone % nPairsInterval == 0)
			{	logPrintf("%d%% ", int(round(nPairsDone*100./nPairsMine)));
				logFlush();
			};
		}
	}
	logPrintf("done.\n"); logFlush();
	
	//Wannierize and output one cell fixed at a time to minimize memory usage:
	string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfHD_"+ds.name, &iSpin);
	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
	FILE* fp = 0;
	if(mpiWorld->isHead())
	{	fp = fopen(fname.c_str(), "wb");
		if(!fp) die_alone("Error opening %s for writing.\n", fname.c_str());
	}
	matrix phase = zeroes(HDtilde.nCols(), prodSup);
	double kPairWeight = 1./(prodSup*prodSup);
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
