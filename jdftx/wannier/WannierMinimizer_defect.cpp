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

//Shared helpers for defect and phonon matrix element Wannierization routines:
void WannierMinimizer::initializeCommensurate(const vector3<int>& sup, int& prodSup, std::vector<int>& ikArr, vector3<> offset)
{
	//Generate list of commensurate k-points (up to offset) in order present in the unit cell calculation
	prodSup = sup[0] * sup[1] * sup[2];
	ikArr.clear();
	for(unsigned ik=0; ik<kMesh.size(); ik++)
	{	vector3<> kSup = (kMesh[ik].point.k - offset) * Diag(sup);
		double roundErr; round(kSup, &roundErr);
		if(roundErr < symmThreshold) //integral => commensurate with supercell
			ikArr.push_back(ik);
	}
	assert(int(ikArr.size()) == prodSup);
}

void WannierMinimizer::initializeCellMaps(
	const matrix3<>& R, const vector3<int>& sup, const vector3<bool>& isTruncated,
	const std::vector<vector3<>>& xCenters, string suffix, int iSpin,
	std::map<vector3<int>,matrix>& cellMap, std::vector<UniqueCell>& uniqueCells)
{
	//Generate and write cell map:
	cellMap = getCellMap(
		R, R * Diag(sup), isTruncated, xCenters, xExpect, wannier.rSmooth,
		wannier.getFilename(Wannier::FilenameDump, "mlwfCellMap"+suffix, &iSpin));
	
	//Corresponding cell weights:
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
	
	//Re-organize by unique cells in supercell:
	uniqueCells.assign(sup[0]*sup[1]*sup[2], UniqueCell());
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


//Make the columns of R for which mask=true orthogonal to those with mask=false
matrix3<> makeOrthogonal(const matrix3<>& R, const vector3<bool>& mask)
{	//Construct projection matrix removing directions with mask=false:
	matrix3<> proj(1, 1, 1); //set to identity
	for(int iDir=0; iDir<3; iDir++)
		if(not mask[iDir])
		{	vector3<> jHat = normalize(proj * R.column(iDir));
			proj -= outer(jHat, jHat);
		}
	//Apply to mask=true directions:
	matrix3<> result;
	for(int iDir=0; iDir<3; iDir++)
	{	vector3<> Ri = R.column(iDir);
		result.set_col(iDir, (mask[iDir] ? (proj * Ri) : Ri));
	}
	return result;
}

//Output defect-related quantities in Wannier basis:
void WannierMinimizer::saveMLWF_defect(int iSpin, DefectSupercell& ds)
{
	//Identify any "perp" periodic directions (for line/plane defects):
	vector3<bool> isTruncated = e.coulombParams.isTruncated();
	vector3<bool> notPerp(false, false, false);
	vector3<int> perp; //perpendicular supercell size
	for(int iDir=0; iDir<3; iDir++)
	{	if((ds.supOut[iDir] == 1) and not isTruncated[iDir])
		{	//Periodic "perp" direction
			perp[iDir] = e.eInfo.kfold[iDir];
			isTruncated[iDir] = true; //make par cellmap orthogonal to this direction
		}
		else
		{	//Translation-symmetry-broken "par" direction
			perp[iDir] = 1;
			notPerp[iDir] = true; //make perp cellmap orthogonal to this direction
		}
	}
	matrix3<> Rperp = makeOrthogonal(e.gInfo.R, notPerp);
	matrix3<> Rpar = makeOrthogonal(e.gInfo.R, isTruncated);
	
	//Get commensurate k-points, cell map and unqiue cells for defect:
	std::vector<vector3<>> xCenter(1, Diag(vector3<>(ds.supIn)) * ds.xCenter); //in unit cell coords
	std::map<vector3<int>,matrix> defectCellMap, defectCellMapPerp;
	std::vector<UniqueCell> uniqueCells, uniqueCellsPerp;
	initializeCellMaps(Rpar, ds.supOut, isTruncated, xCenter, "D_"+ds.name, iSpin, defectCellMap, uniqueCells);
	initializeCellMaps(Rperp, perp, notPerp, xExpect, "PerpD_"+ds.name, iSpin, defectCellMapPerp, uniqueCellsPerp);
	
	//Prepare offsets of the supercell-commensurate k-mesh along periodic directions:
	matrix3<> perpInv = Diag(vector3<>(1.0/perp[0], 1.0/perp[1], 1.0/perp[2]));
	std::vector<vector3<>> offsets;
	for(int iPerp0=0; iPerp0<perp[0]; iPerp0++)
	for(int iPerp1=0; iPerp1<perp[1]; iPerp1++)
	for(int iPerp2=0; iPerp2<perp[2]; iPerp2++)
		offsets.push_back(perpInv * vector3<>(iPerp0, iPerp1, iPerp2));
	int prodPerp = offsets.size();
	
	//Divide supercell-commensurate reciprocal space over MPI:
	int prodSup = ds.supOut[0] * ds.supOut[1] * ds.supOut[2];
	int ikArrStart, ikArrStop; //MPI division for working on commensurate k-points
	TaskDivision ikArrDiv(prodSup, mpiWorld);
	ikArrDiv.myRange(ikArrStart, ikArrStop);
	int iPairStart=ikArrStart*prodSup, iPairStop=ikArrStop*prodSup; //Divide by ik1 and handle all ik2 for each
	int nPairsMine = std::max(1, iPairStop-iPairStart); //avoid zero size matrices below
	int nPairsPerpMine = nPairsMine * prodPerp;
	int nPairsPerpInterval = std::max(1, int(round(nPairsPerpMine/20.))); //interval for reporting progress
	int nPairsPerpDone = 0;
	
	//Compute defect matrix elements in supercell-commensurate reciprocal-space-squared, for each periodic offset:
	logPrintf("Computing matrix elements for defect '%s' ...  ", ds.name.c_str()); logFlush(); 
	std::vector<std::vector<vector3<>>> kArrs(prodPerp);
	std::vector<matrix> HDtildes(prodPerp);
	int nCentersSq = nCenters * nCenters;
	for(int iOffset=0; iOffset<prodPerp; iOffset++)
	{	//Find commensurate k-points
		std::vector<int> ikArr;
		initializeCommensurate(ds.supOut, prodSup, ikArr, offsets[iOffset]);
		for(int ik: ikArr) kArrs[iOffset].push_back(kMesh[ik].point.k);
		
		//Compute defect matrix elements in reciprocal space:
		matrix& HDtilde = HDtildes[iOffset];
		HDtilde = zeroes(nCentersSq, nPairsMine);
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
				nPairsPerpDone++;
				if(nPairsPerpDone % nPairsPerpInterval == 0)
				{	logPrintf("%d%% ", int(round(nPairsPerpDone*100./nPairsPerpMine)));
					logFlush();
				};
			}
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
	//--- initialize Fourier transform phase for perp directions
	matrix phasePerp = zeroes(prodPerp, prodPerp);
	for(int iCellPerp=0; iCellPerp<prodPerp; iCellPerp++)
		for(int iOffset=0; iOffset<prodPerp; iOffset++)
			phasePerp.set(iOffset, iCellPerp, cis(-2*M_PI*dot(offsets[iOffset], uniqueCellsPerp[iCellPerp].iR)));
	matrix phase = zeroes(nPairsMine, prodSup); //parallel-direction fourier transform phase
	double transformWeight = 1./(prodSup*prodSup*prodPerp); //double parallel and single perp Fourier transform weight
	double nrm2totSq = 0., nrm2imSq = 0.;
	for(const UniqueCell& cell1: uniqueCells)
	{	matrix H(nCentersSq*prodSup, prodPerp); //Real-space in parallel, but reciprocal space in perp
		for(int iOffset=0; iOffset<prodPerp; iOffset++)
		{	//Calculate parallel-direction Fourier transform phase (with integration weights):
			int iPairMine = 0;
			const std::vector<vector3<>>& kArr = kArrs[iOffset];
			for(int ikIndex1=ikArrStart; ikIndex1<ikArrStop; ikIndex1++)
			{	const vector3<>& k1 = kArr[ikIndex1];
				for(const vector3<>& k2: kArr)
				{	int iCell2 = 0;
					for(const UniqueCell& cell2: uniqueCells)
						phase.set(iPairMine, iCell2++, 
							transformWeight * cis(2*M_PI*(dot(k1,cell1.iR) - dot(k2,cell2.iR))) );
					iPairMine++;
				}
			}
			//Apply parallel-direction Fourier transform
			matrix Hi = HDtildes[iOffset] * phase;
			mpiWorld->allReduceData(Hi, MPIUtil::ReduceSum);
			callPref(eblas_copy)(H.dataPref() + H.index(0,iOffset), Hi.dataPref(), Hi.nData());
		}
		//Switch to real space in perp direction as well:
		if(prodPerp > 1)
		{	H = H * phasePerp;
			//Swap data order to be cell2, cellPerp, band1, band2
			matrix Hswapped(nCentersSq, prodPerp * prodSup);
			for(int iCell2=0; iCell2<prodSup; iCell2++)
				for(int iCellPerp=0; iCellPerp<prodPerp; iCellPerp++)
					callPref(eblas_copy)(
						Hswapped.dataPref() + nCentersSq*(prodPerp*iCell2+iCellPerp),
						H.dataPref() + nCentersSq*(prodSup*iCellPerp+iCell2), nCentersSq);
			std::swap(H, Hswapped);
		}
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
