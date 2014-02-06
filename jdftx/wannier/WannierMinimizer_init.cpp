/*-------------------------------------------------------------------
Copyright 2014 Ravishankar Sundararaman

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
#include <core/DataIO.h>
#include <electronic/Blip.h>
#include <electronic/operators.h>

//Find a finite difference formula given a list of relative neighbour positions (in cartesian coords)
//[Appendix B of Phys Rev B 56, 12847]
//Returns an empty weight set on failure
std::vector<double> getFDformula(const std::vector< vector3<> >& b)
{	//Group elements of b into shells:
	std::vector<unsigned> shellMax; //cumulative count within each shell
	for(unsigned i=1; i<b.size(); i++)
		if(b[i].length() > b[i-1].length() + symmThreshold)
			shellMax.push_back(i);
	shellMax.push_back(b.size());
	int nShells = shellMax.size();
	//Setup the equations satisfied by the weights:
	const int nEquations = 6;
	matrix Lhs = zeroes(nEquations, nShells);
	complex* LhsData = Lhs.data();
	for(int s=0; s<nShells; s++)
		for(unsigned j = (s ? shellMax[s-1] : 0); j<shellMax[s]; j++)
		{	//Equations from ref.:
			int iEqn = 0;
			//Rank-two sum is identity:
			LhsData[Lhs.index(iEqn++, s)] += b[j][0]*b[j][0];
			LhsData[Lhs.index(iEqn++, s)] += b[j][1]*b[j][1];
			LhsData[Lhs.index(iEqn++, s)] += b[j][2]*b[j][2];
			LhsData[Lhs.index(iEqn++, s)] += b[j][1]*b[j][2];
			LhsData[Lhs.index(iEqn++, s)] += b[j][2]*b[j][0];
			LhsData[Lhs.index(iEqn++, s)] += b[j][0]*b[j][1];
		}
	matrix rhs = zeroes(nEquations, 1);
	for(unsigned i=0; i<3; i++) rhs.data()[i] = 1; //first three components = diagonals of rank-two sum
	//Solve using a singular value decomposition:
	matrix U, Vdag; diagMatrix S;
	Lhs.svd(U, S, Vdag);
	for(double& s: S) s = (s<symmThreshold) ? 0. : 1./s; //invert and zero out small singular values
	matrix wShells = dagger(Vdag) * S * dagger(U(0,nEquations, 0,nShells)) * rhs;
	if(nrm2(Lhs * wShells - rhs) > symmThreshold) //check solution by substitution
		return std::vector<double>(); //Not an exact solution, so quit (and try with more shells)
	//Store the weights in the original indexing:
	complex* wShellsData = wShells.data();
	std::vector<double> w(b.size());
	for(int s=0; s<nShells; s++)
		for(unsigned j = (s ? shellMax[s-1] : 0); j<shellMax[s]; j++)
			w[j] = wShellsData[s].real();
	return w;
}

//Helper function for PeriodicLookup<WannierMinimizer::Kpoint> used in WannierMinimizer::WannierMinimizer
inline vector3<> getCoord(const WannierMinimizer::Kpoint& kpoint) { return kpoint.k; }

WannierMinimizer::WannierMinimizer(const Everything& e, const Wannier& wannier) : e(e), wannier(wannier), sym(e.symm.getMatrices()),
	nCenters(wannier.trialOrbitals.size()), nBands(e.eInfo.nBands),
	nSpins(e.eInfo.spinType==SpinNone ? 1 : 2), qCount(e.eInfo.qnums.size()/nSpins),
	needSuper(wannier.saveWfns || wannier.saveWfnsRealSpace || wannier.numericalOrbitalsFilename.length())
{
	//Create supercell grid:
	logPrintf("\n---------- Initializing supercell grid for Wannier functions ----------\n");
	const Supercell& supercell = *(e.coulombParams.supercell);
	gInfoSuper.R = supercell.Rsuper;
	gInfoSuper.GmaxRho = e.gInfo.Gmax;
	if(needSuper) gInfoSuper.initialize(true, sym);

	//Initialize cell map (for matrix element output):
	WignerSeitz ws(gInfoSuper.R);
	double Rmax = ws.circumRadius();
	vector3<int> iCellMax; //bounding box of supercell Wigner-Seitz circumradius (in unit-cell lattice coordinates)
	for(int l=0; l<3; l++)
		iCellMax[l] = 1 + int(ceil(Rmax * e.gInfo.invR.row(l).length()));
	//--- collect all lattice vectors on or inside ws:
	matrix3<> superInv = inv(matrix3<>(supercell.super));
	vector3<int> iCell;
	iCellMap.clear();
	std::list< vector3<int> > iCellSurface; //list of surface cells
	for(iCell[0]=-iCellMax[0]; iCell[0]<=iCellMax[0]; iCell[0]++)
	for(iCell[1]=-iCellMax[1]; iCell[1]<=iCellMax[1]; iCell[1]++)
	for(iCell[2]=-iCellMax[2]; iCell[2]<=iCellMax[2]; iCell[2]++)
	{	vector3<> xCell = superInv * iCell;
		if(ws.onBoundary(xCell)) iCellSurface.push_back(iCell); //yet to determine multiplicity of surface cells
		else if(ws.boundaryDistance(xCell)>0) iCellMap[iCell] = 1.; //interior cells are unique
	}
	//--- determine multiplicity of surface cells and move them to iCellMap
	for(auto iter=iCellSurface.begin(); iter!=iCellSurface.end(); iter++)
	{	std::vector< vector3<int> > equiv;
		auto iter2=iter; iter2++;
		while(iter2!=iCellSurface.end())
		{	vector3<> dx = superInv * (*iter2 - *iter); //difference in super-lattice coordinates
			bool isEquiv = true;
			for(int l=0; l<3; l++) //each entry of dx must be integral for equivalency
				isEquiv &= (fabs(dx[l]-round(dx[l]))<symmThreshold);
			if(isEquiv)
			{	equiv.push_back(*iter2);
				iter2 = iCellSurface.erase(iter2);
			}
			else iter2++;
		}
		double weight = 1./(1+equiv.size());
		iCellMap[*iter] = weight;
		for(vector3<int> iCell: equiv)
			iCellMap[iCell] = weight;
	}
	//--- check that the weights add up
	{	double weightSum = 0.;
		for(const auto& entry: iCellMap)
			weightSum += entry.second;
		assert(fabs(weightSum / supercell.kmesh.size() - 1.) < symmThreshold);
	}
	//--- write the cell map (for post-processing programs to use)
	if(mpiUtil->isHead())
	{	string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfCellMap");
		logPrintf("Writing '%s' ... ", fname.c_str()); logFlush();
		FILE* fp = fopen(fname.c_str(), "w");
		fprintf(fp, "#i0 i1 i2  x y z  (integer lattice combinations, and cartesian offsets)\n");
		for(const auto& entry: iCellMap)
		{	const vector3<int>& i = entry.first;
			vector3<> r = e.gInfo.R * i;
			fprintf(fp, "%+2d %+2d %+2d  %+11.6lf %+11.6lf %+11.6lf\n", i[0], i[1], i[2], r[0], r[1], r[2]);
		}
		fclose(fp);
		logPrintf("done.\n"); logFlush();
	}
	
	//Determine and output the band ranges:
	for(int iSpin=0; iSpin<nSpins; iSpin++)
	{	std::vector<double> eMin(nBands, DBL_MAX), eMax(nBands, -DBL_MAX);
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++) if(e.eInfo.qnums[q].spin==iSpin)
			for(int b=0; b<nBands; b++)
			{	eMin[b] = std::min(eMin[b], e.eVars.Hsub_eigs[q][b]);
				eMax[b] = std::max(eMax[b], e.eVars.Hsub_eigs[q][b]);
			}
		mpiUtil->allReduce(eMin.data(), eMin.size(), MPIUtil::ReduceMin);
		mpiUtil->allReduce(eMax.data(), eMax.size(), MPIUtil::ReduceMax);
		if(mpiUtil->isHead())
		{	string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfBandRanges", &iSpin);
			logPrintf("Writing '%s' ... ", fname.c_str()); logFlush();
			FILE* fp = fopen(fname.c_str(), "w");
			for(int b=0; b<nBands; b++)
				fprintf(fp, "%+10.5lf %+10.5lf\n", eMin[b], eMax[b]);
			fclose(fp);
			logPrintf("done.\n"); logFlush();
		}
	}
	
	//Determine finite difference formula:
	logPrintf("Setting up finite difference formula on k-mesh ... "); logFlush();
	matrix3<> kbasis = e.gInfo.GT * inv(~matrix3<>(supercell.super)); //basis vectors in reciprocal space for the k-mesh (in cartesian coords)
	std::multimap<double, vector3<> > dkMap; //cartesian offsets from one k-point to other k-points sorted by distance
	vector3<int> ik;
	const int ikBound = 2;
	for(ik[0]=-ikBound; ik[0]<=+ikBound; ik[0]++)
	for(ik[1]=-ikBound; ik[1]<=+ikBound; ik[1]++)
	for(ik[2]=-ikBound; ik[2]<=+ikBound; ik[2]++)
		if(ik.length_squared()) //ignore self
		{	vector3<> dk = kbasis * ik;
			dkMap.insert(std::make_pair<>(dk.length_squared(), dk));
		}
	//--- remove inversion partners (handled by symmetry)
	for(auto iter=dkMap.begin(); iter!=dkMap.end(); iter++)
	{	auto iter2=iter; iter2++;
		while(iter2!=dkMap.end() && iter2->first < iter->first + symmThreshold)
		{	if((iter2->second + iter->second).length_squared() < symmThresholdSq)
				iter2 = dkMap.erase(iter2);
			else iter2++;
		}
	}
	//--- find FD formula shell by shell
	std::vector< vector3<> > b; //list of cartesian offsets to corresponding neighbours
	std::vector<double> wb; //corresponding weights in finite difference formula
	for(auto iter=dkMap.begin(); iter!=dkMap.end(); )
	{	//Add all the neighbours with equivalent distances:
		while(true)
		{	b.push_back(iter->second);
			double prevDist = iter->first;
			iter++;
			if(iter==dkMap.end() || //end of neighbour list (should not be encountered)
				iter->first > prevDist+symmThreshold) //next neighbour is further away beyond tolerance
				break;
		}
		//Check if this list of neighbours is sufficient to get a finite difference formula
		wb = getFDformula(b);
		if(wb.size() == b.size()) break; //success
	}
	if(!wb.size()) die("failed.\n");
	logPrintf("found a %lu neighbour formula.\n", 2*wb.size());
	
	//Create a list of kpoints:
	const std::vector<QuantumNumber>& qnums = e.eInfo.qnums;
	std::vector<Kpoint> kpoints(supercell.kmeshTransform.size());
	for(size_t i=0; i<kpoints.size(); i++)
	{	const Supercell::KmeshTransform& src = supercell.kmeshTransform[i];
		Kpoint& kpoint = kpoints[i];
		//Copy over base class KMeshTransform:
		(Supercell::KmeshTransform&)kpoint = src;
		//Initialize base class QuantumNumber:
		kpoint.k = (~sym[src.iSym]) * (src.invert * qnums[src.iReduced].k) + src.offset;
		kpoint.weight = 1./kpoints.size();
		kpoint.spin = 0;
	}
	PeriodicLookup<WannierMinimizer::Kpoint> plook(kpoints, e.gInfo.GGT); //look-up table for O(1) fuzzy searching
	
	//Determine overall Bloch wavevector of supercell (if any)
	if(needSuper)
	{	qnumSuper.weight = 1.;
		qnumSuper.spin = 0.;
		qnumSuper.k = kpoints[0].k * supercell.super;
		for(int l=0; l<3; l++)
			qnumSuper.k[l] -= floor(0.5+qnumSuper.k[l]);
		if(qnumSuper.k.length_squared()>symmThresholdSq)
		{	logPrintf("WARNING: k-mesh does not contain Gamma point. Orbitals will not be strictly periodic on supercell,\n"
				"\tbecause of an overall Bloch wave-vector: ");
			qnumSuper.k.print(globalLog, " %lf ");
		}
	}
	
	//Determine distribution amongst processes:
	ikStart = (kpoints.size() * mpiUtil->iProcess()) / mpiUtil->nProcesses();
	ikStop = (kpoints.size() * (mpiUtil->iProcess()+1)) / mpiUtil->nProcesses();
	ikStopArr.resize(mpiUtil->nProcesses());
	for(int iProc=0; iProc<mpiUtil->nProcesses(); iProc++)
		ikStopArr[iProc] = (kpoints.size() * (iProc+1)) / mpiUtil->nProcesses();

	//Create the mesh structure with neighbour info:
	kMesh.resize(kpoints.size());
	for(size_t i=0; i<kMesh.size(); i++) //FD formula needed on all nodes for initial matrix generation
	{	//Store the k-point with its FD formula in kMesh
		KmeshEntry& kMeshEntry = kMesh[i];
		kMeshEntry.point = kpoints[i];
		addIndex(kMeshEntry.point);
		for(unsigned j=0; j<wb.size(); j++)
		{	EdgeFD edge;
			edge.wb = wb[j];
			edge.b = b[j];
			//Find neighbour:
			vector3<> kj = kpoints[i].k + inv(e.gInfo.GT) * b[j];
			edge.ik = plook.find(kj);
			edge.point = kpoints[edge.ik];
			for(int l=0; l<3; l++)
				edge.point.offset[l] += int(round(kj[l] - edge.point.k[l])); //extra offset
			edge.point.k = kj;
			addIndex(edge.point);
			kMeshEntry.edge.push_back(edge);
		}
	}
	
	//Create the common reduced basis set (union of all the reduced bases)
	//Collect all referenced full-G indices
	std::set<int> commonSet, commonSuperSet;
	for(auto index: indexMap)
		for(int j=0; j<index.second->nIndices; j++)
		{	commonSet.insert(index.second->data[j]);
			if(needSuper)
				commonSuperSet.insert(index.second->dataSuper[j]);
		}
	//Convert to a Basis object, and create inverse map
	std::vector<int> indexCommon(commonSet.size());
	std::map<int,int> commonInverseMap;
	auto setIter = commonSet.begin();
	for(unsigned j=0; j<indexCommon.size(); j++)
	{	int i = *(setIter++);
		indexCommon[j] = i;
		commonInverseMap[i] = j;
	}
	basis.setup(e.gInfo, e.iInfo, indexCommon);
	//Liekwise for supercell:
	std::vector<int> indexSuperCommon(commonSuperSet.size());
	std::map<int,int> commonSuperInverseMap;
	if(needSuper)
	{	auto superSetIter = commonSuperSet.begin();
		for(unsigned j=0; j<indexSuperCommon.size(); j++)
		{	int i = *(superSetIter++);
			indexSuperCommon[j] = i;
			commonSuperInverseMap[i] = j;
		}
		basisSuper.setup(gInfoSuper, e.iInfo, indexSuperCommon);
	}
	//Update indexMap to point to common reduced basis instead of full G:
	for(auto mapEntry: indexMap)
	{	Index& index = *mapEntry.second;
		for(int j=0; j<index.nIndices; j++)
		{	index.data[j] = commonInverseMap[index.data[j]];
			if(needSuper)
				index.dataSuper[j] = commonSuperInverseMap[index.dataSuper[j]];
		}
		index.set();
	}
	
	//Read numerical trial orbitals if provided:
	if(wannier.numericalOrbitalsFilename.length())
	{	logPrintf("\n--- Initializing numerical trial orbitals ---\n");
		
		//Read the header:
		string fname = wannier.numericalOrbitalsFilename + ".header";
		std::ifstream ifs(fname.c_str());
		if(!ifs.is_open()) die("Could not open '%s' for reading.\n", fname.c_str());
		string comments;
		//--- Line 1:
		int nCols; size_t colLength;
		ifs >> nCols >> colLength;
		getline(ifs, comments);
		logPrintf("%d orbitals with %lu basis elements each.\n", nCols, colLength);
		//--- Line 2:
		matrix3<> GT;
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				ifs >> GT(i,j);
		getline(ifs, comments);
		//--- Basis elements:
		std::vector< vector3<int> > iGarr(colLength);
		vector3<int> hlfSin;
		for(vector3<int>& iG: iGarr)
		{	ifs >> iG[0] >> iG[1] >> iG[2];
			for(int k=0; k<3; k++) hlfSin[k] = std::max(hlfSin[k], abs(iG[k]));
		}
		ifs.close();
		
		//Create input grid:
		GridInfo gInfoIn;
		for(int k=0; k<3; k++)
		{	gInfoIn.S[k] = 2*(2*hlfSin[k] + 1);
			while(!fftSuitable(gInfoIn.S[k])) gInfoIn.S[k] += 2; 
		}
		gInfoIn.R = 2*M_PI * inv(~GT);
		gInfoIn.initialize(true);
		
		//Create input basis:
		Basis basisIn;
		std::vector<int> indexIn(colLength);
		for(size_t i=0; i<colLength; i++)
			indexIn[i] = gInfoIn.fullGindex(iGarr[i]);
		basisIn.setup(gInfoIn, e.iInfo, indexIn);
		
		//Read the wavefunctions:
		fname = wannier.numericalOrbitalsFilename;
		logPrintf("Reading numerical orbitals from '%s' ... ", fname.c_str()); logFlush();
		ColumnBundle Cin(nCols, colLength, &basisIn, &qnumSuper);
		Cin.read(fname.c_str());
		logPrintf("done.\n"); logFlush();
		
		//Convert wavefunctions to supercell basis:
		logPrintf("Converting numerical orbitals to supercell basis ... "); logFlush();
		//--- offset in input
		Cin = translate(Cin, -wannier.numericalOrbitalsOffset);
		//--- resample band-by-band"
		logSuspend();
		BlipResampler resample(gInfoIn, gInfoSuper);
		logResume();
		ColumnBundle C(nCols, basisSuper.nbasis, &basisSuper, &qnumSuper, isGpuEnabled());
		for(int b=0; b<nCols; b++)
			C.setColumn(b, J(resample(Cin.getColumn(b))));
		logPrintf("done.\n"); logFlush();
		
		//Split supercell wavefunction into kpoints:
		logPrintf("Dividing supercell numerical orbitals to k-points ... "); logFlush();
		{	ColumnBundle temp(1, basis.nbasis, &basis, 0, isGpuEnabled());
			for(size_t ik=0; ik<kMesh.size(); ik++) if(isMine_q(ik,0) || isMine_q(ik,1))
			{	const KmeshEntry& ki = kMesh[ik];
				const Index& index = *(indexMap.find(ki.point)->second);
				auto Ck = std::make_shared<ColumnBundle>(nCols, basis.nbasis, &basis, &ki.point, isGpuEnabled());
				Ck->zero();
				for(int b=0; b<nCols; b++)
				{	temp.zero();
					eblas_gather_zdaxpy(index.nIndices, 1., index.dataSuperPref, C.dataPref()+C.index(b,0), temp.dataPref());
					eblas_scatter_zdaxpy(index.nIndices, 1./ki.point.weight, index.dataPref, temp.dataPref(), Ck->dataPref()+Ck->index(b,0));
				}
				numericalOrbitals[ki.point] = Ck;
			}
		}
		logPrintf("done.\n"); logFlush();
	}
}
