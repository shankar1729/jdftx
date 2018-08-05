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

#include <wannier/WannierMinimizerFD.h>

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
	matrix wShells = (nShells < nEquations)
			? ( dagger(Vdag) * S * dagger(U(0,nEquations, 0,nShells)) * rhs )
			: ( dagger(Vdag(0,nEquations, 0,nShells)) * S * dagger(U) * rhs );
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

//Helper function for PeriodicLookup<WannierMinimizer::KmeshEntry> used in WannierMinimizerFD constructor
inline vector3<> getCoord(const WannierMinimizer::KmeshEntry& ke) { return ke.point.k; }

WannierMinimizerFD::WannierMinimizerFD(const Everything& e, const Wannier& wannier)
: WannierMinimizer(e, wannier)
{
	//Determine finite difference formula:
	logPrintf("Setting up finite difference formula on k-mesh ... "); logFlush();
	const Supercell& supercell = *(e.coulombParams.supercell);
	matrix3<> kbasis = e.gInfo.GT * inv(~matrix3<>(supercell.super)); //basis vectors in reciprocal space for the k-mesh (in cartesian coords)
	//--- determine maximum dk among lattice directions:
	double dkMax = 0.;
	for(int iDir=0; iDir<3; iDir++)
		dkMax = std::max(dkMax, kbasis.column(iDir).length());
	//--- determine ik bounding box to contain dkMax sphere (with margins)
	vector3<int> ikBound;
	const matrix3<> kbasisInvT = ~inv(kbasis);
	for(int iDir=0; iDir<3; iDir++)
		ikBound[iDir] = std::max(2, 1+int(ceil(dkMax * kbasisInvT.column(iDir).length())));
	//--- collect non-zero dk's within ik boudning box in ascending order of magnitude:
	std::multimap<double, vector3<> > dkMap; //cartesian offsets from one k-point to other k-points sorted by distance
	vector3<int> ik;
	for(ik[0]=-ikBound[0]; ik[0]<=+ikBound[0]; ik[0]++)
	for(ik[1]=-ikBound[1]; ik[1]<=+ikBound[1]; ik[1]++)
	for(ik[2]=-ikBound[2]; ik[2]<=+ikBound[2]; ik[2]++)
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
	
	//Find edges :
	PeriodicLookup<WannierMinimizer::KmeshEntry> plook(kMesh, e.gInfo.GGT); //look-up table for O(1) fuzzy searching
	edges.resize(kMesh.size(), std::vector<Edge>(wb.size()));
	std::vector<std::set<int>> ranksNeeded(kMesh.size()); //rank of other processes that need each rotation
	for(size_t i=0; i<kMesh.size(); i++)
	{	int iProc = whose(i);
		for(unsigned j=0; j<wb.size(); j++)
		{	Edge& edge = edges[i][j];
			edge.wb = wb[j];
			edge.b = b[j];
			//Find neighbour:
			vector3<> kj = kMesh[i].point.k + inv(e.gInfo.GT) * b[j];
			edge.ik = plook.find(kj);
			edge.point = kMesh[edge.ik].point;
			edge.point.offset += round(kj - edge.point.k); //extra offset
			edge.point.k = kj;
			kpoints.insert(edge.point);
			//Update info about when MPI communication is needed:
			int jProc = whose(edge.ik);
			if(jProc != iProc)
				ranksNeeded[edge.ik].insert(iProc);
		}
	}
	
	//Create MPI communicators for rotations:
	for(size_t i=0; i<kMesh.size(); i++)
		if(ranksNeeded[i].size() //some communication is needed
			&& ( isMine(i) //this  process is the source
				|| ranksNeeded[i].count(mpiWorld->iProcess()) ) ) //this process is a consumer
		{
			std::vector<int> ranks(1, whose(i)); //source rank
			ranks.insert(ranks.end(), ranksNeeded[i].begin(), ranksNeeded[i].end());
			kMesh[i].mpi = std::make_shared<MPIUtil>(mpiWorld, ranks);
		}
	
	//Initialize preconditioner:
	if(wannier.precond)
	{	logPrintf("Initializing preconditioner ... "); logFlush();
		double kappa = M_PI * pow(e.gInfo.detR, 1./3); //inverse screening length (in k-space) set by cell size
		matrix helmholtz = zeroes(kMesh.size(), kMesh.size());
		complex* helmholtzData = helmholtz.data();
		for(size_t ik=ikStart; ik<ikStop; ik++)
		{	double wSum = 0.;
			const double wk = kMesh[ik].point.weight;
			for(const Edge& edge: edges[ik])
			{	wSum += 2*edge.wb;
				helmholtzData[helmholtz.index(ik,edge.ik)] -= wk * edge.wb;
				helmholtzData[helmholtz.index(edge.ik,ik)] -= wk * edge.wb;
			}
			helmholtzData[helmholtz.index(ik,ik)] += wk * (wSum + kappa*kappa);
		}
		mpiWorld->allReduceData(helmholtz, MPIUtil::ReduceSum);
		kHelmholtzInv = dagger_symmetrize(inv(helmholtz))(ikStart,ikStop, 0,kMesh.size()); //invert and split over MPI
		logPrintf("done.\n"); logFlush();
	}
}

void WannierMinimizerFD::initialize(int iSpin)
{
	//Read overlap matrices, if available:
	string fname = wannier.getFilename(Wannier::FilenameDump, "mlwfM0", &iSpin);
	bool M0exists = (fileSize(fname.c_str()) > 0);
	mpiWorld->bcast(M0exists); //Ensure MPI consistency of file check (avoid occassional NFS errors)
	if(wannier.loadRotations && M0exists)
	{	logPrintf("Reading initial overlaps from '%s' ... ", fname.c_str()); logFlush();
		size_t sizePerK = edges[0].size() * nBands*nBands * sizeof(complex);
		MPIUtil::File fp;
		mpiWorld->fopenRead(fp, fname.c_str(), kMesh.size()*sizePerK);
		mpiWorld->fseek(fp, ikStart*sizePerK, SEEK_SET);
		for(size_t ik=ikStart; ik<ikStop; ik++)
			for(Edge& edge: edges[ik])
			{	edge.M0 = zeroes(nBands, nBands);
				mpiWorld->freadData(edge.M0, fp);
			}
		mpiWorld->fclose(fp);
		logPrintf("done.\n"); logFlush();
		return; //read overlaps successfully rom file, so no need to recalculate below
	}
	
	//Compute the overlap matrices for current spin:
	for(int jProcess=0; jProcess<mpiWorld->nProcesses(); jProcess++)
	{	//Send/recv wavefunctions to other processes:
		Cother.assign(e.eInfo.nStates, ColumnBundle());
		VdagCother.clear(); VdagCother.resize(e.eInfo.nStates);
		if(jProcess == mpiWorld->iProcess()) //send
		{	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
			{	mpiWorld->bcastData((ColumnBundle&)e.eVars.C[q], jProcess);
				for(size_t iSp=0; iSp<e.iInfo.species.size(); iSp++)
					if(e.iInfo.species[iSp]->isUltrasoft())
						mpiWorld->bcastData((matrix&)e.eVars.VdagC[q][iSp], jProcess);
			}
		}
		else //recv
		{	for(int q=e.eInfo.qStartOther(jProcess); q<e.eInfo.qStopOther(jProcess); q++)
			{	Cother[q].init(nBands, e.basis[q].nbasis*nSpinor, &e.basis[q], &e.eInfo.qnums[q]);
				mpiWorld->bcastData(Cother[q], jProcess);
				VdagCother[q].resize(e.iInfo.species.size());
				for(size_t iSp=0; iSp<e.iInfo.species.size(); iSp++)
				{	const SpeciesInfo& sp = *(e.iInfo.species[iSp]);
					if(sp.isUltrasoft())
					{	VdagCother[q][iSp].init(sp.nProjectors(), nBands);
						mpiWorld->bcastData(VdagCother[q][iSp], jProcess);
					}
				}
			}
		}
		
		for(size_t ik=0; ik<kMesh.size(); ik++) if(isMine_q(ik,iSpin))
		{	KmeshEntry& ke = kMesh[ik];
			std::vector<matrix> VdagCi, VdagCj;
			ColumnBundle Ci = getWfns(ke.point, iSpin, &VdagCi); //Bloch functions at ik
			//Overlap with neighbours:
			for(Edge& edge: edges[ik])
				if(whose_q(edge.ik,iSpin)==jProcess)
					edge.M0 = overlap(Ci, getWfns(edge.point, iSpin, &VdagCj), &VdagCi, &VdagCj);
		}
	}
	Cother.clear();
	
	//Broadcast and dump the overlap matrices:
	FILE* fp = 0;
	if(mpiWorld->isHead())
	{	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
		fp = fopen(fname.c_str(), "w");
		if(!fp) die_alone("failed to open file for writing.\n");
	}
	for(size_t ik=0; ik<edges.size(); ik++)
		for(Edge& edge: edges[ik])
		{	if(!isMine_q(ik,iSpin)) edge.M0 = zeroes(nBands, nBands);
			mpiWorld->bcastData(edge.M0, whose_q(ik,iSpin));
			if(mpiWorld->isHead()) edge.M0.write(fp);
			if(!isMine(ik)) edge.M0 = matrix(); //not needed any more on this process
		}
	if(mpiWorld->isHead())
	{	fclose(fp);
		logPrintf("done.\n"); logFlush();
	}
}


double WannierMinimizerFD::getOmega(bool grad)
{	static StopWatch watch("WannierMinimizerFD::getOmega"); watch.start();
	
	//Compute the expectation values of r and rSq for each center (split over processes)
	rSqExpect.assign(nCenters, 0.);
	rExpect.assign(nCenters, vector3<>());
	std::vector<matrix> Mcache((ikStop-ikStart)*edges[0].size());
	matrix* McachePtr = Mcache.data();
	for(size_t ik=ikStart; ik<ikStop; ik++)
	{	const KmeshEntry& ki = kMesh[ik];
		for(const Edge& edge: edges[ik])
		{	const KmeshEntry& kj = kMesh[edge.ik];
			matrix& M = *(McachePtr++);
			M = dagger(ki.U) * edge.M0 * kj.U;
			const complex* Mdata = M.data();
			for(int n=0; n<nCenters; n++)
			{	complex Tnn = cis(dot(rPinned[n], edge.b)); //translation phase to rPinned as origin
				complex Mnn = Mdata[M.index(n,n)] * Tnn;
				double argMnn = Mnn.arg();
				rExpect[n] -= (ki.point.weight * edge.wb * argMnn) * edge.b;
				rSqExpect[n] += ki.point.weight * edge.wb * (argMnn*argMnn + 1. - Mnn.norm());
			}
		}
	}
	mpiWorld->allReduceData(rSqExpect, MPIUtil::ReduceSum);
	mpiWorld->allReduceData(rExpect, MPIUtil::ReduceSum);

	//Compute the total variance of the Wannier centers
	double Omega = 0.;
	std::vector< vector3<> > mHalf_Omega_rExpect(nCenters); //-1/2 dOmega/d(rExpect[n])
	for(int n=0; n<nCenters; n++)
	{	const vector3<> r0_n = pinned[n] ? vector3<>() : rExpect[n];
		Omega += (rSqExpect[n] - 2*dot(rExpect[n],r0_n) + r0_n.length_squared());
		mHalf_Omega_rExpect[n] = r0_n;
		//Shift back from rPinned to original origin for stored outputs:
		rSqExpect[n] -= rExpect[n].length_squared();
		rExpect[n] += rPinned[n]; 
		rSqExpect[n] += rExpect[n].length_squared();
	}
	
	//Compute gradients if required:
	if(grad)
	{	const matrix* McachePtr = Mcache.data();
		for(size_t ik=ikStart; ik<ikStop; ik++)
		{	KmeshEntry& ki = kMesh[ik];
			for(Edge& edge: edges[ik])
			{	KmeshEntry& kj = kMesh[edge.ik];
				const matrix& M = *(McachePtr++);
				//Compute dOmega/dM:
				std::vector<complex> Omega_M(nCenters);
				const complex* Mdata = M.data();
				for(int n=0; n<nCenters; n++)
				{	complex Tnn = cis(dot(rPinned[n], edge.b));
					complex Mnn = Mdata[M.index(n,n)] * Tnn;
					double argMnn = Mnn.arg();
					Omega_M[n] = 2. * ki.point.weight * edge.wb * Tnn
						* ((argMnn + dot(mHalf_Omega_rExpect[n],edge.b))*complex(0,-1)/Mnn - Mnn.conj());
				}
				//Propagate Omega_M to Omega_U:
				ki.Omega_UdotU += dagger(M(0,M.nRows(), 0,nCenters) * Omega_M);
				kj.Omega_UdotU += Omega_M * M(0,nCenters, 0,M.nCols());
			}
		}
	}
	watch.stop();
	return Omega;
}

double WannierMinimizerFD::getOmegaI(bool grad)
{	static StopWatch watch("WannierMinimizerFD::getOmegaI"); watch.start();
	double OmegaI = 0.;
	for(size_t ik=ikStart; ik<ikStop; ik++)
	{	KmeshEntry& ki = kMesh[ik];
		for(const Edge& edge: edges[ik])
		{	KmeshEntry& kj = kMesh[edge.ik];
			const matrix M = dagger(ki.U) * edge.M0 * kj.U;
			const auto Msub = M(0,nCenters, 0,nCenters);
			OmegaI += ki.point.weight * edge.wb * (nCenters - trace(Msub * dagger(Msub)).real());
			if(grad)
			{	matrix Omega_M = 2. * ki.point.weight * edge.wb * (-dagger(Msub));
				ki.Omega_UdotU += dagger(M(0,M.nRows(), 0,nCenters) * Omega_M);
				kj.Omega_UdotU += Omega_M * M(0,nCenters, 0,M.nCols());
			}
		}
	}
	mpiWorld->allReduce(OmegaI, MPIUtil::ReduceSum);
	watch.stop();
	return OmegaI;
}

WannierGradient WannierMinimizerFD::precondition(const WannierGradient& grad)
{	static StopWatch watch("WannierMinimizerFD::precond"); watch.start();
	assert(grad.wmin == this);
	WannierGradient Kgrad = grad;
	constrain(Kgrad);
	if(!kHelmholtzInv) //helmholtz preconditioning is disabled
		return Kgrad;
	//Determine range of free bands:
	int nInMax = 0, nFixedMin = nCenters;
	for(size_t ik=ikStart; ik<ikStop; ik++)
	{	nInMax = std::max(nInMax, kMesh[ik].nIn);
		nFixedMin = std::min(nFixedMin, kMesh[ik].nFixed);
	}
	mpiWorld->allReduce(nInMax, MPIUtil::ReduceMax);
	mpiWorld->allReduce(nFixedMin, MPIUtil::ReduceMin);
	int nEntries1 = (nCenters-nFixedMin) * (nInMax-nCenters); //max for B1
	int nEntries2 = (nCenters-nFrozen) * (nCenters-nFrozen); //constant for B2
	int nEntries = nEntries1 + nEntries2;
	//Copy each matrix of gradient into a column of a giant matrix:
	matrix gradMat = zeroes(nEntries, ikStop-ikStart);
    for(size_t ik=ikStart; ik<ikStop; ik++)
	{	const KmeshEntry& ki = kMesh[ik];
		complex* colPtr = gradMat.dataPref() + gradMat.index(0,ik-ikStart);
		if(nCenters > ki.nFixed)
		{	matrix B1 = zeroes(nCenters-nFixedMin, nInMax-nCenters);
			B1.set(ki.nFixed-nFixedMin,nCenters-nFixedMin, 0,ki.nIn-nCenters, Kgrad.B1[ik]);
			callPref(eblas_copy)(colPtr, B1.dataPref(), B1.nData());
		}
		callPref(eblas_copy)(colPtr+nEntries1, Kgrad.B2[ik].dataPref(), nEntries2);
	}
	//Apply preconditioner:
	matrix KgradMat = gradMat * kHelmholtzInv;
	mpiWorld->allReduceData(KgradMat, MPIUtil::ReduceSum); //Note: KgradMat has all k, while gradMat had subset
	//Copy result from each column to a small matrix per k-point:
    for(size_t ik=ikStart; ik<ikStop; ik++)
	{	const KmeshEntry& ki = kMesh[ik];
		const complex* colPtr = KgradMat.dataPref() + KgradMat.index(0,ik);
		if(nCenters > ki.nFixed)
		{	matrix B1 = zeroes(nCenters-nFixedMin, nInMax-nCenters);
			callPref(eblas_copy)(B1.dataPref(), colPtr, B1.nData());
			Kgrad.B1[ik] = B1(ki.nFixed-nFixedMin,nCenters-nFixedMin, 0,ki.nIn-nCenters);
		}
		callPref(eblas_copy)(Kgrad.B2[ik].dataPref(), colPtr+nEntries1, nEntries2);
	}
	constrain(Kgrad);
	watch.stop();
	return Kgrad;
}
