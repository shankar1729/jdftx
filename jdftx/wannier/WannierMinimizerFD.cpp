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
	for(size_t i=0; i<kMesh.size(); i++)
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
		helmholtz.allReduce(MPIUtil::ReduceSum);
		kHelmholtzInv = dagger_symmetrize(inv(helmholtz))(ikStart,ikStop, 0,kMesh.size()); //invert and split over MPI
		logPrintf("done.\n"); logFlush();
	}
}

void WannierMinimizerFD::initialize(int iSpin)
{
	//Compute the overlap matrices for current spin:
	for(int jProcess=0; jProcess<mpiUtil->nProcesses(); jProcess++)
	{	//Send/recv wavefunctions to other processes:
		Cother.assign(e.eInfo.nStates, ColumnBundle());
		if(jProcess == mpiUtil->iProcess()) //send
		{	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
				((ColumnBundle&)e.eVars.C[q]).bcast(jProcess);
		}
		else //recv
		{	for(int q=e.eInfo.qStartOther(jProcess); q<e.eInfo.qStopOther(jProcess); q++)
			{	Cother[q].init(nBands, e.basis[q].nbasis*nSpinor, &e.basis[q], &e.eInfo.qnums[q]);
				Cother[q].bcast(jProcess);
			}
		}
		
		for(size_t ik=0; ik<kMesh.size(); ik++) if(isMine_q(ik,iSpin))
		{	KmeshEntry& ke = kMesh[ik];
			ColumnBundle Ci = getWfns(ke.point, iSpin); //Bloch functions at ik
			//Overlap with neighbours:
			for(Edge& edge: edges[ik])
				if(whose_q(edge.ik,iSpin)==jProcess)
					edge.M0 = overlap(Ci, getWfns(edge.point, iSpin));
		}
	}
	Cother.clear();
	
	//Broadcast the overlap matrices:
	for(size_t ik=0; ik<edges.size(); ik++)
		for(Edge& edge: edges[ik])
		{	if(!isMine_q(ik,iSpin)) edge.M0 = zeroes(nBands, nBands);
			edge.M0.bcast(whose_q(ik,iSpin));
			if(!isMine(ik)) edge.M0 = matrix(); //not needed any more on this process
		}
}


double WannierMinimizerFD::getOmega(bool grad)
{
	//Compute the expectation values of r and rSq for each center (split over processes)
	rSqExpect.assign(nCenters, 0.);
	rExpect.assign(nCenters, vector3<>(0,0,0));
	for(size_t ik=ikStart; ik<ikStop; ik++)
	{	const KmeshEntry& ki = kMesh[ik];
		for(const Edge& edge: edges[ik])
		{	const KmeshEntry& kj = kMesh[edge.ik];
			const matrix M = dagger(ki.U) * edge.M0 * kj.U;
			const complex* Mdata = M.data();
			for(int n=0; n<nCenters; n++)
			{	complex Mnn = Mdata[M.index(n,n)];
				double argMnn = Mnn.arg();
				rExpect[n] -= (ki.point.weight * edge.wb * argMnn) * edge.b;
				rSqExpect[n] += ki.point.weight * edge.wb * (argMnn*argMnn + 1. - Mnn.norm());
			}
		}
	}
	mpiUtil->allReduce(rSqExpect.data(), nCenters, MPIUtil::ReduceSum);
	mpiUtil->allReduce((double*)rExpect.data(), 3*nCenters, MPIUtil::ReduceSum);
	
	//Compute the total variance of the Wannier centers
	double Omega = 0.;
	std::vector< vector3<> > mHalf_Omega_rExpect(nCenters); //-1/2 dOmega/d(rExpect[n])
	for(int n=0; n<nCenters; n++)
	{	const vector3<> r0_n = pinned[n] ? rPinned[n] : rExpect[n];
		Omega += (rSqExpect[n] - 2*dot(rExpect[n],r0_n) + r0_n.length_squared());
		mHalf_Omega_rExpect[n] = r0_n;
	}
	
	//Compute gradients if required:
	if(grad)
	{	for(size_t ik=ikStart; ik<ikStop; ik++)
		{	KmeshEntry& ki = kMesh[ik];
			for(Edge& edge: edges[ik])
			{	KmeshEntry& kj = kMesh[edge.ik];
				const matrix M = dagger(ki.U) * edge.M0 * kj.U;
				//Compute dOmega/dM:
				matrix Omega_M = zeroes(nCenters, nCenters);
				const complex* Mdata = M.data();
				complex* Omega_Mdata = Omega_M.data();
				for(int n=0; n<nCenters; n++)
				{	complex Mnn = Mdata[M.index(n,n)];
					double argMnn = atan2(Mnn.imag(), Mnn.real());
					Omega_Mdata[Omega_M.index(n,n)] =
						2. * ki.point.weight * edge.wb
						* ((argMnn + dot(mHalf_Omega_rExpect[n],edge.b))*complex(0,-1)/Mnn - Mnn.conj());
				}
				//Propagate Omega_M to Omega_U:
				ki.Omega_U += dagger(edge.M0 * kj.U * Omega_M);
				kj.Omega_U += Omega_M * dagger(ki.U) * edge.M0;
			}
		}
	}
	return Omega;
}

double WannierMinimizerFD::getOmegaI(bool grad)
{	double OmegaI = 0.;
	for(size_t ik=ikStart; ik<ikStop; ik++)
	{	KmeshEntry& ki = kMesh[ik];
		for(const Edge& edge: edges[ik])
		{	KmeshEntry& kj = kMesh[edge.ik];
			const matrix M = dagger(ki.U) * edge.M0 * kj.U;
			OmegaI += ki.point.weight * edge.wb * (nCenters - trace(M * dagger(M)).real());
			if(grad)
			{	matrix Omega_M = 2. * ki.point.weight * edge.wb * (-dagger(M));
				ki.Omega_U += dagger(edge.M0 * kj.U * Omega_M);
				kj.Omega_U += Omega_M * dagger(ki.U) * edge.M0;
			}
		}
	}
	mpiUtil->allReduce(OmegaI, MPIUtil::ReduceSum);
	return OmegaI;
}

WannierGradient WannierMinimizerFD::precondition(const WannierGradient& grad)
{	static StopWatch watch("WannierMinimizerFD::precondition"); watch.start();
	assert(grad.size()==kMesh.size());
	WannierGradient Kgrad = grad;
	constrain(Kgrad);
	if(!kHelmholtzInv) //helmholtz preconditioning is disabled
		return Kgrad;
	//Figure out max input bands for any kpoint:
	int nInMax = 0;
	for(size_t ik=ikStart; ik<ikStop; ik++)
		nInMax = std::max(nInMax, kMesh[ik].nIn);
	mpiUtil->allReduce(nInMax, MPIUtil::ReduceMax);
	//Copy each matrix of gradient into a column of a giant matrix:
	matrix gradMat = zeroes(nCenters*nInMax, ikStop-ikStart);
	complex* gradMatData = gradMat.dataPref();
    for(size_t ik=ikStart; ik<ikStop; ik++)
		callPref(eblas_copy)(gradMatData+gradMat.index(0,ik-ikStart), Kgrad[ik].dataPref(), Kgrad[ik].nData());
	//Apply preconditioner:
	matrix KgradMat = gradMat * kHelmholtzInv;
	KgradMat.allReduce(MPIUtil::ReduceSum);
	//Copy result from each column to a small matrix per k-point:
	const complex* KgradMatData = KgradMat.dataPref();
    for(size_t ik=ikStart; ik<ikStop; ik++)
	{	Kgrad[ik].init(nCenters, kMesh[ik].nIn, isGpuEnabled());
		callPref(eblas_copy)(Kgrad[ik].dataPref(), KgradMatData+KgradMat.index(0,ik), Kgrad[ik].nData());
	}
	constrain(Kgrad);
	watch.stop();
	return Kgrad;
}
