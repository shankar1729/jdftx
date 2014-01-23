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
#include <core/LatticeUtils.h>
#include <gsl/gsl_linalg.h>


//Find a finite difference formula given a list of relative neighbour positions (in cartesian coords)
//Generalization of Appendix B from Phys Rev B 56, 12847 to arbitrary k-point meshes
//Returns an empty weight set on failure
std::vector<double> getFDformula(const std::vector< vector3<> >& b)
{	//Group elements of b into shells:
	std::vector<unsigned> shellMax; //cumulative count within each shell
	for(unsigned i=1; i<b.size(); i++)
		if(b[i].length() > b[i-1].length() + symmThreshold)
			shellMax.push_back(i);
	shellMax.push_back(b.size());
	//Setup the equations satisfied by the weights:
	int nEquations = std::max(19, int(shellMax.size())); //pad extra equations to keep nRows>=nCols for SVD
	gsl_matrix* Lhs = gsl_matrix_calloc(nEquations, shellMax.size()); //initializes with 0s
	for(unsigned s=0; s<shellMax.size(); s++)
		for(unsigned j = (s ? shellMax[s-1] : 0); j<shellMax[s]; j++)
		{	//Equations from ref.:
			//Rank-two sum is identity:
			*gsl_matrix_ptr(Lhs, 0, s) += b[j][0]*b[j][0];
			*gsl_matrix_ptr(Lhs, 1, s) += b[j][1]*b[j][1];
			*gsl_matrix_ptr(Lhs, 2, s) += b[j][2]*b[j][2];
			*gsl_matrix_ptr(Lhs, 3, s) += b[j][1]*b[j][2];
			*gsl_matrix_ptr(Lhs, 4, s) += b[j][2]*b[j][0];
			*gsl_matrix_ptr(Lhs, 5, s) += b[j][0]*b[j][1];
			//Additional constraints for arbitrary meshes (always satisfied for Bravais lattices):
			//Rank-one sum vanishes:
			*gsl_matrix_ptr(Lhs, 6, s) += b[j][0];
			*gsl_matrix_ptr(Lhs, 7, s) += b[j][1];
			*gsl_matrix_ptr(Lhs, 8, s) += b[j][2];
			//Rank-three sum vanishes:
			*gsl_matrix_ptr(Lhs,  9, s) += b[j][0]*b[j][0]*b[j][0];
			*gsl_matrix_ptr(Lhs, 10, s) += b[j][1]*b[j][1]*b[j][1];
			*gsl_matrix_ptr(Lhs, 11, s) += b[j][2]*b[j][2]*b[j][2];
			*gsl_matrix_ptr(Lhs, 12, s) += b[j][1]*b[j][1]*b[j][2];
			*gsl_matrix_ptr(Lhs, 13, s) += b[j][2]*b[j][2]*b[j][0];
			*gsl_matrix_ptr(Lhs, 14, s) += b[j][0]*b[j][0]*b[j][1];
			*gsl_matrix_ptr(Lhs, 15, s) += b[j][1]*b[j][2]*b[j][2];
			*gsl_matrix_ptr(Lhs, 16, s) += b[j][2]*b[j][0]*b[j][0];
			*gsl_matrix_ptr(Lhs, 17, s) += b[j][0]*b[j][1]*b[j][1];
			*gsl_matrix_ptr(Lhs, 18, s) += b[j][0]*b[j][1]*b[j][2];
		}
	gsl_vector* rhs = gsl_vector_calloc(nEquations); //initializes with 0s
	for(unsigned i=0; i<3; i++) //first three components = diagonals of rank-two sum
		gsl_vector_set(rhs, i, 1.);
	//Solve using a singular value decomposition:
	gsl_matrix* U = gsl_matrix_alloc(nEquations, shellMax.size());
	gsl_matrix* V = gsl_matrix_alloc(shellMax.size(), shellMax.size());
	gsl_vector* S = gsl_vector_alloc(shellMax.size());
	gsl_vector* work = gsl_vector_alloc(shellMax.size());
	gsl_matrix_memcpy(U, Lhs); //SVD is done in place
	gsl_linalg_SV_decomp(U, V, S, work);
	//Zero out small singular values:
	for(unsigned j=0; j<shellMax.size(); j++)
		if(gsl_vector_get(S,j) < symmThreshold)
			gsl_vector_set(S,j, 0.);
	//Solve for weights:
	gsl_vector* wPairs = gsl_vector_alloc(shellMax.size());
	gsl_linalg_SV_solve(U, V, S, rhs, wPairs);
	gsl_matrix_free(U);
	gsl_matrix_free(V);
	gsl_vector_free(S);
	gsl_vector_free(work);
	//Check solution by substitution:
	cblas_dgemv(CblasRowMajor, CblasNoTrans, nEquations, shellMax.size(), 1., Lhs->data, Lhs->tda,
		wPairs->data, wPairs->stride, -1., rhs->data, rhs->stride); //rhs = Lhs*wPairs - rhs
	if(eblas_dnrm2(nEquations, rhs->data, rhs->stride) > symmThreshold)
		return std::vector<double>(); //Not an exact solution, so quit (and try with more shells)
	gsl_vector_free(rhs);
	gsl_matrix_free(Lhs);
	//Store the weights in the original indexing:
	std::vector<double> w(b.size());
	for(unsigned s=0; s<shellMax.size(); s++)
		for(unsigned j = (s ? shellMax[s-1] : 0); j<shellMax[s]; j++)
			w[j] = gsl_vector_get(wPairs, s);
	return w;
}

//Helper function for PeriodicLookup<WannierMinimizer::Kpoint> used in WannierMinimizer::WannierMinimizer
inline vector3<> getCoord(const WannierMinimizer::Kpoint& kpoint) { return kpoint.k; }

WannierMinimizer::WannierMinimizer(const Everything& e, const Wannier& wannier) : e(e), wannier(wannier), sym(e.symm.getMatrices()),
	nCenters(wannier.centers.size()), nBands(e.eInfo.nBands),
	nSpins(e.eInfo.spinType==SpinNone ? 1 : 2), qCount(e.eInfo.qnums.size()/nSpins)
{
	logPrintf("\n---------- Initializing Wannier Function solver ----------\n");
	logPrintf("Setting up finite difference formula on k-mesh ...\n"); logFlush();
	const std::vector<QuantumNumber>& qnums = e.eInfo.qnums;
	//Create the list of images (closed under the symmetry group):
	const std::vector<int>& invertList = e.symm.getKpointInvertList();
	std::vector<Kpoint> kpoints;
	PeriodicLookup<WannierMinimizer::Kpoint> plook(kpoints, e.gInfo.GGT, qnums.size()*sym.size()); //look-up table for O(1) fuzzy searching
	for(int invert: invertList)
		for(int q=0; q<qCount; q++)
			for(unsigned iRot=0; iRot<sym.size(); iRot++)
			{	vector3<> k = (~sym[iRot]) * qnums[q].k * invert;
				//Find offset that brings it into centered zone
				vector3<int> offset;
				for(int i=0; i<3; i++)
				{	offset[i] = -floor(k[i]+0.5);
					k[i] += offset[i];
				}
				//Add to map if this k-vector has not yet been encountered:
				if(plook.find(k) == string::npos)
				{	Kpoint kpoint = { k, q, iRot, invert, offset };
					plook.addPoint(kpoints.size(), kpoint);
					kpoints.push_back(kpoint);
					addIndex(kpoint);
				}
			}
	
	//Determine distribution amongst processes:
	ikStart = (kpoints.size() * mpiUtil->iProcess()) / mpiUtil->nProcesses();
	ikStop = (kpoints.size() * (mpiUtil->iProcess()+1)) / mpiUtil->nProcesses();
	
	kMesh.resize(kpoints.size());
	for(size_t i=0; i<kMesh.size(); i++) //FD formula needed on all nodes for initial matrix generation
	{	//Create a list of neighbours for FD formula:
		//Collect from 3x3x3 lowest Brillouin zones (for worst case Gamma-only scenario)
		//This could be optimized, but this is unlikely to ever be too expensive
		struct Neighbour
		{	vector3<> dk; //difference from current k-point
			unsigned ik; //index of source k-point
			vector3<int> iG; //brillouin zone index (additional offset)
		};
		std::multimap<double,Neighbour> neighbourMap; //neighbours sorted by distance
		
		vector3<int> iG;
		for(iG[0]=-1; iG[0]<=+1; iG[0]++)
			for(iG[1]=-1; iG[1]<=+1; iG[1]++)
				for(iG[2]=-1; iG[2]<=+1; iG[2]++)
					for(unsigned j=0; j<kpoints.size(); j++)
					{	Neighbour neighbour = { iG + kpoints[j].k - kpoints[i].k, j, iG };
						double dist = sqrt(e.gInfo.GGT.metric_length_squared(neighbour.dk));
						if(dist > symmThreshold) //ignore self
							neighbourMap.insert(std::make_pair(dist, neighbour));
					}
		
		std::vector<Neighbour> neighbours; //list of neighbours chosen so far
		std::vector< vector3<> > b; //list of cartesian offsets to corresponding neighbours
		std::vector<double> wb; //corresponding weights in finite difference formula
		
		for(auto iter=neighbourMap.begin(); iter!=neighbourMap.end(); )
		{	//Add all the neighbours with equivalent distances:
			while(true)
			{	neighbours.push_back(iter->second);
				b.push_back(e.gInfo.GT * iter->second.dk);
				double prevDist = iter->first;
				iter++;
				if(iter==neighbourMap.end() || //end of neighbour list (should not be encountered)
					iter->first > prevDist+symmThreshold) //next neighbour is further away beyond tolerance
					break;
			}
			//Check if this list of neighbours is sufficient to get a finite difference formula
			wb = getFDformula(b);
			if(wb.size() == b.size()) break; //success
		}
		if(!wb.size())
			die("Failed to find a second order finite difference formula around k-point [ %lg %lg %lg ].\n",
				kpoints[i].k[0], kpoints[i].k[1], kpoints[i].k[2]);
		
		//Store the k-point with its FD formula in kMesh
		KmeshEntry& kMeshEntry = kMesh[i];
		kMeshEntry.point = kpoints[i];
		for(unsigned j=0; j<wb.size(); j++)
		{	EdgeFD edge;
			edge.wb = wb[j];
			edge.b = b[j];
			edge.ik = neighbours[j].ik;
			edge.point = kpoints[neighbours[j].ik];
			edge.point.offset += neighbours[j].iG;
			edge.point.k += vector3<>(neighbours[j].iG);
			addIndex(edge.point);
			kMeshEntry.edge.push_back(edge);
		}
		//Find the Brillouin zone integration weight for the k-point:
		int nStabilizer = 0; //size of the stabilizer subgroup
		for(int invert: invertList)
			for(const matrix3<int>& m: sym)
				if(circDistanceSquared(kpoints[i].k, (~m)*qnums[kpoints[i].q].k*invert) < symmThresholdSq)
						nStabilizer++;
		kMeshEntry.wk = nStabilizer * qnums[kpoints[i].q].weight * (0.5*nSpins) / (sym.size() * invertList.size());
	}
	
	//Create the common reduced basis set (union of all the reduced bases)
	//Collect all referenced full-G indices
	std::set<int> commonSet;
	for(auto index: indexMap)
		for(int j=0; j<index.second->nIndices; j++)
			commonSet.insert(index.second->data[j]);
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
	//Update indexMap to point to common reduced basis instead of full G:
	for(auto mapEntry: indexMap)
	{	Index& index = *mapEntry.second;
		for(int j=0; j<index.nIndices; j++)
			index.data[j] = commonInverseMap[index.data[j]];
		index.set();
	}
	
	//Initialize the preconditioner (inverse helhmoltz with kappa chosen to regularize 0 frequency):
	logPrintf("Setting up inverse helmholtz preconditioner on k-mesh ...\n"); logFlush();
	matrix helmholtz(kMesh.size(), kMesh.size());
	double kappa = M_PI * pow(e.gInfo.detR, 1./3); //inverse screening length (in k-space) set by cell size
	helmholtz.zero();
	complex* helmholtzData = helmholtz.data();
	for(size_t ik=ikStart; ik<ikStop; ik++)
	{	double wSum = 0.;
		for(EdgeFD& edge: kMesh[ik].edge)
		{	wSum += edge.wb;
			helmholtzData[helmholtz.index(ik,edge.ik)] -= kMesh[ik].wk * edge.wb;
		}
		helmholtzData[helmholtz.index(ik,ik)] += kMesh[ik].wk * (wSum + kappa*kappa);
	}
	helmholtz.allReduce(MPIUtil::ReduceSum);
	if(nrm2(helmholtz-dagger(helmholtz)) > symmThresholdSq * nrm2(helmholtz))
	{	logPrintf("Laplacian operator on k-mesh not symmetric - using identity preconditioner.\n");
		kHelmholtzInv = eye(kMesh.size());
	}
	else kHelmholtzInv = dagger_symmetrize(inv(helmholtz));
}
