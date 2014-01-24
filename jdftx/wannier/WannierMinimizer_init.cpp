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

	//Create supercell grid:
	logPrintf("\n---------- Initializing supercell grid for Wannier functions ----------\n");
	const Supercell& supercell = *(e.coulombParams.supercell);
	gInfoSuper.R = supercell.Rsuper;
	gInfoSuper.Gmax = e.gInfo.Gmax;
	gInfoSuper.GmaxRho = e.gInfo.GmaxRho;
	gInfoSuper.initialize(true, sym);

	//Determine finite difference formula:
	logPrintf("Setting up finite difference formula on k-mesh ... "); logFlush();
	matrix3<> kbasis = e.gInfo.GT * inv(~matrix3<>(supercell.super)); //basis vectors in reciprocal space for the k-mesh (in cartesian coords)
	std::multimap<double, vector3<> > dkMap; //cartesian offsets from one k-point to other k-points sorted by distance
	vector3<int> ik;
	for(ik[0]=-2; ik[0]<=+2; ik[0]++)
	for(ik[1]=-2; ik[1]<=+2; ik[1]++)
	for(ik[2]=-2; ik[2]<=+2; ik[2]++)
		if(ik.length_squared()) //ignore self
		{	vector3<> dk = kbasis * ik;
			dkMap.insert(std::make_pair<>(dk.length_squared(), dk));
		}
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
	logPrintf("found a %lu neighbour formula.\n", wb.size());
	
	//Create a list of kpoints:
	const std::vector<QuantumNumber>& qnums = e.eInfo.qnums;
	std::vector<Kpoint> kpoints(supercell.kmeshTransform.size());
	for(size_t i=0; i<kpoints.size(); i++)
	{	const Supercell::KmeshTransform& src = supercell.kmeshTransform[i];
		Kpoint& kpoint = kpoints[i];
		kpoint.q = src.iReduced;
		kpoint.iRot = src.iSym;
		kpoint.invert = src.invert;
		kpoint.offset = src.offset;
		kpoint.k = (~sym[src.iSym]) * (src.invert * qnums[src.iReduced].k) + src.offset;
	}
	wk = 1./kpoints.size();
	PeriodicLookup<WannierMinimizer::Kpoint> plook(kpoints, e.gInfo.GGT); //look-up table for O(1) fuzzy searching
	if(plook.find(vector3<>())==string::npos)
		die("k-mesh does not contain Gamma point. Wannier requires uniform Gamma-centered k-mesh.\n");
	
	//Determine distribution amongst processes:
	ikStart = (kpoints.size() * mpiUtil->iProcess()) / mpiUtil->nProcesses();
	ikStop = (kpoints.size() * (mpiUtil->iProcess()+1)) / mpiUtil->nProcesses();
	
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
	auto superSetIter = commonSuperSet.begin();
	for(unsigned j=0; j<indexSuperCommon.size(); j++)
	{	int i = *(superSetIter++);
		indexSuperCommon[j] = i;
		commonSuperInverseMap[i] = j;
	}
	basisSuper.setup(gInfoSuper, e.iInfo, indexSuperCommon);
	//Update indexMap to point to common reduced basis instead of full G:
	for(auto mapEntry: indexMap)
	{	Index& index = *mapEntry.second;
		for(int j=0; j<index.nIndices; j++)
		{	index.data[j] = commonInverseMap[index.data[j]];
			index.dataSuper[j] = commonSuperInverseMap[index.dataSuper[j]];
		}
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
			helmholtzData[helmholtz.index(ik,edge.ik)] -= wk * edge.wb;
		}
		helmholtzData[helmholtz.index(ik,ik)] += wk * (wSum + kappa*kappa);
	}
	helmholtz.allReduce(MPIUtil::ReduceSum);
	if(nrm2(helmholtz-dagger(helmholtz)) > symmThresholdSq * nrm2(helmholtz))
	{	logPrintf("Laplacian operator on k-mesh not symmetric - using identity preconditioner.\n");
		kHelmholtzInv = eye(kMesh.size());
	}
	else kHelmholtzInv = dagger_symmetrize(inv(helmholtz));
}
