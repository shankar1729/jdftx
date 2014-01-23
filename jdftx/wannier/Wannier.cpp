/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

#include <wannier/Wannier.h>
#include <electronic/operators.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/SpeciesInfo_internal.h>
#include <core/LatticeUtils.h>
#include <core/BlasExtra.h>
#include <core/Random.h>
#include <gsl/gsl_linalg.h>

typedef std::vector<matrix> WannierGradient;

static WannierGradient clone(const WannierGradient& grad) { return grad; }
static double dot(const WannierGradient& x, const WannierGradient& y)
{	assert(x.size()==y.size());
	double result = 0.;
	for(unsigned i=0; i<x.size(); i++)
		result += dotc(x[i], y[i]).real();
	return result;
}
static WannierGradient& operator*=(WannierGradient& x, double alpha)
{	for(unsigned i=0; i<x.size(); i++) x[i] *= alpha;
	return x;
}
static void axpy(double alpha, const WannierGradient& x, WannierGradient& y)
{	assert(x.size()==y.size());
	for(unsigned i=0; i<x.size(); i++) axpy(alpha, x[i], y[i]);
}

static matrix randomMatrix(int nRows, int nCols)
{	matrix ret(nRows, nCols, false);
	complex* retData = ret.data();
	for(unsigned j=0; j<ret.nData(); j++)
		retData[j] = Random::normalComplex();
	return ret;
}
static void randomize(WannierGradient& x)
{	for(unsigned i=0; i<x.size(); i++)
		x[i] = dagger_symmetrize(randomMatrix(x[i].nRows(), x[i].nCols()));
}

class WannierEval : public Minimizable<WannierGradient>
{
public:
	WannierEval(const Everything& e, const Wannier& wannier);
	virtual ~WannierEval() {}
	
	void saveMLWF(const std::vector<Wannier::Center>& centers); //save wannier functions for all spins
	void saveMLWF(const std::vector<Wannier::Center>& centers, int iSpin); //save for specified spin
	
	//Interface for minimize:
	void step(const WannierGradient& grad, double alpha);
	double compute(WannierGradient* grad);
	WannierGradient precondition(const WannierGradient& grad);
	double sync(double x) const { mpiUtil->bcast(x); return x; } //!< All processes minimize together; make sure scalars are in sync to round-off error

	//! Entries in the k-point mesh
	struct Kpoint
	{	vector3<> k; //!< k-point in reciprocal lattice coordinates
		int q; //!< source state
		unsigned iRot; //!< symmetry matrix index from state q to this kpoint
		int invert; //!< sign of transformation matrix (if inversion symmetry of k-mesh is being used)
		vector3<int> offset; //!< translation after rotation to bring to k
		
		bool operator<(const Kpoint& other) const;
		bool operator==(const Kpoint& other) const;
	};
	
private:
	const Everything& e;
	const Wannier& wannier;
	const std::vector< matrix3<int> >& sym;
	int nSpins, qCount; //!< number of spins, and number of states per spin
	std::vector<double> rSqExpect; //!< Expectation values for r^2 per center in current group
	std::vector< vector3<> > rExpect; //!< Expectation values for r per center in current group
	matrix kHelmholtzInv; //!< Inverse of Hemlholtz equation in k-space (preconditioner)
	
	//!An edge of the k-mesh involved in the finite difference formula
	struct EdgeFD
	{	double wb; //!< weight of neighbour
		vector3<> b; //!< displacement to neighbour
		unsigned ik; //!< index of neighbour in kMesh
		Kpoint point; //!< description of neighbour (source state, rotation, translation etc.)
		matrix M0; //!< initial overlap matrix for this pair (after applying the initial unitary rotations)
	};
	
	//!Entries in the k-point mesh with FD formula
	struct KmeshEntry
	{	Kpoint point; //!< point in the mesh
		double wk; //!< weight of k-point for Brillouin zone integration
		std::vector<EdgeFD> edge; //!< neighbours with weights defining finite difference formula
		//State of system for Wannier minimize:
		matrix B; //!< Independent variable for minimization
		matrix Bevecs; //!< Eigenvectors of B
		diagMatrix Beigs; //!< Eigenvalues of B
		matrix U0; //!< Initial unitary rotation (from trial wave functions)
		matrix V; //!< Subsequent unitary rotation of k-point given by e^iB (net unitary rotation = U0 * V)
	};
	
	//!Indices from reduced basis to full G-space or a union of all reduced bases
	struct Index
	{	const int nIndices;
		int* data;
		#ifdef GPU_ENABLED
		int* dataGpu;
		#endif
		int* dataPref;
		
		Index(int nIndices=0); //!< allocate space for indices
		~Index();
		void set(); //!< update gpu-copies
		
		//Non-copyable:
		Index(const Index&)=delete;
		Index& operator=(const Index&)=delete;
	};
	
	std::vector<KmeshEntry> kMesh; //!< k-point mesh with FD formula
	std::map<Kpoint, std::shared_ptr<Index> > indexMap; //!< wave-function indexing from each k-point to the common union basis
	Basis basis; //!< common basis (with indexing into full G-space)
	
	//k-mesh MPI division:
	size_t ikStart, ikStop;
	bool isMine(size_t ik) const { return ik>=ikStart && ik<ikStop; }

	//state MPI division (wrappers to ElecInfo)
	bool isMine_q(int ik, int iSpin) const { return e.eInfo.isMine(kMesh[ik].point.q + iSpin*qCount); }
	int whose_q(int ik, int iSpin) const { return e.eInfo.whose(kMesh[ik].point.q + iSpin*qCount); }
	
	void addIndex(const Kpoint& kpoint); //!< Add index for a given kpoint to indexMap, with indices pointing to full G-space
	
	//! Get the wavefunctions for a particular k-point for bands involved in current group
	//! The wavefunctions are returned in the common basis
	ColumnBundle getWfns(const Kpoint& kpoint, const std::vector<Wannier::Center>& centers, int iSpin) const;
	std::vector<ColumnBundle> Cother; //wavefunctions from another process
	
	//! Get the trial wavefunctions (gaussians) for the group of centers in the common basis
	ColumnBundle trialWfns(const Kpoint& kpoint, const std::vector<Wannier::Center>& centers) const;
};

//--------------------- class Wannier implementation ----------------------

Wannier::Wannier() : eval(0)
{
}

Wannier::~Wannier()
{	if(eval) delete eval;
}


void Wannier::setup(const Everything& everything)
{	e = &everything;
	//Initialize minimization parameters:
	minParams.fpLog = globalLog;
	minParams.linePrefix = "WannierMinimize: ";
	minParams.energyLabel = "rVariance";
	//Initialize evaluator:
	eval = new WannierEval(*e, *this);
	Citations::add("Maximally-localized Wannier functions",
		"N. Marzari and D. Vanderbilt, Phys. Rev. B 56, 12847 (1997)");
}

void Wannier::saveMLWF()
{	for(auto g: group)
		eval->saveMLWF(g);
}

//--------------------- class WannierEval implementation ---------------------

inline double toUnitInterval(double x) { return x-floor(x); }
inline double toCenteredUnitInterval(double x)
{	double xUnit = toUnitInterval(x);
	if(xUnit>0.5) xUnit-=1.;
	return xUnit;
}

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

//Helper function for PeriodicLookup<WannierEval::Kpoint> used in WannierEval::WannierEval
inline vector3<> getCoord(const WannierEval::Kpoint& kpoint) { return kpoint.k; }

WannierEval::WannierEval(const Everything& e, const Wannier& wannier) : e(e), wannier(wannier), sym(e.symm.getMatrices()),
	nSpins(e.eInfo.spinType==SpinNone ? 1 : 2), qCount(e.eInfo.qnums.size()/nSpins)
{
	logPrintf("\n---------- Initializing Wannier Function solver ----------\n");
	logPrintf("Setting up finite difference formula on k-mesh ...\n"); logFlush();
	const std::vector<QuantumNumber>& qnums = e.eInfo.qnums;
	//Create the list of images (closed under the symmetry group):
	const std::vector<int>& invertList = e.symm.getKpointInvertList();
	std::vector<Kpoint> kpoints;
	PeriodicLookup<WannierEval::Kpoint> plook(kpoints, e.gInfo.GGT, qnums.size()*sym.size()); //look-up table for O(1) fuzzy searching
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

void WannierEval::saveMLWF(const std::vector<Wannier::Center>& centers)
{	for(int iSpin=0; iSpin<nSpins; iSpin++)
		saveMLWF(centers, iSpin);
}

void WannierEval::saveMLWF(const std::vector<Wannier::Center>& centers, int iSpin)
{	logPrintf("Computing wannier functions for bands (");
	for(auto center: centers) logPrintf(" %d", center.band);
	if(nSpins==1) logPrintf(" )\n"); else logPrintf(" ) spin %s\n", iSpin==0 ? "up" : "dn");
	logFlush();
	
	//Check for initial state:
	string UvarName;
	{	ostringstream ossUvarName;
		int nMin = centers.front().band;
		for(const auto& c: centers) nMin = std::min(nMin, c.band);
		ossUvarName << nMin << ".Umlwf";
		if(nSpins==2) ossUvarName << (iSpin==0 ? "Up" : "Dn");
		UvarName = ossUvarName.str();
	}
	bool readInitialMats = false;
	if(e.eVars.wfnsFilename.length())
	{	//Replace wfns with above avriable name to locate file:
		string filename = e.eVars.wfnsFilename;
		filename.replace(filename.find("wfns"),4, UvarName);
		//Check if file exists:
		{	FILE* fp = fopen(filename.c_str(), "r");
			if(fp)
			{	logPrintf("Reading initial matrices from %s (ignoring trial projections).\n", filename.c_str());
				for(size_t ik=0; ik<kMesh.size(); ik++)
				{	kMesh[ik].U0.init(centers.size(), centers.size());
					kMesh[ik].U0.read(fp);
				}
				fclose(fp);
				readInitialMats = true;
			}
		}
	}
	
	//Compute the overlap matrices and initial rotations for current group of centers:
	for(int jProcess=0; jProcess<mpiUtil->nProcesses(); jProcess++)
	{	//Send/recv wavefunctions to other processes:
		Cother.assign(e.eInfo.nStates, ColumnBundle());
		if(jProcess == mpiUtil->iProcess()) //send
		{	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
				((ColumnBundle&)e.eVars.C[q]).bcast(jProcess);
		}
		else //recv
		{	for(int q=e.eInfo.qStartOther(jProcess); q<e.eInfo.qStopOther(jProcess); q++)
			{	Cother[q].init(e.eInfo.nBands, e.basis[q].nbasis, &e.basis[q], &e.eInfo.qnums[q]);
				Cother[q].bcast(jProcess);
			}
		}
		
		for(size_t ik=0; ik<kMesh.size(); ik++) if(isMine_q(ik,iSpin))
		{	ColumnBundle Ci = getWfns(kMesh[ik].point, centers, iSpin); //Bloch functions at ik
			ColumnBundle OCi = O(Ci);
			//Overlap with neighbours:
			for(EdgeFD& edge: kMesh[ik].edge)
				if(whose_q(edge.ik,iSpin)==jProcess)
				{	//Pick up result from reverse edge if it has already been computed:
					bool foundReverse = false;
					if(jProcess==mpiUtil->iProcess() && edge.ik<ik)
					{	auto neighbourEdges = kMesh[edge.ik].edge;
						for(const EdgeFD& reverseEdge: neighbourEdges)
							if(reverseEdge.ik==ik)
							{	edge.M0 = dagger(reverseEdge.M0);
								foundReverse = true;
								break;
							}
					}
					//Compute overlap if reverse edge not yet computed:
					if(!foundReverse)
						edge.M0 = OCi ^ getWfns(edge.point, centers, iSpin);
				}
			if(!jProcess) //Do only once (will get here multiple times for local wfns)
			{	//Initial rotation:
				if(!readInitialMats)
				{	matrix CdagOg = OCi ^ trialWfns(kMesh[ik].point, centers);
					kMesh[ik].U0 = CdagOg * invsqrt(dagger(CdagOg) * CdagOg);
				}
			}
		}
	}
	
	//Broadcast overlaps and initial rotations:
	for(size_t ik=0; ik<kMesh.size(); ik++)
	{	for(EdgeFD& edge: kMesh[ik].edge)
		{	if(!isMine_q(ik,iSpin)) edge.M0 = zeroes(centers.size(), centers.size());
			edge.M0.bcast(whose_q(ik,iSpin));
			if(!isMine(ik)) edge.M0 = matrix(); //not needed any more on this process
		}
		if(!readInitialMats)
		{	if(!isMine_q(ik,iSpin)) kMesh[ik].U0 = zeroes(centers.size(), centers.size());
			kMesh[ik].U0.bcast(whose_q(ik,iSpin));
		}
		kMesh[ik].B = zeroes(centers.size(), centers.size());
	}
	
	//Apply initial rotations to the overlap matrices:
	for(size_t ik=ikStart; ik<ikStop; ik++)
		for(EdgeFD& edge: kMesh[ik].edge)
			edge.M0 = dagger(kMesh[ik].U0) * edge.M0 * kMesh[edge.ik].U0;
	
	//Minimize:
	minimize(wannier.minParams);
	
	//Save the matrices:
	{	string fname = e.dump.getFilename(UvarName);
		logPrintf("Dumping '%s' ... ", fname.c_str());
		if(mpiUtil->isHead())
		{	FILE* fp = fopen(fname.c_str(), "w");
			for(const auto& kMeshEntry: kMesh)
				(kMeshEntry.U0 * kMeshEntry.V).write(fp);
			fclose(fp);
		}
		logPrintf("done.\n"); logFlush();
	}
	
	//Save supercell wavefunctions:
	const vector3<int>& nSuper = wannier.supercell;
	vector3<int> offsetSuper(-nSuper[0]/2, -nSuper[1]/2, -nSuper[2]/2); //supercell offset
	const vector3<int>& S = e.gInfo.S;
	size_t nrSuper = e.gInfo.nr * nSuper[0]*nSuper[1]*nSuper[2];
	complex *psiSuper = new complex[nrSuper], *phaseSuper = new complex[nSuper[2]];
	for(unsigned n=0; n<centers.size(); n++)
	{	//Generate filename
		ostringstream varName;
		varName << centers[n].band << ".mlwf";
		if(nSpins==2) varName << (iSpin==0 ? "Up" : "Dn");
		string fname = e.dump.getFilename(varName.str());
		logPrintf("Dumping '%s':\n", fname.c_str());
		//Print stats for function:
		vector3<> rCoords = e.iInfo.coordsType==CoordsCartesian
			? rExpect[n] : e.gInfo.invR * rExpect[n]; //r in coordinate system of choice
		logPrintf("\tCenter: [ %lg %lg %lg ] (%s coords)\n", rCoords[0], rCoords[1], rCoords[2],
			e.iInfo.coordsType==CoordsCartesian ? "cartesian" : "lattice");
		logPrintf("\tSpread: %lg bohrs\n", sqrt(rSqExpect[n] - rExpect[n].length_squared()));
		logFlush();
		
		//Generate supercell function:
		eblas_zero(nrSuper, psiSuper);
		for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
		{	complexDataRptr psi = 
				I( (getWfns(kMesh[i].point, centers, iSpin) //original wavefunctions transformed to common basis
					* (kMesh[i].U0 * kMesh[i].V)(0,centers.size(), n,n+1) //combined to n'th localized function
					).getColumn(0) ); //expand to full-G space and then put in real space
			multiplyBlochPhase(psi, kMesh[i].point.k); //multiply by exp(i k.r) (for r in base unit cell)
			//Accumulate with appropriate phases in each unit cell
			complex* psiData = psi->data();
			vector3<int> iSuper, iCell; //index of supercell and index within cell
			for(iCell[0]=0; iCell[0]<S[0]; iCell[0]++)
			for(iCell[1]=0; iCell[1]<S[1]; iCell[1]++)
				for(iSuper[0]=0; iSuper[0]<nSuper[0]; iSuper[0]++)
				for(iSuper[1]=0; iSuper[1]<nSuper[1]; iSuper[1]++)
				{	size_t inOffset = S[2]*(iCell[1] + S[1]*iCell[0]);
					size_t outOffset = S[2]*nSuper[2]*
						(iCell[1] + S[1]*(iSuper[1] + nSuper[1]*
							(iCell[0] + S[0]*iSuper[0]) ) );
					//Compute the supercell phases for each iSuper[2]:
					for(iSuper[2]=0; iSuper[2]<nSuper[2]; iSuper[2]++)
						phaseSuper[iSuper[2]] = cis(2*M_PI*dot(kMesh[i].point.k, offsetSuper+iSuper));
					//Accumulate wavefunction for all iSuper[2] and iCell[2] using BLAS2:
					complex alpha = kMesh[i].wk;
					cblas_zgeru(CblasColMajor, S[2], nSuper[2], &alpha,
						psiData+inOffset, 1, phaseSuper, 1,
						psiSuper+outOffset, S[2]);
				}
		}
		mpiUtil->allReduce((double*)psiSuper, nrSuper*2, MPIUtil::ReduceSum);
		FILE* fp = 0;
		if(mpiUtil->isHead())
		{	fp= fopen(fname.c_str(), "wb");
			if(!fp) die("Failed to open file '%s' for binary write.\n", fname.c_str());
		}
		if(wannier.convertReal)
		{	//Convert to a real wavefunction:
			double meanPhase, sigmaPhase, rmsImagErr;
			removePhase(nrSuper, psiSuper, meanPhase, sigmaPhase, rmsImagErr);
			logPrintf("\tPhase = %lf +/- %lf\n", meanPhase, sigmaPhase); logFlush();
			logPrintf("\tRMS imaginary part = %le (after phase removal)\n", rmsImagErr);
			logFlush();
			//Write real part of supercell wavefunction to file:
			if(fp)
				for(size_t i=0; i<nrSuper; i++)
					fwrite(psiSuper+i, sizeof(double), 1, fp);
		}
		else
		{	//Write complex function as is:
			if(fp) fwrite(psiSuper, sizeof(complex), nrSuper, fp);
		}
		if(fp) fclose(fp);
	}
	delete[] psiSuper;
	delete[] phaseSuper;
	
	return; //TODO: Add Hsub read-in so that the following works
	//Save Hamiltonian in Wannier basis:
	std::vector<matrix> Hwannier(nSuper[0]*nSuper[1]*nSuper[2]);
	for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
	{	//Fetch Hamiltonian for subset of bands in center:
		matrix Hsub(centers.size(), centers.size());
		{	complex* HsubData = Hsub.data();
			const matrix& HsubFull = e.eVars.Hsub[kMesh[i].point.q + iSpin*qCount];
			const complex* HsubFullData = HsubFull.data();
			for(size_t c1=0; c1<centers.size(); c1++)
				for(size_t c2=0; c2<centers.size(); c2++)
					HsubData[Hsub.index(c1,c2)] = HsubFullData[HsubFull.index(centers[c1].band,centers[c2].band)];
		}
		//Apply MLWF-optimizd rotation:
		matrix U = kMesh[i].U0 * kMesh[i].V;
		Hsub = dagger(U) * Hsub * U;
		//Accumulate with each requested Bloch phase
		vector3<int> iSuper;
		std::vector<matrix>::iterator HwannierIter = Hwannier.begin();
		for(iSuper[0]=0; iSuper[0]<nSuper[0]; iSuper[0]++)
		for(iSuper[1]=0; iSuper[1]<nSuper[1]; iSuper[1]++)
		for(iSuper[2]=0; iSuper[2]<nSuper[2]; iSuper[2]++)
			*(HwannierIter++) += (kMesh[i].wk * cis(2*M_PI*dot(kMesh[i].point.k, offsetSuper+iSuper))) * Hsub;
	}
	for(matrix& H: Hwannier) H.allReduce(MPIUtil::ReduceSum);
	//-- save to file
	{	string Hvarname = UvarName;
		Hvarname.replace(Hvarname.find("Umlwf"), 5, "Hmlwf");
		string fname = e.dump.getFilename(Hvarname);
		logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
		if(mpiUtil->isHead())
		{	FILE* fp = fopen(fname.c_str(), "wb");
			if(!fp) die("Failed to open file '%s' for binary write.\n", fname.c_str());
			for(matrix& H: Hwannier)
			{	if(wannier.convertReal)
					H.write_real(fp);
				else
					H.write(fp);
			}
			fclose(fp);
		}
		logPrintf("done.\n"); logFlush();
	}
}


void WannierEval::step(const WannierGradient& grad, double alpha)
{	assert(grad.size()==kMesh.size());
	for(unsigned i=0; i<kMesh.size(); i++)
		axpy(alpha, grad[i], kMesh[i].B);
}

double WannierEval::compute(WannierGradient* grad)
{	int nCenters = kMesh[0].B.nRows();
	//Compute the unitary matrices:
	for(size_t i=0; i<kMesh.size(); i++)
		kMesh[i].V = cis(kMesh[i].B, &kMesh[i].Bevecs, &kMesh[i].Beigs);
	
	//Compute the expectation values of r and rSq for each center (split over processes)
	rSqExpect.assign(nCenters, 0.);
	rExpect.assign(nCenters, vector3<>(0,0,0));
	for(size_t i=ikStart; i<ikStop; i++)
		for(EdgeFD& edge: kMesh[i].edge)
		{	unsigned j = edge.ik;
			const matrix M = dagger(kMesh[i].V) * edge.M0 * kMesh[j].V;
			const complex* Mdata = M.data();
			for(int n=0; n<nCenters; n++)
			{	complex Mnn = Mdata[M.index(n,n)];
				double argMnn = atan2(Mnn.imag(), Mnn.real());
				rExpect[n] -= (kMesh[i].wk * edge.wb * argMnn) * edge.b;
				rSqExpect[n] += kMesh[i].wk * edge.wb * (argMnn*argMnn + 1. - Mnn.norm());
			}
		}
	mpiUtil->allReduce(rSqExpect.data(), nCenters, MPIUtil::ReduceSum);
	mpiUtil->allReduce((double*)rExpect.data(), 3*nCenters, MPIUtil::ReduceSum);
	
	//Compute the mean variance of the Wannier centers
	double rVariance = 0.;
	for(int n=0; n<nCenters; n++)
		rVariance += (1./nCenters) * (rSqExpect[n] - rExpect[n].length_squared());
	
	//Compute the gradients of the mean variance (if required)
	if(grad)
	{	//Allocate and initialize all gradients to zero:
		grad->resize(kMesh.size());
		for(size_t i=0; i<kMesh.size(); i++)
			(*grad)[i] = zeroes(nCenters, nCenters);
		//Accumulate gradients from each edge (split over processes):
		for(size_t i=ikStart; i<ikStop; i++)
			for(EdgeFD& edge: kMesh[i].edge)
			{	unsigned j = edge.ik;
				const matrix M = dagger(kMesh[i].V) * edge.M0 * kMesh[j].V;
				//Compute drVariance/dM:
				matrix rVariance_M = zeroes(nCenters, nCenters);
				const complex* Mdata = M.data();
				complex* rVariance_Mdata = rVariance_M.data();
				for(int n=0; n<nCenters; n++)
				{	complex Mnn = Mdata[M.index(n,n)];
					double argMnn = atan2(Mnn.imag(), Mnn.real());
					rVariance_Mdata[rVariance_M.index(n,n)] =
						(2./nCenters) * kMesh[i].wk * edge.wb
						* ((argMnn + dot(rExpect[n],edge.b))*complex(0,-1)/Mnn - Mnn.conj());
				}
				//Propagate to drVariance/dBi and drVariance/dBj:
				matrix F0 = kMesh[j].V * rVariance_M * dagger(kMesh[i].V);
				(*grad)[i] -= dagger_symmetrize(cis_grad(edge.M0 * F0, kMesh[i].Bevecs, kMesh[i].Beigs));
				(*grad)[j] += dagger_symmetrize(cis_grad(F0 * edge.M0, kMesh[j].Bevecs, kMesh[j].Beigs));
			}
		for(size_t i=0; i<kMesh.size(); i++) (*grad)[i].allReduce(MPIUtil::ReduceSum);
	}
	return rVariance;
}

WannierGradient WannierEval::precondition(const WannierGradient& grad)
{	assert(grad.size()==kMesh.size());
	int nCenters = grad[0].nRows();
	//Copy each matrix of gradient into a column of a giant matrix:
	matrix gradMat(nCenters*nCenters, kMesh.size());
	complex* gradMatData = gradMat.dataPref();
    for(unsigned i=0; i<kMesh.size(); i++)
		callPref(eblas_copy)(gradMatData+gradMat.index(0,i), grad[i].dataPref(), grad[i].nData());
	//Apply preconditioner:
	const matrix KgradMat = gradMat * kHelmholtzInv;
	//Copy result from each column to a small matrix per k-point:
	WannierGradient Kgrad(grad.size());
	const complex* KgradMatData = KgradMat.dataPref();
    for(unsigned i=0; i<kMesh.size(); i++)
	{	Kgrad[i].init(nCenters, nCenters, isGpuEnabled());
		callPref(eblas_copy)(Kgrad[i].dataPref(), KgradMatData+KgradMat.index(0,i), Kgrad[i].nData());
	}
	return Kgrad;
}


bool WannierEval::Kpoint::operator<(const WannierEval::Kpoint& other) const
{	if(q!=other.q) return q<other.q;
	if(iRot!=other.iRot) return iRot<other.iRot;
	if(invert!=other.invert) return invert<other.invert;
	if(!(offset==other.offset)) return offset<other.offset;
	return false; //all equal
}

bool WannierEval::Kpoint::operator==(const WannierEval::Kpoint& other) const
{	if(q!=other.q) return false;
	if(iRot!=other.iRot) return false;
	if(invert!=other.invert) return false;
	if(!(offset==other.offset)) return false;
	return true;
}


WannierEval::Index::Index(int nIndices) : nIndices(nIndices), dataPref(0)
{	data = new int[nIndices];
	#ifdef GPU_ENABLED
	dataGpu = 0;
	#endif
}
WannierEval::Index::~Index()
{	delete[] data;
	#ifdef GPU_ENABLED
	if(dataGpu) cudaFree(dataGpu);
	#endif
}
void WannierEval::Index::set()
{
	#ifdef GPU_ENABLED
	cudaMalloc(&dataGpu, sizeof(int)*nIndices); gpuErrorCheck();
	cudaMemcpy(dataGpu, data, sizeof(int)*nIndices, cudaMemcpyHostToDevice); gpuErrorCheck();
	dataPref = dataGpu;
	#else
	dataPref = data;
	#endif
}

void WannierEval::addIndex(const WannierEval::Kpoint& kpoint)
{	if(indexMap.find(kpoint)!=indexMap.end()) return; //previously computed
	//Compute transformed index array (mapping to full G-space)
	const Basis& basis = e.basis[kpoint.q];
	std::shared_ptr<Index> index(new Index(basis.nbasis));
	const matrix3<int> mRot = (~sym[kpoint.iRot]) * kpoint.invert;
	for(int j=0; j<index->nIndices; j++)
		index->data[j] = e.gInfo.fullGindex(mRot * basis.iGarr[j] - kpoint.offset);
	//Save to map:
	indexMap[kpoint] = index;
}

ColumnBundle WannierEval::getWfns(const WannierEval::Kpoint& kpoint, const std::vector<Wannier::Center>& centers, int iSpin) const
{	const Index& index = *(indexMap.find(kpoint)->second);
	ColumnBundle ret(centers.size(), basis.nbasis, &basis, 0, isGpuEnabled());
	ret.zero();
	//Pick required bands, and scatter from reduced basis to common basis with transformations:
	int q = kpoint.q + iSpin*qCount;
	const ColumnBundle& C = e.eInfo.isMine(q) ? e.eVars.C[q] : Cother[q];
	assert(C);
	for(unsigned c=0; c<centers.size(); c++)
		callPref(eblas_scatter_zdaxpy)(index.nIndices, 1., index.dataPref,
			C.dataPref()+C.index(centers[c].band,0), ret.dataPref()+ret.index(c,0));
	//Complex conjugate if inversion symmetry employed:
	if(kpoint.invert < 0)
		callPref(eblas_dscal)(ret.nData(), -1., ((double*)ret.dataPref())+1, 2); //negate the imaginary parts
	return ret;
}

//Fourier transform of hydrogenic orbitals
inline double hydrogenicTilde(double G, double a, int nIn, int l, double normPrefac)
{	int n = nIn+1 + l; //conventional principal quantum number
	double nG = n*G*a/(l+1), nGsq = nG*nG;
	double prefac = normPrefac / pow(1.+nGsq, n+1);
	switch(l)
	{	case 0:
			switch(n)
			{	case 1: return prefac;
				case 2: return prefac*8.*(-1.+nGsq);
				case 3: return prefac*9.*(3.+nGsq*(-10.+nGsq*3.));
				case 4: return prefac*64.*(-1.+nGsq*(7.+nGsq*(-7.+nGsq)));
			}
		case 1:
			switch(n)
			{	case 2: return prefac*16.*nG;
				case 3: return prefac*144.*nG*(-1.+nGsq);
				case 4: return prefac*128.*nG*(5.+nGsq*(-14.+nGsq*5.));
			}
		case 2:
			switch(n)
			{	case 3: return prefac*288.*nGsq;
				case 4: return prefac*3072.*nGsq*(-1.+nGsq);
			}
		case 3:
			switch(n)
			{	case 4: return prefac*6144.*nG*nGsq;
			}
	}
	return 0.;
}

ColumnBundle WannierEval::trialWfns(const WannierEval::Kpoint& kpoint, const std::vector<Wannier::Center>& centers) const
{	ColumnBundle ret(centers.size(), basis.nbasis, &basis, 0, isGpuEnabled());
	#ifdef GPU_ENABLED
	vector3<>* pos; cudaMalloc(&pos, sizeof(vector3<>));
	#endif
	complex* retData = ret.dataPref();
	for(auto c: centers)
	{	const DOS::Weight::OrbitalDesc& od = c.orbitalDesc;
		//--- Copy the center to GPU if necessary:
		#ifdef GPU_ENABLED
		cudaMemcpy(pos, &c.r, sizeof(vector3<>), cudaMemcpyHostToDevice);
		#else
		const vector3<>* pos = &c.r;
		#endif
		//--- Create the radial part:
		RadialFunctionG atRadial;
		double normPrefac = pow((od.l+1)/c.a,3);
		for(unsigned p=od.n+1; p<=od.n+1+2*od.l; p++)
			normPrefac *= p;
		normPrefac = 16*M_PI/(e.gInfo.detR * sqrt(normPrefac));
		atRadial.init(od.l, 0.02, e.gInfo.GmaxSphere, hydrogenicTilde, c.a, od.n, od.l, normPrefac);
		//--- Initialize the projector:
		callPref(Vnl)(basis.nbasis, basis.nbasis, 1, od.l, od.m, kpoint.k, basis.iGarrPref, e.gInfo.G, pos, atRadial, retData);
		callPref(eblas_zscal)(basis.nbasis, cis(0.5*M_PI*od.l), retData,1); //ensures odd l projectors are real
		retData += basis.nbasis;
	}
	#ifdef GPU_ENABLED
	cudaFree(pos);
	#endif
	return ret;
}
