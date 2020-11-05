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

#include <electronic/ExactExchange.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/ColumnBundleTransform.h>
#include <core/Util.h>
#include <core/GpuUtil.h>
#include <core/Operators.h>
#include <core/LatticeUtils.h>
#include <core/Random.h>
#include <list>

//! Internal computation object for ExactExchange
class ExactExchangeEval
{
public:
	ExactExchangeEval(const Everything& e);
	
	struct KpairEntry
	{	vector3<> k; //transformed k
		SpaceGroupOp sym; //symmetry operation
		int invert; //whether inversion involved (+/-1)
		double weight;
		std::shared_ptr<Basis> basis;
		std::shared_ptr<ColumnBundleTransform> transform; //wavefunction transformation from reduced set
		
		void setup(const Everything& e, int q, bool needTransform=true)
		{	basis = std::make_shared<Basis>();
			logSuspend();
			basis->setup((e.gInfoWfns ? *(e.gInfoWfns) : e.gInfo), e.iInfo, e.cntrl.Ecut, k);
			logResume();
			if(needTransform)
				transform = std::make_shared<ColumnBundleTransform>(e.eInfo.qnums[q].k, e.basis[q], k, *basis, e.eInfo.spinorLength(), sym, invert);
		}
	};
	std::vector<std::vector<std::vector<KpairEntry>>> kpairs; //list of transformations (inner index) for each pair of untransformed q (middle index) and transformed k (outer index)
	
	//! Calculate exchange operator and return energy
	double compute(double aXX, double omega,
		const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C,
		std::vector<ColumnBundle>* HC=0, matrix3<>* EXX_RRT=0) const;
	
	//! Calculate exchange contributions for one pair of transformed ik and untransformed iq.
	//! Gradients are only accumulated to untransformed iq, taking advantage of Hermitian symmetry of exchange operator.
	double computePair(int ikReduced, int iqReduced, size_t& progress, size_t& progressTarget, double aXX, double omega,
		const diagMatrix& Fk, const ColumnBundle& CkRed, const diagMatrix& Fq, const ColumnBundle& Cq,
		ColumnBundle* HCq, matrix3<>* EXX_RRT=0) const;
	
private:
	friend class ExactExchange;
	const Everything& e;
	const std::vector<SpaceGroupOp>& sym; //!< symmetry matrices in lattice coordinates
	const std::vector<int>& invertList; //!< whether to add inversion explicitly over the symmetry group
	const int nSpins, nSpinor, qCount; //!< number of spin channels, spinor components and states per spin channel
	const int blockSize; //!< number of bands FFT'd together
	double omegaACE; //!< omega for which ACE has been initialized (NAN if none)
	std::vector<ColumnBundle> psiACE; //!< projectors for ACE representation of exchange Hamiltonian

	//Local chunks of untransformed q-state wavefunctions used for re-organizing q-states for load balancing
	struct LocalState
	{	int iqReduced; //index into reduced q mesh (without spin)
		int bStart, bStop; //start and stop indices of current chunk within original wavefunctions
		ColumnBundle Cq, HCq; //chunk of wavefunctions and gradients to be processed
		diagMatrix Fq; //corresponding occupations
	};
	std::vector<std::vector<LocalState>> localStates; //local state descriptions on each process
	std::vector<LocalState>& localStatesMine; //reference to local state on this process
	size_t progressMax, progressInterval; //for progress reporting in compute / computePair
};


ExactExchange::ExactExchange(const Everything& e) : e(e)
{
	logPrintf("\n---------- Setting up exact exchange ----------\n");
	eval = new ExactExchangeEval(e);
}

ExactExchange::~ExactExchange()
{
	delete eval;
}

double ExactExchange::operator()(double aXX, double omega, 
	const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C,
	std::vector<ColumnBundle>* HC, matrix3<>* EXX_RRTptr) const
{
	if((omega == eval->omegaACE) and (not EXX_RRTptr))
	{	//Use previously initialized ACE representation
		double EXX = 0.;
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		{	ColumnBundle HCtmp; 
			ColumnBundle& HCq = HC ? HC->at(q) : HCtmp;
			EXX += applyHamiltonian(aXX, omega, q, F[q], C[q], HCq);
		}
		mpiWorld->allReduce(EXX, MPIUtil::ReduceSum);
		return EXX;
	}
	else
	{	//Compute full operator (if no ACE ready, or need lattice gradient):
		logPrintf("Computing exact exchange ... "); logFlush();
		double EXX = eval->compute(aXX, omega, F, C, HC, EXX_RRTptr);
		logPrintf("done.\n");
		return EXX;
	}
}

//Initialize the ACE (Adiabatic Compression of Exchange) representation in preparation for applyHamiltonian
void ExactExchange::prepareHamiltonian(double omega, const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C)
{	logPrintf("Constructing ACE exchange operator ... "); logFlush();
	//Compute the full exchange operator evaluated on each orbital of C in W
	std::vector<ColumnBundle> W(e.eInfo.nStates);
	eval->compute(-1., omega, F, C, &W); //compute with aXX = 1; applyHamiltonian handles the overall scale factor aXX
	//Construct ACE representation:
	eval->psiACE.assign(e.eInfo.nStates, ColumnBundle()); //clear any previous results
	bool isSingularAny = false;
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{	matrix M = C[q] ^ W[q]; //positive semi-definite matrix = dagger(C) * (-Vxx) * C (because aXX = -1 used above)
		bool isSingular = false;
		eval->psiACE[q] = W[q] * invsqrt(M, 0, 0, &isSingular); //Same effect as Cholesky factorization, but done via diagonalization
		W[q].free(); //clear memory
		isSingularAny = isSingularAny or isSingular;
	}
	logPrintf("done.\n");
	//Check and report any singular inversions:
	mpiWorld->allReduce(isSingularAny, MPIUtil::ReduceLOr);
	if(isSingularAny) logPrintf("WARNING: singularity encountered in constructing ACE representation.\n");
	//Mark ACE ready at specified omega:
	eval->omegaACE = omega;
}

//Apply Hamiltonian using ACE representation initialized previously
double ExactExchange::applyHamiltonian(double aXX, double omega, int q, const diagMatrix& Fq, const ColumnBundle& Cq, ColumnBundle& HCq) const
{	assert(omega == eval->omegaACE); //Confirm that ACE representation is ready at required omega
	const ColumnBundle& psi = eval->psiACE[q]; //ACE projectors for current q
	matrix psiDagC = psi ^ Cq;
	if(HCq) HCq -= aXX * (psi * psiDagC);
	return (-0.5 * aXX * e.eInfo.qnums[q].weight) * trace(psiDagC * Fq * dagger(psiDagC)).real();
}

//Construct dense version of EXX Hamiltonian
void ExactExchange::addHamiltonian(double aXX, double omega, int q, matrix& H,
	const std::vector<int>& iRowsMine, const std::vector<int>& iColsMine) const
{
	logPrintf("\t\tInitializing EXX:"); logFlush();
	assert(eval->nSpinor == 1); //currently implemented only for non-spinorial calculations
	
	//Setup:
	int iSpin = q / eval->qCount;
	int iqReduced = q - iSpin*eval->qCount;
	int nColOrbitals = 0;
	for(int ikReduced=0; ikReduced<eval->qCount; ikReduced++)
		nColOrbitals += eval->kpairs[ikReduced][iqReduced].size() * e.eInfo.nBands * iColsMine.size();
	int iColOrbInterval = std::max(1, int(round(nColOrbitals/20.))); //interval for reporting progress
	const double Fcut = 1e-8; //threshold for determining occupied states
	const QuantumNumber& qnum_q = e.eInfo.qnums[q];
	const Basis& basis_q = e.basis[q];
	const vector3<int>* iGqArr = basis_q.iGarr.data();
	const GridInfo& gInfoWfns = e.gInfoWfns ? *(e.gInfoWfns) : e.gInfo;
	
	//Loop over k orbitals:
	int iColOrbital = 0;
	for(int ikReduced=0; ikReduced<eval->qCount; ikReduced++)
	{	int ikSrc = ikReduced + iSpin*eval->qCount;
		for(ExactExchangeEval::KpairEntry kpair: eval->kpairs[ikReduced][iqReduced])
		{	//Set up qnum and basis on all processes (and column transform on ikSrc owner):
			QuantumNumber qnum_k = e.eInfo.qnums[ikSrc]; qnum_k.k =  kpair.k;
			kpair.setup(e, ikReduced, e.eInfo.isMine(ikSrc));
			const Basis& basis_k = *(kpair.basis);
			const vector3<int>* iGkArr = basis_k.iGarr.data();
			ColumnBundle Ckb(1, basis_k.nbasis, &basis_k, &qnum_k, isGpuEnabled());
			//Broadcast fillings:
			diagMatrix Fk(e.eInfo.nBands);
			if(e.eInfo.isMine(ikSrc))
				Fk = e.eVars.F[ikSrc];
			mpiWorld->bcastData(Fk, e.eInfo.whose(ikSrc));
			//Loop over occupied k bands:
			for(int b=0; b<e.eInfo.nBands; b++) 
			{	if(Fk[b] > Fcut)
				{
					//Broadcast k orbital:
					if(e.eInfo.isMine(ikSrc))
					{	Ckb.zero();
						kpair.transform->scatterAxpy(1., e.eVars.C[ikSrc],b, Ckb,0);
					}
					mpiWorld->bcastData(Ckb, e.eInfo.whose(ikSrc));
					const complex* CkbData = Ckb.data();
					const double prefac = -aXX * Fk[b] * kpair.weight / qnum_q.weight;
					
					//Loop over local columns of iG:
					complexScalarFieldTilde psi_kbConj;
					nullToZero(psi_kbConj, gInfoWfns);
					complex* Hdata = H.data();
					for(int j: iColsMine)
					{
						//Expand psi_kb with Gj offset:
						psi_kbConj->zero();
						complex* psiData = psi_kbConj->data();
						for(size_t n=0; n<basis_k.nbasis; n++)
							psiData[gInfoWfns.fullGindex(iGqArr[j] - iGkArr[n])] = CkbData[n].conj();
						
						//Apply Coulomb operator:
						complexScalarFieldTilde Kpsi = (*e.coulombWfns)(psi_kbConj, qnum_q.k-qnum_k.k, omega);
						const complex* KpsiData = Kpsi->data();
						
						//Collect results for each row:
						for(int i: iRowsMine)
						{	complex result;
							for(size_t n=0; n<basis_k.nbasis; n++)
								result += CkbData[n] * KpsiData[gInfoWfns.fullGindex(iGqArr[i] - iGkArr[n])];
							*(Hdata++) += prefac * result;
						}
						
						//Report progress:
						iColOrbital++;
						if(iColOrbital % iColOrbInterval == 0)
						{	logPrintf(" %d%%", int(round(iColOrbital*100./nColOrbitals)));
							logFlush();
						}
					}
				}
				else iColOrbital += iColsMine.size();
			}
		}
	}
	logPrintf(" done.\n"); logFlush();
}

//--------------- class ExactExchangeEval implementation ----------------------


ExactExchangeEval::ExactExchangeEval(const Everything& e)
: e(e),
	sym(e.symm.getMatrices()),
	invertList(e.symm.getKpointInvertList()),
	nSpins(e.eInfo.nSpins()),
	nSpinor(e.eInfo.spinorLength()),
	qCount(e.eInfo.nStates/nSpins),
	blockSize(e.cntrl.exxBlockSize),
	omegaACE(NAN),
	localStates(mpiWorld->nProcesses()),
	localStatesMine(localStates[mpiWorld->iProcess()])
{
	//Find all symmtries relating each kmesh point to corresponding reduced point:
	const Supercell& supercell = *(e.coulombParams.supercell);
	struct Ktransform { SpaceGroupOp sym; int invert; int multiplicity;
		inline bool operator==(const Ktransform& kt) const
		{	if(invert != kt.invert) return false;
			if(sym.rot != kt.sym.rot) return false;
			return circDistanceSquared(sym.a, kt.sym.a) < symmThresholdSq;
		}
	};
	struct KmeshEntry { int iReduced; std::vector<Ktransform> transform; };
	struct Kmesh : public std::vector<KmeshEntry>
	{	int qCount;
		Kmesh(size_t n, int qCount) : std::vector<KmeshEntry>(n), qCount(qCount) {}
		//Function to generate unique pair transforms for each transform choice (and thereby score this choice)
		int score(const std::vector<int> choice, std::vector<std::vector<std::vector<Ktransform>>>& transforms)
		{	transforms.assign(qCount, std::vector<std::vector<Ktransform>>(qCount));
			int nTransforms = 0;
			for(size_t ik=0; ik<size(); ik++)
			{	const int iq = (*this)[ik].iReduced;
				const Ktransform& kti = (*this)[ik].transform[choice[ik]];
				for(size_t jk=0; jk<size(); jk++)
				{	const int jq = (*this)[jk].iReduced;
					const Ktransform& ktj = (*this)[jk].transform[choice[jk]];
					//Create net transform:
					Ktransform kt;
					kt.sym = kti.sym * ktj.sym.inv();
					kt.invert = kti.invert * ktj.invert;
					kt.multiplicity = 1;
					//Add if unique:
					bool found = false;
					for(Ktransform& prev: transforms[iq][jq])
						if(prev == kt)
						{	found = true;
							prev.multiplicity++;
							break;
						}
					if(not found)
					{	transforms[iq][jq].push_back(kt);
						nTransforms++;
					}
				}
			}
			return nTransforms;
		}
	}
	kmesh(supercell.kmesh.size(), qCount);
	
	for(unsigned ik=0; ik<supercell.kmesh.size(); ik++)
	{	KmeshEntry& ki = kmesh[ik];
		ki.iReduced = supercell.kmeshTransform[ik].iReduced;
		for(int invert: invertList)
			for(const SpaceGroupOp& op: sym)
			{	vector3<> k = invert * op.applyRecip(e.eInfo.qnums[ki.iReduced].k);
				if(circDistanceSquared(k, supercell.kmesh[ik]) < symmThresholdSq)
				{	Ktransform kt = { op, invert };
					ki.transform.push_back(kt);
				}
			}
	}
	
	//Find best choice for transforms:
	logPrintf("Optimizing transforms to minimize k-point pairs ... "); logFlush();
	std::vector<std::vector<std::vector<Ktransform>>> transforms;
	std::vector<int> choices(kmesh.size(), 0), bestChoices = choices;
	int nSteps = 0;
	for(size_t ik=0; ik<kmesh.size(); ik++)
	{	int nChoices = kmesh[ik].transform.size();
		choices[ik] = Random::uniformInt(nChoices);
		nSteps += 10*(nChoices-1);
	}
	int score = kmesh.score(choices, transforms), bestScore = score;
	for(int step=0; step<nSteps;)
	{	//Random step:
		std::vector<int> newChoices = choices;
		int ikStep = Random::uniformInt(kmesh.size());
		newChoices[ikStep] = Random::uniformInt(kmesh[ikStep].transform.size());
		if(newChoices[ikStep]==choices[ikStep]) continue; //null step; don't count
		int newScore = kmesh.score(newChoices, transforms);
		//Update best encountered case:
		if(newScore < bestScore)
		{	bestScore = newScore;
			bestChoices = newChoices;
		}
		//Metropolis accept-reject:
		if(Random::uniform() < exp(score-newScore))
		{	score = newScore;
			choices = newChoices;
		}
		step++;
	}
	int bestProc = mpiWorld->iProcess();
	mpiWorld->allReduce(bestScore, bestProc, MPIUtil::ReduceMin);
	mpiWorld->bcastData(bestChoices, bestProc);
	kmesh.score(bestChoices, transforms);
	logPrintf("done (%d steps).\n", nSteps); logFlush();
	
	//Set up selected transforms:
	kpairs.assign(qCount, std::vector<std::vector<KpairEntry>>(qCount));
	size_t nTransformsMin = transforms[0][0].size(), nTransformsMax = 0;
	std::vector<size_t> jCost(qCount); //estimated relative cost of a band in each jq
	for(int iq=0; iq<qCount; iq++)
	for(int jq=0; jq<qCount; jq++)
	{	for(const Ktransform& kt: transforms[iq][jq])
		{	KpairEntry kpair;
			kpair.sym = kt.sym;
			kpair.invert = kt.invert;
			kpair.k = kpair.sym.applyRecip(e.eInfo.qnums[iq].k) * kpair.invert; 
			kpair.weight = e.eInfo.spinWeight * pow(kmesh.size(),-2) * kt.multiplicity;
			kpairs[iq][jq].push_back(kpair); //note that kpair setup is run below after determining load balancing
		}
		nTransformsMin = std::min(nTransformsMin, transforms[iq][jq].size());
		nTransformsMax = std::max(nTransformsMax, transforms[iq][jq].size());
		jCost[jq] += transforms[iq][jq].size();
	}
	size_t nkPairs = bestScore;
	logPrintf("Reduced %lu k-pairs to %lu under symmetries.\n", kmesh.size()*kmesh.size(), nkPairs);
	logPrintf("Transforms per reduced k-pair: %lu min, %lu max, %.1lf mean.\n",
		nTransformsMin, nTransformsMax, double(nkPairs)/(qCount*qCount));
	
	//Determine band/state division for load balancing:
	std::vector<int> iqStartProc(mpiWorld->nProcesses()+1, 0);
	std::vector<int> bStartProc(mpiWorld->nProcesses()+1, 0);
	if(mpiWorld->isHead())
	{	//Determine cumulative cost estimate before a given band
		std::vector<size_t> jCostPrev(qCount*e.eInfo.nBands+1);
		int jIndex = 0;
		for(int jq=0; jq<qCount; jq++)
			for(int bq=0; bq<e.eInfo.nBands; bq++)
			{	jCostPrev[jIndex+1] = jCostPrev[jIndex] + jCost[jq];
				jIndex++;
			}
		//Determine start jIndex for each process:
		for(int jProc=0; jProc<=mpiWorld->nProcesses(); jProc++)
		{	size_t costTarget = (jProc * jCostPrev.back())/mpiWorld->nProcesses(); //for equal distribution
			jIndex = std::lower_bound(jCostPrev.begin(), jCostPrev.end(), costTarget) - jCostPrev.begin();
			//Split to reduced-state and band indices:
			iqStartProc[jProc] = jIndex / e.eInfo.nBands;
			bStartProc[jProc] = jIndex - iqStartProc[jProc] * e.eInfo.nBands;
		}
	}
	mpiWorld->bcastData(iqStartProc);
	mpiWorld->bcastData(bStartProc);
	for(int jProc=0; jProc<mpiWorld->nProcesses(); jProc++)
	{	int nLocal = iqStartProc[jProc+1] - iqStartProc[jProc] + (bStartProc[jProc+1] ? 1 : 0);
		localStates[jProc].resize(nLocal);
		for(int iLocal=0; iLocal<nLocal; iLocal++)
		{	LocalState& ls = localStates[jProc][iLocal];
			ls.iqReduced = iqStartProc[jProc] + iLocal;
			ls.bStart = (ls.iqReduced==iqStartProc[jProc]) ? bStartProc[jProc] : 0;
			ls.bStop = (ls.iqReduced==iqStartProc[jProc+1]) ? bStartProc[jProc+1] : e.eInfo.nBands;
		}
	}
	progressMax = 0;
	for(const LocalState& ls: localStatesMine)
	{	for(int iq=0; iq<qCount; iq++)
			for(KpairEntry& kpair: kpairs[iq][ls.iqReduced])
			{	kpair.setup(e, iq);
				progressMax += nSpins*(ls.bStop-ls.bStart);
			}
	}
	progressInterval = std::max(1, int(round(progressMax/50.))); //interval for reporting progress
	
	//Print cost estimate to give the user some idea of how long it might take!
	double costFFT = e.eInfo.nStates * e.eInfo.nBands * 9.*e.gInfo.nr*log(e.gInfo.nr);
	double costBLAS3 = e.eInfo.nStates * pow(e.eInfo.nBands,2) * e.basis[0].nbasis;
	double costSemiLocal = (8 * costBLAS3 + 3 * costFFT) * (qCount * e.eInfo.nBands); //rough estimate
	double costEXX = costFFT * ( nkPairs * std::pow(e.eInfo.nBands,2) );
	double relativeCost = 1+costEXX/costSemiLocal;
	double relativeCostOrder = pow(10, floor(log(relativeCost)/log(10)));
	relativeCost = round(relativeCost/relativeCostOrder) * relativeCostOrder; //eliminate extra sigfigs
	logPrintf("Per-iteration cost relative to semi-local calculation ~ %lg\n", relativeCost);
}


double ExactExchangeEval::compute(double aXX, double omega, 
	const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C,
	std::vector<ColumnBundle>* HC, matrix3<>* EXX_RRTptr) const
{
	static StopWatch watch("ExactExchange"); watch.start();
	
	double EXX = 0.;
	matrix3<> EXX_RRT; //computed only if EXX_RRTptr is non-null
	size_t progress=0, progressTarget = progressInterval; //Current progress and next threshold before reporting
	
	for(int iSpin=0; iSpin<nSpins; iSpin++)
	{
		//Redistribute chunks of q-state wavefunctions:
		std::vector<MPIUtil::Request> requests; requests.reserve(2*(localStates.size()+mpiWorld->nProcesses()));
		for(int jProc=0; jProc<mpiWorld->nProcesses(); jProc++)
		{	for(LocalState& ls: (std::vector<LocalState>&)localStates[jProc])
			{	int iqSrc = ls.iqReduced + iSpin*qCount; //source state number
				if(jProc == mpiWorld->iProcess())
				{	//Allocate:
					ls.Cq.init(ls.bStop-ls.bStart, e.basis[iqSrc].nbasis*nSpinor, &(e.basis[iqSrc]), &(e.eInfo.qnums[iqSrc]), isGpuEnabled());
					if(HC) { ls.HCq = ls.Cq.similar(); ls.HCq.zero(); }
					//Post async recv's or copy local data:
					if(e.eInfo.isMine(iqSrc))
					{	callPref(eblas_copy)(ls.Cq.dataPref(), e.eVars.C[iqSrc].dataPref() + ls.Cq.colLength()*ls.bStart, ls.Cq.nData());
						ls.Fq = e.eVars.F[iqSrc](ls.bStart, ls.bStop);
					}
					else
					{	//Recv wavefunctions:
						int tagC = iqSrc*e.eInfo.nBands+ls.bStart; requests.push_back(MPIUtil::Request());
						mpiWorld->recvData(ls.Cq, e.eInfo.whose(iqSrc), tagC, &requests.back());
						//Recv occupations:
						ls.Fq.resize(ls.bStop-ls.bStart);
						int tagF = tagC + e.eInfo.nBands*e.eInfo.nStates; requests.push_back(MPIUtil::Request());
						mpiWorld->recvData(ls.Fq, e.eInfo.whose(iqSrc), tagF, &requests.back());
					}
				}
				else if(e.eInfo.isMine(iqSrc))
				{	//Send wavefunctions:
					const ColumnBundle& Cq = e.eVars.C[iqSrc];
					int tagC = iqSrc*e.eInfo.nBands+ls.bStart; requests.push_back(MPIUtil::Request());
					mpiWorld->send(Cq.data()+Cq.index(ls.bStart,0), Cq.colLength()*(ls.bStop-ls.bStart), jProc, tagC, &requests.back());
					//Send occupations:
					int tagF = tagC + e.eInfo.nBands*e.eInfo.nStates;
					mpiWorld->send(e.eVars.F[iqSrc].data()+ls.bStart, ls.bStop-ls.bStart, jProc, tagF, &requests.back());
				}
			}
		}
		mpiWorld->waitAll(requests); requests.clear();
		
		//Compute exchange for this spin channel:
		for(int ikReduced=0; ikReduced<qCount; ikReduced++)
		{
			//Prepare (reduced) ik state on all processes:
			int ikSrc = ikReduced + iSpin*qCount; //source state number
			ColumnBundle CkTmp;
			diagMatrix Fk(e.eInfo.nBands);
			if(e.eInfo.isMine(ikSrc))
				Fk = F[ikSrc];
			else
				CkTmp.init(e.eInfo.nBands, e.basis[ikSrc].nbasis*nSpinor, &(e.basis[ikSrc]), &(e.eInfo.qnums[ikSrc]), isGpuEnabled());
			ColumnBundle& CkRed = e.eInfo.isMine(ikSrc) ? (ColumnBundle&)(C[ikSrc]) : CkTmp; 
			mpiWorld->bcastData(CkRed, e.eInfo.whose(ikSrc));
			mpiWorld->bcastData(Fk, e.eInfo.whose(ikSrc));
			
			//Calculate energy (and gradient):
			for(LocalState& ls: localStatesMine)
				EXX += computePair(ikReduced, ls.iqReduced, progress, progressTarget, aXX, omega,
					Fk, CkRed, ls.Fq, ls.Cq, HC ? &(ls.HCq) : 0, EXX_RRTptr ? &EXX_RRT : 0);
		}
		
		//Free local wavefunction chunks:
		for(LocalState& ls: localStatesMine) ls.Cq.free();
		
		//Send back gradient chunks if needed:
		if(HC)
		{	//Move ls.HCq to process where HC[q] is local:
			std::vector<MPIUtil::Request> requests; requests.reserve(localStates.size()+mpiWorld->nProcesses());
			for(int jProc=0; jProc<mpiWorld->nProcesses(); jProc++)
			{	for(LocalState& ls: (std::vector<LocalState>&)localStates[jProc])
				{	int iqSrc = ls.iqReduced + iSpin*qCount; //source state number
					if(jProc == mpiWorld->iProcess())
					{	//Post async sends (local accumulate dealt with below):
						if(not e.eInfo.isMine(iqSrc))
						{	int tagC = iqSrc*e.eInfo.nBands+ls.bStart; requests.push_back(MPIUtil::Request());
							mpiWorld->sendData(ls.HCq, e.eInfo.whose(iqSrc), tagC, &requests.back());
						}
					}
					else if(e.eInfo.isMine(iqSrc))
					{	//Post async recv to temporary buffer in ls.HCq:
						ls.HCq.init(ls.bStop-ls.bStart, e.basis[iqSrc].nbasis*nSpinor, &(e.basis[iqSrc]), &(e.eInfo.qnums[iqSrc]), isGpuEnabled());
						int tagC = iqSrc*e.eInfo.nBands+ls.bStart; requests.push_back(MPIUtil::Request());
						mpiWorld->recvData(ls.HCq, jProc, tagC, &requests.back());
					}
				}
			}
			mpiWorld->waitAll(requests);
			//Accumulate HC contributions from ls.HCq to HC:
			for(int jProc=0; jProc<mpiWorld->nProcesses(); jProc++)
				for(LocalState& ls: (std::vector<LocalState>&)localStates[jProc])
				{	int iqSrc = ls.iqReduced + iSpin*qCount; //source state number
					ColumnBundle& HCq = (*HC)[iqSrc];
					if(e.eInfo.isMine(iqSrc))
					{	if(!HCq) { HCq = C[iqSrc].similar(); HCq.zero(); } //Allocate to zero if null
						callPref(eblas_zaxpy)(ls.HCq.nData(), 1., ls.HCq.dataPref(),1, HCq.dataPref()+HCq.index(ls.bStart,0),1);
						ls.HCq.free();
					}
				}
		}
	}
	mpiWorld->allReduce(EXX, MPIUtil::ReduceSum, true);
	if(EXX_RRTptr)
	{	mpiWorld->allReduce(EXX_RRT, MPIUtil::ReduceSum, true);
		*EXX_RRTptr += EXX_RRT;
	}
	watch.stop();
	return EXX;
}

double ExactExchangeEval::computePair(int ikReduced, int iqReduced, size_t& progress, size_t& progressTarget, double aXX, double omega,
	const diagMatrix& Fk, const ColumnBundle& CkRed, const diagMatrix& Fq, const ColumnBundle& Cq,
	ColumnBundle* HCq, matrix3<>* EXX_RRT) const
{
	const QuantumNumber& qnum_q = *(Cq.qnum);
	if(CkRed.qnum->spin != qnum_q.spin) return 0.;
	int nBlocks = ceildiv(Cq.nCols(), blockSize);
	double EXX = 0.;
	//Loop over blocks:
	int bqStart = 0;
	for(int iBlock=0; iBlock<nBlocks; iBlock++)
	{	int bqStop = std::min(bqStart+blockSize, Cq.nCols());
		//Prepare q-states in real space:
		std::vector<std::vector<complexScalarField>> Ipsiq(bqStop-bqStart), grad_Ipsiq;
		if(HCq) grad_Ipsiq.assign(bqStop-bqStart, std::vector<complexScalarField>(nSpinor));
		for(int bq=bqStart; bq<bqStop; bq++)
		{	Ipsiq[bq-bqStart].resize(nSpinor);
			for(int s=0; s<nSpinor; s++)
				Ipsiq[bq-bqStart][s] = I(Cq.getColumn(bq,s));
		}
		//Loop over k-states:
		for(int bk=0; bk<CkRed.nCols(); bk++)
		{	//Loop over symmetry transformations of this k-state:
			for(const KpairEntry& kpair: kpairs[ikReduced][iqReduced])
			{	//Prepare transformed k-state in reciprocal space:
				const Basis& basis_k = *(kpair.basis);
				QuantumNumber qnum_k = *(CkRed.qnum); qnum_k.k =  kpair.k;
				ColumnBundle Ck(1, basis_k.nbasis*nSpinor, &basis_k, &qnum_k, isGpuEnabled());
				Ck.zero();
				kpair.transform->scatterAxpy(1., CkRed,bk, Ck,0);
				const double prefac = -0.5*aXX * kpair.weight / (qnum_k.weight * qnum_q.weight);
				//Put this state in real space:
				std::vector<complexScalarField> Ipsik(nSpinor);
				for(int s=0; s<nSpinor; s++)
					Ipsik[s] = I(Ck.getColumn(0,s));
				double wFk = qnum_k.weight * Fk[bk];
				//Loop over q-bands within block:
				for(int bq=bqStart; bq<bqStop; bq++)
				{	double wFq = qnum_q.weight * Fq[bq];
					if(!wFk && !wFq) continue; //at least one of the orbitals must be occupied
					complexScalarField In; //state pair density
					for(int s=0; s<nSpinor; s++)
						In += conj(Ipsik[s]) * Ipsiq[bq-bqStart][s];
					complexScalarFieldTilde n = J(In);
					complexScalarFieldTilde Kn = O((*e.coulombWfns)(n, qnum_q.k-qnum_k.k, omega)); //Electrostatic potential due to n
					EXX += (prefac*wFk*wFq) * dot(n,Kn).real();
					if(HCq)
					{	complexScalarField E_In = Jdag(Kn);
						for(int s=0; s<nSpinor; s++)
							grad_Ipsiq[bq-bqStart][s] += (2.*prefac*wFk) * E_In * Ipsik[s]; //factor of 2 to count grad_Ipsik using Hermitian symmetry
					}
					if(EXX_RRT) *EXX_RRT += (prefac*wFk*wFq) * e.coulombWfns->latticeGradient(n, qnum_q.k-qnum_k.k, omega); //Stress contribution
				}
			}
		}
		//Convert q-state gradients back to reciprocal space (if needed):
		if(HCq)
		{	for(int bq=bqStart; bq<bqStop; bq++)
				for(int s=0; s<nSpinor; s++)
					if(grad_Ipsiq[bq-bqStart][s])
						HCq->accumColumn(bq,s, Idag(grad_Ipsiq[bq-bqStart][s]));
		}
		//Report progress:
		progress += (bqStop-bqStart)*kpairs[ikReduced][iqReduced].size();
		if(progress >= progressTarget)
		{	logPrintf("%d%% ", int(round(progress*100./progressMax))); logFlush();
			progressTarget = std::min(progressTarget+progressInterval, progressMax); //next target for reporting
		}
		bqStart = bqStop;
	}
	return EXX;
}
