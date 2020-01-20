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
	
	//! Calculate for one pair of transformed ik and untransformed iq
	double calc(int ikReduced, int iqReduced, double aXX, double omega,
		const diagMatrix& Fk, const ColumnBundle& CkRed, ColumnBundle* HCkRed,
		const diagMatrix& Fq, const ColumnBundle& Cq, ColumnBundle* HCq, matrix3<>* EXX_RRT=0) const;
	
private:
	friend class ExactExchange;
	const Everything& e;
	const std::vector<SpaceGroupOp>& sym; //!< symmetry matrices in lattice coordinates
	const std::vector<int>& invertList; //!< whether to add inversion explicitly over the symmetry group
	const int nSpins, nSpinor, qCount; //!< number of spin channels, spinor components and states per spin channel
	const int blockSize; //!< number of bands FFT'd together
	
	//Occupied state properties for fixed-H mode:
	std::vector<diagMatrix> Focc;
	std::vector<ColumnBundle> Cocc;
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
{	static StopWatch watch("ExactExchange"); watch.start();

	//Prepare outputs:
	if(HC)
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
			if(!(*HC)[q])
			{	(*HC)[q] = C[q].similar();
				(*HC)[q].zero();
			}
	
	//Calculate:
	double EXX = 0.;
	matrix3<> EXX_RRT; //computed only if EXX_RRTptr is non-null
	for(int iSpin=0; iSpin<eval->nSpins; iSpin++)
		for(int ikReduced=0; ikReduced<eval->qCount; ikReduced++)
		{
			//Prepare (reduced) ik state and gradient on all processes:
			int ikSrc = ikReduced + iSpin*eval->qCount; //source state number
			ColumnBundle CkTmp, HCkTmp;
			diagMatrix Fk(e.eInfo.nBands);
			if(e.eInfo.isMine(ikSrc))
				Fk = F[ikSrc];
			else
			{	CkTmp.init(e.eInfo.nBands, e.basis[ikSrc].nbasis*eval->nSpinor, &(e.basis[ikSrc]), &(e.eInfo.qnums[ikSrc]), isGpuEnabled());
				if(HC) { HCkTmp = CkTmp.similar(); HCkTmp.zero(); }
			}
			ColumnBundle& CkRed = e.eInfo.isMine(ikSrc) ? (ColumnBundle&)(C[ikSrc]) : CkTmp; 
			mpiWorld->bcastData(CkRed, e.eInfo.whose(ikSrc));
			mpiWorld->bcastData(Fk, e.eInfo.whose(ikSrc));
			ColumnBundle& HCkRed = (HC && e.eInfo.isMine(ikSrc)) ? (*HC)[ikSrc] : HCkTmp;
			
			//Calculate energy (and gradient):
			for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
				EXX += eval->calc(ikReduced, q-iSpin*eval->qCount, aXX, omega,
					Fk, CkRed, HC ? &HCkRed : 0,
					F[q], C[q], HC ? &(HC->at(q)) : 0,
					EXX_RRTptr ? &EXX_RRT : 0);
			
			//Move ik state gradient back to host process (if necessary):
			if(HC) mpiWorld->reduceData(HCkRed, MPIUtil::ReduceSum, e.eInfo.whose(ikSrc));
		}
	mpiWorld->allReduce(EXX, MPIUtil::ReduceSum, true);
	if(EXX_RRTptr)
	{	mpiWorld->allReduce(EXX_RRT, MPIUtil::ReduceSum, true);
		*EXX_RRTptr += EXX_RRT;
		
		logPrintf("\nEXX_RRT:\n"); EXX_RRT.print(globalLog, "%12lg ");
	}
	watch.stop();
	return EXX;
}


void ExactExchange::setOccupied(const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C)
{	const double Fcut = 1e-8; //threshold for determining occupied states
	eval->Focc.resize(e.eInfo.nStates);
	eval->Cocc.resize(e.eInfo.nStates);
	
	//Make occupied fillings and wavefunctions available on all processes:
	for(int q=0; q<e.eInfo.nStates; q++)
	{	//Determine number of occupied bands:
		int nOccupied = 0;
		if(e.eInfo.isMine(q))
		{	nOccupied = F[q].nRows();
			while(nOccupied>1 and F[q][nOccupied-1]<Fcut)
				nOccupied--;
		}
		mpiWorld->bcast(nOccupied, e.eInfo.whose(q));
		//Broadcast F and C:
		if(e.eInfo.isMine(q))
		{	eval->Cocc[q] = C[q].getSub(0,nOccupied);
			eval->Focc[q] = F[q](0,nOccupied);
		}
		else
		{	eval->Cocc[q].init(nOccupied, e.basis[q].nbasis*eval->nSpinor, &(e.basis[q]), &(e.eInfo.qnums[q]), isGpuEnabled());
			eval->Focc[q].resize(nOccupied);
		}
		mpiWorld->bcastData(eval->Cocc[q], e.eInfo.whose(q));
		mpiWorld->bcastData(eval->Focc[q], e.eInfo.whose(q));
	}
}


void ExactExchange::applyHamiltonian(double aXX, double omega, int q, const ColumnBundle& Cq, ColumnBundle& HCq)
{	int iSpin = q / eval->qCount;
	diagMatrix Fq = eval->Focc[q]; Fq.resize(Cq.nCols());
	for(int ikReduced=0; ikReduced<eval->qCount; ikReduced++)
	{	int ikSrc = ikReduced + iSpin*eval->qCount;
		eval->calc(ikReduced, q-iSpin*eval->qCount, 2.*aXX, omega, //factor of 2 to account for k <-> q symmetry
				eval->Focc[ikSrc], eval->Cocc[ikSrc], 0,
				Fq, Cq, &HCq);
	}
}


//--------------- class ExactExchangeEval implementation ----------------------


ExactExchangeEval::ExactExchangeEval(const Everything& e)
: e(e),
	sym(e.symm.getMatrices()),
	invertList(e.symm.getKpointInvertList()),
	nSpins(e.eInfo.nSpins()),
	nSpinor(e.eInfo.spinorLength()),
	qCount(e.eInfo.nStates/nSpins),
	blockSize(e.cntrl.exxBlockSize)
{
	//Find all symmtries relating each kmesh point to corresponding reduced point:
	const Supercell& supercell = *(e.coulombParams.supercell);
	struct Ktransform { SpaceGroupOp sym; int invert; int multiplicity;
		inline bool operator==(const Ktransform& kt) const
		{	if(invert != kt.invert) return false;
			if(sym.rot != kt.sym.rot) return false;
			return (sym.a - kt.sym.a).length_squared() < symmThresholdSq;
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
	for(int iq=0; iq<qCount; iq++)
	for(int jq=0; jq<qCount; jq++)
		for(const Ktransform& kt: transforms[iq][jq])
		{	KpairEntry kpair;
			kpair.sym = kt.sym;
			kpair.invert = kt.invert;
			kpair.k = kpair.sym.applyRecip(e.eInfo.qnums[iq].k) * kpair.invert; 
			kpair.weight = e.eInfo.spinWeight * pow(kmesh.size(),-2) * kt.multiplicity;
			if(e.eInfo.isMine(jq) || e.eInfo.isMine(jq + qCount))
				kpair.setup(e, iq);
			kpairs[iq][jq].push_back(kpair);
		}
	
	//Print cost estimate to give the user some idea of how long it might take!
	size_t nkPairs = bestScore;
	logPrintf("Reduced %lu k-pairs to %lu under symmetries.\n", kmesh.size()*kmesh.size(), nkPairs);
	double costFFT = e.eInfo.nStates * e.eInfo.nBands * 9.*e.gInfo.nr*log(e.gInfo.nr);
	double costBLAS3 = e.eInfo.nStates * pow(e.eInfo.nBands,2) * e.basis[0].nbasis;
	double costSemiLocal = (8 * costBLAS3 + 3 * costFFT) * (qCount * e.eInfo.nBands); //rough estimate
	double costEXX = costFFT * ( nkPairs * std::pow(e.eInfo.nBands,2) );
	double relativeCost = 1+costEXX/costSemiLocal;
	double relativeCostOrder = pow(10, floor(log(relativeCost)/log(10)));
	relativeCost = round(relativeCost/relativeCostOrder) * relativeCostOrder; //eliminate extra sigfigs
	logPrintf("Per-iteration cost relative to semi-local calculation ~ %lg\n", relativeCost);
}


double ExactExchangeEval::calc(int ikReduced, int iqReduced, double aXX, double omega,
	const diagMatrix& Fk, const ColumnBundle& CkRed, ColumnBundle* HCkRed,
	const diagMatrix& Fq, const ColumnBundle& Cq, ColumnBundle* HCq, matrix3<>* EXX_RRT) const
{
	const QuantumNumber& qnum_q = *(Cq.qnum);
	if(CkRed.qnum->spin != qnum_q.spin) return 0.;
	int nBlocks = ceildiv(Cq.nCols(), blockSize);
	double EXX = 0.;
	for(const KpairEntry& kpair: kpairs[ikReduced][iqReduced])
	{	//Prepare ik state and gradient:
		const Basis& basis_k = *(kpair.basis);
		QuantumNumber qnum_k = *(CkRed.qnum); qnum_k.k =  kpair.k;
		ColumnBundle Ck(CkRed.nCols(), basis_k.nbasis*nSpinor, &basis_k, &qnum_k, isGpuEnabled()), HCk;
		if(HCkRed) { HCk = Ck.similar(); HCk.zero(); }
		Ck.zero();
		kpair.transform->scatterAxpy(1., CkRed, Ck,0,1);
		const double prefac = -0.5*aXX * kpair.weight / (qnum_k.weight * qnum_q.weight);
		
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
			for(int bk=0; bk<Ck.nCols(); bk++)
			{	//Put this state in real space:
				std::vector<complexScalarField> Ipsik(nSpinor), grad_Ipsik(nSpinor);
				for(int s=0; s<nSpinor; s++)
					Ipsik[s] = I(Ck.getColumn(bk,s));
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
					if(HCq or HCk)
					{	complexScalarField E_In = Jdag(Kn);
						for(int s=0; s<nSpinor; s++)
						{	if(HCk) grad_Ipsik[s] += (prefac*wFq) * conj(E_In) * Ipsiq[bq-bqStart][s];
							if(HCq) grad_Ipsiq[bq-bqStart][s] += (prefac*wFk) * E_In * Ipsik[s];
						}
					}
					if(EXX_RRT) *EXX_RRT += (prefac*wFk*wFq) * e.coulombWfns->latticeGradient(n, qnum_q.k-qnum_k.k, omega); //Stress contribution
				}
				//Convert k-state gradients back to reciprocal space (if needed):
				if(HCk)
				{	for(int s=0; s<nSpinor; s++)
						if(grad_Ipsik[s])
							HCk.accumColumn(bk,s, Idag(grad_Ipsik[s]));
				}
			}
			//Convert q-state gradients back to reciprocal space (if needed):
			if(HCq)
			{	for(int bq=bqStart; bq<bqStop; bq++)
					for(int s=0; s<nSpinor; s++)
						if(grad_Ipsiq[bq-bqStart][s])
							HCq->accumColumn(bq,s, Idag(grad_Ipsiq[bq-bqStart][s]));
			}
			bqStart = bqStop;
		}
		
		//Transform HCk to HCkRed (if needed)
		if(HCk) kpair.transform->gatherAxpy(1., HCk,0,1, *HCkRed);
	}
	return EXX;
}

void ExactExchange::addHamiltonian(double aXX, double omega, int q, matrix& H,
	const std::vector<int>& iRowsMine, const std::vector<int>& iColsMine)
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
