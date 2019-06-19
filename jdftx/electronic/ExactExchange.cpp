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
#include <list>

//! Internal computation object for ExactExchange
class ExactExchangeEval
{
public:
	ExactExchangeEval(const Everything& e);
	
	//! Calculate for one entry of the k-mesh at a particular spin:
	double calc(int iSpin, unsigned iReduced, unsigned iInvert, unsigned iSym, 
		double aXX, double omega, const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C, std::vector<ColumnBundle>* HC) const;
	
private:
	friend class ExactExchange;
	const Everything& e;
	const std::vector<SpaceGroupOp>& sym; //!< symmetry matrices in lattice coordinates
	const std::vector<int>& invertList; //!< whether to add inversion explicitly over the symmetry group
	int nSpins;
	int nSpinor;
	int qCount; //!< number of states of each spin
	//Symmetry rotation map:
	struct KmapEntry
	{	vector3<> k;
		Basis basis;
		std::shared_ptr<ColumnBundleTransform> transform; //wavefunction transformation from reduced set
	};
	std::vector<KmapEntry> kmap;
	inline int kmapIndex(int iReduced, int iInvert, int iSym) const { return (iReduced*invertList.size() + iInvert)*sym.size() + iSym; }
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
	std::vector<ColumnBundle>* HC) const
{	static StopWatch watch("ExactExchange"); watch.start();

	//prepare outputs
	if(HC)
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
			if(!(*HC)[q])
			{	(*HC)[q] = C[q].similar();
				(*HC)[q].zero();
			}
	
	//Calculate:
	double EXX = 0.0;
	for(int iSpin=0; iSpin<eval->nSpins; iSpin++)
		for(int iReduced=0; iReduced<eval->qCount; iReduced++)
		for(unsigned iInvert=0; iInvert<eval->invertList.size(); iInvert++)
		for(unsigned iSym=0; iSym<eval->sym.size(); iSym++)
			EXX += eval->calc(iSpin, iReduced, iInvert, iSym, aXX, omega, F, C, HC);
	watch.stop();
	return EXX;
}

//--------------- class ExactExchangeEval implementation ----------------------



ExactExchangeEval::ExactExchangeEval(const Everything& e)
: e(e),
	sym(e.symm.getMatrices()),
	invertList(e.symm.getKpointInvertList()),
	nSpins(e.eInfo.nSpins()),
	nSpinor(e.eInfo.spinorLength()),
	qCount(e.eInfo.nStates/nSpins),
	kmap(qCount * invertList.size() * sym.size())
{
	//Print cost estimate to give the user some idea of how long it might take!
	double costFFT = e.eInfo.nStates * e.eInfo.nBands * 9.*e.gInfo.nr*log(e.gInfo.nr);
	double costBLAS3 = e.eInfo.nStates * pow(e.eInfo.nBands,2) * e.basis[0].nbasis;
	double costSemiLocal = 8 * costBLAS3 + 3 * costFFT; //very rough estimates of course!
	double costEXX = costFFT * ( sym.size() * invertList.size() * qCount * e.eInfo.nBands );
	double relativeCost = 1+costEXX/costSemiLocal;
	double relativeCostOrder = pow(10, floor(log(relativeCost)/log(10)));
	relativeCost = round(relativeCost/relativeCostOrder) * relativeCostOrder; //eliminate extra sigfigs
	logPrintf("Per-iteration cost relative to semi-local calculation ~ %lg\n", relativeCost);
	if(qCount==1 && sym.size()>1)
		logPrintf("HINT: For gamma-point only calculations, turn off symmetries to speed up exact exchange.\n");
	
	//Initialize kmap:
	logSuspend();
	for(int iReduced=0; iReduced<qCount; iReduced++)
	for(unsigned iInvert=0; iInvert<invertList.size(); iInvert++)
	for(unsigned iSym=0; iSym<sym.size(); iSym++)
	{	KmapEntry& ki = kmap[kmapIndex(iReduced, iInvert, iSym)];
		ki.k = e.eInfo.qnums[iReduced].k * sym[iSym].rot * invertList[iInvert];
		ki.basis.setup(e.gInfo, e.iInfo, e.cntrl.Ecut, ki.k);
		if(e.eInfo.isMine(iReduced) || e.eInfo.isMine(iReduced + qCount))
			ki.transform = std::make_shared<ColumnBundleTransform>(e.eInfo.qnums[iReduced].k, e.basis[iReduced],
				ki.k, ki.basis, nSpinor, sym[iSym], invertList[iInvert]);
	}
	logResume();
}

double ExactExchangeEval::calc(int iSpin, unsigned iReduced, unsigned iInvert, unsigned iSym,
	double aXX, double omega, const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C, std::vector<ColumnBundle>* HC) const
{
	//Prepare ik state and gradient on all processes:
	const KmapEntry& ki = kmap[kmapIndex(iReduced, iInvert, iSym)];
	int ikSrc = iReduced + iSpin*qCount; //source state number
	const Basis& basis_k = ki.basis;
	QuantumNumber qnum_k = e.eInfo.qnums[ikSrc]; qnum_k.k =  ki.k;
	ColumnBundle Ck(e.eInfo.nBands, basis_k.nbasis*nSpinor, &basis_k, &qnum_k, isGpuEnabled()), HCk;
	diagMatrix Fk(e.eInfo.nBands);
	if(e.eInfo.isMine(ikSrc))
	{	Ck.zero();
		ki.transform->scatterAxpy(1., C[ikSrc], Ck,0,1);
		Fk = F[ikSrc];
	}
	mpiWorld->bcastData(Ck, e.eInfo.whose(ikSrc));
	mpiWorld->bcastData(Fk, e.eInfo.whose(ikSrc));
	if(HC) { HCk = Ck.similar(); HCk.zero(); }
	
	//Calculate energy (and gradient):
	const double prefac = -0.5*aXX / (sym.size()*invertList.size()*e.eInfo.spinWeight);
	double EXX = 0.;
	for(int bk=0; bk<e.eInfo.nBands; bk++)
	{	//Put this state in real space:
		std::vector<complexScalarField> Ipsik(nSpinor), grad_Ipsik(nSpinor);
		for(int s=0; s<nSpinor; s++)
			Ipsik[s] = I(Ck.getColumn(bk,s));
		double wFk = qnum_k.weight * Fk[bk];
		
		//Loop over states of same spin belonging to this MPI process:
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		{	const QuantumNumber& qnum_q = e.eInfo.qnums[q];
			if(qnum_k.spin != qnum_q.spin) continue;
			for(int bq=0; bq<e.eInfo.nBands; bq++)
			{	double wFq = qnum_q.weight * F[q][bq];
				if(!wFk && !wFq) continue; //at least one of the orbitals must be occupied
				
				std::vector<complexScalarField> Ipsiq(nSpinor);
				complexScalarField In; //state pair density
				for(int s=0; s<nSpinor; s++)
				{	Ipsiq[s] = I(C[q].getColumn(bq,s));
					In += conj(Ipsik[s]) * Ipsiq[s];
				}
				complexScalarFieldTilde n = J(In);
				complexScalarFieldTilde Kn = O((*e.coulomb)(n, qnum_q.k-qnum_k.k, omega)); //Electrostatic potential due to n
				EXX += (prefac*wFk*wFq) * dot(n,Kn).real();
				
				if(HC)
				{	complexScalarField E_In = Jdag(Kn);
					for(int s=0; s<nSpinor; s++)
					{	grad_Ipsik[s] += (prefac*wFq) * conj(E_In) * Ipsiq[s];
						(*HC)[q].accumColumn(bq,s, Idag((prefac*wFk) * E_In * Ipsik[s]));
					}
				}
			}
		}
		if(HC)
		{	for(int s=0; s<nSpinor; s++)
				if(grad_Ipsik[s])
					HCk.accumColumn(bk,s, Idag(grad_Ipsik[s]));
		}
	}
	mpiWorld->allReduce(EXX, MPIUtil::ReduceSum, true);
	
	//Move ik state gradient back to host process (if necessary):
	if(HC)
	{	mpiWorld->allReduceData(HCk, MPIUtil::ReduceSum);
		if(e.eInfo.isMine(ikSrc))
			ki.transform->gatherAxpy(1., HCk,0,1, (*HC)[ikSrc]);
	}
	return EXX;
}
