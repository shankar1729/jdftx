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
#include <electronic/operators.h>
#include <core/Util.h>
#include <core/GpuUtil.h>
#include <core/Operators.h>
#include <list>

//! Internal computation object for ExactExchange
class ExactExchangeEval
{
public:
	ExactExchangeEval(const Everything& e);
	
	//! Calculate exchange energy and accumulate gradients for a fixed left-orbital
	double calc_sub(int q1, int b1, std::mutex* lock, double aXX, double omega,
		const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C, std::vector<ColumnBundle>* HC);
	
	//! Run calc_sub over all possible s, k1, b1 in threads with load-balancing:
	//! iDone points to an integer used for job management, the result is accumulated to EXX
	//! The targets of iDone and EXX should be set to zero before this call
	static void calc_thread(int iThread, int nThreads,
		ExactExchangeEval* eval, double aXX, double omega, double* EXX, int* iDone, std::mutex* lock,
		const std::vector<diagMatrix>* F, const std::vector<ColumnBundle>* C, std::vector<ColumnBundle>* HC);
	
private:
	friend class ExactExchange;
	const Everything& e;
	const std::vector< matrix3<int> >& sym; //!< symmetry matrices in lattice coordinates
	const std::vector< matrix3<int> >& symMesh; //!< symmetry matrices in mesh coordinates
	const std::vector<int>& invertList; //!< whether to add inversion explicitly over the symmetry group
	int nSpins;
	int qCount; //!< number of states of each spin
	double wSpinless; //!< factor multiplying state weights to get to spinless weights = nSpins/2
};


ExactExchange::ExactExchange(const Everything& e) : e(e)
{
	logPrintf("\n---------- Setting up exact exchange ----------\n");
	if(mpiUtil->nProcesses()>1) die("Exact exchange not yet implemented in MPI mode.\n");
	if(e.eInfo.isNoncollinear()) die("Exact exchange not yet implemented for noncollinear spins.\n");
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
		for(int q=0; q<e.eInfo.nStates; q++)
			if(!(*HC)[q])
			{	(*HC)[q] = C[q].similar();
				(*HC)[q].zero();
			}
	double EXX = 0.0;
	int iDone = 0; //integer used for job management
	std::mutex lock; //lock used for thread synchronization
	
	threadLaunch(isGpuEnabled() ? 1 : 0, ExactExchangeEval::calc_thread, 0,
		eval, aXX, omega, &EXX, &iDone, &lock, &F, &C, HC);
	
	watch.stop();
	return EXX;
}

//--------------- class ExactExchangeEval implementation ----------------------



ExactExchangeEval::ExactExchangeEval(const Everything& e)
: e(e),
	sym(e.symm.getMatrices()),
	symMesh(e.symm.getMeshMatrices()),
	invertList(e.symm.getKpointInvertList()),
	nSpins(e.eInfo.spinType==SpinNone ? 1 : 2), 
	qCount(e.eInfo.nStates/nSpins),
	wSpinless(0.5*nSpins)
{
	//Print cost estimate to give the user some idea of how long it might take!
	double costFFT = e.eInfo.nStates * e.eInfo.nBands * 9.*e.gInfo.nr*log(e.gInfo.nr);
	double costBLAS3 = e.eInfo.nStates * pow(e.eInfo.nBands,2) * e.basis[0].nbasis;
	double costSemiLocal = 8 * costBLAS3 + 3 * costFFT; //very rough estimates of course!
	double costEXX = costFFT * sym.size()
		* ( (e.eInfo.spinType==SpinNone ? 1 : 2) * e.eInfo.nBands
			* pow(qCount,2) * 0.5/e.eInfo.nStates );
	double relativeCost = 1+costEXX/costSemiLocal;
	double relativeCostOrder = pow(10, floor(log(relativeCost)/log(10)));
	relativeCost = round(relativeCost/relativeCostOrder) * relativeCostOrder; //eliminate extra sigfigs
	logPrintf("Per-iteration cost relative to semi-local calculation ~ %lg\n", relativeCost);
	if(qCount==1 && sym.size()>1)
		logPrintf("HINT: For gamma-point only calculations, turn off symmetries to speed up exact exchange.\n");
}

double ExactExchangeEval::calc_sub(int q1, int b1, std::mutex* lock, double aXX, double omega,
	const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C, std::vector<ColumnBundle>* HC)
{
	double scale = -aXX * wSpinless / (sym.size() * invertList.size()); 
	//includes negative sign for exchange, empirical scale, weight for the restricted case and symmetries

	const std::vector<QuantumNumber>& qnums = e.eInfo.qnums;
	complexScalarField Ipsi1 = I(C[q1].getColumn(b1,0)), grad_Ipsi1;
	double wF1 = F[q1][b1] * qnums[q1].weight;
	vector3<> k1 = qnums[q1].k;
	double EXX = 0.;
	
	for(int q2 = q1<qCount ? 0 : qCount; q2<q1+1; q2++)
		for(int b2=0; b2<(q1==q2 ? b1+1 : e.eInfo.nBands); b2++)
		{	
			if(F[q1][b1]<=e.cntrl.occupiedThreshold && F[q2][b2]<=e.cntrl.occupiedThreshold)
				continue; //at least one of the orbitals must be occupied
			
			bool diag = q1==q2 && b1==b2; // whether this is a diagonal (self) term
			complexScalarField Ipsi2 = diag ? Ipsi1 : I(C[q2].getColumn(b2,0)), grad_Ipsi2;
			double wF2 = F[q2][b2] * qnums[q2].weight;
			for(int invert: invertList)
				for(unsigned iSym=0; iSym<sym.size(); iSym++)
				{	vector3<> k2 = (~sym[iSym]) * qnums[q2].k * invert;
					complexScalarField RIpsi2 = pointGroupGather(Ipsi2, symMesh[iSym]);
					if(invert==-1) RIpsi2 = conj(RIpsi2); //complex-conjugate for explicitly inverted k-point
					complexScalarFieldTilde n = J(conj(Ipsi1) * RIpsi2); //Compute state-pair 'density'
					complexScalarFieldTilde Kn = O((*e.coulomb)(n, k2-k1, omega)); //Electrostatic potential due to n
					EXX += (diag ? 0.5 : 1.) * scale*wF1*wF2 * dot(n,Kn).real();
					if(HC)
					{	complexScalarField EXX_In = Jdag(Kn);
						grad_Ipsi1 += (scale*wF2) * conj(EXX_In) * RIpsi2;
						if(!diag)
						{	complexScalarField grad_RIpsi2 = (scale*wF1) * EXX_In * Ipsi1;
							if(invert==-1) grad_RIpsi2 = conj(grad_RIpsi2); //complex-conjugate for explicitly inverted k-point
							grad_Ipsi2 += pointGroupScatter(grad_RIpsi2, symMesh[iSym]);
						}
					}
				}
			if(HC && grad_Ipsi2)
			{	complexScalarFieldTilde grad_psi2 = Idag(grad_Ipsi2);
				lock->lock();
				(*HC)[q2].accumColumn(b2,0, grad_psi2);
				lock->unlock();
			}
		}
	if(HC && grad_Ipsi1)
	{	complexScalarFieldTilde grad_psi1 = Idag(grad_Ipsi1);
		lock->lock();
		(*HC)[q1].accumColumn(b1,0, grad_psi1);
		lock->unlock();
	}
	return EXX;
}

void ExactExchangeEval::calc_thread(int iThread, int nThreads,
	ExactExchangeEval* eval, double aXX, double omega, double* EXX, int* iDone, std::mutex* lock,
	const std::vector<diagMatrix>* F, const std::vector<ColumnBundle>* C, std::vector<ColumnBundle>* HC)
{
	const Everything& e = eval->e;
	int nJobs = e.eInfo.nStates * e.eInfo.nBands;
	
	while(true)
	{	//Each thread picks next available job when done with previous one
		lock->lock();
		int iJob = (*iDone)++;
		lock->unlock();
		if(iJob >= nJobs) break; //no jobs left
		//Process current job:
		int q1 = iJob / e.eInfo.nBands;
		int b1 = iJob - q1*e.eInfo.nBands;
		double EXX_sub = eval->calc_sub(q1, b1, lock, aXX, omega, *F, *C, HC);
		lock->lock(); *EXX += EXX_sub; lock->unlock();
	}
}
