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
#include <electronic/ExactExchange_internal.h>
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
	ExactExchangeEval(const Everything& e, double aXX, double omega);
	
	//! Calculate exchange energy and accumulate gradients for a fixed left-orbital
	double calc_sub(int q1, int b1, std::mutex* lock,
		const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C, std::vector<ColumnBundle>* HC);
	
	//! Run calc_sub over all possible s, k1, b1 in threads with load-balancing:
	//! iDone points to an integer used for job management, the result is accumulated to EXX
	//! The targets of iDone and EXX should be set to zero before this call
	static void calc_thread(int iThread, int nThreads,
		ExactExchangeEval* eval, double* EXX, int* iDone, std::mutex* lock,
		const std::vector<diagMatrix>* F, const std::vector<ColumnBundle>* C, std::vector<ColumnBundle>* HC);
	
private:
	friend class ExactExchange;
	const Everything& e;
	double aXX; //!< Exact exchange scale factor (=1 for Hartree-Fock)
	double omegaSq; //!< Square of exact exchange range parameter (=0 for long-range i.e. unscreened)
	const std::vector< matrix3<int> >& sym; //!< symmetry matrices in lattice coordinates
	const std::vector< matrix3<int> >& symMesh; //!< symmetry matrices in mesh coordinates
	std::vector<double> weightedVzero; //!< singular correction for each state (scaled suitably for direct use in K())
	int nSpins;
	int qCount; //!< number of states of each spin
	double wSpinless; //!< factor multiplying state weights to get to spinless weights = nSpins/2
	
	//! Singular-corrected, optionally screened, coulomb operator
	//! kDiff = k2 - k1 when applying to the Bloch state product Conj(u1)*u2
	//! and q1 is the first state index (used for singular corrections, which only contribute for q1==q2)
	complexDataGptr K(const complexDataGptr& in, vector3<> kDiff, int q1) const;
};


ExactExchange::ExactExchange(const Everything& e, double aXX, double omega) : e(e)
{
	logPrintf("\n---------- Setting up exact exchange ----------\n");
	if(e.cntrl.fixed_n && !e.cntrl.fixOccupied)
		die("Exact exchange band structure calculations need fixed occupied states (command fix-occupied)\n");
	logPrintf("Exact exchange scale factor = %lf\n", aXX);
	logPrintf("Range parameter omega = %lf bohr^-1\n", omega);
	
	eval = new ExactExchangeEval(e, aXX, omega);
}

ExactExchange::~ExactExchange()
{
	delete eval;
}

double ExactExchange::operator()(const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C,
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
	
	threadOperators = false;
	threadLaunch(isGpuEnabled() ? 1 : 0, ExactExchangeEval::calc_thread, 0,
		eval, &EXX, &iDone, &lock, &F, &C, HC);
	threadOperators = true;
	
	watch.stop();
	return EXX;
}

//--------------- class ExactExchangeEval implementation ----------------------


// Exchange singular correction based on:
// [ P. Carrier, S. Rohra and A. Gorling, Phys. Rev. B 75, 205126 (2007) ]

//Reciprocal space function (with correct periodicity) singular at G=0 as 1/G^2
//g is a dimensionless coordinate (coefficients to reciprocal lattice vectors)
//  When omega!=0, the G=0 behavior is (1-exp(-omega^2 G^2/4))/G^2, which is not singular,
//  but still benefits from this correction for brillouin zones larger than omega
inline double fSingular(const vector3<>& g, const matrix3<>& GGT, double omegaSq)
{	vector3<> sinPi; for(int i=0; i<3; i++) sinPi[i] = sin(M_PI*g[i]);
	vector3<> sin2Pi; for(int i=0; i<3; i++) sin2Pi[i] = sin(2*M_PI*g[i]);
	double effGsq =  (1./(M_PI*M_PI)) * //Periodic function from ref. which is >0 and ~ G^2 near G=0
		( sinPi[0]*sinPi[0]*GGT(0,0) + sinPi[1]*sinPi[1]*GGT(1,1) + sinPi[2]*sinPi[2]*GGT(2,2)
		+ 0.5*(sin2Pi[0]*sin2Pi[1]*GGT(0,1) + sin2Pi[1]*sin2Pi[2]*GGT(1,2) + sin2Pi[2]*sin2Pi[0]*GGT(2,0)));
	return erfcTilde(effGsq, omegaSq);
}

//Integrate fSingular on a box of size scale and center gCenter (in reciprocal lattice coordinates)
double fSingularIntegralBox(double scale, const vector3<> gCenter, const matrix3<>& GGT, double omegaSq)
{	//Weights and abscissae of the 15-point gauss quadrature:
	const int N = 7;
	static const double w[N+1] =
	{	0.030753241996117268354628393577204, 0.070366047488108124709267416450667,
		0.107159220467171935011869546685869, 0.139570677926154314447804794511028,
		0.166269205816993933553200860481209, 0.186161000015562211026800561866423,
		0.198431485327111576456118326443839, 0.202578241925561272880620199967519
	};
	static const double x[N+1] =
	{	0.987992518020485428489565718586613, 0.937273392400705904307758947710209,
		0.848206583410427216200648320774217, 0.724417731360170047416186054613938,
		0.570972172608538847537226737253911, 0.394151347077563369897207370981045,
		0.201194093997434522300628303394596, 0.000000000000000000000000000000000
	};
	double ret = 0.0;
	double h = 0.5 * scale; 
	vector3<> g;
	for(int i0=-N; i0<=N; i0++)
	{	g[0] = gCenter[0] + h*x[N-abs(i0)]*(i0>0?1:-1);
		double w0 = w[N-abs(i0)];
		for(int i1=-N; i1<=N; i1++)
		{	g[1] = gCenter[1] + h*x[N-abs(i1)]*(i1>0?1:-1);
			double w01 = w0 * w[N-abs(i1)];
			for(int i2=-N; i2<=N; i2++)
			{	g[2] = gCenter[2] + h*x[N-abs(i2)]*(i2>0?1:-1);
				double w012 = w01 * w[N-abs(i2)];
				ret += w012 * fSingular(g, GGT, omegaSq);
			}
		}
	}
	return h*h*h * ret;
}

//Integrate fSingular between the box of size scale and scale/3 centered at th eorigin (in reciprocal lattice coordinates)
double fSingularIntegralBoxDiff(double scale, const matrix3<>& GGT, double omegaSq)
{	double scaleBy3 = scale/3.;
	double ret = 0.0;
	vector3<int> ig;
	for(ig[0]=-1; ig[0]<=1; ig[0]++)
	for(ig[1]=-1; ig[1]<=1; ig[1]++)
	for(ig[2]=-1; ig[2]<=1; ig[2]++)
		if(ig.length_squared()) //except the center box
			ret += fSingularIntegralBox(scaleBy3, scaleBy3*ig, GGT, omegaSq);
	return ret;
}

ExactExchangeEval::ExactExchangeEval(const Everything& e, double aXX, double omega)
: e(e), aXX(aXX), omegaSq(omega*omega),
	sym(e.symm.getMatrices()),
	symMesh(e.symm.getMeshMatrices()),
	nSpins(e.eInfo.spinType==SpinNone ? 1 : 2), 
	qCount(e.eInfo.nStates/nSpins),
	wSpinless(0.5*nSpins)
{
	//Compute the integral of fSingular over the brillouin zone
	double fSingularIntegral = 0.0;
	for(double scale=1.0; scale>1e-16; scale/=3.)
		fSingularIntegral += (4*M_PI/e.gInfo.detR) * fSingularIntegralBoxDiff(scale, e.gInfo.GGT, omegaSq);

	//Compute the integral of fSingular as sampled by the discrete k-point sum (for each state):
	weightedVzero.resize(e.eInfo.nStates);
	double VzeroSum=0., VzeroSqSum=0., VzeroDen=0.;
	for(int q1=0; q1<qCount; q1++)
	{	double w1 = wSpinless * e.eInfo.qnums[q1].weight;
		vector3<> k1 = e.eInfo.qnums[q1].k;
		double fSingularSum = 0.0;
		int countGzero = 0;
		for(int q2=0; q2<qCount; q2++)
		{	double w2mean = wSpinless * e.eInfo.qnums[q2].weight / sym.size();
			for(const matrix3<int>& m: sym)
			{	vector3<> k2 = (~m) * e.eInfo.qnums[q2].k;
				if(circDistanceSquared(k1, k2) > GzeroTOL)
					fSingularSum += (4*M_PI/e.gInfo.detR) * w2mean * fSingular(k2-k1, e.gInfo.GGT, omegaSq);
				else countGzero++;
			}
		}
		double Vzero = fSingularIntegral - fSingularSum;
		//Collect statistics
		VzeroDen += w1;
		VzeroSum += w1 * Vzero;
		VzeroSqSum += w1 * Vzero*Vzero;
		//Save with a scale factor appropriate for direct use by K()
		weightedVzero[q1] = Vzero * (sym.size()/(w1*countGzero)) * (e.gInfo.detR/(4*M_PI));
		if(q1+qCount < e.eInfo.nStates)
			weightedVzero[q1+qCount] = weightedVzero[q1];
	}
	//Print statistics of Vzero:
	double VzeroMean = VzeroSum/VzeroDen;
	logPrintf("Singular Vxx(G=0) correction = %le", VzeroMean);
	double VzeroVariance = VzeroSqSum/VzeroDen - VzeroMean*VzeroMean;
	if(VzeroVariance > 1e-15*VzeroMean*VzeroMean)
		logPrintf(" (on average, with std dev %le)", sqrt(VzeroVariance));
	logPrintf("\n");
	
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
}

double ExactExchangeEval::calc_sub(int q1, int b1, std::mutex* lock,
	const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C, std::vector<ColumnBundle>* HC)
{
	double scale = -aXX * wSpinless / sym.size(); 
	//includes negative sign for exchange, empirical scale, weight for the restricted case and symmetries

	const std::vector<QuantumNumber>& qnums = e.eInfo.qnums;
	complexDataRptr Ipsi1 = I(C[q1].getColumn(b1)), grad_Ipsi1;
	double wF1 = F[q1][b1] * qnums[q1].weight;
	vector3<> k1 = qnums[q1].k;
	double EXX = 0.;
	
	for(int q2 = q1<qCount ? 0 : qCount; q2<q1+1; q2++)
		for(int b2=0; b2<(q1==q2 ? b1+1 : e.eInfo.nBands); b2++)
		{	
			if(F[q1][b1]<=e.cntrl.occupiedThreshold && F[q2][b2]<=e.cntrl.occupiedThreshold)
				continue; //at least one of the orbitals must be occupied
			
			bool diag = q1==q2 && b1==b2; // whether this is a diagonal (self) term
			complexDataRptr Ipsi2 = diag ? Ipsi1 : I(C[q2].getColumn(b2)), grad_Ipsi2;
			double wF2 = F[q2][b2] * qnums[q2].weight;
			for(unsigned iSym=0; iSym<sym.size(); iSym++)
			{	vector3<> k2 = (~sym[iSym]) * qnums[q2].k;
				complexDataRptr RIpsi2 = pointGroupGather(Ipsi2, symMesh[iSym]);
				complexDataGptr n = J(conj(Ipsi1) * RIpsi2); //Compute state-pair 'density'
				complexDataGptr Kn = K(n, k2-k1, q1); //Electrostatic potential due to n
				EXX += (diag ? 0.5 : 1.) * scale*wF1*wF2 * dot(n,Kn).real();
				if(HC)
				{	complexDataRptr EXX_In = Jdag(Kn);
					grad_Ipsi1 += (scale*wF2) * conj(EXX_In) * RIpsi2;
					if(!diag)
						grad_Ipsi2 += pointGroupScatter((scale*wF1) * EXX_In * Ipsi1, symMesh[iSym]);
				}
			}
			if(HC && grad_Ipsi2)
			{	complexDataGptr grad_psi2 = Idag(grad_Ipsi2);
				lock->lock();
				(*HC)[q2].accumColumn(b2, grad_psi2);
				lock->unlock();
			}
		}
	if(HC && grad_Ipsi1)
	{	complexDataGptr grad_psi1 = Idag(grad_Ipsi1);
		lock->lock();
		(*HC)[q1].accumColumn(b1, grad_psi1);
		lock->unlock();
	}
	return EXX;
}

void ExactExchangeEval::calc_thread(int iThread, int nThreads,
	ExactExchangeEval* eval, double* EXX, int* iDone, std::mutex* lock,
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
		double EXX_sub = eval->calc_sub(q1, b1, lock, *F, *C, HC);
		lock->lock(); *EXX += EXX_sub; lock->unlock();
	}
}

//Compute the screened singular-corrected coulomb potential for a complex density
void screenedCoulombK_sub(int iStart, int iStop,
	const vector3<int>& S, const matrix3<>& GGT, const complex* in, complex* out,
	const vector3<>& kDiff, double weightedVzero, double omegaSq)
{	THREAD_fullGspaceLoop( out[i] = in[i] * screenedCoulombK_calc(iG, GGT, kDiff, weightedVzero, omegaSq); )
}
#ifdef GPU_ENABLED //implemented in Operators.cu
void screenedCoulombK_gpu(const vector3<int>& S, const matrix3<>& GGT, const complex* in, complex* out,
	const vector3<>& kDiff, double weightedVzero, double omegaSq);
#endif
complexDataGptr ExactExchangeEval::K(const complexDataGptr& in, vector3<> kDiff, int q1) const
{	complexDataGptr out(complexDataG::alloc(e.gInfo, isGpuEnabled()));
	out->scale = (4*M_PI*e.gInfo.detR) * in->scale;
	#ifdef GPU_ENABLED
	screenedCoulombK_gpu(
		e.gInfo.S, e.gInfo.GGT, in->dataGpu(false), out->dataGpu(false),
		kDiff, weightedVzero[q1], omegaSq);
	#else
	threadLaunch(threadOperators ? 0 : 1, screenedCoulombK_sub, e.gInfo.nr,
		e.gInfo.S, e.gInfo.GGT, in->data(false), out->data(false),
		kDiff, weightedVzero[q1], omegaSq);
	#endif
	return out;
}
