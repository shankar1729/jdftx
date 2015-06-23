/*-------------------------------------------------------------------
Copyright 2013 Deniz Gunceler, Ravishankar Sundararaman

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

#include <electronic/SCF.h>
#include <electronic/ElecMinimizer.h>
#include <core/ScalarFieldIO.h>
#include <queue>

inline void setKernels(int i, double Gsq, bool mixDensity, double mixFraction, double qKerkerSq, double qMetricSq, double* kerkerMix, double* diisMetric)
{	if(mixDensity)
	{	kerkerMix[i] = mixFraction * (qKerkerSq ? Gsq/(Gsq + qKerkerSq) : 1.);
		diisMetric[i] = Gsq ? (Gsq + qMetricSq)/Gsq : 0.;
	}
	else
	{	kerkerMix[i] = mixFraction;
		diisMetric[i] = qMetricSq ? Gsq / (qMetricSq + Gsq) : 1.;
	}
}

inline ScalarFieldArray operator*(const RealKernel& K, const ScalarFieldArray& x)
{	ScalarFieldArray Kx(x.size());
	for(size_t i=0; i<x.size(); i++) Kx[i] = I(K * J(x[i]));
	return Kx;
}

SCF::SCF(Everything& e): Pulay<SCFvariable>(e.scfParams), e(e), kerkerMix(e.gInfo), diisMetric(e.gInfo)
{	SCFparams& sp = e.scfParams;
	eigenShiftInit();
	mixTau = e.exCorr.needsKEdensity();
	
	//Initialize the preconditioner and metric kernels:
	applyFuncGsq(e.gInfo, setKernels, sp.mixedVariable==SCFparams::MV_Density,
		sp.mixFraction, pow(sp.qKerker,2), pow(sp.qMetric,2), kerkerMix.data, diisMetric.data);
	kerkerMix.set(); diisMetric.set();
	
	//Load history if available:
	if(sp.historyFilename.length())
	{	loadState(sp.historyFilename.c_str());
		sp.historyFilename.clear(); //make sure it doesn't get loaded again on subsequent SCFs (eg. in an ionic loop)
	}
}

void SCF::minimize()
{	
	ElecInfo& eInfo = e.eInfo;
	ElecVars& eVars = e.eVars;
	SCFparams& sp = e.scfParams;
	
	logPrintf("Will mix electronic %s%s at each iteration.\n",
		(mixTau ? "and kinetic " : ""),
		(sp.mixedVariable==SCFparams::MV_Density ? "density" : "potential"));
	
	string Elabel = relevantFreeEnergyName(e);
	if(!e.exCorr.hasEnergy())
	{	e.scfParams.energyDiffThreshold = 0.;
		logPrintf("Turning off total energy convergence threshold for potential-only functionals.\n");
		Elabel += "~"; //remind that the total energy is at best an estimate
	}
	sp.energyLabel = Elabel.c_str();
	sp.linePrefix = "SCF: ";
	sp.energyFormat = "%+.15lf";
	sp.fpLog = globalLog;
	
	bool subspaceRotation=false;
	std::swap(eInfo.subspaceRotation, subspaceRotation); //Switch off subspace rotation for SCF
	
	//Backup electronic minimize params that are modified below:
	double eMinThreshold = e.elecMinParams.energyDiffThreshold;
	int eMinIterations = e.elecMinParams.nIterations;

	//Compute energy for the initial guess
	double E = eVars.elecEnergyAndGrad(e.ener, 0, 0, true); mpiUtil->bcast(E); //Compute energy (and ensure consistency to machine precision)
	
	//Optimize using Pulay mixer:
	std::vector<string> extraNames(1, "deigs");
	std::vector<double> extraThresh(1, sp.eigDiffThreshold);
	Pulay<SCFvariable>::minimize(E, extraNames, extraThresh);
	e.iInfo.augmentDensityGridGrad(e.eVars.Vscloc); //to make sure grid projections are compatible with final Vscloc
	
	std::swap(eInfo.subspaceRotation, subspaceRotation); //Restore subspaceRotation to its original state

	//Restore electronic minimize params that were modified above:
	e.elecMinParams.energyDiffThreshold = eMinThreshold;
	e.elecMinParams.nIterations = eMinIterations;
	
	//Update gradient of energy w.r.t overlap matrix (important for ultrasoft forces)
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		eVars.grad_CdagOC[q] = -(eVars.Hsub_eigs[q] * eVars.F[q]);
	
	//Set auxiliary Hamiltonian to subspace Hamiltonian
	// (for general compatibility with minimizer; also important for lattice gradient in metals)
	if(e.eInfo.fillingsUpdate == ElecInfo::FermiFillingsAux)
	{	eVars.B = eVars.Hsub;
		eVars.B_eigs = eVars.Hsub_eigs;
		eVars.B_evecs = eVars.Hsub_evecs;
	}
}

double SCF::sync(double x) const
{	mpiUtil->bcast(x);
	return x;
}

double SCF::cycle(double dEprev, std::vector<double>& extraValues)
{	const SCFparams& sp = e.scfParams;
	
	//Cache required quantities:
	std::vector<diagMatrix> eigsPrev = e.eVars.Hsub_eigs;
	if(e.eInfo.fillingsUpdate == ElecInfo::MaximumOverlapMethod)
		Cold = e.eVars.C;
	
	//Band-structure minimize:
	if(not sp.verbose) { logSuspend(); e.elecMinParams.fpLog = nullLog; } // Silence eigensolver output
	e.elecMinParams.energyDiffThreshold = std::min(1e-6, 0.1*fabs(dEprev));
	if(sp.nEigSteps) e.elecMinParams.nIterations = sp.nEigSteps;
	bandMinimize(e);
	if(not sp.verbose) { logResume(); e.elecMinParams.fpLog = globalLog; }  // Resume output

	//Compute new density and energy
	e.ener.Eband = 0.; //only affects printing (if non-zero Energies::print assumes band structure calc)
	if(e.eInfo.fillingsUpdate != ElecInfo::ConstantFillings) updateFillings(); // Update fillings
		
	double E = e.eVars.elecEnergyAndGrad(e.ener); mpiUtil->bcast(E); //ensure consistency to machine precision

	extraValues[0] = eigDiffRMS(eigsPrev, e.eVars.Hsub_eigs);
	return E;
}

void SCF::report(int iter)
{
	if(e.cntrl.shouldPrintEigsFillings) print_Hsub_eigs(e);
	if(e.cntrl.shouldPrintEcomponents) { logPrintf("\n"); e.ener.print(); logPrintf("\n"); }
	logFlush();

	e.dump(DumpFreq_Electronic, iter);
	//--- write SCF history if dumping state:
	if(e.dump.count(std::make_pair(DumpFreq_Electronic,DumpState)) && e.dump.checkInterval(DumpFreq_Electronic,iter))
	{	string fname = e.dump.getFilename("scfHistory");
		logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
		saveState(fname.c_str());
		logPrintf("done\n"); logFlush();
	}
}


void SCF::axpy(double alpha, const SCFvariable& X, SCFvariable& Y) const
{	//Density:
	Y.n.resize(e.eVars.n.size());
	::axpy(alpha, X.n, Y.n);
	//KE density:
	if(mixTau)
	{	Y.tau.resize(e.eVars.n.size());
		::axpy(alpha, X.tau, Y.tau);
	}
	//Atomic density matrices:
	if(e.eInfo.hasU)
	{	if(!Y.rhoAtom.size()) e.iInfo.rhoAtom_initZero(Y.rhoAtom);
		for(size_t i=0; i<X.rhoAtom.size(); i++)
			::axpy(alpha, X.rhoAtom[i], Y.rhoAtom[i]);
	}
}

double SCF::dot(const SCFvariable& X, const SCFvariable& Y) const
{	double ret = 0.;
	//Density:
	ret += e.gInfo.dV * ::dot(X.n, Y.n);
	//KE density:
	if(mixTau)
		ret += e.gInfo.dV * ::dot(X.tau, Y.tau);
	//Atomic density matrices:
	if(e.eInfo.hasU)
	{	for(size_t i=0; i<X.rhoAtom.size(); i++)
			ret += dotc(X.rhoAtom[i],Y.rhoAtom[i]).real();
	}
	return ret;
}

size_t SCF::variableSize() const
{	size_t nDoubles = e.gInfo.nr * e.eVars.n.size() * (mixTau ? 2 : 1); //n and optionally tau
	if(e.eInfo.hasU)
	{	std::vector<matrix> rhoAtom;
		e.iInfo.rhoAtom_initZero(rhoAtom);
		for(const matrix& m: rhoAtom) nDoubles += 2*m.nData();
	}
	return nDoubles * sizeof(double);
}

void SCF::readVariable(SCFvariable& v, FILE* fp) const
{	//Density:
	nullToZero(v.n, e.gInfo, e.eVars.n.size());
	for(ScalarField& X: v.n) loadRawBinary(X, fp);
	//KE density:
	if(mixTau)
	{	nullToZero(v.tau, e.gInfo, e.eVars.n.size());
		for(ScalarField& X: v.tau) loadRawBinary(X, fp);
	}
	//Atomic density matrices:
	if(e.eInfo.hasU)
	{	e.iInfo.rhoAtom_initZero(v.rhoAtom);
		for(matrix& m: v.rhoAtom) m.read(fp);
	}
}

void SCF::writeVariable(const SCFvariable& v, FILE* fp) const
{	//Density:
	for(const ScalarField& X: v.n) saveRawBinary(X, fp);
	//KE density:
	if(mixTau)
	{	for(const ScalarField& X: v.tau) saveRawBinary(X, fp);
	}
	//Atomic density matrices:
	if(e.eInfo.hasU)
	{	for(const matrix& m: v.rhoAtom) m.write(fp);
	}
}

namespace Magnetization
{	//Conversions between spin-density(-matrix) and magnetization
	ScalarFieldArray fromSpinDensity(const ScalarFieldArray& n)
	{	ScalarFieldArray m(n.size());
		if(n.size()==1)
		{	m[0] = clone(n[0]);
		}
		else
		{	m[0] = n[0] + n[1]; //total density
			m[1] = n[0] - n[1]; //z-magnetization
			if(n.size()==4)
			{	m[2] = (+2)*n[2]; //x-magnetization
				m[3] = (-2)*n[3]; //y-magnetization
			}
		}
		return m;
	}
	ScalarFieldArray toSpinDensity(const ScalarFieldArray& m)
	{	ScalarFieldArray n(m.size());
		if(m.size()==1)
		{	n[0] = clone(m[0]);
		}
		else
		{	n[0] = 0.5*(m[0] + m[1]); //up-up
			n[1] = 0.5*(m[0] - m[1]); //dn-dn
			if(m.size()==4)
			{	n[2] = (+0.5)*m[2]; //Re(up-dn)
				n[3] = (-0.5)*m[3]; //Im(up-dn)
			}
		}
		return n;
	}
}

SCFvariable SCF::getVariable() const
{	bool mixDensity = (e.scfParams.mixedVariable==SCFparams::MV_Density);
	SCFvariable v;
	//Density:
	v.n = Magnetization::fromSpinDensity(mixDensity ? e.eVars.n : e.eVars.Vscloc);
	//KE density:
	if(mixTau)
		v.tau = Magnetization::fromSpinDensity(mixDensity ? e.eVars.tau : e.eVars.Vtau);
	//Atomic density matrices:
	if(e.eInfo.hasU)
		v.rhoAtom = (mixDensity ? e.eVars.rhoAtom : e.eVars.U_rhoAtom);
	return v;
}

void SCF::setVariable(const SCFvariable& v)
{	bool mixDensity = (e.scfParams.mixedVariable==SCFparams::MV_Density);
	//Density:
	(mixDensity ? e.eVars.n : e.eVars.Vscloc) = Magnetization::toSpinDensity(v.n);
	//KE density:
	if(mixTau)
		(mixDensity ? e.eVars.tau : e.eVars.Vtau) = Magnetization::toSpinDensity(v.tau);
	//Atomic density matrices:
	if(e.eInfo.hasU)
		(mixDensity ? e.eVars.rhoAtom : e.eVars.U_rhoAtom) = v.rhoAtom;
	//Update precomputed quantities of one-particle Hamiltonian:
	if(mixDensity) e.eVars.EdensityAndVscloc(e.ener); //Recompute Vscloc (Vtau) if mixing density
	e.iInfo.augmentDensityGridGrad(e.eVars.Vscloc); //update Vscloc atom projections for ultrasoft psp's 
}

SCFvariable SCF::precondition(const SCFvariable& v) const
{	SCFvariable vOut;
	double magEnhance = e.scfParams.mixFractionMag / e.scfParams.mixFraction;
	//Density:
	vOut.n = kerkerMix * v.n;
	for(size_t s=1; s<vOut.n.size(); s++)
		vOut.n[s] *= magEnhance;
	//KE density:
	if(mixTau)
	{	vOut.tau = kerkerMix * v.tau;
		for(size_t s=1; s<vOut.tau.size(); s++)
			vOut.tau[s] *= magEnhance;
	}
	//Atomic density matrices:
	if(e.eInfo.hasU)
	{	vOut.rhoAtom = v.rhoAtom;
		for(matrix& m: vOut.rhoAtom)
			m *= e.scfParams.mixFraction;
	}
	return vOut;
}

SCFvariable SCF::applyMetric(const SCFvariable& v) const
{	SCFvariable vOut;
	//Density:
	vOut.n = diisMetric * v.n;
	//KE density:
	if(mixTau)
		vOut.tau = diisMetric * v.tau;
	//Atomic density matrices:
	if(e.eInfo.hasU)
		vOut.rhoAtom = v.rhoAtom;
	return vOut;
}

void SCF::updateFillings()
{
	// The maximum-overlap-method is used, if enabled
	if(e.eInfo.fillingsUpdate == ElecInfo::MaximumOverlapMethod) 
	{	updateMOM(); 
		return;
	}
  
	ElecInfo& eInfo = e.eInfo; 
	ElecVars& eVars = e.eVars;
	
	//Update nElectrons from mu, or mu from nElectrons as appropriate:
	eigenShiftApply(false);
	double mu; // Electron chemical potential
	if(std::isnan(eInfo.mu)) mu = eInfo.findMu(eVars.Hsub_eigs, eInfo.nElectrons);
	else
	{	mu = eInfo.mu; 
		((ElecInfo&)eInfo).nElectrons = eInfo.nElectronsFermi(mu, eVars.Hsub_eigs); 
	}
	//Compute fillings from aux hamiltonian eigenvalues:
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		eVars.F[q] = eInfo.fermi(mu, eVars.Hsub_eigs[q]);
	//Update TS and muN:
	eInfo.updateFillingsEnergies(e.eVars.F, e.ener);
	eigenShiftApply(true);
	
	// Print filling update information
	eInfo.printFermi("Update", &mu);
}

void SCF::updateMOM()
{
	std::vector<matrix> overlap_matrix(e.eInfo.nStates);
	for(size_t s=0; s<overlap_matrix.size(); s++)
	{	
		// Overlap of the old and new wavefunctions
		overlap_matrix[s] = Cold[s]^O(e.eVars.C[s]);
		// Total no electrons in this channel
		int nElec = 0.;
		for(int j=0; j<e.eInfo.nBands; j++) nElec +=e.eVars.F[s][j];
		// For each new wavefunction, calculate the total overlap with the span of the old ones
		std::vector<double> span_overlap(e.eInfo.nBands);		
		for(int i=0; i<e.eInfo.nBands; i++)  // Loop over new wavefunctions
		{	span_overlap[i] = 0.;
			for(int j=0; j<e.eInfo.nBands; j++)  // Loop over old wavefunctions
				span_overlap[i] += e.eVars.F[s][j] * abs(overlap_matrix[s](i,j));
			//logPrintf("%.2f ", span_overlap[i]);
		}
		//logPrintf("\n");
		
		// Construct the new fillings matrix
		std::priority_queue<std::pair<double, int>> q;
		for (size_t i = 0; i < span_overlap.size(); ++i) {
			q.push(std::pair<double, int>(span_overlap[i], i));
		}
		diagMatrix Fnew = 0.*e.eVars.F[s];
		for(int k=0; k<nElec; k++)
		{	int index = q.top().second;
			Fnew[index] = 1.;
			q.pop();
		}
		//Fnew.print(globalLog);
		
		e.eVars.F[s] = Fnew;
	}
	logPrintf("\tRecomputed orbital fillings using the maximum overlap method.\n");

}

double SCF::eigDiffRMS(const std::vector<diagMatrix>& eigs1, const std::vector<diagMatrix>& eigs2, const Everything& e)
{	double rmsNum=0., rmsDen=0.;
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{	double wq = e.eInfo.qnums[q].weight;
		for(int b=0; b<e.eInfo.nBands; b++)
		{	double de = eigs1[q][b] - eigs2[q][b];
			rmsNum += wq * de*de;
			rmsDen += wq;
		}
	}
	mpiUtil->allReduce(rmsNum, MPIUtil::ReduceSum);
	mpiUtil->allReduce(rmsDen, MPIUtil::ReduceSum);
	return sqrt(rmsNum/rmsDen);
}

double SCF::eigDiffRMS(const std::vector<diagMatrix>& eigs1, const std::vector<diagMatrix>& eigs2) const
{	return eigDiffRMS(eigs1, eigs2, e);
}

void SCF::eigenShiftInit()
{	for(SCFparams::EigenShift& es: e.scfParams.eigenShifts)
	{	// Correct for the 0
		if(es.fromHOMO)
		{	es.n += e.eInfo.findHOMO(es.q);
			es.fromHOMO = false; // Needed if SCF is called multiple times, e.g. from an ionic minimize
		}
		// Check for a meaningful q and n
		if(es.q < 0) die("Eigenshift quantum number (q) must be greater than 0!\n");
		if(es.q > e.eInfo.nStates) die("Eigenshift quantum number (q) must be less than nStates=%i!\n", e.eInfo.nStates);
		if(es.n < 0) die("Eigenshift band index (n) must be greater than 0!\n");
		if(es.n > e.eInfo.nBands) die("Eigenshift band index (n) must be less than nBands=%i!\n", e.eInfo.nBands);		
	}
}

void SCF::eigenShiftApply(bool reverse)
{	int sign = reverse ? -1 : +1;
	for(const SCFparams::EigenShift& es: e.scfParams.eigenShifts)
		if(e.eInfo.isMine(es.q))
			e.eVars.Hsub_eigs[es.q][es.n] += sign * es.shift;
}

#define cutoff 1e-6
void varInverse_kernel(int i, vector3<> r, double* varInverse, double* var)
{	varInverse[i] = (var[i]>cutoff ? 1./var[i] : 0.);
}

void spness_kernel(int i, vector3<> r, double* tauW, double* tau, double* spness)
{	spness[i] = (tau[i]>cutoff ? tauW[i]/tau[i] : 1.);
}
void tauW_kernel(int i, vector3<> r, double* tauW, double* grad_n_sq, double* n)
{	tauW[i] = (n[i]>cutoff ? 0.125*grad_n_sq[i]/n[i] : 0.);
}

void tauW_n_kernel(int i, vector3<> r, double* tauW_n, double* grad_n_sq, double* n, complex* lap_n)
{	tauW_n[i] = (n[i]>cutoff ? 0.125*grad_n_sq[i]/pow(n[i],2) - 0.25*real(lap_n[i])/n[i] : 0.);
}

void Vtau_kernel(int i, vector3<> r, double* Vtau, double* fz_z, double* tauW, double* tau, double* Ex, double* Eh, double* n, double N)
{	Vtau[i] = (n[i]>cutoff ? fz_z[i]*tauW[i]*(Ex[i] + Eh[i]/N)/pow(tau[i], 2) : 0.);
}