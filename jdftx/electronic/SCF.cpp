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

inline DataRptrCollection operator*(const RealKernel& K, const DataRptrCollection& x)
{	DataRptrCollection Kx(x.size());
	for(size_t i=0; i<x.size(); i++) Kx[i] = I(K * J(x[i]));
	return Kx;
}

SCF::SCF(Everything& e): e(e), kerkerMix(e.gInfo), diisMetric(e.gInfo)
{	//Set up the caching size:
	SCFparams& sp = e.scfParams;
	if(sp.vectorExtrapolation == SCFparams::VE_Plain) //No history needed for plain mixing
		sp.history = 1;
	overlap.init(sp.history, sp.history);
	
	eigenShiftInit();
	
	//Initialize the preconditioner and metric kernels:
	double Gmin = sqrt(std::min(e.gInfo.GGT(0,0), std::min(e.gInfo.GGT(1,1), e.gInfo.GGT(2,2))));
	if(sp.qKerker<0) { sp.qKerker=Gmin; logPrintf("Setting default qKerker = %lg a_0^{-1} (shortest reciprocal lattice vector)\n", sp.qKerker); }
	if(sp.qMetric<0) { sp.qMetric=Gmin; logPrintf("Setting default qMetric = %lg a_0^{-1} (shortest reciprocal lattice vector)\n", sp.qMetric); }
	applyFuncGsq(e.gInfo, setKernels, sp.mixedVariable==SCFparams::MV_Density,
		sp.mixFraction, pow(sp.qKerker,2), pow(sp.qMetric,2), kerkerMix.data, diisMetric.data);
	kerkerMix.set(); diisMetric.set();
}

void SCF::minimize()
{	
	ElecInfo& eInfo = e.eInfo;
	ElecVars& eVars = e.eVars;
	SCFparams sp = e.scfParams;
	
	logPrintf("Will mix electronic %s%s at each iteration.\n",
		(e.exCorr.needsKEdensity() ? "and kinetic " : ""),
		(sp.mixedVariable==SCFparams::MV_Density ? "density" : "potential"));
	
	if(e.exCorr.orbitalDep)
	{	e.scfParams.energyDiffThreshold = 0.;
		logPrintf("Turning off total energy convergence threshold for orbital-dependent functional.\n");
	}
	
	bool subspaceRotation=false;
	std::swap(eInfo.subspaceRotation, subspaceRotation); //Switch off subspace rotation for SCF
	
	//Compute energy for the initial guess
	if(e.eInfo.fillingsUpdate != ElecInfo::ConstantFillings) //Compute Hsub and update fillings
	{	eVars.elecEnergyAndGrad(e.ener, 0, 0, true);
		updateFillings();
	}
	double E = eVars.elecEnergyAndGrad(e.ener, 0, 0, true); mpiUtil->bcast(E); //Compute energy (and ensure consistency to machine precision)
	
	double Eprev = 0.;
	double dE = E;
	std::vector<diagMatrix> eigsPrev;
	for(int scfCounter=0; scfCounter<e.scfParams.nIterations; scfCounter++)
	{
		//If history is full, remove oldest member
		if(((int)pastResiduals.size() >= sp.history) or ((int)pastVariables.size() >= sp.history))
		{	size_t ndim = pastResiduals.size();
			if(ndim>1) overlap.set(0,ndim-1, 0,ndim-1, overlap(1,ndim, 1,ndim));
			pastVariables.erase(pastVariables.begin());
			pastResiduals.erase(pastResiduals.begin());
		}
		
		//Cache the old energy and variables
		Eprev = E;
		eigsPrev = e.eVars.Hsub_eigs;
		pastVariables.push_back(getVariable());
		
		//Band-structure minimize:
		if(not sp.verbose) { logSuspend(); e.elecMinParams.fpLog = nullLog; } // Silence eigensolver output
		e.elecMinParams.energyDiffThreshold = std::min(1e-6, std::abs(dE)/10.);
		bandMinimize(e);
		if(not sp.verbose) { logResume(); e.elecMinParams.fpLog = globalLog; }  // Resume output
		
		//Compute new density and energy
		if(e.eInfo.fillingsUpdate != ElecInfo::ConstantFillings) // Update fillings
			updateFillings();
		E = eVars.elecEnergyAndGrad(e.ener); mpiUtil->bcast(E); //ensure consistency to machine precision
		string Elabel(relevantFreeEnergyName(e));
		if(e.exCorr.orbitalDep) Elabel += "~"; //remind that the total energy is at best an estimate
		
		//Debug output:
		if(e.cntrl.shouldPrintEigsFillings) print_Hsub_eigs(e);
		if(e.cntrl.shouldPrintEcomponents) { logPrintf("\n"); e.ener.print(); logPrintf("\n"); }

		//Calculate and cache residual:
		double residualNorm = 0.;
		{	DataRptrCollection residual = getVariable() - pastVariables.back();
			pastResiduals.push_back(residual);
			residualNorm = e.gInfo.dV * sqrt(dot(residual,residual)); mpiUtil->bcast(residualNorm); //ensure consistency to machine precision
		}
		double deigs = eigDiffRMS(eigsPrev, eVars.Hsub_eigs); mpiUtil->bcast(deigs); //ensure consistency to machine precision
		dE = E-Eprev; //change in energy is also used to determine energyDiffThreshold of the bandminimizer
		logPrintf("SCF: Cycle: %2i   %s: %.15e   dE: %.3e   |deigs|: %.3e   |Residual|: %.3e\n",
			scfCounter, Elabel.c_str(), E, dE, deigs, residualNorm);
		
		//Check for convergence and update variable:
		if(fabs(E - Eprev) < sp.energyDiffThreshold) { logPrintf("SCF: Converged (|Delta E|<%le).\n\n", sp.energyDiffThreshold); break; }
		else if(deigs < sp.eigDiffThreshold)         { logPrintf("SCF: Converged (|deigs|<%le).\n\n", sp.eigDiffThreshold); break; }
		else if(residualNorm < sp.residualThreshold) { logPrintf("SCF: Converged (|Residual|<%le).\n\n", sp.residualThreshold); break; }
		else
		{	switch(sp.vectorExtrapolation)
			{	case SCFparams::VE_Plain: mixPlain(); break;
				case SCFparams::VE_DIIS: mixDIIS(); break;
			}
		}
		logFlush();
	}
	
	std::swap(eInfo.subspaceRotation, subspaceRotation); //Restore subspaceRotation to its original state
}

void SCF::mixPlain()
{	setVariable(pastVariables.back() + kerkerMix*pastResiduals.back());
}

void SCF::mixDIIS()
{	size_t ndim = pastResiduals.size();

	//Apply the metric to the latest residual
	DataRptrCollection M_lastResidual = diisMetric * pastResiduals.back();
	
	//Update the overlap matrix
	for(size_t j=0; j<ndim; j++)
	{	double thisOverlap = dot(pastResiduals[j], M_lastResidual);
		overlap.set(j, ndim-1, thisOverlap);
		overlap.set(ndim-1, j, thisOverlap);
	}
	
	//Plain mixing on the first few steps:
	if(ndim<3) { mixPlain(); return; }
	
	//Invert the residual overlap matrix to get the minimum of residual
	matrix cOverlap(ndim+1, ndim+1); // Add row and column to enforce normalization constraint
	cOverlap.set(0, ndim, 0, ndim, overlap(0, ndim, 0, ndim));
	for(size_t j=0; j<ndim; j++)
	{	cOverlap.set(j, ndim, 1);
		cOverlap.set(ndim, j, 1);
	}
	cOverlap.set(ndim, ndim, 0);
	matrix cOverlap_inv = inv(cOverlap);
	
	//Update variable:
	complex* coefs = cOverlap_inv.data();
	DataRptrCollection variable(pastVariables.back().size()); //zero
	for(size_t j=0; j<ndim; j++)
		variable += coefs[cOverlap_inv.index(j, ndim)].real() * (pastVariables[j] + kerkerMix*pastResiduals[j]);
	setVariable(variable);
}

DataRptrCollection SCF::getVariable() const
{	bool mixDensity = (e.scfParams.mixedVariable==SCFparams::MV_Density);
	DataRptrCollection variable = mixDensity ? e.eVars.n : e.eVars.Vscloc;
	if(e.exCorr.needsKEdensity()) //append the relevant KE variable:
		for(DataRptr v: (mixDensity ? e.eVars.tau : e.eVars.Vtau)) variable.push_back(v);
	return clone(variable);
}

void SCF::setVariable(DataRptrCollection variable)
{	bool mixDensity = (e.scfParams.mixedVariable==SCFparams::MV_Density);
	size_t iVar=0;
	for(DataRptr& v: (mixDensity ? e.eVars.n : e.eVars.Vscloc)) v = variable[iVar++];
	if(e.exCorr.needsKEdensity())
		for(DataRptr& v: (mixDensity ? e.eVars.tau : e.eVars.Vtau)) v = variable[iVar++];
	assert(iVar == variable.size());
	if(mixDensity) e.eVars.EdensityAndVscloc(e.ener); //Recompute Vscloc if mixing density
	e.iInfo.augmentDensityGridGrad(e.eVars.Vscloc); //update Vscloc atom projections for ultrasoft psp's 
}


void SCF::updateFillings()
{
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
	if(e.eInfo.fillingsUpdate)
	{	logPrintf("\tFillingsUpdate:  mu: %.15le  nElectrons: %.15le", mu, e.eInfo.nElectrons);
		if(e.eInfo.spinType == SpinZ)
		{	double spinPol = integral(e.eVars.n[0] - e.eVars.n[1]);
			logPrintf("  magneticMoment: %.5f", spinPol);
		}
		logPrintf("\n"); logFlush();
	}
}

double SCF::eigDiffRMS(const std::vector<diagMatrix>& eigs1, const std::vector<diagMatrix>& eigs2) const
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
