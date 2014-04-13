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
#include <core/DataIO.h>

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

SCF::SCF(Everything& e): e(e), kerkerMix(e.gInfo), diisMetric(e.gInfo), skipInitialFillings(false)
{	SCFparams& sp = e.scfParams;
	overlap.init(sp.history, sp.history);
	eigenShiftInit();
	
	//Initialize the preconditioner and metric kernels:
	applyFuncGsq(e.gInfo, setKernels, sp.mixedVariable==SCFparams::MV_Density,
		sp.mixFraction, pow(sp.qKerker,2), pow(sp.qMetric,2), kerkerMix.data, diisMetric.data);
	kerkerMix.set(); diisMetric.set();
	
	//Load history if available:
	if(sp.historyFilename.length())
	{	size_t nBytesCycle = 2 * variableSize(); //number of bytes per history entry
		size_t nBytesFile = fileSize(sp.historyFilename.c_str());
		size_t ndim = nBytesFile / nBytesCycle;
		size_t dimOffset = 0;
		if(int(ndim) > sp.history)
		{	dimOffset = ndim - sp.history;
			ndim = sp.history;
		}
		if(nBytesFile % nBytesCycle != 0)
			die("SCF history file '%s' does not contain an integral multiple of the mixed variables and residuals.\n", sp.historyFilename.c_str());
		logPrintf("Reading %lu past variables and residuals from '%s' ... ", ndim, sp.historyFilename.c_str()); logFlush();
		pastVariables.resize(ndim);
		pastResiduals.resize(ndim);
		FILE* fp = fopen(sp.historyFilename.c_str(), "r");
		if(dimOffset) fseek(fp, dimOffset*nBytesCycle, SEEK_SET);
		for(size_t idim=0; idim<ndim; idim++)
		{	readVariable(pastVariables[idim], fp);
			readVariable(pastResiduals[idim], fp);
		}
		fclose(fp);
		logPrintf("done.\n"); logFlush();
		sp.historyFilename.clear(); //make sure it doesn't get loaded again on subsequent SCFs (eg. in an ionic loop)
		//Compute overlaps of loaded history:
		for(size_t i=0; i<ndim; i++)
		{	Variable Mresidual_i = applyMetric(pastResiduals[i]);
			for(size_t j=0; j<=i; j++)
			{	double thisOverlap = dot(pastResiduals[j], Mresidual_i);
				overlap.set(i,j, thisOverlap);
				overlap.set(j,i, thisOverlap);
			}
		}
		skipInitialFillings = true;
	}
}

void SCF::minimize()
{	
	ElecInfo& eInfo = e.eInfo;
	ElecVars& eVars = e.eVars;
	SCFparams& sp = e.scfParams;
	
	logPrintf("Will mix electronic %s%s at each iteration.\n",
		(e.exCorr.needsKEdensity() ? "and kinetic " : ""),
		(sp.mixedVariable==SCFparams::MV_Density ? "density" : "potential"));
	
	if(!e.exCorr.hasEnergy())
	{	e.scfParams.energyDiffThreshold = 0.;
		logPrintf("Turning off total energy convergence threshold for potential-only functionals.\n");
	}
	
	bool subspaceRotation=false;
	std::swap(eInfo.subspaceRotation, subspaceRotation); //Switch off subspace rotation for SCF
	
	//Compute energy for the initial guess
	if(e.eInfo.fillingsUpdate!=ElecInfo::ConstantFillings && !skipInitialFillings && !eVars.Hsub[eInfo.qStart]) //Compute Hsub and update fillings
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
		if(sp.nEigSteps) e.elecMinParams.nIterations = sp.nEigSteps;
		bandMinimize(e);
		if(not sp.verbose) { logResume(); e.elecMinParams.fpLog = globalLog; }  // Resume output
		
		//Compute new density and energy
		e.ener.Eband = 0.; //only affects printing (if non-zero Energies::print assumes band structure calc)
		if(e.eInfo.fillingsUpdate != ElecInfo::ConstantFillings) // Update fillings
			updateFillings();
		E = eVars.elecEnergyAndGrad(e.ener); mpiUtil->bcast(E); //ensure consistency to machine precision
		if(e.scfParams.sp_constraint)
			single_particle_constraint(e.scfParams.sp_constraint); //ensure the single particle limit of exchange 
		string Elabel(relevantFreeEnergyName(e));
		if(!e.exCorr.hasEnergy()) Elabel += "~"; //remind that the total energy is at best an estimate
		
		//Debug output:
		if(e.cntrl.shouldPrintEigsFillings) print_Hsub_eigs(e);
		if(e.cntrl.shouldPrintEcomponents) { logPrintf("\n"); e.ener.print(); logPrintf("\n"); }
		
		//Calculate and cache residual:
		double residualNorm = 0.;
		{	Variable residual = getVariable(); axpy(-1., pastVariables.back(), residual);
			pastResiduals.push_back(residual);
			residualNorm = sqrt(dot(residual,residual)); mpiUtil->bcast(residualNorm); //ensure consistency to machine precision
		}
		double deigs = eigDiffRMS(eigsPrev, eVars.Hsub_eigs); mpiUtil->bcast(deigs); //ensure consistency to machine precision
		dE = E-Eprev; //change in energy is also used to determine energyDiffThreshold of the bandminimizer
		logPrintf("SCF: Cycle: %2i   %s: %.15e   dE: %.3e   |deigs|: %.3e   |Residual|: %.3e\n",
			scfCounter, Elabel.c_str(), E, dE, deigs, residualNorm);
		
		//Per-cycle dumps if any (must do before the next mix step / convergence exit):
		e.dump(DumpFreq_Gummel, scfCounter);
		//--- write SCF history if dumping state:
		if(e.dump.count(std::make_pair(DumpFreq_Gummel,DumpState)) && e.dump.checkInterval(DumpFreq_Gummel,scfCounter))
		{	string fname = e.dump.getFilename("scfHistory");
			logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
			if(mpiUtil->isHead())
			{	FILE* fp = fopen(fname.c_str(), "w");
				for(size_t idim=0; idim<pastVariables.size(); idim++)
				{	writeVariable(pastVariables[idim], fp);
					writeVariable(pastResiduals[idim], fp);
				}
				fclose(fp);
			}
			logPrintf("done\n"); logFlush();
		}
		
		//Check for convergence and update variable:
		if(fabs(E - Eprev) < sp.energyDiffThreshold) { logPrintf("SCF: Converged (|Delta E|<%le).\n\n", sp.energyDiffThreshold); break; }
		else if(deigs < sp.eigDiffThreshold)         { logPrintf("SCF: Converged (|deigs|<%le).\n\n", sp.eigDiffThreshold); break; }
		else if(residualNorm < sp.residualThreshold) { logPrintf("SCF: Converged (|Residual|<%le).\n\n", sp.residualThreshold); break; }
		else mixDIIS();
		logFlush();
	}
	
	std::swap(eInfo.subspaceRotation, subspaceRotation); //Restore subspaceRotation to its original state
	
	//Update gradient of energy w.r.t overlap matrix (important for ultrasoft forces)
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		eVars.grad_CdagOC[q] = -(eVars.Hsub_eigs[q] * eVars.F[q]);
}

void SCF::axpy(double alpha, const SCF::Variable& X, SCF::Variable& Y) const
{	//Density:
	Y.n.resize(e.eVars.n.size());
	::axpy(alpha, X.n, Y.n);
	//KE density:
	if(e.exCorr.needsKEdensity())
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

double SCF::dot(const SCF::Variable& X, const SCF::Variable& Y) const
{	double ret = 0.;
	//Density:
	ret += e.gInfo.dV * ::dot(X.n, Y.n);
	//KE density:
	if(e.exCorr.needsKEdensity())
		ret += e.gInfo.dV * ::dot(X.tau, Y.tau);
	//Atomic density matrices:
	if(e.eInfo.hasU)
	{	for(size_t i=0; i<X.rhoAtom.size(); i++)
			ret += dotc(X.rhoAtom[i],Y.rhoAtom[i]).real();
	}
	return ret;
}

size_t SCF::variableSize() const
{	size_t nDoubles = e.gInfo.nr * e.eVars.n.size() * (e.exCorr.needsKEdensity() ? 2 : 1); //n and optionally tau
	if(e.eInfo.hasU)
	{	std::vector<matrix> rhoAtom;
		e.iInfo.rhoAtom_initZero(rhoAtom);
		for(const matrix& m: rhoAtom) nDoubles += 2*m.nData();
	}
	return nDoubles * sizeof(double);
}

void SCF::readVariable(SCF::Variable& v, FILE* fp) const
{	//Density:
	nullToZero(v.n, e.gInfo, e.eVars.n.size());
	for(DataRptr& X: v.n) loadRawBinary(X, fp);
	//KE density:
	if(e.exCorr.needsKEdensity())
	{	nullToZero(v.tau, e.gInfo, e.eVars.n.size());
		for(DataRptr& X: v.tau) loadRawBinary(X, fp);
	}
	//Atomic density matrices:
	if(e.eInfo.hasU)
	{	e.iInfo.rhoAtom_initZero(v.rhoAtom);
		for(matrix& m: v.rhoAtom) m.read(fp);
	}
}

void SCF::writeVariable(const SCF::Variable& v, FILE* fp) const
{	//Density:
	for(const DataRptr& X: v.n) saveRawBinary(X, fp);
	//KE density:
	if(e.exCorr.needsKEdensity())
	{	for(const DataRptr& X: v.tau) saveRawBinary(X, fp);
	}
	//Atomic density matrices:
	if(e.eInfo.hasU)
	{	for(const matrix& m: v.rhoAtom) m.write(fp);
	}
}

SCF::Variable SCF::getVariable() const
{	bool mixDensity = (e.scfParams.mixedVariable==SCFparams::MV_Density);
	Variable v;
	//Density:
	v.n = clone(mixDensity ? e.eVars.n : e.eVars.Vscloc);
	//KE density:
	if(e.exCorr.needsKEdensity())
		v.tau = clone(mixDensity ? e.eVars.tau : e.eVars.Vtau);
	//Atomic density matrices:
	if(e.eInfo.hasU)
		v.rhoAtom = (mixDensity ? e.eVars.rhoAtom : e.eVars.U_rhoAtom);
	return v;
}

void SCF::setVariable(const SCF::Variable& v)
{	bool mixDensity = (e.scfParams.mixedVariable==SCFparams::MV_Density);
	//Density:
	(mixDensity ? e.eVars.n : e.eVars.Vscloc) = v.n;
	//KE density:
	if(e.exCorr.needsKEdensity())
		(mixDensity ? e.eVars.tau : e.eVars.Vtau) = v.tau;
	//Atomic density matrices:
	if(e.eInfo.hasU)
		(mixDensity ? e.eVars.rhoAtom : e.eVars.U_rhoAtom) = v.rhoAtom;
	//Update precomputed quantities of one-particle Hamiltonian:
	if(mixDensity) e.eVars.EdensityAndVscloc(e.ener); //Recompute Vscloc (Vtau) if mixing density
	e.iInfo.augmentDensityGridGrad(e.eVars.Vscloc); //update Vscloc atom projections for ultrasoft psp's 
}

SCF::Variable SCF::applyKerker(const SCF::Variable& v) const
{	Variable vOut;
	//Density:
	vOut.n = kerkerMix * v.n;
	//KE density:
	if(e.exCorr.needsKEdensity())
		vOut.tau = kerkerMix * v.tau;
	//Atomic density matrices:
	if(e.eInfo.hasU)
	{	vOut.rhoAtom = v.rhoAtom;
		for(matrix& m: vOut.rhoAtom)
			m *= e.scfParams.mixFraction;
	}
	return vOut;
}

SCF::Variable SCF::applyMetric(const SCF::Variable& v) const
{	Variable vOut;
	//Density:
	vOut.n = diisMetric * v.n;
	//KE density:
	if(e.exCorr.needsKEdensity())
		vOut.tau = diisMetric * v.tau;
	//Atomic density matrices:
	if(e.eInfo.hasU)
		vOut.rhoAtom = v.rhoAtom;
	return vOut;
}


void SCF::mixDIIS()
{	//Update the overlap matrix
	size_t ndim = pastResiduals.size();
	Variable MlastResidual = applyMetric(pastResiduals.back());
	for(size_t j=0; j<ndim; j++)
	{	double thisOverlap = dot(pastResiduals[j], MlastResidual);
		overlap.set(j, ndim-1, thisOverlap);
		overlap.set(ndim-1, j, thisOverlap);
	}
	
	//Invert the residual overlap matrix to get the minimum of residual
	matrix cOverlap(ndim+1, ndim+1); //Add row and column to enforce normalization constraint
	cOverlap.set(0, ndim, 0, ndim, overlap(0, ndim, 0, ndim));
	for(size_t j=0; j<ndim; j++)
	{	cOverlap.set(j, ndim, 1);
		cOverlap.set(ndim, j, 1);
	}
	cOverlap.set(ndim, ndim, 0);
	matrix cOverlap_inv = inv(cOverlap);
	
	//Update variable:
	complex* coefs = cOverlap_inv.data();
	Variable v;
	for(size_t j=0; j<ndim; j++)
	{	double alpha = coefs[cOverlap_inv.index(j, ndim)].real();
		axpy(alpha, pastVariables[j], v);
		axpy(alpha, applyKerker(pastResiduals[j]), v);
	}
	setVariable(v);
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

void SCF::single_particle_constraint(double sp_constraint)
{	
	static StopWatch watch("single particle constraint"); watch.start(); 
	
	const auto& n = e.eVars.n;
	const auto& tau = e.eVars.KEdensity(); //KE density
	
	assert(n.size() == 2);
	
	// Calculate 1/n and 1/tau
	DataRptrCollection nInverse(n.size()), tauInverse(n.size());
	std::vector<double> N(n.size()), NInverse(n.size());
	for(size_t j=0; j<n.size(); j++)
	{	nullToZero(nInverse[j], e.gInfo); applyFunc_r(e.gInfo, varInverse_kernel, nInverse[j]->data(), n[j]->data());
		nullToZero(tauInverse[j], e.gInfo); applyFunc_r(e.gInfo, varInverse_kernel, tauInverse[j]->data(), tau[j]->data());
		N[j] = integral(n[j]);
		NInverse[j] = (N[j]>cutoff ? 1./N[j] : 0.);
	}
	
	// Single particle KE density and the single-particleness
	DataRptrCollection tauW(n.size());
	DataRptrCollection spness(n.size());
	DataRptrCollection fz(n.size()), fz_z(n.size()); // Interpolation cofficient and its gradient wrt z
	DataRptrCollection grad_n_sq(n.size()); // Square of grad_n
	for(size_t j=0; j<n.size(); j++)
	{	grad_n_sq[j] = lengthSquared(gradient(n[j]));	
		nullToZero(tauW[j], e.gInfo); applyFunc_r(e.gInfo, tauW_kernel, tauW[j]->data(), grad_n_sq[j]->data(), n[j]->data());
		nullToZero(spness[j], e.gInfo); applyFunc_r(e.gInfo, spness_kernel, tauW[j]->data(), tau[j]->data(), spness[j]->data());
		fz[j] = pow(spness[j], sp_constraint);
		fz_z[j] = sp_constraint*pow(spness[j], sp_constraint-1);
	}

	// Hartree energy and potential of a single channel
	DataRptrCollection Vhartree(n.size()), Ehartree(n.size());
	DataGptrCollection nTilde = J(n);
	for(size_t j=0; j<n.size(); j++)
	{	Vhartree[j] = I((*(e.coulomb))(nTilde[j]));
		Vhartree[j] = Vhartree[j];
		Ehartree[j] = 0.5*n[j]*Vhartree[j];
	}

	// Gradient of tauW wrt n
	DataRptrCollection tauW_n(n.size());
	for(size_t j=0; j<n.size(); j++)
	{	DataRptr lap_n = I(L(J(n[j])));
		nullToZero(tauW_n[j], e.gInfo);
		tauW_n[j] = 0.125*grad_n_sq[j]*pow(nInverse[j],2) - 0.25*lap_n*nInverse[j];
		//applyFunc_r(e.gInfo, tauW_n_kernel, tauW_n[j]->data(), grad_n_sq[j]->data(), n[j]->data(), lap_n->data());
	}

	// Slater exchange
	DataRptrCollection Vslater(n.size());
	DataRptrCollection Eslater(n.size());
	for(size_t j=0; j<n.size(); j++)
	{	double coef = pow(3/M_PI, 1./3.);
		Vslater[j] = - coef * pow(2.*n[j], 1./3.);
		Vslater[j] = Vslater[j];
		Eslater[j] = -(3./4.)*coef*pow(2.*n[j], 4./3.) * 0.5;  // x2 and x0.5 are because of spin polarization
		Eslater[j] = Eslater[j];
	}
	
	for(size_t j=0; j<n.size(); j++)
	{
		// Update Vtau
		e.eVars.Vtau[j] = fz_z[j]*tauW[j]*pow(tauInverse[j],2)*(Eslater[j] + Ehartree[j]*NInverse[j]);
		//applyFunc_r(e.gInfo, Vtau_kernel, e.eVars.Vtau[j]->data(), fz_z[j]->data(), tauW[j]->data(), tau[j]->data(), Eslater[j]->data(),
		//Ehartree[j]->data(), n[j]->data(), integral(n[j]));
		e.eVars.Vtau[j] = JdagOJ(e.eVars.Vtau[j]);
		
		// Update Vscloc
		e.eVars.Vscloc[j] += JdagOJ(- fz[j] * (Vslater[j] + Vhartree[j]*NInverse[j]));
		e.eVars.Vscloc[j] += JdagOJ(- tauInverse[j] * fz_z[j] * tauW_n[j] * (Eslater[j] + Ehartree[j]*NInverse[j]));
	}
	
	// Calculate energy
	e.ener.E["Econstraint"] = 0.;
	for(size_t j=0; j<n.size(); j++)
		e.ener.E["Econstraint"] += integral(-fz[j]*(Eslater[j] + (1./integral(n[j]))*Ehartree[j]));

	watch.stop();
}
