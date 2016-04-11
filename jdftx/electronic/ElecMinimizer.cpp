/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

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

#include <electronic/ElecMinimizer.h>
#include <electronic/BandMinimizer.h>
#include <electronic/BandDavidson.h>
#include <electronic/Everything.h>
#include <electronic/operators.h>
#include <electronic/Dump.h>
#include <fluid/FluidSolver.h>
#include <core/Random.h>
#include <core/ScalarField.h>
#include <ctime>
#include <electronic/SCF.h>

void ElecGradient::init(const Everything& e)
{	eInfo = &e.eInfo;
	Y.resize(eInfo->nStates);
	B.resize(eInfo->nStates);
}

ElecGradient& ElecGradient::operator*=(double alpha)
{	for(int q=eInfo->qStart; q<eInfo->qStop; q++)
	{	if(Y[q]) Y[q] *= alpha;
		if(B[q]) B[q] *= alpha;
	}
	return *this;
}

void axpy(double alpha, const ElecGradient& x, ElecGradient& y)
{	assert(x.eInfo == y.eInfo);
	for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++)
	{	if(x.Y[q]) { if(y.Y[q]) axpy(alpha, x.Y[q], y.Y[q]); else y.Y[q] = alpha*x.Y[q]; }
		if(x.B[q]) { if(y.B[q]) axpy(alpha, x.B[q], y.B[q]); else y.B[q] = alpha*x.B[q]; }
	}
}

double dot(const ElecGradient& x, const ElecGradient& y)
{	assert(x.eInfo == y.eInfo);
	complex result(0,0);
	for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++)
	{	if(x.Y[q] && y.Y[q]) result += dotc(x.Y[q], y.Y[q])*2.0;
		if(x.B[q] && y.B[q]) result += dotc(x.B[q], y.B[q]);
	}
	mpiUtil->allReduce(result.real(), MPIUtil::ReduceSum);
	return result.real();
}

ElecGradient clone(const ElecGradient& x)
{	return x;
}

void randomize(ElecGradient& x)
{	randomize(x.Y, *x.eInfo);
	for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++)
		if(x.B[q])
		{	randomize(x.B[q]);
			x.B[q] = dagger_symmetrize(x.B[q]); //make hermitian
		}
}



ElecMinimizer::ElecMinimizer(Everything& e)
: e(e), eVars(e.eVars), eInfo(e.eInfo), Knorm(0.)
{	Kgrad.init(e);
}

void ElecMinimizer::step(const ElecGradient& dir, double alpha)
{	assert(dir.eInfo == &eInfo);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
	{	if(dir.Y[q]) axpy(alpha, dir.Y[q], eVars.Y[q]);
		if(dir.B[q]) axpy(alpha, dir.B[q], eVars.B[q]);
	}
}

double ElecMinimizer::compute(ElecGradient* grad)
{	if(grad) grad->init(e);
	double ener = e.eVars.elecEnergyAndGrad(e.ener, grad, grad ? &Kgrad : 0);
	if(grad) Knorm = sync(dot(*grad, Kgrad));
	return ener;
}

ElecGradient ElecMinimizer::precondition(const ElecGradient& grad)
{	return Kgrad;
}

bool ElecMinimizer::report(int iter)
{
	if(e.cntrl.shouldPrintEcomponents)
	{	//Print the iteration header
		time_t timenow = time(0);
		logPrintf("------------------------------------------------------\n");
		logPrintf("Iteration %d   %s\n",iter, ctime(&timenow));
		//Print the energies
		e.ener.print(); logPrintf("\n"); logFlush();
	}
	if(e.cntrl.shouldPrintEigsFillings)
	{	//Print the eigenvalues if requested
		print_Hsub_eigs(e);
		logPrintf("\n"); logFlush();
	}
	
	//Auxiliary hamiltonian fillings printout:
	if(eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux)
		eInfo.printFermi("Aux");
	
	//Fillings mix update:
	bool stateModified = false;
	if(eInfo.fillingsUpdate==ElecInfo::FermiFillingsMix
		&& iter % eInfo.mixInterval==0)
	{	eInfo.mixFillings(eVars.F, e.ener);
		stateModified = true;
	}
	
	//Dump:
	e.dump(DumpFreq_Electronic, iter);
	
	//Handle indefinite preconditioner issue for FermiFillingsAux:
	if(Knorm < 0.)
	{	logPrintf("%s\tPreconditioner indefiniteness detected (grad_K will become NAN): ", e.elecMinParams.linePrefix);
		convergeEmptyStates(e);
		return true;
	}
	
	//Re-orthogonalize wavefunctions if necessary:
	if(e.cntrl.overlapCheckInterval
		&& (iter % e.cntrl.overlapCheckInterval == 0)
		&& (eVars.overlapCondition > e.cntrl.overlapConditionThreshold) )
	{
		logPrintf("%s\tCondition number of orbital overlap matrix (%lg) exceeds threshold (%lg): ",
			e.elecMinParams.linePrefix, eVars.overlapCondition, e.cntrl.overlapConditionThreshold);
		eVars.setEigenvectors();
		return true;
	}
	
	return stateModified;
}

void ElecMinimizer::constrain(ElecGradient& dir)
{	if(e.cntrl.fixOccupied)
	{	//Project out occupied directions:
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{	int nOcc = eVars.nOccupiedBands(q);
			if(nOcc)
				callPref(eblas_zero)(dir.Y[q].colLength()*nOcc, dir.Y[q].dataPref());
		}
	}
}

double ElecMinimizer::sync(double x) const
{	mpiUtil->bcast(x);
	return x;
}

void bandMinimize(Everything& e)
{	bool fixed_H = true; std::swap(fixed_H, e.cntrl.fixed_H); //remember fixed_H flag and temporarily set it to true
	logPrintf("Minimization will be done independently for each quantum number.\n");
	e.ener.Eband = 0.;
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{	logPrintf("\n---- Minimization of quantum number: "); e.eInfo.kpointPrint(globalLog, q, true); logPrintf(" ----\n");
		switch(e.cntrl.elecEigenAlgo)
		{	case ElecEigenCG: { BandMinimizer(e, q, true).minimize(e.elecMinParams); break; }
			case ElecEigenDavidson: { BandDavidson(e, q).minimize(); break; }
		}
		e.ener.Eband += e.eInfo.qnums[q].weight * trace(e.eVars.Hsub_eigs[q]);
	}
	mpiUtil->allReduce(e.ener.Eband, MPIUtil::ReduceSum);
	if(e.cntrl.shouldPrintEigsFillings)
	{	//Print the eigenvalues if requested
		print_Hsub_eigs(e);
		logPrintf("\n"); logFlush();
	}
	std::swap(fixed_H, e.cntrl.fixed_H); //restore fixed_H flag
	e.eVars.setEigenvectors();
}

void elecMinimize(Everything& e)
{	
	if(e.cntrl.scf)
	{	SCF scf(e);
		scf.minimize();
	}
	else if((not e.cntrl.fixed_H) or e.exCorr.exxFactor())
	{	ElecMinimizer emin(e);
		emin.minimize(e.elecMinParams);
		if (!e.ionDynamicsParams.tMax) e.eVars.setEigenvectors(); //Don't spend time with this if running MD
	}
	else
	{	bandMinimize(e);
	}
	e.eVars.isRandom = false; //wavefunctions are no longer random
	//Converge empty states if necessary:
	if(e.cntrl.convergeEmptyStates and (not e.cntrl.fixed_H))
		convergeEmptyStates(e);
}

void elecFluidMinimize(Everything &e)
{	Control &cntrl = e.cntrl;
	ElecVars &eVars = e.eVars;
	ElecInfo& eInfo = e.eInfo;
	IonInfo& iInfo = e.iInfo;
	Energies &ener = e.ener;

	if(!eVars.HauxInitialized && eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux)
	{	if(std::isnan(eInfo.mu)) //Constant nElectrons mode
		{	logPrintf("\nSetting the auxilliary hamiltonian equal to the subspace hamiltonian.\n");
			//calculate Hsub at current fillings:
			eInfo.fillingsUpdate=ElecInfo::ConstantFillings; eInfo.subspaceRotation=false;
			eVars.elecEnergyAndGrad(e.ener, 0, 0, true);
			eInfo.fillingsUpdate=ElecInfo::FermiFillingsAux; eInfo.subspaceRotation=true;
			//Update B:
			for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++) eVars.B[q] = eVars.Hsub[q];
			eVars.HauxInitialized = true;
		}
		else //constant mu mode
			die("Constant chemical potential auxilliary hamiltonian fillings requires the\n"
				"auxilliary hamiltonian to be read in using Haux-initial (or initial-state)\n");
	}
	
	//Prevent change in mu from abruptly changing electron count:
	if(eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux && !std::isnan(eInfo.mu))
	{	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
			eVars.B[q].diagonalize(eVars.B_evecs[q], eVars.B_eigs[q]);
		double Bz, mu = eInfo.findMu(eVars.B_eigs, eInfo.nElectrons, Bz);
		logPrintf("Shifting auxilliary hamiltonian by %lf to set nElectrons=%lf\n", eInfo.mu-mu, eInfo.nElectrons);
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
			eVars.B[q] += eye(eInfo.nBands)*(eInfo.mu-mu);
	}
	
	double Evac0 = NAN;
	if(eVars.isRandom && eVars.fluidParams.fluidType!=FluidNone)
	{	logPrintf("Fluid solver invoked on fresh (partially random / LCAO) wavefunctions\n");
		logPrintf("Running a vacuum solve first:\n");
		FluidType origType = eVars.fluidParams.fluidType;
		eVars.fluidParams.fluidType = FluidNone; //temporarily disable the fluid
		double muOrig = eInfo.mu;
		eInfo.mu = NAN; //temporarily disable fixed mu (if present)
		logPrintf("\n-------- Initial electronic minimization -----------\n"); logFlush();
		elecMinimize(e); //minimize without fluid
		Evac0 = relevantFreeEnergy(e);
		logPrintf("Vacuum energy after initial minimize, %s = %+.15f\n\n", relevantFreeEnergyName(e), Evac0);
		eVars.fluidParams.fluidType = origType; //restore fluid flag
		eInfo.mu = muOrig; //restore mu target (if any)
	}
	eVars.elecEnergyAndGrad(ener);
	
	//First electronic minimization (with fluid if present) in most cases
	if(!( eVars.fluidParams.fluidType!=FluidNone
		&& eVars.fluidSolver->useGummel()
		&& (fabs(eInfo.nElectrons-iInfo.getZtot())>=1e-12) )) //skip only if non-neutral Gummel-requiring fluid
	{
		logPrintf("\n-------- Electronic minimization -----------\n"); logFlush();
		elecMinimize(e); //non-gummel fluid will be minimized each EdensityAndVscloc() [see ElecVars.cpp]
	}
	
	if(eVars.fluidParams.fluidType!=FluidNone && eVars.fluidSolver->useGummel())
	{	//gummel loop
		logPrintf("\n-------- Electron <-> Fluid self-consistency loop -----------\n"); logFlush();
		double dAtyp = 1.;
		bool converged = false;
		for(int iGummel=0; iGummel<cntrl.fluidGummel_nIterations && !killFlag; iGummel++)
		{
			//Fluid-side:
			logPrintf("\n---------------------- Fluid Minimization # %d -----------------------\n", iGummel+1); logFlush();
			double A_diel_prev = ener.E["A_diel"];
			e.fluidMinParams.energyDiffThreshold = std::min(1e-5, 0.01*dAtyp);
			eVars.fluidSolver->minimizeFluid();
			ener.E["A_diel"] = eVars.fluidSolver->get_Adiel_and_grad(&eVars.d_fluid, &eVars.V_cavity);
			double dAfluid = ener.E["A_diel"] - A_diel_prev;
			logPrintf("\nFluid minimization # %d changed total free energy by %le at t[s]: %9.2lf\n", iGummel+1, dAfluid, clock_sec());

			//Electron-side:
			logPrintf("\n-------------------- Electronic Minimization # %d ---------------------\n", iGummel+1); logFlush();
			double A_JDFT_prev = relevantFreeEnergy(e);
			e.elecMinParams.energyDiffThreshold = std::min(1e-5, 0.01*dAtyp);
			elecMinimize(e);
			double dAelec = relevantFreeEnergy(e) - A_JDFT_prev;
			logPrintf("\nElectronic minimization # %d changed total free energy by %le at t[s]: %9.2lf\n", iGummel+1, dAelec, clock_sec());
			
			//Dump:
			e.dump(DumpFreq_Gummel, iGummel);

			//Check self-consistency:
			dAtyp = std::max(fabs(dAfluid), fabs(dAelec));
			if(dAtyp<cntrl.fluidGummel_Atol)
			{	logPrintf("\nFluid<-->Electron self-consistency loop converged to %le hartrees after %d minimization pairs at t[s]: %9.2lf.\n",
					cntrl.fluidGummel_Atol, iGummel+1, clock_sec());
				converged = true;
				break;
			}
		}
		if(!converged)
			logPrintf("\nFluid<-->Electron self-consistency loop not yet converged to %le hartrees after %d minimization pairs at t[s]: %9.2lf.\n",
				cntrl.fluidGummel_Atol, cntrl.fluidGummel_nIterations, clock_sec());
	}
	
	if(!std::isnan(Evac0))
		logPrintf("Single-point solvation energy estimate, Delta%s = %+.15f\n", relevantFreeEnergyName(e), relevantFreeEnergy(e)-Evac0);
}

void convergeEmptyStates(Everything& e)
{	logPrintf("Converging empty states (this may take a while): "); logFlush();
	std::vector<diagMatrix> eigsPrev = e.eVars.Hsub_eigs;
	logSuspend(); e.elecMinParams.fpLog = nullLog;
	bandMinimize(e); //this will also set state to eigenvectors
	logResume(); e.elecMinParams.fpLog = globalLog;
	e.ener.Eband = 0.; //only affects printing (if non-zero Energies::print assumes band structure calc)
	logPrintf("|deigs|: %.3e\n", SCF::eigDiffRMS(e.eVars.Hsub_eigs, eigsPrev, e)); logFlush();
}

