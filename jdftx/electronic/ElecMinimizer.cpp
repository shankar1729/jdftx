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
#include <electronic/Everything.h>
#include <electronic/operators.h>
#include <electronic/Dump.h>
#include <electronic/FluidSolver.h>
#include <core/Random.h>
#include <core/Data.h>
#include <ctime>
#include <electronic/SCF.h>

void ElecGradient::init(Everything& e)
{	init(e.eInfo.nStates);
}

void ElecGradient::init(int nStates)
{	Y.resize(nStates);
	B.resize(nStates);
}

ElecGradient& ElecGradient::operator*=(double alpha)
{	for(unsigned s=0; s<Y.size(); s++)
	{	if(Y[s]) Y[s] *= alpha;
		if(B[s]) B[s] *= alpha;
	}
	return *this;
}

void axpy(double alpha, const ElecGradient& x, ElecGradient& y)
{	assert(x.Y.size() == y.Y.size());
	for(unsigned s=0; s<x.Y.size(); s++)
	{	if(x.Y[s]) { if(y.Y[s]) axpy(alpha, x.Y[s], y.Y[s]); else y.Y[s] = alpha*x.Y[s]; }
		if(x.B[s]) { if(y.B[s]) axpy(alpha, x.B[s], y.B[s]); else y.B[s] = alpha*x.B[s]; }
	}
}

double dot(const ElecGradient& x, const ElecGradient& y)
{	assert(x.Y.size() == y.Y.size());
	complex result(0,0);
	for(unsigned s=0; s<x.Y.size(); s++)
	{	if(x.Y[s] && y.Y[s]) result += dotc(x.Y[s], y.Y[s])*2.0;
		if(x.B[s] && y.B[s]) result += dotc(x.B[s], y.B[s]);
	}
	return result.real();
}

ElecGradient clone(const ElecGradient& x)
{	return x;
}

void randomize(ElecGradient& x)
{	randomize(x.Y);
	for(unsigned s=0; s<x.Y.size(); s++)
		if(x.B[s])
		{	complex* Bdata = x.B[s].data();
			for(unsigned i=0; i<x.B[s].nData(); i++)
				Bdata[i] = Random::normalComplex();
			x.B[s] = dagger_symmetrize(x.B[s]); //make hermitian
		}
}



ElecMinimizer::ElecMinimizer(Everything& e, bool precond)
: e(e), eVars(e.eVars), eInfo(e.eInfo), precond(precond)
{	if(precond) Kgrad.init(e);
}

void ElecMinimizer::step(const ElecGradient& dir, double alpha)
{	assert(dir.Y.size() == eVars.Y.size());
	for(unsigned s=0; s<dir.Y.size(); s++)
	{	if(dir.Y[s]) axpy(alpha, dir.Y[s], eVars.Y[s]);
		if(dir.B[s]) axpy(alpha, dir.B[s], eVars.B[s]);
	}
}

double ElecMinimizer::compute(ElecGradient* grad)
{
	if(grad) grad->init(e);
	return e.eVars.elecEnergyAndGrad(e.ener, grad, precond ? &Kgrad : 0);
}

ElecGradient ElecMinimizer::precondition(const ElecGradient& grad)
{	return precond ? Kgrad : grad;
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
	{	double mu = std::isnan(eInfo.mu)
			? eInfo.findMu(eVars.Hsub_eigs, eInfo.nElectrons) //get mu from eigenvalues
			: eInfo.mu; //report the target mu
		logPrintf("FillingsAux:  mu: %.15le  nElectrons: %.15le\n", mu, eInfo.nElectrons);
	}
	
	//Fillings mix update:
	bool stateModified = false;
	if(eInfo.fillingsUpdate==ElecInfo::FermiFillingsMix
		&& iter % eInfo.mixInterval==0)
	{	eInfo.mixFillings(eVars.F, e.ener);
		stateModified = true;
	}
	
	//Dump:
	e.dump(DumpFreq_Electronic, iter);
	
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
		for(unsigned q=0; q<dir.Y.size(); q++) if(dir.Y[q])
		{	int nOcc = eVars.nOccupiedBands(q);
			if(nOcc)
				callPref(eblas_zero)(dir.Y[q].colLength()*nOcc, dir.Y[q].dataPref());
		}
	}
}


void elecMinimize(Everything& e)
{	
	if(e.cntrl.minimisingResidual)
	{	SCF scf(e);
		scf.minimize();
	}
	else if((not e.cntrl.fixed_n) or e.exCorr.exxFactor() or e.eInfo.hasU)
	{	ElecMinimizer emin(e, true);
		emin.minimize(e.elecMinParams);
	}
	else
	{	logPrintf("Minimization will be done independently for each quantum number.\n");
		for(int q = 0; q < e.eInfo.nStates; q++)
		{	logPrintf("\n---- Minimization of quantum number: %i  kpoint: [%1.6f, %1.6f, %1.6f]  weight: %1.6f",
						q, e.eInfo.qnums[q].k[0], e.eInfo.qnums[q].k[1], e.eInfo.qnums[q].k[2], e.eInfo.qnums[q].weight);
			if(e.eInfo.qnums[q].spin)
				logPrintf("  spin: %i", e.eInfo.qnums[q].spin);
			logPrintf(" ----\n");
			BandMinimizer bmin(e, q, true);
			bmin.minimize(e.elecMinParams);
		}
		// Recompute energy to get the same Eband as the simultaneous minimization
		e.eVars.activeSubspace = -1;
		e.eVars.elecEnergyAndGrad(e.ener);
	}
	e.eVars.setEigenvectors();
	e.eVars.isRandom = false; //wavefunctions are no longer random
}

void elecFluidMinimize(Everything &e)
{	Control &cntrl = e.cntrl;
	ElecVars &eVars = e.eVars;
	ElecInfo& eInfo = e.eInfo;
	Energies &ener = e.ener;

	if(!eVars.HauxInitialized && eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux)
	{	if(std::isnan(eInfo.mu)) //Constant nElectrons mode
		{	logPrintf("\nSetting the auxilliary hamiltonian equal to the subspace hamiltonian.\n");
			//calculate Hsub at current fillings:
			eInfo.fillingsUpdate=ElecInfo::ConstantFillings; eInfo.subspaceRotation=false;
			eVars.elecEnergyAndGrad(e.ener, 0, 0, true);
			eInfo.fillingsUpdate=ElecInfo::FermiFillingsAux; eInfo.subspaceRotation=true;
			//Update B:
			for(int q=0; q<eInfo.nStates; q++) eVars.B[q] = eVars.Hsub[q];
			eVars.HauxInitialized = true;
		}
		else //constant mu mode
			die("Constant chemical potential auxilliary hamiltonian fillings requires the\n"
				"auxilliary hamiltonian to be read in using Haux-initial (or initial-state)\n");
	}
	
	//Prevent change in mu from abruptly changing electron count:
	if(eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux && !std::isnan(eInfo.mu))
	{	for(int q=0; q<eInfo.nStates; q++)
			eVars.B[q].diagonalize(eVars.B_evecs[q], eVars.B_eigs[q]);
		double mu = eInfo.findMu(eVars.B_eigs, eInfo.nElectrons);
		logPrintf("Shifting auxilliary hamiltonian by %lf to set nElectrons=%lf\n", eInfo.mu-mu, eInfo.nElectrons);
		for(int q=0; q<eInfo.nStates; q++)
			eVars.B[q] += eye(eInfo.nBands)*(eInfo.mu-mu);
	}
	
	if(eVars.isRandom && eVars.fluidParams.fluidType!=FluidNone)
	{	logPrintf("Fluid solver invoked on fresh (partially random / LCAO) wavefunctions\n");
		logPrintf("Running a vacuum solve first:\n");
		FluidType origType = eVars.fluidParams.fluidType;
		eVars.fluidParams.fluidType = FluidNone; //temporarily disable the fluid
		logPrintf("\n-------- Initial electronic minimization -----------\n"); logFlush();
		elecMinimize(e); //minimize without fluid
		eVars.fluidParams.fluidType = origType; //restore fluid flag
	}
	eVars.elecEnergyAndGrad(ener);
	
	if(eVars.fluidParams.fluidType==FluidNone || !eVars.fluidSolver->needsGummel())
	{	//no fluid, or fluid doesn't require gummel
		logPrintf("\n-------- Electronic minimization -----------\n"); logFlush();
		elecMinimize(e); //non-gummel fluid will be minimized each EdensityAndVscloc() [see EnergyAndGradient.cpp]
	}
	else
	{	//gummel loop
		logPrintf("\n-------- Electron <-> Fluid self-consistency loop -----------\n"); logFlush();
		for(int iGummel=0; iGummel<cntrl.fluidGummel_nIterations && !killFlag; iGummel++)
		{
			//Fluid-side:
			logPrintf("\n---------------------- Fluid Minimization # %d -----------------------\n", iGummel+1); logFlush();
			double A_diel_prev = ener.E["A_diel"];
			eVars.fluidSolver->minimizeFluid();
			ener.E["A_diel"] = eVars.fluidSolver->get_Adiel_and_grad(eVars.d_fluid, eVars.V_cavity, eVars.fluidForces)
						+ (eInfo.nElectrons - e.iInfo.getZtot())*e.iInfo.ionWidthMuCorrection() ;
			double dAfluid = ener.E["A_diel"] - A_diel_prev;
			logPrintf("\nFluid minimization # %d changed total free energy by %le\n", iGummel+1, dAfluid);

			//Electron-side:
			logPrintf("\n-------------------- Electronic Minimization # %d ---------------------\n", iGummel+1); logFlush();
			double A_JDFT_prev = relevantFreeEnergy(e) + dAfluid; //add dAfluid because change in Adiel not yet included in Etot
			elecMinimize(e); 
			double dAelec = relevantFreeEnergy(e) - A_JDFT_prev;
			logPrintf("\nElectronic minimization # %d changed total free energy by %le\n", iGummel+1, dAelec);
			
			//Dump:
			e.dump(DumpFreq_Gummel, iGummel);

			//Check self-consistency:
			if(fabs(dAelec)<cntrl.fluidGummel_Atol && fabs(dAfluid)<cntrl.fluidGummel_Atol)
			{
				logPrintf("\nFluid<-->Electron self-consistency loop converged to %le hartrees after %d minimization pairs.\n",
					cntrl.fluidGummel_Atol, iGummel+1);
				return;
			}
		}
		logPrintf("\nFluid<-->Electron self-consistency loop not yet converged to %le hartrees after %d minimization pairs.\n",
			cntrl.fluidGummel_Atol, cntrl.fluidGummel_nIterations);
	}
}
