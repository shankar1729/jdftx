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
#include <electronic/ColumnBundle.h>
#include <electronic/Everything.h>
#include <electronic/Dump.h>
#include <fluid/FluidSolver.h>
#include <core/Random.h>
#include <core/ScalarField.h>
#include <ctime>
#include <electronic/SCF.h>

void ElecGradient::init(const Everything& e)
{	eInfo = &e.eInfo;
	C.resize(eInfo->nStates);
	Haux.resize(eInfo->nStates);
}

ElecGradient& ElecGradient::operator*=(double alpha)
{	for(int q=eInfo->qStart; q<eInfo->qStop; q++)
	{	if(C[q]) C[q] *= alpha;
		if(Haux[q]) Haux[q] *= alpha;
	}
	return *this;
}

void axpy(double alpha, const ElecGradient& x, ElecGradient& y)
{	assert(x.eInfo == y.eInfo);
	for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++)
	{	if(x.C[q]) { if(y.C[q]) axpy(alpha, x.C[q], y.C[q]); else y.C[q] = alpha*x.C[q]; }
		if(x.Haux[q]) { if(y.Haux[q]) axpy(alpha, x.Haux[q], y.Haux[q]); else y.Haux[q] = alpha*x.Haux[q]; }
	}
}

double dot(const ElecGradient& x, const ElecGradient& y, double* auxContrib)
{	assert(x.eInfo == y.eInfo);
	std::vector<double> result(2, 0.); //calculate wavefunction and auxiliary contributions separately
	for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++)
	{	if(x.C[q] && y.C[q]) result[0] += dotc(x.C[q], y.C[q]).real()*2.0;
		if(x.Haux[q] && y.Haux[q]) result[1] += dotc(x.Haux[q], y.Haux[q]).real();
	}
	mpiUtil->allReduce(result.data(), 2, MPIUtil::ReduceSum);
	if(auxContrib) *auxContrib=result[1]; //store auxiliary contribution, if requested
	return result[0]+result[1]; //return total
}

ElecGradient clone(const ElecGradient& x)
{	return x;
}

void randomize(ElecGradient& x)
{	randomize(x.C, *x.eInfo);
	for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++)
		if(x.Haux[q])
		{	randomize(x.Haux[q]);
			x.Haux[q] = dagger_symmetrize(x.Haux[q]); //make hermitian
		}
}

struct SubspaceRotationAdjust
{
	Everything& e;
	bool adjust; //whether adjustment is active (otherwise only check for Knorm indefiniteness)
	std::vector<matrix> KgPrevHaux; //preconditioned auxiliary gradient at previous step
	double gDotKgPrevHaux; //overlap of current auxiliary gradient with KgPrevHaux
	double cumulatedScale; //net change in subspace rotation factor since last direction reset
	double KnormTot, KnormAux; //latest preconditioned gradient norms: total and auxiliary alone
	
	SubspaceRotationAdjust(Everything& e)
	: e(e), adjust(e.cntrl.subspaceRotationAdjust), gDotKgPrevHaux(0.), cumulatedScale(1.)
	{
	}
	
	void cacheGradientOverlaps(const ElecGradient& grad, const ElecGradient& Kgrad)
	{	//Calculate overlaps of current gradient:
		KnormTot = dot(grad, Kgrad, &KnormAux);
		mpiUtil->bcast(KnormTot);
		mpiUtil->bcast(KnormAux);
		
		//Compute auxiliary overlap with previous preconditioned gradient:
		if(KgPrevHaux.size())
		{	gDotKgPrevHaux = 0.;
			for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
				gDotKgPrevHaux += dotc(grad.Haux[q], KgPrevHaux[q]).real();
			mpiUtil->allReduce(gDotKgPrevHaux, MPIUtil::ReduceSum);
			mpiUtil->bcast(gDotKgPrevHaux);
		}
	}
	
	//Handle indefinite Knorm, adjust subspace rotation factor if necessary and report changes.
	//Returns true if CG needs to be reset
	bool report(const std::vector<matrix>& KgradHaux)
	{
		//Handle indefinite preconditioner issues:
		if(KnormTot <= 0.)
		{	logPrintf("%s\tPreconditioner indefiniteness detected (grad_K will become NAN): ", e.elecMinParams.linePrefix);
			convergeEmptyStates(e);
			return true; //resets CG
		}
		//Adjust subspace rotation factor if necessary:
		if(adjust)
		{	if(gDotKgPrevHaux)
			{	double& kappa = e.cntrl.subspaceRotationFactor;
				double Emin_kappa = gDotKgPrevHaux / KnormTot; //dimensionless version of dEmin/d(log kappa)
				double kappaMax = 0.5*kappa*fabs(KnormTot)/fabs(KnormAux); //aux gradient magnitude should not exceed half that of total
				//Heuristic for adjusting subspace rotation factor (called kappa here):
				double scaleFactor = exp(Emin_kappa/hypot(1,Emin_kappa)); //saturated to range (1/e,e)
				scaleFactor = std::min(scaleFactor, kappaMax/kappa); //cap scale factor at maximum
				kappa *= scaleFactor;
				cumulatedScale *= scaleFactor;
				logPrintf("\tSubspaceRotationAdjust: set factor to %.3lg\n", kappa);
				if(fabs(log(cumulatedScale)) > 2.) //cumulated scale adjustment > e^2
				{	logPrintf("\tSubspaceRotationAdjust: resetting CG because factor has changed by %lg\n", cumulatedScale);
					cumulatedScale = 1.;
					return true; //resets CG
				}
			}
			KgPrevHaux = KgradHaux; //Cached preconditioned gradient for next step
		}
		return false;
	}
};


ElecMinimizer::ElecMinimizer(Everything& e)
: e(e), eVars(e.eVars), eInfo(e.eInfo), rotPrev(eInfo.nStates), rotPrevC(eInfo.nStates), rotPrevCinv(eInfo.nStates)
{
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
	{	rotPrev[q] = eye(eInfo.nBands);
		rotPrevC[q] = eye(eInfo.nBands);
		rotPrevCinv[q] = eye(eInfo.nBands);
	}
	rotExists = false; //rotation is identity
	
	//Initialize subspace rotation adjuster if required:
	if(e.cntrl.subspaceRotationAdjust && ( eInfo.fillingsUpdate==ElecInfo::FillingsHsub || !eInfo.scalarFillings) )
		sra = std::make_shared<SubspaceRotationAdjust>(e);
}

void ElecMinimizer::step(const ElecGradient& dir, double alpha)
{	assert(dir.eInfo == &eInfo);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
	{	axpy(alpha, rotExists ? dir.C[q]*rotPrevC[q] : dir.C[q], eVars.C[q]);
		if(eInfo.fillingsUpdate==ElecInfo::FillingsConst && eInfo.scalarFillings)
		{	//Constant scalar fillings: no rotations required
			eVars.orthonormalize(q);
		}
		else
		{	//Haux or non-scalar fillings: rotations required
			assert(dir.Haux[q]);
			matrix rot;
			if(eInfo.fillingsUpdate == ElecInfo::FillingsHsub)
			{	//Haux fillings:
				matrix Haux = eVars.Haux_eigs[q];
				axpy(alpha, rotExists ? dagger(rotPrev[q])*dir.Haux[q]*rotPrev[q] : dir.Haux[q], Haux);
				Haux.diagonalize(rot, eVars.Haux_eigs[q]); //rotation chosen to diagonalize auxiliary matrix
			}
			else
			{	//Non-scalar fillings:
				assert(!eInfo.scalarFillings);
				rot = cis(alpha * dir.Haux[q]); //auxiliary matrix directly generates rotations
			}
			matrix rotC = rot;
			eVars.orthonormalize(q, &rotC);
			rotPrev[q] = rotPrev[q] * rot;
			rotPrevC[q] = rotPrevC[q] * rotC;
			rotPrevCinv[q] = inv(rotC) * rotPrevCinv[q];
			rotExists = true; //rotation is no longer identity
		}
	}
}

double ElecMinimizer::compute(ElecGradient* grad, ElecGradient* Kgrad)
{	if(grad) grad->init(e);
	if(Kgrad) Kgrad->init(e);
	double ener = e.eVars.elecEnergyAndGrad(e.ener, grad, Kgrad);
	if(grad)
	{	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{	//Rotate wavefunction gradients if necessary:
			if(rotExists)
			{	grad->C[q] = grad->C[q] * rotPrevCinv[q];
				Kgrad->C[q] = Kgrad->C[q] * rotPrevCinv[q];
			}
			//Subspace gradient handling depends on mode:
			if(eInfo.fillingsUpdate == ElecInfo::FillingsHsub)
			{	//Haux fillings: rotate gradient computed by ElecVars if necessary
				if(rotExists)
				{	grad->Haux[q] = rotPrev[q] * grad->Haux[q] * dagger(rotPrev[q]);
					Kgrad->Haux[q] = rotPrev[q] * Kgrad->Haux[q] * dagger(rotPrev[q]);
				}
			}
			else if(!eInfo.scalarFillings)
			{	//Non-scalar fillings:
				grad->Haux[q] = dagger_symmetrize(complex(0,1) * (eVars.F[q]*eVars.Hsub[q] - eVars.Hsub[q]*eVars.F[q]));
				Kgrad->Haux[q] = e.cntrl.subspaceRotationFactor * grad->Haux[q];
			}
			//else: constant scalar fillings (no subspace gradient)
		}
		
		//Cache gradient overlaps, if needed, for subspace rotation handling:
		if(sra) sra->cacheGradientOverlaps(*grad, *Kgrad);
		KgradHaux = Kgrad->Haux;
	}
	return ener;
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
	
	//Fillings update report:
	if(eInfo.fillingsUpdate==ElecInfo::FillingsHsub)
		eInfo.smearReport();
	
	//Dump:
	e.dump(DumpFreq_Electronic, iter);
	
	//Re-unitarize rotations:
	if(rotExists)
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{	rotPrevC[q] = rotPrev[q];
			rotPrevCinv[q] = dagger(rotPrev[q]);
		}
	
	//Subspace rotation preconditioner handling:
	if(sra) return sra->report(KgradHaux);
	return false;
}

void ElecMinimizer::constrain(ElecGradient& dir)
{	assert(dir.eInfo == &eInfo);
	//Project component of search direction along current wavefunctions:
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		dir.C[q] -= eVars.C[q] * (eVars.C[q]^O(dir.C[q]));
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
		{	case ElecEigenCG: { BandMinimizer(e, q).minimize(e.elecMinParams); break; }
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


void muOuterLoop(Everything& e)
{	logPrintf("\n-------- mu target loop -----------\n"); logFlush();
	assert(!std::isnan(e.eInfo.mu));
	double muTarget = NAN;
	std::swap(muTarget, e.eInfo.mu); //set ElecInfo::mu to NAN (fixed charge mode)

	#define RunFixedCharge \
		{	e.ener.E["minusMuN"] = -muTarget * Ncur; \
			e.eInfo.nElectrons = Ncur; \
			elecMinimize(e); \
			double Bz; \
			muCur = e.eInfo.findMu(e.eVars.Haux_eigs, Ncur, Bz); \
			Gcur = relevantFreeEnergy(e); \
		}
	
	//Initial calculation:
	double Ncur = e.eInfo.nElectrons, muCur, Gcur;
	RunFixedCharge
	
	//Secant-method:
	double Nprev = NAN, muPrev = NAN, Gprev = NAN; //previous point for secant method
	double EdiffThreshold = e.elecMinParams.energyDiffThreshold;
	double muThreshold = EdiffThreshold;
	double dN0 = 0.1; //initial step size
	for(int iter=0; iter<e.elecMinParams.nIterations; iter++)
	{	//Print and check convergence:
		logPrintf("MuTargetLoop: Iter: %3d  G: %+.15lf  mu: %+.9lf  nElectrons: %.6lf\n", iter, Gcur, muCur, Ncur);
		if(fabs(muCur-muTarget) < muThreshold)
		{	logPrintf("MuTargetLoop: Converged (|mu-muTarget|<%le).\n", muThreshold);
			break;
		}
		if(!std::isnan(Gprev) && (Gcur<Gprev) && fabs(Gcur-Gprev)<EdiffThreshold)
		{	logPrintf("MuTargetLoop: Converged (|Delta G|<%le)\n", EdiffThreshold);
			break;
		}
		//Determine next point:
		double Nnext;
		if(std::isnan(Nprev)) //No previous point, so take a fixed step:
			Nnext = Ncur + copysign(dN0, muTarget-muCur);
		else //Secant step:
			Nnext = Ncur + (Ncur - Nprev)*(muTarget - muCur)/(muCur - muPrev);
		//Switch to next point and calculate:
		Nprev = Ncur; muPrev = muCur; Gprev = Gcur;
		Ncur = Nnext;
		RunFixedCharge
		#undef RunFixedCharge
	}
	std::swap(muTarget, e.eInfo.mu); //restore finite ElecInfo::mu (fixed mu mode)
	e.ener.E["minusMuN"] = 0.; e.ener.muN = e.eInfo.mu*e.eInfo.nElectrons; //restore normal storage of muN component
}


void elecMinimize(Everything& e)
{	
	if(!std::isnan(e.eInfo.mu) && e.eInfo.muLoop)
	{	muOuterLoop(e); //Run a loop over fixed charge calculations to target mu
	}
	else if(e.cntrl.scf)
	{	SCF scf(e);
		scf.minimize();
	}
	else if(e.cntrl.fixed_H)
	{	bandMinimize(e);
	}
	else
	{	ElecMinimizer emin(e);
		emin.minimize(e.elecMinParams);
		if (!e.ionDynamicsParams.tMax) e.eVars.setEigenvectors(); //Don't spend time with this if running MD
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

	if(!eVars.HauxInitialized && eInfo.fillingsUpdate==ElecInfo::FillingsHsub)
	{	if(std::isnan(eInfo.mu)) //Constant nElectrons mode
		{	logPrintf("\nSetting the auxilliary hamiltonian equal to the subspace hamiltonian.\n");
			//calculate Hsub at current fillings:
			eInfo.fillingsUpdate=ElecInfo::FillingsConst;
			eVars.elecEnergyAndGrad(e.ener, 0, 0, true);
			eInfo.fillingsUpdate=ElecInfo::FillingsHsub;
			//Update B:
			for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++) eVars.Haux_eigs[q] = eVars.Hsub_eigs[q];
			eVars.HauxInitialized = true;
		}
		else //constant mu mode
			die("Constant chemical potential auxilliary hamiltonian fillings requires the band\n"
				"eigenvalues to either be read in using one of the two commands initial-state\n"
				"or elec-initial-eigenvals, or be automatically initialized during LCAO.");
	}
	
	double Evac0 = NAN;
	if(eVars.isRandom && eVars.fluidParams.fluidType!=FluidNone)
	{	logPrintf("Fluid solver invoked on fresh (partially random / LCAO) wavefunctions\n");
		logPrintf("Running a vacuum solve first:\n");
		FluidType origType = eVars.fluidParams.fluidType;
		eVars.fluidParams.fluidType = FluidNone; //temporarily disable the fluid
		double muOrig = eInfo.mu;
		eInfo.mu = NAN; //temporarily disable fixed mu (if present)
		if(!std::isnan(muOrig)) e.ener.E["minusMuN"] = -muOrig*e.eInfo.nElectrons; //but add constant to make G printout correct
		logPrintf("\n-------- Initial electronic minimization -----------\n"); logFlush();
		elecMinimize(e); //minimize without fluid
		Evac0 = relevantFreeEnergy(e);
		logPrintf("Vacuum energy after initial minimize, %s = %+.15f\n\n", relevantFreeEnergyName(e), Evac0);
		eVars.fluidParams.fluidType = origType; //restore fluid flag
		eInfo.mu = muOrig; //restore mu target (if any)
		e.ener.E["minusMuN"] = 0.; //remove fixed-mu energy correction (if any)
	}
	
	//Prevent change in mu from abruptly changing electron count:
	if(eInfo.fillingsUpdate==ElecInfo::FillingsHsub && !std::isnan(eInfo.mu))
	{	double Bz, mu = eInfo.findMu(eVars.Haux_eigs, eInfo.nElectrons, Bz);
		logPrintf("Shifting auxilliary hamiltonian by %lf to set nElectrons=%lf\n", eInfo.mu-mu, eInfo.nElectrons);
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
			eVars.Haux_eigs[q] += eye(eInfo.nBands)*(eInfo.mu-mu);
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

