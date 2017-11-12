/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver, Deniz Gunceler

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

#include <electronic/Everything.h>
#include <fluid/NonlinearPCM.h>
#include <fluid/LinearPCM.h>
#include <fluid/PCM_internal.h>
#include <core/ScalarFieldIO.h>
#include <core/Util.h>

//Utility functions to extract/set the members of a MuEps
inline ScalarField& getMuPlus(ScalarFieldMuEps& X) { return X[0]; }
inline const ScalarField& getMuPlus(const ScalarFieldMuEps& X) { return X[0]; }
inline ScalarField& getMuMinus(ScalarFieldMuEps& X) { return X[1]; }
inline const ScalarField& getMuMinus(const ScalarFieldMuEps& X) { return X[1]; }
inline VectorField getEps(ScalarFieldMuEps& X) { return VectorField(&X[2]); }
inline const VectorField getEps(const ScalarFieldMuEps& X) { return VectorField(&X[2]); }
inline void setMuEps(ScalarFieldMuEps& mueps, ScalarField muPlus, ScalarField muMinus, VectorField eps) { mueps[0]=muPlus; mueps[1]=muMinus; for(int k=0; k<3; k++) mueps[k+2]=eps[k]; }


inline double setPreconditioner(double G, double kappaSqByEpsilon, double muByEpsSq)
{	double G2 = G*G;
	double den = G2 + kappaSqByEpsilon;
	return den ? muByEpsSq*G2/(den*den) : 0.;
}

inline void setMetric(int i, double Gsq, double qMetricSq, double* metric)
{	metric[i] = qMetricSq ? Gsq / (qMetricSq + Gsq) : 1.;
}

NonlinearPCM::NonlinearPCM(const Everything& e, const FluidSolverParams& fsp)
: PCM(e, fsp), Pulay(fsp.scfParams), pMol(0.), ionNbulk(0.), ionZ(0.), screeningEval(0), dielectricEval(0)
{
	const auto& solvent = fsp.solvents[0];
	pMol = solvent->pMol ? solvent->pMol : solvent->molecule.getDipole().length();
	
	//Initialize dielectric evaluation class:
	dielectricEval = new NonlinearPCMeval::Dielectric(fsp.linearDielectric,
		fsp.T, solvent->Nbulk, pMol, epsBulk, solvent->epsInf);

	//Check and setup ionic screening:
	if(fsp.cations.size() > 1) die("NonlinearPCM currently only supports a single cationic component.\n");
	if(fsp.anions.size() > 1) die("NonlinearPCM currently only supports a single anionic component.\n");
	assert(fsp.anions.size() == fsp.cations.size()); //this should be ensured by charge neutrality check in FluidSolver constructor
	if(fsp.cations.size())
	{	//Ensure charge balanced:
		if(fabs(fsp.cations[0]->molecule.getCharge() + fsp.anions[0]->molecule.getCharge())>1e-12)
			die("NonlinearPCM currently only supports charge-balanced (Z:Z) electrolytes.\n");
		ionNbulk = fsp.cations[0]->Nbulk;
		ionZ = fsp.anions[0]->molecule.getCharge();
		double VhsCation = fsp.cations[0]->molecule.getVhs(); if(!VhsCation) VhsCation = (4*M_PI/3)*pow(fsp.cations[0]->Rvdw,3);
		double VhsAnion = fsp.anions[0]->molecule.getVhs(); if(!VhsAnion) VhsAnion = (4*M_PI/3)*pow(fsp.anions[0]->Rvdw,3);
		screeningEval = new NonlinearPCMeval::Screening(fsp.linearScreening, fsp.T, ionNbulk, ionZ, VhsCation, VhsAnion, epsBulk);
	}
	else
	{	ionNbulk = 0.;
		screeningEval = 0;
	}
	
	if(fsp.nonlinearSCF)
	{	//Initialize lookup tables for SCF version:
		//--- Dielectric lookup table
		double dxMapped = 1./512;
		std::vector<double> samples;
		for(double xMapped=0.; xMapped<=1.; xMapped+=dxMapped)
		{	double x;
			if(xMapped==0.) x = 1e-12;
			else if(xMapped==1.) x = 1e+12;
			else x = xMapped/(1.-xMapped); //inverse of xMapped = x / (1 + x)
			samples.push_back(dielectricEval->eps_from_x(x)/x);
		}
		gLookup.init(0, samples, dxMapped);
		//--- Screening lookup table
		if(screeningEval)
		{
			double dVmapped = 1./512;
			std::vector<double> samples;
			for(double Vmapped=-1.; Vmapped<=1.; Vmapped+=dVmapped)
			{	if(fabs(Vmapped)==1)
					samples.push_back(0.);
				else
				{	double V = std::pow(Vmapped/(1.-Vmapped*Vmapped), 3.); //inverse of Vmapped = copysign(2cbrt(V) / (1 + sqrt(1 + (2cbrt(V))^2)), V)
					double x = screeningEval->x_from_V(V);
					double xMapped = 1./(1.+x); //maps [0,infty) -> (0,1]
					samples.push_back(xMapped);
				}
			}
			xLookup.init(1, samples, dVmapped);
		}
		//--- Pulay metric
		metric = std::make_shared<RealKernel>(gInfo);
		applyFuncGsq(e.gInfo, setMetric, std::pow(fsp.scfParams.qMetric,2), metric->data());
	}
	else
	{	//Initialize preconditioner (for mu channel):
		double muByEps = (ionZ/pMol) * (1.-dielectricEval->alpha/3); //relative scale between mu and eps
		preconditioner.init(0, 0.02, gInfo.GmaxGrid, setPreconditioner, k2factor/epsBulk, muByEps*muByEps);
	}
}

NonlinearPCM::~NonlinearPCM()
{	delete dielectricEval;
	if(screeningEval) delete screeningEval;
	if(preconditioner) preconditioner.free();
	if(gLookup) gLookup.free();
	if(xLookup) xLookup.free();
}


void NonlinearPCM::set_internal(const ScalarFieldTilde& rhoExplicitTilde, const ScalarFieldTilde& nCavityTilde)
{	
	bool setPhiFromState = false; //whether to set linearPCM phi from state (first time in SCF version when state has been read in)
	
	if(fsp.nonlinearSCF || (!state))
	{	if(!linearPCM)
		{	logSuspend();
			linearPCM = std::make_shared<LinearPCM>(e, fsp);
			logResume();
			if(state) setPhiFromState = true; //set phi from state which has already been loaded (only in SCF mode)
		}
		logSuspend();
		linearPCM->atpos = atpos;
		linearPCM->set_internal(rhoExplicitTilde, nCavityTilde);
		logResume();
	}
	
	//Initialize state if required:
	if(!state)
	{	logPrintf("Initializing state of NonlinearPCM using a similar LinearPCM: "); logFlush();
		FILE*& fpLog = ((MinimizeParams&)e.fluidMinParams).fpLog;
		fpLog = nullLog; //disable iteration log from LinearPCM
		linearPCM->minimizeFluid();
		fpLog = globalLog; //restore usual iteration log
		logFlush();
		//Guess nonlinear states based on the electrostatic potential of the linear version:
		//mu:
		ScalarField mu;
		if(screeningEval && screeningEval->linear)
		{	mu = (-ionZ/fsp.T) * I(linearPCM->state);
			mu -= integral(mu)/gInfo.detR; //project out G=0
		}
		else initZero(mu, gInfo); //initialization logic does not work well with hard sphere limit
		//eps:
		VectorField eps = (-pMol/fsp.T) * I(gradient(linearPCM->state));
		ScalarField E = sqrt(eps[0]*eps[0] + eps[1]*eps[1] + eps[2]*eps[2]);
		ScalarField Ecomb = 0.5*((dielectricEval->alpha-3.) + E);
		ScalarField epsByE = inv(E) * (Ecomb + sqrt(Ecomb*Ecomb + 3.*E));
		eps *= epsByE; //enhancement due to correlations
		//collect:
		setMuEps(state, mu, clone(mu), eps);
	}
	
	this->rhoExplicitTilde = rhoExplicitTilde; zeroNyquist(this->rhoExplicitTilde);
	this->nCavity = I(nCavityTilde + getFullCore());
	
	updateCavity();
	
	if(setPhiFromState)
	{	ScalarFieldMuEps gradUnused;
		ScalarFieldTilde phiFluidTilde;
		(*this)(state, gradUnused, &phiFluidTilde);
		linearPCM->state = phiFluidTilde + coulomb(rhoExplicitTilde);
	}
}

double NonlinearPCM::operator()(const ScalarFieldMuEps& state, ScalarFieldMuEps& Adiel_state,
	ScalarFieldTilde* Adiel_rhoExplicitTilde, ScalarFieldTilde* Adiel_nCavityTilde, IonicGradient* forces) const
{
	EnergyComponents& Adiel = ((NonlinearPCM*)this)->Adiel;
	ScalarFieldArray Adiel_shape; if(Adiel_nCavityTilde) nullToZero(Adiel_shape, gInfo, shape.size());
	
	ScalarFieldTilde rhoFluidTilde;
	ScalarField muPlus, muMinus, Adiel_muPlus, Adiel_muMinus;
	initZero(Adiel_muPlus, gInfo); initZero(Adiel_muMinus, gInfo);
	double mu0 = 0., Qexp = 0., Adiel_Qexp = 0.;
	if(screeningEval)
	{	//Get neutrality Lagrange multiplier:
		muPlus = getMuPlus(state);
		muMinus = getMuMinus(state);
		Qexp = integral(rhoExplicitTilde);
		mu0 = screeningEval->neutralityConstraint(muPlus, muMinus, shape.back(), Qexp);
		//Compute ionic free energy and bound charge
		ScalarField Aout, rhoIon;
		initZero(Aout, gInfo);
		initZero(rhoIon, gInfo);
		callPref(screeningEval->freeEnergy)(gInfo.nr, mu0, muPlus->dataPref(), muMinus->dataPref(), shape.back()->dataPref(),
			rhoIon->dataPref(), Aout->dataPref(), Adiel_muPlus->dataPref(), Adiel_muMinus->dataPref(), Adiel_shape.size() ? Adiel_shape.back()->dataPref() : 0);
		Adiel["Akappa"] = integral(Aout);
		rhoFluidTilde += J(rhoIon); //include bound charge due to ions
	}
	
	//Compute the dielectric free energy and bound charge:
	VectorField eps = getEps(state), Adiel_eps;
	{	ScalarField Aout; VectorField p;
		initZero(Aout, gInfo);
		nullToZero(p, gInfo);
		nullToZero(Adiel_eps, gInfo);
		callPref(dielectricEval->freeEnergy)(gInfo.nr, eps.const_dataPref(), shape[0]->dataPref(),
			p.dataPref(), Aout->dataPref(), Adiel_eps.dataPref(), Adiel_shape.size() ? Adiel_shape[0]->dataPref() : 0);
		Adiel["Aeps"] = integral(Aout);
		rhoFluidTilde -= divergence(J(p)); //include bound charge due to dielectric
	} //scoped to automatically deallocate temporaries
	
	//Compute the electrostatic terms:
	ScalarFieldTilde phiFluidTilde = coulomb(rhoFluidTilde);
	ScalarFieldTilde phiExplicitTilde = coulomb(rhoExplicitTilde);
	Adiel["Coulomb"] = dot(rhoFluidTilde, O(0.5*phiFluidTilde + phiExplicitTilde));
	
	if(screeningEval)
	{	//Propagate gradients from rhoIon to mu, shape
		ScalarField Adiel_rhoIon = I(phiFluidTilde+phiExplicitTilde);
		callPref(screeningEval->convertDerivative)(gInfo.nr, mu0, muPlus->dataPref(), muMinus->dataPref(), shape.back()->dataPref(),
			Adiel_rhoIon->dataPref(), Adiel_muPlus->dataPref(), Adiel_muMinus->dataPref(), Adiel_shape.size() ? Adiel_shape.back()->dataPref() : 0);
		//Propagate gradients from mu0 to mu, shape, Qexp:
		double Adiel_mu0 = integral(Adiel_muPlus) + integral(Adiel_muMinus), mu0_Qexp;
		ScalarField mu0_muPlus, mu0_muMinus, mu0_shape;
		screeningEval->neutralityConstraint(muPlus, muMinus, shape.back(), Qexp, &mu0_muPlus, &mu0_muMinus, &mu0_shape, &mu0_Qexp);
		Adiel_muPlus += Adiel_mu0 * mu0_muPlus;
		Adiel_muMinus += Adiel_mu0 * mu0_muMinus;
		if(Adiel_shape.size()) Adiel_shape.back() += Adiel_mu0 * mu0_shape;
		Adiel_Qexp = Adiel_mu0 * mu0_Qexp;
	}
	
	//Propagate gradients from p to eps, shape
	{	VectorField Adiel_p = I(gradient(phiFluidTilde+phiExplicitTilde)); //Because dagger(-divergence) = gradient
		callPref(dielectricEval->convertDerivative)(gInfo.nr, eps.const_dataPref(), shape[0]->dataPref(),
			Adiel_p.const_dataPref(), Adiel_eps.dataPref(), Adiel_shape.size() ? Adiel_shape[0]->dataPref() : 0);
	}
	
	//Optional outputs:
	if(Adiel_rhoExplicitTilde)
	{	if(fsp.nonlinearSCF && !useGummel())
		{	//Non-variational version (from inner solve):
			*Adiel_rhoExplicitTilde = linearPCM->state - phiExplicitTilde;
		}
		else
		{	//Variational version (from bound charge with neutrality constraint):
			*Adiel_rhoExplicitTilde = phiFluidTilde;
			(*Adiel_rhoExplicitTilde)->setGzero(Adiel_Qexp);
		}
	}
	if(Adiel_nCavityTilde)
	{	ScalarField Adiel_nCavity;
		propagateCavityGradients(Adiel_shape, Adiel_nCavity, *Adiel_rhoExplicitTilde, forces);
		*Adiel_nCavityTilde = J(Adiel_nCavity);
	}
	
	//Collect energy and gradient pieces:
	setMuEps(Adiel_state, Adiel_muPlus, Adiel_muMinus, Adiel_eps);
	Adiel_state *= gInfo.dV; //converts variational derivative to total derivative
	return Adiel;
}

void NonlinearPCM::minimizeFluid()
{	if(fsp.nonlinearSCF)
	{	clearState();
		Pulay<ScalarFieldTilde>::minimize(compute(0,0));
	}
	else
		Minimizable<ScalarFieldMuEps>::minimize(e.fluidMinParams);
}

void NonlinearPCM::loadState(const char* filename)
{	nullToZero(state, gInfo);
	state.loadFromFile(filename);
}

void NonlinearPCM::saveState(const char* filename) const
{	if(mpiWorld->isHead()) state.saveToFile(filename);
}

double NonlinearPCM::get_Adiel_and_grad_internal(ScalarFieldTilde& Adiel_rhoExplicitTilde, ScalarFieldTilde& Adiel_nCavityTilde, IonicGradient* extraForces) const
{	ScalarFieldMuEps Adiel_state;
	double A = (*this)(state, Adiel_state, &Adiel_rhoExplicitTilde, &Adiel_nCavityTilde, extraForces);
	accumExtraForces(extraForces, Adiel_nCavityTilde);
	return A;
}

void NonlinearPCM::step(const ScalarFieldMuEps& dir, double alpha)
{	::axpy(alpha, dir, state);
}

double NonlinearPCM::compute(ScalarFieldMuEps* grad, ScalarFieldMuEps* Kgrad)
{	ScalarFieldMuEps gradUnused;
	double E = (*this)(state, grad ? *grad : gradUnused);
	//Compute preconditioned gradient:
	if(Kgrad)
	{	const ScalarFieldMuEps& in = grad ? *grad : gradUnused;
		double dielPrefac = 1./(gInfo.dV * dielectricEval->NT);
		double ionsPrefac = screeningEval ? 1./(gInfo.dV * screeningEval->NT) : 0.;
		setMuEps(*Kgrad,
			ionsPrefac * I(preconditioner*J(getMuPlus(in))),
			ionsPrefac * I(preconditioner*J(getMuMinus(in))),
			dielPrefac * getEps(in));
	}
	return E;
}


void NonlinearPCM::dumpDensities(const char* filenamePattern) const
{	PCM::dumpDensities(filenamePattern);

	//Output dielectric bound charge:
	string filename;
	{	ScalarField Aout; initZero(Aout, gInfo);
		VectorField p; nullToZero(p, gInfo);
		VectorField Adiel_eps; nullToZero(Adiel_eps, gInfo);
		callPref(dielectricEval->freeEnergy)(gInfo.nr, getEps(state).const_dataPref(), shape[0]->dataPref(),
			p.dataPref(), Aout->dataPref(), Adiel_eps.dataPref(), 0);
		ScalarField rhoDiel = -divergence(p); //include bound charge due to dielectric
		FLUID_DUMP(rhoDiel, "RhoDiel");
	}
	
	//Output ionic bound charge (if any):
	if(screeningEval)
	{	ScalarField Nplus, Nminus;
		{	ScalarField muPlus = getMuPlus(state);
			ScalarField muMinus = getMuMinus(state);
			double Qexp = integral(rhoExplicitTilde);
			double mu0 = screeningEval->neutralityConstraint(muPlus, muMinus, shape.back(), Qexp);
			Nplus = ionNbulk * shape.back() * (fsp.linearScreening ? 1.+(mu0+muPlus) : exp(mu0+muPlus));
			Nminus = ionNbulk * shape.back() * (fsp.linearScreening ? 1.-(mu0+muMinus) : exp(-(mu0+muMinus)));
		}
		FLUID_DUMP(Nplus, "N+");
		FLUID_DUMP(Nminus, "N-");
		FLUID_DUMP(ionZ*(Nplus-Nminus), "RhoIon");
	}
}

//--------- Interface for Pulay<ScalarFieldTilde> ---------

double NonlinearPCM::cycle(double dEprev, std::vector<double>& extraValues)
{	//Update epsilon / kappaSq based on current phi:
	phiToState(false);
	//Inner linear solve
	FILE*& fpLog = ((MinimizeParams&)e.fluidMinParams).fpLog;
	fpLog = nullLog; //disable iteration log from LinearPCM
	linearPCM->minimizeFluid();
	fpLog = globalLog; //restore usual iteration log
	//Update state from new phi:
	phiToState(true);
	return compute(0,0);
}


void NonlinearPCM::readVariable(ScalarFieldTilde& X, FILE* fp) const
{	nullToZero(X, gInfo);
	loadRawBinary(X, fp);
}

void NonlinearPCM::writeVariable(const ScalarFieldTilde& X, FILE* fp) const
{	saveRawBinary(X, fp);
}

ScalarFieldTilde NonlinearPCM::getVariable() const
{	return clone(linearPCM->state);
}

void NonlinearPCM::setVariable(const ScalarFieldTilde& X)
{	linearPCM->state = clone(X);
}

ScalarFieldTilde NonlinearPCM::precondition(const ScalarFieldTilde& X) const
{	return fsp.scfParams.mixFraction * X;	
}

ScalarFieldTilde NonlinearPCM::applyMetric(const ScalarFieldTilde& X) const
{	return (*metric) * X;
}

void NonlinearPCM::phiToState(bool setState)
{	//Initialize inputs:
	const ScalarField phi = I(linearPCM->state);
	const VectorField Dphi = I(gradient(linearPCM->state));
	//Prepare outputs:
	ScalarField epsilon, kappaSq;
	if(!setState)
	{	nullToZero(epsilon, gInfo);
		if(screeningEval)
			nullToZero(kappaSq, gInfo);
	}
	VectorField eps = getEps(state);
	ScalarField& muPlus = getMuPlus(state);
	ScalarField& muMinus = getMuMinus(state);
	//Calculate eps/mu or epsilon/kappaSq as needed:
	vector3<double*> vecDataUnused(0,0,0); double* dataUnused=0;
	callPref(dielectricEval->phiToState)(gInfo.nr, Dphi.dataPref(), shape[0]->dataPref(), gLookup, setState,
		setState ? eps.dataPref() : vecDataUnused,
		setState ? dataUnused : epsilon->dataPref() );
	if(screeningEval)
		callPref(screeningEval->phiToState)(gInfo.nr, phi->dataPref(), shape.back()->dataPref(), xLookup, setState,
			setState ? muPlus->dataPref() : dataUnused,
			setState ? muMinus->dataPref() : dataUnused, 
			setState ? dataUnused : kappaSq->dataPref() );
	//Save to global state or linearPCM as required:
	if(setState)
		setMuEps(state, muPlus, muMinus, eps);
	else
		linearPCM->override(epsilon, kappaSq);
}
