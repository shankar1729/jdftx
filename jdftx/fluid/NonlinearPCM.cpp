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
#include <electronic/operators.h>
#include <core/DataIO.h>
#include <core/Util.h>

//Utility functions to extract/set the members of a MuEps
inline DataRptr& getMuPlus(DataRMuEps& X) { return X[0]; }
inline const DataRptr& getMuPlus(const DataRMuEps& X) { return X[0]; }
inline DataRptr& getMuMinus(DataRMuEps& X) { return X[1]; }
inline const DataRptr& getMuMinus(const DataRMuEps& X) { return X[1]; }
inline DataRptrVec getEps(DataRMuEps& X) { return DataRptrVec(X.component+2); }
inline const DataRptrVec getEps(const DataRMuEps& X) { return DataRptrVec(X.component+2); }
inline void setMuEps(DataRMuEps& mueps, DataRptr muPlus, DataRptr muMinus, DataRptrVec eps) { mueps[0]=muPlus; mueps[1]=muMinus; for(int k=0; k<3; k++) mueps[k+2]=eps[k]; }


inline double setPreconditioner(double G, double kappaSqByEpsilon, double muByEpsSq)
{	double G2 = G*G;
	double den = G2 + kappaSqByEpsilon;
	return den ? muByEpsSq*G2/(den*den) : 0.;
}

NonlinearPCM::NonlinearPCM(const Everything& e, const FluidSolverParams& fsp)
: PCM(e, fsp), pMol(0.), ionNbulk(0.), ionZ(0.), screeningEval(0), dielectricEval(0)
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
	
	//Initialize preconditioner (for mu channel):
	double muByEps = (ionZ/pMol) * (1.-dielectricEval->alpha/3); //relative scale between mu and eps
	preconditioner.init(0, 0.02, gInfo.GmaxGrid, setPreconditioner, k2factor/epsBulk, muByEps*muByEps);
}

NonlinearPCM::~NonlinearPCM()
{	preconditioner.free();
	delete dielectricEval;
	if(screeningEval)
		delete screeningEval;
}


void NonlinearPCM::set_internal(const DataGptr& rhoExplicitTilde, const DataGptr& nCavityTilde)
{	
	//Initialize state if required:
	if(!state)
	{	logPrintf("Initializing state of NonlinearPCM using a similar LinearPCM:\n");
		FILE*& fpLog = ((MinimizeParams&)e.fluidMinParams).fpLog;
		fpLog = fopen("/dev/null", "w"); //disable iteration log from LinearPCM
		LinearPCM linearPCM(e, fsp);
		linearPCM.set_internal(rhoExplicitTilde, nCavityTilde);
		linearPCM.minimizeFluid();
		fclose(fpLog);
		fpLog = globalLog; //retsore usual iteration log
		//Guess nonlinear states based on the electrostatic potential of the linear version:
		//mu:
		DataRptr mu;
		if(screeningEval && screeningEval->linear)
		{	mu = (-ionZ/fsp.T) * I(linearPCM.state);
			mu -= integral(mu)/gInfo.detR; //project out G=0
		}
		else initZero(mu, gInfo); //initialization logic does not work well with hard sphere limit
		//eps:
		DataRptrVec eps = (-pMol/fsp.T) * I(gradient(linearPCM.state));
		DataRptr E = sqrt(eps[0]*eps[0] + eps[1]*eps[1] + eps[2]*eps[2]);
		DataRptr Ecomb = 0.5*((dielectricEval->alpha-3.) + E);
		DataRptr epsByE = inv(E) * (Ecomb + sqrt(Ecomb*Ecomb + 3.*E));
		eps *= epsByE; //enhancement due to correlations
		//collect:
		setMuEps(state, mu, clone(mu), eps);
	}
	
	this->rhoExplicitTilde = rhoExplicitTilde; zeroNyquist(this->rhoExplicitTilde);
	this->nCavity = I(nCavityTilde);

	updateCavity();
}

double NonlinearPCM::operator()(const DataRMuEps& state, DataRMuEps& Adiel_state, DataGptr* Adiel_rhoExplicitTilde, DataGptr* Adiel_nCavityTilde) const
{
	EnergyComponents& Adiel = ((NonlinearPCM*)this)->Adiel;
	DataRptr Adiel_shape; if(Adiel_nCavityTilde) nullToZero(Adiel_shape, gInfo);
	
	DataGptr rhoFluidTilde;
	DataRptr muPlus, muMinus, Adiel_muPlus, Adiel_muMinus;
	initZero(Adiel_muPlus, gInfo); initZero(Adiel_muMinus, gInfo);
	double mu0 = 0., Qexp = 0., Adiel_Qexp = 0.;
	if(screeningEval)
	{	//Get neutrality Lagrange multiplier:
		muPlus = getMuPlus(state);
		muMinus = getMuMinus(state);
		Qexp = integral(rhoExplicitTilde);
		mu0 = screeningEval->neutralityConstraint(muPlus, muMinus, shape, Qexp);
		//Compute ionic free energy and bound charge
		DataRptr Aout, rhoIon;
		initZero(Aout, gInfo);
		initZero(rhoIon, gInfo);
		callPref(screeningEval->freeEnergy)(gInfo.nr, mu0, muPlus->dataPref(), muMinus->dataPref(), shape->dataPref(),
			rhoIon->dataPref(), Aout->dataPref(), Adiel_muPlus->dataPref(), Adiel_muMinus->dataPref(), Adiel_shape ? Adiel_shape->dataPref() : 0);
		Adiel["Akappa"] = integral(Aout);
		rhoFluidTilde += J(rhoIon); //include bound charge due to ions
	}
	
	//Compute the dielectric free energy and bound charge:
	DataRptrVec eps = getEps(state), Adiel_eps;
	{	DataRptr Aout; DataRptrVec p;
		initZero(Aout, gInfo);
		nullToZero(p, gInfo);
		nullToZero(Adiel_eps, gInfo);
		callPref(dielectricEval->freeEnergy)(gInfo.nr, eps.const_dataPref(), shape->dataPref(),
			p.dataPref(), Aout->dataPref(), Adiel_eps.dataPref(), Adiel_shape ? Adiel_shape->dataPref() : 0);
		Adiel["Aeps"] = integral(Aout);
		rhoFluidTilde -= divergence(J(p)); //include bound charge due to dielectric
	} //scoped to automatically deallocate temporaries
	
	//Compute the electrostatic terms:
	DataGptr phiFluidTilde = coulomb(rhoFluidTilde);
	DataGptr phiExplicitTilde = coulomb(rhoExplicitTilde);
	Adiel["Coulomb"] = dot(rhoFluidTilde, O(0.5*phiFluidTilde + phiExplicitTilde));
	
	if(screeningEval)
	{	//Propagate gradients from rhoIon to mu, shape
		DataRptr Adiel_rhoIon = I(phiFluidTilde+phiExplicitTilde);
		callPref(screeningEval->convertDerivative)(gInfo.nr, mu0, muPlus->dataPref(), muMinus->dataPref(), shape->dataPref(),
			Adiel_rhoIon->dataPref(), Adiel_muPlus->dataPref(), Adiel_muMinus->dataPref(), Adiel_shape ? Adiel_shape->dataPref() : 0);
		//Propagate gradients from mu0 to mu, shape, Qexp:
		double Adiel_mu0 = integral(Adiel_muPlus) + integral(Adiel_muMinus), mu0_Qexp;
		DataRptr mu0_muPlus, mu0_muMinus, mu0_shape;
		screeningEval->neutralityConstraint(muPlus, muMinus, shape, Qexp, &mu0_muPlus, &mu0_muMinus, &mu0_shape, &mu0_Qexp);
		Adiel_muPlus += Adiel_mu0 * mu0_muPlus;
		Adiel_muMinus += Adiel_mu0 * mu0_muMinus;
		if(Adiel_shape) Adiel_shape += Adiel_mu0 * mu0_shape;
		Adiel_Qexp = Adiel_mu0 * mu0_Qexp;
	}
	
	//Propagate gradients from p to eps, shape
	{	DataRptrVec Adiel_p = I(gradient(phiFluidTilde+phiExplicitTilde), true); //Because dagger(-divergence) = gradient
		callPref(dielectricEval->convertDerivative)(gInfo.nr, eps.const_dataPref(), shape->dataPref(),
			Adiel_p.const_dataPref(), Adiel_eps.dataPref(), Adiel_shape ? Adiel_shape->dataPref() : 0);
	}
	
	//Optional outputs:
	if(Adiel_rhoExplicitTilde)
	{	*Adiel_rhoExplicitTilde = phiFluidTilde;
		(*Adiel_rhoExplicitTilde)->setGzero(Adiel_Qexp);
	}
	if(Adiel_nCavityTilde)
	{	DataRptr Adiel_nCavity;
		propagateCavityGradients(Adiel_shape, Adiel_nCavity, *Adiel_rhoExplicitTilde);
		*Adiel_nCavityTilde = J(Adiel_nCavity);
	}
	
	//Collect energy and gradient pieces:
	setMuEps(Adiel_state, Adiel_muPlus, Adiel_muMinus, Adiel_eps);
	Adiel_state *= gInfo.dV; //converts variational derivative to total derivative
	return Adiel;
}

void NonlinearPCM::minimizeFluid()
{	minimize(e.fluidMinParams);
}

void NonlinearPCM::loadState(const char* filename)
{	nullToZero(state, gInfo);
	state.loadFromFile(filename);
}

void NonlinearPCM::saveState(const char* filename) const
{	if(mpiUtil->isHead()) state.saveToFile(filename);
}

double NonlinearPCM::get_Adiel_and_grad_internal(DataGptr& Adiel_rhoExplicitTilde, DataGptr& Adiel_nCavityTilde, IonicGradient* extraForces) const
{	DataRMuEps Adiel_state;
	double A = (*this)(state, Adiel_state, &Adiel_rhoExplicitTilde, &Adiel_nCavityTilde);
	setExtraForces(extraForces, Adiel_nCavityTilde);
	return A;
}

void NonlinearPCM::step(const DataRMuEps& dir, double alpha)
{	axpy(alpha, dir, state);
}

double NonlinearPCM::compute(DataRMuEps* grad)
{	DataRMuEps gradUnused;
	return (*this)(state, grad ? *grad : gradUnused);
}

DataRMuEps NonlinearPCM::precondition(const DataRMuEps& in)
{	DataRMuEps out;
	double dielPrefac = 1./(gInfo.dV * dielectricEval->NT);
	double ionsPrefac = screeningEval ? 1./(gInfo.dV * screeningEval->NT) : 0.;
	setMuEps(out,
		ionsPrefac * I(preconditioner*J(getMuPlus(in))),
		ionsPrefac * I(preconditioner*J(getMuMinus(in))),
		dielPrefac * getEps(in));
	return out;
}


void NonlinearPCM::dumpDensities(const char* filenamePattern) const
{	PCM::dumpDensities(filenamePattern);

	if(screeningEval)
	{
		DataRptr Nplus, Nminus;
		{	DataRptr muPlus = getMuPlus(state);
			DataRptr muMinus = getMuMinus(state);
			double Qexp = integral(rhoExplicitTilde);
			double mu0 = screeningEval->neutralityConstraint(muPlus, muMinus, shape, Qexp);
			Nplus = ionNbulk * shape * (fsp.linearScreening ? 1.+(mu0+muPlus) : exp(mu0+muPlus));
			Nminus = ionNbulk * shape * (fsp.linearScreening ? 1.-(mu0+muMinus) : exp(-(mu0+muMinus)));
		}
		string filename;
		FLUID_DUMP(Nplus, "N+");
		FLUID_DUMP(Nminus, "N-");
	}
}
