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
#include <electronic/NonlinearPCM.h>
#include <electronic/PCM_internal.h>
#include <core/Units.h>
#include <core/DataIO.h>
#include <core/Operators.h>
#include <electronic/operators.h>
#include <core/Util.h>
#include <core/DataMultiplet.h>

//Utility functions to extract/set the members of a MuEps
inline DataRptr& getMu(DataRMuEps& X) { return X[0]; }
inline const DataRptr& getMu(const DataRMuEps& X) { return X[0]; }
inline DataRptrVec getEps(DataRMuEps& X) { return DataRptrVec(X.component+1); }
inline const DataRptrVec getEps(const DataRMuEps& X) { return DataRptrVec(X.component+1); }
inline void setMuEps(DataRMuEps& mueps, DataRptr mu, DataRptrVec eps) { mueps[0]=mu; for(int k=0; k<3; k++) mueps[k+1]=eps[k]; }


inline void setPreconditioner(int i, double G2, double* preconditioner, double epsilon, double kappaSq)
{	preconditioner[i] = (i || kappaSq) ? (4*M_PI)/(epsilon*G2 + kappaSq) : 0.;
}


NonlinearPCM::NonlinearPCM(const Everything& e, const FluidSolverParams& fsp)
: FluidSolver(e), params(fsp), preconditioner(e.gInfo)
{
	//Initialize dielectric evaluation class:
	dielectricEval = new NonlinearPCMeval::Dielectric(params.linearDielectric,
		params.T, params.Nbulk, params.pMol, params.epsilonBulk, params.epsInf);
	
	//Optionally initialize screening evaluation class:
	screeningEval = 0;
	if(params.ionicConcentration)
		screeningEval = new NonlinearPCMeval::Screening(params.linearScreening,
			params.T, params.ionicConcentration, params.ionicZelectrolyte, params.epsilonBulk);
		
	//Initialize preconditioner:
	double bulkKappaSq = (8*M_PI/params.T) * params.ionicConcentration * pow(params.ionicZelectrolyte,2);
	applyFuncGsq(e.gInfo, setPreconditioner, preconditioner.data, params.epsilonBulk, bulkKappaSq);
	preconditioner.set();
}

NonlinearPCM::~NonlinearPCM()
{
	delete dielectricEval;
	if(screeningEval)
		delete screeningEval;
}


void NonlinearPCM::set(const DataGptr& rhoExplicitTilde, const DataGptr& nCavityTilde)
{	this->rhoExplicitTilde = rhoExplicitTilde;
	this->nCavity = I(nCavityTilde);
	//Initialize point-dipole density
	nullToZero(shape,e.gInfo);
	pcmShapeFunc(nCavity, shape, params.nc, params.sigma);

	//Compute the cavitation energy and gradient
	Acavity = cavitationEnergyAndGrad(shape, Acavity_shape, params.cavityTension, params.cavityPressure);
	
	//Initialize state if required:
	if(!state) nullToZero(state, e.gInfo);
}

double NonlinearPCM::operator()(const DataRMuEps& state, DataRMuEps& Adiel_state, DataGptr* Adiel_rhoExplicitTilde, DataGptr* Adiel_nCavityTilde) const
{
	DataRptr Adiel_shape; if(Adiel_nCavityTilde) nullToZero(Adiel_shape, e.gInfo);
	
	DataGptr rhoFluidTilde;
	DataRptr mu, Adiel_mu; initZero(Adiel_mu, e.gInfo);
	double Akappa = 0., mu0 = 0., Qexp = 0., Adiel_Qexp = 0.;
	if(screeningEval)
	{	//Get neutrality Lagrange multiplier:
		mu = getMu(state);
		Qexp = integral(rhoExplicitTilde);
		mu0 = screeningEval->neutralityConstraint(mu, shape, Qexp);
		//Compute ionic free energy and bound charge
		DataRptr Aout, rhoIon;
		initZero(Aout, e.gInfo);
		initZero(rhoIon, e.gInfo);
		callPref(screeningEval->freeEnergy)(e.gInfo.nr, mu0, mu->dataPref(), shape->dataPref(),
			rhoIon->dataPref(), Aout->dataPref(), Adiel_mu->dataPref(), Adiel_shape ? Adiel_shape->dataPref() : 0);
		Akappa = integral(Aout);
		rhoFluidTilde += J(rhoIon); //include bound charge due to ions
	}
	
	//Compute the dielectric free energy and bound charge:
	DataRptrVec eps = getEps(state), Adiel_eps; double Aeps = 0.;
	{	DataRptr Aout; DataRptrVec p;
		initZero(Aout, e.gInfo);
		nullToZero(p, e.gInfo);
		nullToZero(Adiel_eps, e.gInfo);
		callPref(dielectricEval->freeEnergy)(e.gInfo.nr, eps.const_dataPref(), shape->dataPref(),
			p.dataPref(), Aout->dataPref(), Adiel_eps.dataPref(), Adiel_shape ? Adiel_shape->dataPref() : 0);
		Aeps = integral(Aout);
		rhoFluidTilde -= divergence(J(p)); //include bound charge due to dielectric
	} //scoped to automatically deallocate temporaries
	
	//Compute the electrostatic terms:
	DataGptr phiFluidTilde = (*e.coulomb)(rhoFluidTilde);
	DataGptr phiExplicitTilde = (*e.coulomb)(rhoExplicitTilde);
	double U = dot(rhoFluidTilde, O(0.5*phiFluidTilde + phiExplicitTilde));
	
	if(screeningEval)
	{	//Propagate gradients from rhoIon to mu, shape
		DataRptr Adiel_rhoIon = I(phiFluidTilde+phiExplicitTilde);
		callPref(screeningEval->convertDerivative)(e.gInfo.nr, mu0, mu->dataPref(), shape->dataPref(),
			Adiel_rhoIon->dataPref(), Adiel_mu->dataPref(), Adiel_shape ? Adiel_shape->dataPref() : 0);
		//Propagate gradients from mu0 to mu, shape, Qexp:
		double Adiel_mu0 = integral(Adiel_mu), mu0_Qexp; DataRptr mu0_mu, mu0_shape;
		screeningEval->neutralityConstraint(mu, shape, Qexp, &mu0_mu, &mu0_shape, &mu0_Qexp);
		Adiel_mu += Adiel_mu0 * mu0_mu;
		if(Adiel_shape) Adiel_shape += Adiel_mu0 * mu0_shape;
		Adiel_Qexp = Adiel_mu0 * mu0_Qexp;
	}
	
	//Propagate gradients from p to eps, shape
	{	DataRptrVec Adiel_p = I(gradient(phiFluidTilde+phiExplicitTilde), true); //Because dagger(-divergence) = gradient
		callPref(dielectricEval->convertDerivative)(e.gInfo.nr, eps.const_dataPref(), shape->dataPref(),
			Adiel_p.const_dataPref(), Adiel_eps.dataPref(), Adiel_shape ? Adiel_shape->dataPref() : 0);
	}
	
	//Optional outputs:
	if(Adiel_rhoExplicitTilde)
	{	*Adiel_rhoExplicitTilde = phiFluidTilde;
		(*Adiel_rhoExplicitTilde)->setGzero(Adiel_Qexp);
	}
	if(Adiel_nCavityTilde)
	{	Adiel_shape  += Acavity_shape; // Add the cavitation contribution to Adiel_shape
		DataRptr Adiel_nCavity(DataR::alloc(e.gInfo));
		pcmShapeFunc_grad(nCavity, Adiel_shape, Adiel_nCavity, params.nc, params.sigma);
		*Adiel_nCavityTilde = J(Adiel_nCavity);
	}
	
	//Collect energy and gradient pieces:
	setMuEps(Adiel_state, Adiel_mu, Adiel_eps);
	Adiel_state *= e.gInfo.dV; //converts variational derivative to total derivative
	return Akappa + Aeps + U + Acavity;
}

void NonlinearPCM::minimizeFluid()
{	minimize(e.fluidMinParams);
}

void NonlinearPCM::loadState(const char* filename)
{	nullToZero(state, e.gInfo);
	state.loadFromFile(filename);
}

void NonlinearPCM::saveState(const char* filename) const
{	state.saveToFile(filename);
}

double NonlinearPCM::get_Adiel_and_grad(DataGptr& Adiel_rhoExplicitTilde, DataGptr& Adiel_nCavityTilde, IonicGradient& extraForces)
{	DataRMuEps Adiel_state;
	return (*this)(state, Adiel_state, &Adiel_rhoExplicitTilde, &Adiel_nCavityTilde);
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
	setMuEps(out, I(preconditioner*J(getMu(in))), getEps(in));
	return out;
}

void NonlinearPCM::dumpDensities(const char* filenamePattern) const
{	
	string filename(filenamePattern);
	filename.replace(filename.find("%s"), 2, "Shape");
	logPrintf("Dumping '%s'... ", filename.c_str());  logFlush();
	saveRawBinary(shape, filename.c_str());
	logPrintf("done.\n"); logFlush();
	
	/*
	if (params.ionicConcentration)
	{
		double nPosGzero = 0.0, nNegGzero = 0.0, rhoExplicitGzero = 0.0;
		complex* rhoPtr = rhoExplicitTilde->data();
		rhoExplicitGzero = rhoPtr[0].real()*e.gInfo.detR;		
		
		const DataRptr& mu = getMu(state);	
		DataRptr muEff;	 nullToZero(muEff,e.gInfo);
		DataRptr nPos;  nullToZero(nPos,e.gInfo);
		DataRptr nNeg;  nullToZero(nNeg,e.gInfo);
		
		if (params.linearScreening)
		{	nNegGzero=params.ionicConcentration*dot(shape,1.0-mu+0.5*mu*mu)*e.gInfo.dV;
			nPosGzero=params.ionicConcentration*dot(shape,1.0+mu+0.5*mu*mu)*e.gInfo.dV;
		}
		else
		{	nNegGzero=params.ionicConcentration*dot(shape,exp(-mu))*e.gInfo.dV;
			nPosGzero=params.ionicConcentration*dot(shape,exp(mu))*e.gInfo.dV;
		}

		double fmu = log((-rhoExplicitGzero+sqrt(rhoExplicitGzero*rhoExplicitGzero+4.0*nNegGzero*nPosGzero))/(2.0*nNegGzero));
		muEff = mu-fmu;
		
		threadedLoop(calc_Nion, e.gInfo.nr, muEff->data(), shape->data(), nNeg->data(), nPos->data(), &params);
				
		filename = filenamePattern;
		filename.replace(filename.find("%s"), 2, "N+");
		logPrintf("Dumping '%s'... ", filename.c_str());  logFlush();
		saveRawBinary(nPos, filename.c_str());
		logPrintf("done.\n"); logFlush();
		
		filename = filenamePattern;
		filename.replace(filename.find("%s"), 2, "N-");
		logPrintf("Dumping '%s'... ", filename.c_str());  logFlush();
		saveRawBinary(nNeg, filename.c_str());
		logPrintf("done.\n"); logFlush();
	}
	*/
}
