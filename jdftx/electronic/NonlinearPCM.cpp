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

void randomize(DataGptr& X)
{	DataRptr tmp(DataR::alloc(X->gInfo));
	initRandom(tmp);
	X = J(tmp);
}


NonlinearPCM::NonlinearPCM(const Everything& e, const FluidSolverParams& fsp)
: FluidSolver(e), params(fsp), preconditioner(e.gInfo)
{
	//Initialize dielectric evaluation class:
	dielectricEval = new NonlinearPCMeval::Dielectric(params.linearDielectric,
		params.T, params.Nbulk, params.pMol, params.epsilonBulk, params.epsInf);
	dielectricEval->initLookupTables();
	
	//Optionally initialize screening evaluation class:
	screeningEval = 0;
	if(params.ionicConcentration)
		screeningEval = new NonlinearPCMeval::Screening(params.linearScreening,
			params.T, params.ionicConcentration, params.ionicZelectrolyte);
}

NonlinearPCM::~NonlinearPCM()
{
	dielectricEval->freeLookupTables();
	delete dielectricEval;
	if(screeningEval)
		delete screeningEval;
}


inline void setPreconditioner(int i, double G2, double* preconditioner, double epsilon, double kappaSq)
{	preconditioner[i] = (i || kappaSq) ? (4*M_PI)/(epsilon*G2 + kappaSq) : 0.;
}

// inline void setPreconditioner(int i, double G2, double* preconditioner, double epsilon, double kappaSq)
// {	preconditioner[i] = (i || kappaSq) ? sqrt((4*M_PI)/((epsilon-1)*G2 + kappaSq)) : 0.;
// }

void NonlinearPCM::set(const DataGptr& rhoExplicitTilde, const DataGptr& nCavityTilde)
{	this->rhoExplicitTilde = rhoExplicitTilde;
	this->nCavity = I(nCavityTilde);
	//Initialize point-dipole density
	nullToZero(shape,e.gInfo);
	pcmShapeFunc(nCavity, shape, params.nc, params.sigma);

	// Compute the cavitation energy and gradient
	Acavity = cavitationEnergyAndGrad(shape, Acavity_shape, params.cavityTension, params.cavityPressure);
	
	//Update preconditioner:
	double bulkKappaSq = (8*M_PI/params.T) * params.ionicConcentration * pow(params.ionicZelectrolyte,2);
// 	double fluidFraction = integral(shape)/e.gInfo.detR; //fraction of cell occupied by fluid
// 	double meanEpsilon = 1. + (params.epsilonBulk-1.)*fluidFraction;
// 	double meanKappaSq = fluidFraction * bulkKappaSq;
// 	logPrintf("meanEpsilon = %lf    meanKappaSq = %lf\n", meanEpsilon, meanKappaSq);
// 	applyFuncGsq(e.gInfo, setPreconditioner, preconditioner.data, meanEpsilon, meanKappaSq);
	applyFuncGsq(e.gInfo, setPreconditioner, preconditioner.data, params.epsilonBulk, bulkKappaSq);
	preconditioner.set();

	if(!state) initZero(state, e.gInfo); // = preconditioner * rhoExplicitTilde; //initialize using bulk linear-response approx
}

double NonlinearPCM::operator()(const DataGptr& phiTilde, DataGptr& Adiel_phiTilde, DataGptr* Adiel_rhoExplicitTilde, DataGptr* Adiel_nCavityTilde) const
{
	DataRptr Adiel_shape; if(Adiel_nCavityTilde) nullToZero(Adiel_shape, e.gInfo);
	
	//Compute the ionic free energy and bound charge:
	DataGptr rhoFluidTilde;
	DataRptr phi = I(phiTilde), Adiel_phi; double Akappa = 0.;
	initZero(Adiel_phi, e.gInfo);
	if(screeningEval)
	{	DataRptr Aout, rhoIon;
		initZero(Aout, e.gInfo);
		initZero(rhoIon, e.gInfo);
		callPref(screeningEval->freeEnergy)(e.gInfo.nr, phi->dataPref(), shape->dataPref(),
			rhoIon->dataPref(), Aout->dataPref(), Adiel_phi->dataPref(), Adiel_shape ? Adiel_shape->dataPref() : 0);
		Akappa = integral(Aout);
		rhoFluidTilde += J(rhoIon); //include bound charge due to ions
	}
	
	//Compute the dielectric free energy and bound charge:
	DataRptrVec gradPhi = I(gradient(phiTilde), true), Adiel_gradPhi; double Aeps = 0.;
	{	DataRptr Aout; DataRptrVec p;
		initZero(Aout, e.gInfo);
		nullToZero(p, e.gInfo);
		nullToZero(Adiel_gradPhi, e.gInfo);
		callPref(dielectricEval->freeEnergy)(e.gInfo.nr, gradPhi.const_dataPref(), shape->dataPref(),
			p.dataPref(), Aout->dataPref(), Adiel_gradPhi.dataPref(), Adiel_shape ? Adiel_shape->dataPref() : 0);
		Aeps = integral(Aout);
		rhoFluidTilde -= divergence(J(p)); //include bound charge due to dielectric
	} //scoped to automatically deallocate temporaries
	
	//Compute the electrostatic terms:
	DataGptr phiFluidTilde = (*e.coulomb)(rhoFluidTilde);
	DataGptr phiExplicitTilde = (*e.coulomb)(rhoExplicitTilde);
	double U = dot(rhoFluidTilde, O(0.5*phiFluidTilde + phiExplicitTilde));
	
	//Propagate gradients from rhoIon to phi, shape
	if(screeningEval)
	{	DataRptr Adiel_rhoIon = I(phiFluidTilde+phiExplicitTilde);
		callPref(screeningEval->convertDerivative)(e.gInfo.nr, phi->dataPref(), shape->dataPref(),
			Adiel_rhoIon->dataPref(), Adiel_phi->dataPref(), Adiel_shape ? Adiel_shape->dataPref() : 0);
	}
	
	//Propagate gradients from p to gradPhi, shape
	{	DataRptrVec Adiel_p = I(gradient(phiFluidTilde+phiExplicitTilde), true); //Because dagger(-divergence) = gradient
		callPref(dielectricEval->convertDerivative)(e.gInfo.nr, gradPhi.const_dataPref(), shape->dataPref(),
			Adiel_p.const_dataPref(), Adiel_gradPhi.dataPref(), Adiel_shape ? Adiel_shape->dataPref() : 0);
	}
	
	//Optional outputs:
	if(Adiel_rhoExplicitTilde)
	{	*Adiel_rhoExplicitTilde = phiFluidTilde;
	}
	if(Adiel_nCavityTilde)
	{	Adiel_shape  += Acavity_shape; // Add the cavitation contribution to Adiel_shape
		DataRptr Adiel_nCavity(DataR::alloc(e.gInfo));
		pcmShapeFunc_grad(nCavity, Adiel_shape, Adiel_nCavity, params.nc, params.sigma);
		*Adiel_nCavityTilde = J(Adiel_nCavity);
	}
	
	//Collect energy and gradient pieces:
	Adiel_phiTilde = e.gInfo.dV * (Idag(Adiel_phi) - divergence(Idag(Adiel_gradPhi))); //Because dagger(gradient) = -divergence
	return Akappa + Aeps + U + Acavity;
}

void NonlinearPCM::minimizeFluid()
{	minimize(e.fluidMinParams);
}

void NonlinearPCM::loadState(const char* filename)
{	DataRptr Istate(DataR::alloc(e.gInfo));
	loadRawBinary(Istate, filename);
	state = J(Istate);
}

void NonlinearPCM::saveState(const char* filename) const
{	saveRawBinary(I(state), filename);
}

double NonlinearPCM::get_Adiel_and_grad(DataGptr& Adiel_rhoExplicitTilde, DataGptr& Adiel_nCavityTilde, IonicGradient& extraForces)
{	DataGptr Adiel_state;
	return (*this)(state, Adiel_state, &Adiel_rhoExplicitTilde, &Adiel_nCavityTilde);
}

void NonlinearPCM::step(const DataGptr& dir, double alpha)
{	axpy(alpha, dir, state);
}

double NonlinearPCM::compute(DataGptr* grad)
{	DataGptr gradUnused;
	return (*this)(state, grad ? *grad : gradUnused);
}

DataGptr NonlinearPCM::precondition(const DataGptr& in)
{	return preconditioner * in;
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
