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

inline DataRptr& getMu(DataRMuEps& X) { return X[0]; }
inline const DataRptr& getMu(const DataRMuEps& X) { return X[0]; }
inline DataRptrVec getEps(DataRMuEps& X) { return DataRptrVec(X.component+1); }
inline const DataRptrVec getEps(const DataRMuEps& X) { return DataRptrVec(X.component+1); }
inline void setMuEps(DataRMuEps& mueps, DataRptr mu, DataRptrVec eps) { mueps[0]=mu; for(int k=0; k<3; k++) mueps[k+1]=eps[k]; }

//! Returns the cavitation energy contribution (volume and surface terms)
inline double cavitationEnergyAndGrad(const DataRptr& shape, DataRptr& grad, double cavityTension, double cavityPressure)
{
	DataRptrVec shape_x = gradient(shape);
	DataRptr surfaceDensity = sqrt(shape_x[0]*shape_x[0] + shape_x[1]*shape_x[1] + shape_x[2]*shape_x[2]);
	double surfaceArea = integral(surfaceDensity);
	double volume = integral(1.-shape);

	if(grad)
	{	DataRptr invSurfaceDensity = inv(surfaceDensity);
		grad += -cavityTension*divergence(shape_x*invSurfaceDensity); // Surface term
		grad += -cavityPressure;
	}
	
	return surfaceArea*cavityTension + volume*cavityPressure;
}

//Ratio of the fractional polarization (p/shape*pMol) to eps
//(with regularization around 0)
inline double polFrac(double x, double chiFactor, bool linearDielectric)
{	if(linearDielectric) return 1.0/3;
	else
	{	register double xSq=x*x;
		if(x<1e-1) return 1.0/3 + chiFactor + xSq*(-1.0/45 + xSq*(2.0/945 + xSq*(-1.0/4725)));
		return (x/tanh(x)-1)/xSq + chiFactor;
	}
}

inline double polFrac_x(double x, bool linearDielectric)
{	if(linearDielectric) return 0.0;
	else
	{	register double xSq=x*x;
		if(x<1e-1) return x*(-2.0/45 + xSq*(8.0/945 + xSq*(-6.0/4725)));
		return (2 - x/tanh(x) - pow(x/sinh(x),2))/(x*xSq);
	}
}

inline double polFrac_x_over_x(double x, bool linearDielectric)
{	if(linearDielectric) return 0.0;
	else
	{	register double xSq=x*x;
		if(x<1e-1) return -2.0/45 + xSq*(8.0/945 + xSq*(-6.0/4725));
		return (2 - x/tanh(x) - pow(x/sinh(x),2))/(xSq*xSq);
	}
}

inline double logsinch(double x, bool linearDielectric)
{	if(linearDielectric)return x*x/6;
	else
	{	register double xSq=x*x;
		if(x<1e-1) return xSq*(1.0/6 + xSq*(-1.0/180 + xSq*(1.0/2835)));
		return log(sinh(x)/x);
	}
}

//! Switch between nonlinear and linear screening
struct ScreeningFunctions
{	bool linearScreening; //! Whether screening is linearized
	ScreeningFunctions(bool linearScreening) : linearScreening(linearScreening) {}
	
	inline double exp(double x) { return linearScreening ? 1.+x : ::exp(x); }
	inline double sinh(double x) { return linearScreening ? x : ::sinh(x); }
	inline double cosh(double x) { return linearScreening ? 1.+0.5*x*x : ::cosh(x); }
	inline double xcosh(double x) { return linearScreening ? x : x*::cosh(x); }
};


//calculates f[mu] (the G=0 component of mu required to allow charge balance between ions and explicit system)
inline double calc_grad_fmu(size_t i, double rhoExplicitGzero, double nNegGzero, double nPosGzero,
	const double* shapeData, const double* muData, const NonlinearPCMparams* params)
{
	const double &shape=shapeData[i];
	const double &mu=muData[i];
	ScreeningFunctions screenFunc(params->linearScreening);
	double alphaPos = shape * params->ionicConcentration*screenFunc.exp(+mu) / nPosGzero;
	double alphaNeg = shape * params->ionicConcentration*screenFunc.exp(-mu) / nNegGzero;
	double rhoFac = rhoExplicitGzero / sqrt(4*nPosGzero*nNegGzero + rhoExplicitGzero*rhoExplicitGzero);
	return 0.5 * (alphaPos*(1+rhoFac) + alphaNeg*(1-rhoFac));
}


//calculates derivative of f[mu] with respect to the shape function
inline double calc_derivfmu_shape(size_t i, double rhoExplicitGzero, double nNegGzero, double nPosGzero,
	const double* muData, const NonlinearPCMparams* params)
{
	const double &mu=muData[i];
	ScreeningFunctions screenFunc(params->linearScreening);
	double betaPos = +params->ionicConcentration*screenFunc.exp(+mu) / nPosGzero;
	double betaNeg = -params->ionicConcentration*screenFunc.exp(-mu) / nNegGzero;
	double rhoFac = rhoExplicitGzero / sqrt(4*nPosGzero*nNegGzero + rhoExplicitGzero*rhoExplicitGzero);
	return 0.5 * (betaPos*(1+rhoFac) + betaNeg*(1-rhoFac));
}

//Compute polarization density, entropy and dipole correlation terms, and their gradient
inline double dipoleCorrEntropy(size_t i, vector3<const double*> epsData, const double* shapeData,
		vector3<double*> pData, vector3<double*> grad_epsData, double* grad_shapeData,
		const NonlinearPCMparams* params)
{	const double &shape=shapeData[i];
	vector3<> epsVec = loadVector(epsData, i);
	register double eps = epsVec.length();
	register double frac = polFrac(eps, params->chiFactor, params->linearDielectric);
	register double fracPrime = polFrac_x(eps, params->linearDielectric);

	vector3<> pVec = epsVec * params->Nbulk * shape * params->pMol * frac;
	register double grad_shape = params->Nbulk*(params->T*(eps*eps*frac - logsinch(eps, params->linearDielectric)) - 0.5*params->Kdip*pow(eps*frac,2));
	register double A = shape*grad_shape;
	register double grad_eps_by_eps = shape*params->Nbulk*(params->T*(frac + eps*fracPrime) - params->Kdip*frac*(frac+eps*fracPrime));
	vector3<> grad_epsVec = grad_eps_by_eps*epsVec;

	storeVector(pVec, pData, i);
	storeVector(grad_epsVec, grad_epsData, i);
	if(grad_shapeData) grad_shapeData[i] += grad_shape;
	return A;
}

//Compute polarization density, entropy and dipole correlation terms, (ionic screening terms if present) and their gradient
inline double IonEnergy(size_t i,const double* muEffData, const double* shapeData, double* grad_muEffData,
						double* grad_shapeData, const NonlinearPCMparams* params)
{
	const double &shape=shapeData[i];
	const double &muEff=muEffData[i];
	ScreeningFunctions screenFunc(params->linearScreening);
	
	register double grad_shape = 2.0*params->ionicConcentration*params->T*(muEff*screenFunc.sinh(muEff)-screenFunc.cosh(muEff));
	register double A = shape*grad_shape;

	register double grad_muEff =  2.0*params->ionicConcentration*params->T*shape*screenFunc.xcosh(muEff);
	grad_muEffData[i] += grad_muEff;
	if(grad_shapeData) grad_shapeData[i] += grad_shape;
	return A;
}

//Accumulate contributions from grad_p into grad_eps:
inline void convert_grad_p_eps(size_t i, vector3<const double*> epsData, const double* shapeData,
		vector3<const double*> grad_pData, vector3<double*> grad_epsData, double* grad_shapeData,
		const NonlinearPCMparams* params)
{	const double &shape=shapeData[i];
	vector3<> epsVec = loadVector(epsData, i);
	vector3<> grad_pVec = loadVector(grad_pData, i);

	register double eps = epsVec.length();
	register double frac=polFrac(eps, params->chiFactor, params->linearDielectric);
	register double fracPrime_by_eps=polFrac_x_over_x(eps, params->linearDielectric);

	vector3<> grad_epsVec = shape*params->Nbulk*params->pMol*(grad_pVec*frac + fracPrime_by_eps*epsVec*dot(epsVec,grad_pVec));
	accumVector(grad_epsVec, grad_epsData, i);

	if(grad_shapeData) grad_shapeData[i] += dot(epsVec,grad_pVec)*params->pMol*frac*params->Nbulk;
}



//Accumulate contributions from grad_muEff into grad_mu:
inline void convert_grad_muEff_mu(size_t i, const double* muData, const double* shapeData,
		const double* grad_muEffData, double* grad_muData, double* grad_shapeData, const double rhoExplicitGzero,
		const double nNegGzero, const double nPosGzero, const double sumGradMuEff, const NonlinearPCMparams* params)
{
	double grad_muEff=grad_muEffData[i];
	double grad_fmu=calc_grad_fmu(i, rhoExplicitGzero, nNegGzero, nPosGzero, shapeData, muData, params);
	grad_muData[i] = grad_muEff-sumGradMuEff*grad_fmu;

	if(grad_shapeData)
	{	double deriv_fmu_shape = calc_derivfmu_shape(i, rhoExplicitGzero, nNegGzero, nPosGzero, muData, params);
		grad_shapeData[i] -= sumGradMuEff*deriv_fmu_shape;
	}
}

//Accumulate contributions from grad_muEff into grad_mu:
inline void calc_rhoIon(size_t i, const double* muEffData, const double* shapeData,
		double* rhoIonData, double* grad_rhoIon_muEffData, double* deriv_rhoIon_shapeData,const NonlinearPCMparams* params)
{
	register double muEff, shape;
	muEff = muEffData[i];
	shape = shapeData[i];
	ScreeningFunctions screenFunc(params->linearScreening);
	rhoIonData[i]=2.0*params->ionicConcentration*shape*screenFunc.sinh(muEff);
	grad_rhoIon_muEffData[i]=2.0*params->ionicConcentration*shape*screenFunc.xcosh(muEff)/muEff;

	if(deriv_rhoIon_shapeData)
		deriv_rhoIon_shapeData[i] = 2.0*params->ionicConcentration*screenFunc.sinh(muEff);
}

NonlinearPCM::NonlinearPCM(const Everything& e, const FluidSolverParams& fsp)
: FluidSolver(e), params(fsp), MuKernel(e.gInfo)
{	
	//Initialize extra parameters:
	params.Kdip = (3.0 * params.T * (params.epsilonBulk - params.epsInf) - 4.*M_PI*params.Nbulk*pow(params.pMol, 2))/(params.epsilonBulk - 1.);
	params.k2factor= (8*M_PI/params.T) * params.ionicConcentration * pow(params.ionicZelectrolyte,2);
	params.chiFactor = (params.T/pow(params.pMol, 2))*(params.epsInf - 1.)/(4*M_PI*params.Nbulk);
	
	// Check whether the linear response makes sense
	if(params.Kdip < 0.)
		die("\nCurrent Nonlinear PCM parameters imply negative correlation for dipole rotations.\n"
			"\tHINT: Reduce solvent molecule dipole or epsInf to fix this");
}

inline void setMuKernel(int i, double G2, double* MuKernel, double meanKsq)
{	if(i==0) MuKernel[i] = meanKsq ? 1.0/meanKsq : 0.0;
	else MuKernel[i] = 1.0/(G2 + meanKsq);
}

void NonlinearPCM::set(const DataGptr& rhoExplicitTilde, const DataGptr& nCavityTilde)
{	this->rhoExplicitTilde = rhoExplicitTilde;
	this->nCavity = I(nCavityTilde);
	//Initialize point-dipole density
	nullToZero(shape,e.gInfo);
	pcmShapeFunc(nCavity, shape, params.nc, params.sigma);

	DataRptr epsilon = 1 + (params.epsilonBulk-1)*shape;
	DataRptr kappaSq = params.ionicConcentration ? params.k2factor*shape : 0;

	//Update preconditioner:
	double meanKsq = (kappaSq ? sum(kappaSq)/sum(epsilon) : 0.0);
	applyFuncGsq(e.gInfo, setMuKernel, MuKernel.data, meanKsq);
	MuKernel.set();

	if(!state)
	{	//Initialize eps using a linear response approx:
		RealKernel gauss(e.gInfo); initGaussianKernel(gauss, 2.0);
		nullToZero(state, e.gInfo); initRandomFlat(state); state=e.gInfo.dV*I(gauss*J(state));
		//setMuEps(state, mu, eps);
	}
}

double NonlinearPCM::operator()(const DataRMuEps& state, DataRMuEps& grad_state, DataGptr* grad_rhoExplicitTilde, DataGptr* grad_nCavityTilde) const
{
	DataRptr grad_shape; if(grad_nCavityTilde) nullToZero(grad_shape, e.gInfo);

	//Calculate the energy contribution/gradient from the dipole entropy
	DataRptrVec p(e.gInfo);
	const DataRptrVec& eps = getEps(state);
	const DataRptr& mu = getMu(state);
	DataRptr grad_mu; DataRptrVec grad_eps; nullToZero(grad_mu, e.gInfo); nullToZero(grad_eps, e.gInfo);
	DataRptr rhoIon; nullToZero(rhoIon, e.gInfo);	//initialize because needed in bound charge

	double AcorrNI = e.gInfo.dV*threadedAccumulate(dipoleCorrEntropy, e.gInfo.nr, eps.data(), shape->data(),
		p.data(), grad_eps.data(), grad_shape ? grad_shape->data() : 0, &params);

	//declare those parameters needed outside the loop
	double nPosGzero = 0.0, nNegGzero = 0.0,rhoExplicitGzero = 0.0;
	DataRptr muEff, grad_muEff, grad_rhoIon_muEff, grad_rhoIon_shape;

	//! The cavitation energy
	double Acavity = cavitationEnergyAndGrad(shape, grad_shape, params.cavityTension, params.cavityPressure);
	
	if(params.ionicConcentration)
	{
		nullToZero(muEff,e.gInfo);
		nullToZero(grad_muEff,e.gInfo);

		complex* rhoPtr = rhoExplicitTilde->data();
		rhoExplicitGzero = rhoPtr[0].real();
		//rhoExplicitGzero=0.0;

		if (params.linearScreening)
		{	nNegGzero=params.ionicConcentration*dot(shape,1.0-mu+0.5*mu*mu)*e.gInfo.dV;
			nPosGzero=params.ionicConcentration*dot(shape,1.0+mu+0.5*mu*mu)*e.gInfo.dV;
		}
		else
		{	nNegGzero=params.ionicConcentration*dot(shape,exp(-mu))*e.gInfo.dV;
			nPosGzero=params.ionicConcentration*dot(shape,exp(mu))*e.gInfo.dV;
		}

		double fmu = log((-rhoExplicitGzero+sqrt(rhoExplicitGzero*rhoExplicitGzero+4.0*nNegGzero*nPosGzero))/(2.0*nNegGzero));
		muEff=mu-fmu;

		//compute rhoIon from muEff=mu+f[mu]
		nullToZero(grad_rhoIon_muEff, e.gInfo);
		if(grad_shape) nullToZero(grad_rhoIon_shape, e.gInfo);

		threadedLoop(calc_rhoIon, e.gInfo.nr, muEff->data(), shape->data(), rhoIon->data(), grad_rhoIon_muEff->data(),
		         (grad_rhoIon_shape ? grad_rhoIon_shape->data() : 0), &params);
	}

	DataGptr nBound = divergence(J(p))+J(rhoIon); //adds rhoIon here if we have ionic screening

	DataGptr phiBound = -4*M_PI*Linv(O(nBound)); //d_fluid
	DataGptr phiExplicit = -4*M_PI*Linv(O(rhoExplicitTilde)); //d_vac
	double Acoulomb = 0.5*dot(nBound,O(phiBound)) + dot(rhoExplicitTilde,O(phiBound));

	DataRptrVec grad_p = -I(gradient(phiBound+phiExplicit));

	threadedLoop(convert_grad_p_eps, e.gInfo.nr, eps.const_data(), shape->data(), grad_p.const_data(),
		         grad_eps.data(),(grad_shape ? grad_shape->data() : 0), &params);

    double Aion = 0.0;
	if(params.ionicConcentration)
	{
		Aion = e.gInfo.dV*threadedAccumulate(IonEnergy, e.gInfo.nr, muEff->data(), shape->data(), grad_muEff->data(),
				(grad_shape ? grad_shape->data() : 0), &params);
		//calculate the gradient with respect to muEff & shape if necessary of the coulomb piece here
		grad_muEff += grad_rhoIon_muEff*I(phiBound+phiExplicit);
		if (grad_shape) grad_shape += grad_rhoIon_shape*I(phiBound+phiExplicit);
		//convert mu_eff to mu here.
		double sumGradMuEff = sum(grad_muEff)*e.gInfo.dV;
		threadedLoop(convert_grad_muEff_mu, e.gInfo.nr, mu->data(), shape->data(),
			grad_muEff->data(), grad_mu->data(),(grad_shape ? grad_shape->data() : 0),
			rhoExplicitGzero, nNegGzero, nPosGzero, sumGradMuEff, &params);
	}

	if(grad_rhoExplicitTilde)
	{	*grad_rhoExplicitTilde = phiBound;
	}
	if(grad_nCavityTilde)
	{	DataRptr grad_nCavity(DataR::alloc(e.gInfo));
		pcmShapeFunc_grad(nCavity, grad_shape, grad_nCavity, params.nc, params.sigma);
		*grad_nCavityTilde = J(grad_nCavity);
	}
	setMuEps(grad_state, grad_mu, grad_eps);
	grad_state *= e.gInfo.dV; //converts variational derivative to total derivative
	return AcorrNI + Acoulomb + Aion + Acavity;
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

double NonlinearPCM::get_Adiel_and_grad(DataGptr& grad_rhoExplicitTilde, DataGptr& grad_nCavityTilde)
{	DataRMuEps grad_state;//change to grad_state
	return (*this)(state, grad_state, &grad_rhoExplicitTilde, &grad_nCavityTilde);
}

void NonlinearPCM::step(const DataRMuEps& dir, double alpha)
{	axpy(alpha, dir, state);
}

double NonlinearPCM::compute(DataRMuEps* grad)
{	DataRMuEps gradUnused;
	return (*this)(state, grad ? *grad : gradUnused);
}

DataRMuEps NonlinearPCM::precondition(const DataRMuEps& in) const
{	DataRMuEps PreconIn;
	setMuEps(PreconIn, I(MuKernel*J(getMu(in))), getEps(in));
	return PreconIn;
}
