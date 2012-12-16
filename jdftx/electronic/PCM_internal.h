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

#ifndef JDFTX_ELECTRONIC_PCM_INTERNAL_H
#define JDFTX_ELECTRONIC_PCM_INTERNAL_H

#include <core/vector3.h>

//----------- Common PCM functions (top level interface not seen by .cu files) ------------
#ifndef __in_a_cu_file__

#include <core/Operators.h>

//! Compute the shape function (0 to 1) given the cavity-determining electron density
void pcmShapeFunc(const DataRptr& nCavity, DataRptr& shape, const double nc, const double sigma);

//!Compute derivative with respect to cavity-determining electron density, given derivative with respect to shape function
void pcmShapeFunc_grad(const DataRptr& nCavity, const DataRptr& grad_shape, DataRptr& grad_nCavity, const double nc, const double sigma);

//! Returns the cavitation energy contribution (volume and surface terms) and sets the gradient to grad
double cavitationEnergyAndGrad(const DataRptr& shape, DataRptr& grad, double cavityTension, double cavityPressure);

#endif


//--------- Compute kernels (shared by CPU and GPU implementations) --------

//Cavity shape function and gradient
__hostanddev__ void pcmShapeFunc_calc(int i, const double* nCavity, double* shape, const double nc, const double sigma)
{	shape[i] = erfc(sqrt(0.5)*log(fabs(nCavity[i])/nc)/sigma)*0.5;
}
__hostanddev__ void pcmShapeFunc_grad_calc(int i, const double* nCavity, const double* grad_shape, double* grad_nCavity, const double nc, const double sigma)
{	grad_nCavity[i] = (-1.0/(nc*sigma*sqrt(2*M_PI))) * grad_shape[i]
		* exp(0.5*(pow(sigma,2) - pow(log(fabs(nCavity[i])/nc)/sigma + sigma, 2)));
}


//------------- Helper classes for NonlinearPCM  -------------
namespace NonlinearPCMeval
{
	//!Helper class for ionic screening portion of NonlinearPCM
	struct Screening
	{
		bool linear; //!< whether ionic screening is linearized
		double NT, NZ; //!< NT and NZ, where T=temperature, N=bulk ionic concentration and Z=charge (assumed balanced)
		
		Screening(bool linear, double T, double Nion, double Zion, double epsBulk); //epsBulk is used only for printing screening length
		
		#ifndef __in_a_cu_file__
		//! Compute the neutrality Lagrange multiplier mu0 and optionally its derivatives
		inline double neutralityConstraint(const DataRptr& mu, const DataRptr& shape, double Qexp,
			DataRptr* mu0_mu=0, DataRptr* mu0_shape=0, double* mu0_Qexp=0)
		{
			if(linear)
			{	double Qsum = 2.*NZ * integral(shape);
				double Qdiff = 2.*NZ * integral(shape * mu);
				//Compute the constraint function and its derivatives w.r.t above moments:
				double mu0 = -(Qexp + Qdiff)/Qsum;
				double mu0_Qdiff = -1./Qsum;
				double mu0_Qsum = (Qexp + Qdiff)/(Qsum*Qsum);
				//Collect reuslt and optional gradients:
				if(mu0_mu) *mu0_mu = (mu0_Qdiff * 2.*NZ) * shape;
				if(mu0_shape) *mu0_shape = 2.*NZ * (mu0_Qdiff * mu + mu0_Qsum);
				if(mu0_Qexp) *mu0_Qexp = -1./Qsum;
				return mu0;
			}
			else
			{	DataRptr etaPlus  = exp(mu);
				DataRptr etaMinus = exp(-mu);
				double Qplus  = +NZ * integral(shape * etaPlus);
				double Qminus = -NZ * integral(shape * etaMinus);
				//Compute the constraint function and its derivatives w.r.t above moments:
				double mu0, mu0_Qplus, mu0_Qminus;
				double disc = sqrt(Qexp*Qexp - 4.*Qplus*Qminus); //discriminant for quadratic
				//Pick the numerically stable path (to avoid roundoff problems when |Qplus*Qminus| << Qexp^2):
				if(Qexp<0) 
				{	mu0 = log((disc-Qexp)/(2.*Qplus));
					mu0_Qplus  = -2.*Qminus/(disc*(disc-Qexp)) - 1./Qplus;
					mu0_Qminus = -2.*Qplus/(disc*(disc-Qexp));
				}
				else
				{	mu0 = log(-2.*Qminus/(disc+Qexp));
					mu0_Qplus  = 2.*Qminus/(disc*(disc+Qexp));
					mu0_Qminus = 2.*Qplus/(disc*(disc+Qexp)) + 1./Qminus;
				}
				//Collect result and optional gradients:
				if(mu0_mu) *mu0_mu = NZ * shape * (mu0_Qplus * etaPlus  + mu0_Qminus * etaMinus);
				if(mu0_shape) *mu0_shape = NZ * (mu0_Qplus * etaPlus - mu0_Qminus * etaMinus);
				if(mu0_Qexp) *mu0_Qexp = -1./disc;
				return mu0;
			}
		}
		#endif
		
		//! Compute the nonlinear functions in the free energy and charge density prior to scaling by shape function
		//! Note that x here is muEff = mu(r) + mu0, i.e. after imposing charge neutrality constraint
		__hostanddev__ void compute(double x, double& F, double& F_x, double& Rho, double& Rho_x) const
		{	if(linear)
			{	F = 2.*NT * (0.5*x*x);
				F_x = 2.*NT * x;
				Rho = 2.*NZ * x;
				Rho_x = 2.*NZ;
			}
			else
			{	double s = sinh(x), c = cosh(x);
				F = 2.*NT * (x*s + 1. - c);
				F_x = 2.*NT * (x*c);
				Rho = 2.*NZ * s;
				Rho_x = 2.*NZ * c;
			}
		}
		
		//! Given shape function s and potential mu, compute induced charge rho, free energy density A and accumulate its derivatives
		__hostanddev__ void freeEnergy_calc(size_t i, double mu0, const double* mu, const double* s, double* rho, double* A, double* A_mu, double* A_s) const
		{	double F, F_mu, Rho, Rho_mu;
			compute(mu[i]+mu0, F, F_mu, Rho, Rho_mu);
			A[i] = s[i] * F;
			A_mu[i] += s[i] * F_mu;
			if(A_s) A_s[i] += F;
			rho[i] = s[i] * Rho;
		}
		void freeEnergy(size_t N, double mu0, const double* mu, const double* s, double* rho, double* A, double* A_mu, double* A_s) const;
		#ifdef GPU_ENABLED
		void freeEnergy_gpu(size_t N, double mu0, const double* mu, const double* s, double* rho, double* A, double* A_mu, double* A_s) const;
		#endif
		
		//! Propagate derivative A_rho and accumulate to A_mu and A_s
		__hostanddev__ void convertDerivative_calc(size_t i, double mu0, const double* mu, const double* s, const double* A_rho, double* A_mu, double* A_s) const
		{	double F, F_mu, Rho, Rho_mu;
			compute(mu[i]+mu0, F, F_mu, Rho, Rho_mu);
			A_mu[i] += s[i] * Rho_mu * A_rho[i];
			if(A_s) A_s[i] += Rho * A_rho[i];
		}
		void convertDerivative(size_t N, double mu0, const double* mu, const double* s, const double* A_rho, double* A_mu, double* A_s) const;
		#ifdef GPU_ENABLED
		void convertDerivative_gpu(size_t N, double mu0, const double* mu, const double* s, const double* A_rho, double* A_mu, double* A_s) const;
		#endif
	};
	
	//!Helper class for dielectric portion of NonlinearPCM
	struct Dielectric
	{
		bool linear; //!< whether dielectric is linearized
		double Np, NT; //!< N*p, p/T, N*T and N*chi where N is molecular density, p is molecular dipole and chi is molecular polarizability
		double alpha, X; //!< dipole correlation factor and chi*T/p^2 where chi is the molecular susceptibility
		
		Dielectric(bool linear, double T, double Nmol, double pMol, double epsBulk, double epsInf);
		
		//! Compute the nonlinear functions in the free energy and effective susceptibility (p/eps) prior to scaling by shape function
		__hostanddev__ void compute(double epsSqHlf, double& F, double& F_epsSqHlf, double& ChiEff, double& ChiEff_epsSqHlf) const
		{	double epsSq = 2.*epsSqHlf, eps = sqrt(epsSq);
			//----- Nonlinear functions of eps
			double frac, frac_epsSqHlf, logsinch;
			if(linear)
			{	frac = 1.0/3;
				frac_epsSqHlf = 0.;
				logsinch = epsSq*(1.0/6);
			}
			else
			{	if(eps < 1e-1) //Use series expansions
				{	frac = 1.0/3 + epsSq*(-1.0/45 + epsSq*(2.0/945 + epsSq*(-1.0/4725)));
					frac_epsSqHlf = -2.0/45 + epsSq*(8.0/945 + epsSq*(-6.0/4725));
					logsinch = epsSq*(1.0/6 + epsSq*(-1.0/180 + epsSq*(1.0/2835)));
				}
				else
				{	frac = (eps/tanh(eps)-1)/epsSq;
					frac_epsSqHlf = (2 - eps/tanh(eps) - pow(eps/sinh(eps),2))/(epsSq*epsSq);
					logsinch = eps<20. ? log(sinh(eps)/eps) : eps - log(2.*eps);
				}
			}
			//----- Free energy and derivative
			double screen = 1 - alpha*frac; //correlation screening factor = (pE/T) / eps where E is the real electric field
			F = NT * (epsSq*(frac - 0.5*alpha*frac*frac + 0.5*X*screen*screen) - logsinch);
			F_epsSqHlf = NT * (frac + X*screen + epsSq*frac_epsSqHlf*(1.-X*alpha)) * screen; //Simplified using logsinch_epsSqHlf = frac
			//----- Effective susceptibility and derivative
			ChiEff = Np * (frac + X*screen);
			ChiEff_epsSqHlf = Np * frac_epsSqHlf * (1.-X*alpha);
		}
		
		//! Given shape function s and gradient of phi eps, compute polarization p, free energy density A and accumulate its derivatives
		__hostanddev__ void freeEnergy_calc(size_t i, vector3<const double*> eps, const double* s, vector3<double*> p, double* A, vector3<double*> A_eps, double* A_s) const
		{	vector3<> epsVec = loadVector(eps, i);
			double epsSqHlf = 0.5*epsVec.length_squared();
			double F, F_epsSqHlf, ChiEff, ChiEff_epsSqHlf;
			compute(epsSqHlf, F, F_epsSqHlf, ChiEff, ChiEff_epsSqHlf);
			A[i] = F * s[i];
			accumVector((F_epsSqHlf * s[i]) * epsVec, A_eps, i);
			if(A_s) A_s[i] += F;
			storeVector((ChiEff * s[i]) * epsVec, p, i);
		}
		void freeEnergy(size_t N, vector3<const double*> eps, const double* s, vector3<double*> p, double* A, vector3<double*> A_eps, double* A_s) const;
		#ifdef GPU_ENABLED
		void freeEnergy_gpu(size_t N, vector3<const double*> eps, const double* s, vector3<double*> p, double* A, vector3<double*> A_eps, double* A_s) const;
		#endif
		
		//! Propagate derivative A_p and accumulate to A_eps and A_s
		__hostanddev__ void convertDerivative_calc(size_t i, vector3<const double*> eps, const double* s, vector3<const double*> A_p, vector3<double*> A_eps, double* A_s) const
		{	vector3<> epsVec = loadVector(eps, i);
			double epsSqHlf = 0.5*epsVec.length_squared();
			double F, F_epsSqHlf, ChiEff, ChiEff_epsSqHlf;
			compute(epsSqHlf, F, F_epsSqHlf, ChiEff, ChiEff_epsSqHlf);
			//Propagate derivatives:
			vector3<> A_pVec = loadVector(A_p, i);
			accumVector(s[i]*( ChiEff*A_pVec + ChiEff_epsSqHlf*epsVec*dot(A_pVec, epsVec)), A_eps, i);
			if(A_s) A_s[i] += ChiEff * dot(A_pVec, epsVec);
		}
		void convertDerivative(size_t N, vector3<const double*> eps, const double* s, vector3<const double*> A_p, vector3<double*> A_eps, double* A_s) const;
		#ifdef GPU_ENABLED
		void convertDerivative_gpu(size_t N, vector3<const double*> eps, const double* s, vector3<const double*> A_p, vector3<double*> A_eps, double* A_s) const;
		#endif
	};
}

#endif // JDFTX_ELECTRONIC_PCM_INTERNAL_H
