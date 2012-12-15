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

#include <core/Data.h>

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
		double NT2, ZbyT, m2NZ; //!< 2NT, Z/T and -2NZ, where T=temperature, N=bulk ionic concentration and Z=charge (assumed balanced)
		
		Screening(bool linear, double T, double Nion, double Zion);
		
		//! Compute the nonlinear functions in the free energy and charge density prior to scaling by shape function
		__hostanddev__ void compute(double phi, double& F, double& F_phi, double& Rho, double& Rho_phi) const
		{	double Phi = ZbyT * phi; //dimensionless combination
			if(linear)
			{	F = NT2 * (0.5*Phi*Phi);
				F_phi = NT2 * Phi * ZbyT;
				Rho = m2NZ * Phi;
				Rho_phi = m2NZ * ZbyT;
			}
			else
			{	double s = sinh(Phi), c = cosh(Phi);
				F = NT2 * (Phi*s + 1. - c);
				F_phi = NT2 * (Phi*c) * ZbyT;
				Rho = m2NZ * s;
				Rho_phi = m2NZ * c * ZbyT;
			}
		}
		
		//! Given shape function s and potential phi, compute induced charge rho, free energy density A and accumulate its derivatives
		__hostanddev__ void freeEnergy_calc(size_t i, const double* phi, const double* s, double* rho, double* A, double* A_phi, double* A_s) const
		{	double F, F_phi, Rho, Rho_phi;
			compute(phi[i], F, F_phi, Rho, Rho_phi);
			A[i] = s[i] * F;
			A_phi[i] += s[i] * F_phi;
			if(A_s) A_s[i] += F;
			rho[i] = s[i] * Rho;
		}
		void freeEnergy(size_t N, const double* phi, const double* s, double* rho, double* A, double* A_phi, double* A_s) const;
		#ifdef GPU_ENABLED
		void freeEnergy_gpu(size_t N, const double* phi, const double* s, double* rho, double* A, double* A_phi, double* A_s) const;
		#endif
		
		//! Propagate derivative A_rho and accumulate to A_phi and A_s
		__hostanddev__ void convertDerivative_calc(size_t i, const double* phi, const double* s, const double* A_rho, double* A_phi, double* A_s) const
		{	double F, F_phi, Rho, Rho_phi;
			compute(phi[i], F, F_phi, Rho, Rho_phi);
			A_phi[i] += s[i] * Rho_phi * A_rho[i];
			if(A_s) A_s[i] += Rho * A_rho[i];
		}
		void convertDerivative(size_t N, const double* phi, const double* s, const double* A_rho, double* A_phi, double* A_s) const;
		#ifdef GPU_ENABLED
		void convertDerivative_gpu(size_t N, const double* phi, const double* s, const double* A_rho, double* A_phi, double* A_s) const;
		#endif
	};
	
	//!Helper class for dielectric portion of NonlinearPCM
	struct Dielectric
	{
		bool linear; //!< whether dielectric is linearized
		double Np, pByT, NT, Nchi; //!< N*p, p/T, N*T and N*chi where N is molecular density, p is molceular dipole and chi is molecular polarizability
		double alpha; //!< dipole correlation factor
		
		Dielectric(bool linear, double T, double Nmol, double pMol, double epsBulk, double epsInf);
		void initLookupTables(); //!< Initialize fBar lookup tables
		void freeLookupTables(); //!< Free fBar lookup tables
		
		//! Compute the screening factor due to correlations for the electric field squared
		__hostanddev__ double get_fBarSq(double Esq, double& fBarSq_Esq) const
		{	//if(linear)
			{	double fBar = 3./(3-alpha);
				fBarSq_Esq = 0.;
				return fBar * fBar;
			}
		}
		
		//! Compute the nonlinear functions in the free energy and effective susceptibility prior to scaling by shape function
		__hostanddev__ void compute(double gradPhiSqHlf, double& F, double& F_gradPhiSqHlf, double& ChiEff, double& ChiEff_gradPhiSqHlf) const
		{	//Contribution from molecular polarizability:
			F = Nchi * gradPhiSqHlf;
			F_gradPhiSqHlf = Nchi;
			ChiEff = Nchi;
			ChiEff_gradPhiSqHlf = 0.;
			//Contribution from rotations:
			double Esq = pByT*pByT * gradPhiSqHlf*2; //dimensionless field magnitude
			double fBarSq_Esq, fBarSq = get_fBarSq(Esq, fBarSq_Esq); //local field screening factor
			double epsSq = fBarSq * Esq, eps = sqrt(epsSq);
			double frac, frac_epsSqHlf, logsinch;
			//----- Nonlinear functions in th escreened field eps
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
			//----- Rotational free energy and derivative
			F += NT * (epsSq*frac*(1.-0.5*alpha*frac) - logsinch);
			double F_epsSqHlf = NT * (frac + epsSq*frac_epsSqHlf) * (1.-alpha*frac); //Simplified using logsinch_epsSqHlf = frac
			double epsSq_Esq = fBarSq + Esq*fBarSq_Esq;
			F_gradPhiSqHlf += F_epsSqHlf * epsSq_Esq * pByT*pByT;
			//----- Effective susceptibility and derivative
			double fBar = sqrt(fBarSq);
			ChiEff += (Np*pByT) * frac * fBar;
			ChiEff_gradPhiSqHlf += (Np*pByT) * (frac_epsSqHlf*epsSq_Esq*fBar + frac*fBarSq_Esq/fBar) * pByT*pByT;
		}
		
		//! Given shape function s and gradient of phi gradPhi, compute polarization p, free energy density A and accumulate its derivatives
		__hostanddev__ void freeEnergy_calc(size_t i, vector3<const double*> gradPhi, const double* s, vector3<double*> p, double* A, vector3<double*> A_gradPhi, double* A_s) const
		{	vector3<> gradPhiVec = loadVector(gradPhi, i);
			double gradPhiSqHlf = 0.5*gradPhiVec.length_squared();
			double F, F_gradPhiSqHlf, ChiEff, ChiEff_gradPhiSqHlf;
			compute(gradPhiSqHlf, F, F_gradPhiSqHlf, ChiEff, ChiEff_gradPhiSqHlf);
			A[i] = F * s[i];
			accumVector((F_gradPhiSqHlf * s[i]) * gradPhiVec, A_gradPhi, i);
			if(A_s) A_s[i] += F;
			storeVector((-ChiEff * s[i]) * gradPhiVec, p, i);
		}
		void freeEnergy(size_t N, vector3<const double*> gradPhi, const double* s, vector3<double*> p, double* A, vector3<double*> A_gradPhi, double* A_s) const;
		#ifdef GPU_ENABLED
		void freeEnergy_gpu(size_t N, vector3<const double*> gradPhi, const double* s, vector3<double*> p, double* A, vector3<double*> A_gradPhi, double* A_s) const;
		#endif
		
		//! Propagate derivative A_p and accumulate to A_gradPhi and A_s
		__hostanddev__ void convertDerivative_calc(size_t i, vector3<const double*> gradPhi, const double* s, vector3<const double*> A_p, vector3<double*> A_gradPhi, double* A_s) const
		{	vector3<> gradPhiVec = loadVector(gradPhi, i);
			double gradPhiSqHlf = 0.5*gradPhiVec.length_squared();
			double F, F_gradPhiSqHlf, ChiEff, ChiEff_gradPhiSqHlf;
			compute(gradPhiSqHlf, F, F_gradPhiSqHlf, ChiEff, ChiEff_gradPhiSqHlf);
			//Propagate derivatives:
			vector3<> A_pVec = loadVector(A_p, i);
			accumVector(-s[i]*( ChiEff*A_pVec + ChiEff_gradPhiSqHlf*gradPhiVec*dot(A_pVec, gradPhiVec)), A_gradPhi, i);
			if(A_s) A_s[i] += (-ChiEff) * dot(A_pVec, gradPhiVec);
		}
		void convertDerivative(size_t N, vector3<const double*> gradPhi, const double* s, vector3<const double*> A_p, vector3<double*> A_gradPhi, double* A_s) const;
		#ifdef GPU_ENABLED
		void convertDerivative_gpu(size_t N, vector3<const double*> gradPhi, const double* s, vector3<const double*> A_p, vector3<double*> A_gradPhi, double* A_s) const;
		#endif
	};
}

#endif // JDFTX_ELECTRONIC_PCM_INTERNAL_H
