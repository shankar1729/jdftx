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

#ifndef JDFTX_FLUID_FEX_SCALAREOS_INTERNAL_H
#define JDFTX_FLUID_FEX_SCALAREOS_INTERNAL_H

#include <core/Units.h>
#include <core/scalar.h>

struct ScalarEOS_eval
{	double T, b; //!< temperature and exclusion volume
	
	double vdwRadius() const
	{	return pow(3.*b/(16.*M_PI), 1./3);
	}
	
	__hostanddev__ double getAhs(double N, double& Ahs_N, double Vhs) const
	{	double n3 = Vhs*N; if(n3 >= 1.) { Ahs_N = NAN; return NAN; }
		double den = 1./(1-n3);
		Ahs_N = T*Vhs * (den*den*den)*2*(2-n3);
		return T * (den*den)*n3*(4-3*n3); //corresponds to Carnahan-Starling EOS
	}
};

struct JeffereyAustinEOS_eval : public ScalarEOS_eval
{
	double alpha, prefacHB, prefacVW1, prefacVW2; //!< prefactors to the HB and VW terms
	double lambda, C1, nHB, dnHB;
	double nc, VPzi; //!< vapor-pressure correction temperature dependent prefactor

	JeffereyAustinEOS_eval(double T)
	{	this->T = T;
		
		const double TB = 1408.4 * Kelvin; //Boyle temperature
		const double vB = 4.1782e-5 * pow(meter,3)/mol; //Boyle volume
		const double Tf = 273.16 * Kelvin; //Triple point temperature

		//Some of the random Jeff-Austin EOS parameters (see their paper for details)
		alpha = 2.145*vB;
		const double bStar = 1.0823*vB;
		const double aVW = 0.5542 * Joule*pow(meter,3)/pow(mol,2);
		lambda = 0.3241;
		C1 = 0.7140;
		nHB = (0.8447e6/18.01528) * mol/pow(meter,3);
		dnHB = 0.1687*nHB;

		const double epsHB = -11.49 * KJoule/mol;
		const double S0 = -61.468 * Joule/(mol*Kelvin), omega0 = exp(-S0);
		const double SHB = -5.128 * Joule/(mol*Kelvin), omegaHB = exp(-SHB);
		const double b1 = 0.25081;
		const double b2 = 0.99859;
		
		b = vB * (0.2*exp(-21.4*pow(T/TB+0.0445,3)) - b1*exp(1.016*T/TB) + b2);
		prefacHB = -2*T * log((omega0+omegaHB*exp(-epsHB/T))/(omega0+omegaHB)) * exp(-0.18*pow(T/Tf,8)) * (1+C1);
		prefacVW1 = -T*alpha/(lambda*b);
		prefacVW2 = bStar*T + aVW;
		
		//vapor pressure correction prefactor:
		const double A1 = -2.9293;
		const double A2 = 0.1572;
		const double A5 = 1.873;
		const double kappa = 0.8921;
		const double Tc = 647.096*Kelvin;
		nc = (322/18.01528) * mol/liter;
		VPzi = A1*exp(-A5*pow(T/Tc,6)) * (pow(T/Tc - kappa, 2) + A2);
	}

	__hostanddev__ double getVPphiInt(double n, double& VPphiInt_n) const
	{	const double A4 = 0.0610;
		const double d0 = 1.917;
		const double d1 = 26.01;
		const double beta = 3.24;
		double x = n/nc, xPowBetaM1 = pow(x, beta-1.), xPow57 = pow(x,5.7);
		double den = 1./(d0*(1 + d0*x*(1 + d0*x)) + x*d1*xPowBetaM1);
		double den_x = -den*den*(d0*d0*(1 + d0*x*2) + beta*d1*xPowBetaM1);
		double expTerm = exp(-A4*xPow57*x);
		VPphiInt_n = expTerm*(A4*xPow57*6.7*den-den_x);
		return -nc*expTerm*den;
	}
	
	//Compute the per-particle free energies at each grid point, and the gradient w.r.t the weighted density
	__hostanddev__ void operator()(int i, const double* Nbar, double* Aex, double* Aex_Nbar, double Vhs) const
	{
		if(Nbar[i]<0.)
		{	Aex[i] = 0.;
			Aex_Nbar[i] = 0.;
			return;
		}
		//HB part:
		double gaussHB = exp(pow((Nbar[i]-nHB)/dnHB,2));
		double fHBden = C1 + gaussHB;
		double fHBdenPrime = gaussHB * 2*(Nbar[i]-nHB)/pow(dnHB,2);
		double AHB = prefacHB / fHBden;
		double AHB_Nbar = -AHB * fHBdenPrime/fHBden;
		//VW part:
		double Ginv = 1 - lambda*b*Nbar[i]; if(Ginv<=0.0) { Aex_Nbar[i] = NAN; Aex[i]=NAN; return; }
		double VPphiInt_Nbar, VPphiInt = getVPphiInt(Nbar[i], VPphiInt_Nbar);
		double AVW = prefacVW1*log(Ginv) + (VPzi*VPphiInt -  Nbar[i])*prefacVW2;
		double AVW_Nbar = T*alpha/Ginv + (VPzi*VPphiInt_Nbar - 1.) * prefacVW2;
		//FMT part:
		double AFMT_Nbar, AFMT = getAhs(Nbar[i], AFMT_Nbar, Vhs);
		//Total
		Aex_Nbar[i] = AHB_Nbar + AVW_Nbar - AFMT_Nbar;
		Aex[i] = AHB + AVW - AFMT;
	}
};

//!Tao-Mason equation of state [F. Tao and E. A. Mason, J. Chem. Phys. 100, 9075 (1994)]
struct TaoMasonEOS_eval : public ScalarEOS_eval
{
	double lambda, prefacQuad, prefacVap, prefacPole;
	
	//! Construct the equation of state for temperature T, given critical point and acentricity (all in atomic units)
	TaoMasonEOS_eval(double T, double Tc, double Pc, double omega)
	{	this->T = T;
		
		//Constants in vapor pressure correction:
		double kappa = 1.093 + 0.26*(sqrt(0.002+omega) + 4.5*(0.002+omega));
		double A1 = 0.143;
		double A2 = 1.64 + 2.65*(exp(kappa-1.093)-1);
		//Boyle temperature and volume:
		double TB = Tc * (2.6455 - 1.1941*omega);
		double vB = (Tc/Pc) * (0.1646 + 0.1014*omega);
		//Temperature dependent functions:
		const double a1 = -0.0648, a2 = 1.8067, c1 = 2.6038, c2 = 0.9726;
		double alpha = vB*( a1*exp(-c1*T/TB) + a2*(1 - exp(-c2*pow(TB/T,0.25))) );
		double B = (Tc/Pc) * (0.1445+omega*0.0637 + (Tc/T)*(-0.330 + (Tc/T)*(-0.1385+0.331*omega + (Tc/T)*(-0.0121-0.423*omega + pow(Tc/T,5)*(-0.000607-0.008*omega)))));
		b = vB*( a1*(1-c1*T/TB)*exp(-c1*T/TB) + a2*(1 - (1+0.25*c2*pow(TB/T,0.25))*exp(-c2*pow(TB/T,0.25))) );
		lambda = 0.4324 - 0.3331*omega;
		//Coalesce prefactors for each free energy term:
		prefacQuad = -T*(alpha - B);
		prefacVap = prefacQuad * (-A1 * (exp(kappa*Tc/T)-A2)) / (2*sqrt(1.8)*b);
		prefacPole = alpha*T/(lambda*b);
	}
	
	//Compute the per-particle free energies at each grid point, and the gradient w.r.t the weighted density
	__hostanddev__ void operator()(int i, const double* Nbar, double* Aex, double* Aex_Nbar, double Vhs) const
	{
		if(Nbar[i]<0.)
		{	Aex[i] = 0.;
			Aex_Nbar[i] = 0.;
			return;
		}
		//VW part:
		double Ginv = 1 - lambda*b*Nbar[i]; if(Ginv<=0.0) { Aex_Nbar[i] = NAN; Aex[i]=NAN; return; }
		double AVW = Nbar[i]*prefacQuad + prefacPole*(-log(Ginv));
		double AVW_Nbar = prefacQuad + prefacPole*(lambda*b/Ginv);
		//Vapor pressure correction:
		double b2term = sqrt(1.8)*b*b;
		double bn2term = b2term*Nbar[i]*Nbar[i];
		double Avap = prefacVap * atan(bn2term);
		double Avap_Nbar = prefacVap * b2term*Nbar[i]*2. / (1 + bn2term*bn2term);
		//FMT part:
		double AFMT_Nbar, AFMT = getAhs(Nbar[i], AFMT_Nbar, Vhs);
		//Total
		Aex_Nbar[i] = AVW_Nbar + Avap_Nbar - AFMT_Nbar;
		Aex[i] = AVW + Avap - AFMT;
	}
};

#endif // JDFTX_FLUID_FEX_SCALAREOS_INTERNAL_H
