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

#ifndef JDFTX_FLUID_FEX_H2O_SCALAREOS_INTERNAL_H
#define JDFTX_FLUID_FEX_H2O_SCALAREOS_INTERNAL_H

#include <core/Units.h>
#include <core/scalar.h>
#include <fluid/ErfFMTweight.h>

//High freq cutoff on the coulomb kernel expressed as a site charge kernel
inline void setCoulombCutoffKernel(int i, double G2, double* siteChargeKernel)
{	const double Gc = 0.33;
	double G = sqrt(G2);
	siteChargeKernel[i] = 1./sqrt(1 + pow(G/Gc,4));
}

struct ScalarEOS_eval
{
	double T, alpha;
	double b; //!< exclusion volume
	double prefacHB, prefacVW1, prefacVW2; //!< prefactors to the HB and VW terms
	
	double sphereRadius, zi3; //!< Hard sphere radius and volume (= scalar weight function w3 at G=0)
	
	ScalarEOS_eval(double T) : T(T)
	{
		const double TB = 1408.4 * Kelvin; //Boyle temperature
		const double vB = 4.1782e-5 * pow(meter,3)/mol; //Boyle volume
		const double Tf = 273.16 * Kelvin; //Triple point temperatue

		//Some of the random Jeff-Austin EOS parameters (see their paper for details)
		alpha = 2.145*vB;
		const double bStar = 1.0823*vB;
		const double aVW = 0.5542 * Joule*pow(meter,3)/pow(mol,2);
		const double C1 = 0.7140;
		const double lambda = 0.3241;

		const double epsHB = -11.49 * KJoule/mol;
		const double S0 = -61.468 * Joule/(mol*Kelvin), omega0 = exp(-S0);
		const double SHB = -5.128 * Joule/(mol*Kelvin), omegaHB = exp(-SHB);
		const double b1 = 0.25081;
		const double b2 = 0.99859;
		
		b = vB * (0.2*exp(-21.4*pow(T/TB+0.0445,3)) - b1*exp(1.016*T/TB) + b2);
		prefacHB = -2*T * log((omega0+omegaHB*exp(-epsHB/T))/(omega0+omegaHB)) * exp(-0.18*pow(T/Tf,8)) * (1+C1);
		prefacVW1 = -T*alpha/(lambda*b);
		prefacVW2 = bStar*T + aVW;
		
		//Hard sphere properties:
		sphereRadius = 1.36 * Angstrom;
		zi3 = (4*M_PI/3) * pow(sphereRadius,3);
	}

	//Compute the per-particle free energies at each grid point, and the gradient w.r.t the weighted density
	__hostanddev__ void operator()(int i, const double* Nbar, double* Aex, double* grad_Nbar) const
	{
		if(Nbar[i]<0.)
		{	Aex[i] = 0.;
			grad_Nbar[i] = 0.;
			return;
		}
		//More Jeff-Austin parameters:
		const double nHB = (0.8447e6/18.01528) * mol/pow(meter,3);
		const double dnHB = 0.1687*nHB;
		const double C1 = 0.7140;
		const double lambda = 0.3241;
		const double NaN = 0.0/0.0; //returned by FMT/VW when beyond pole
		//HB part:
		double gaussHB = exp(pow((Nbar[i]-nHB)/dnHB,2));
		double fHBden = C1 + gaussHB;
		double fHBdenPrime = gaussHB * 2*(Nbar[i]-nHB)/pow(dnHB,2);
		double AHB = prefacHB / fHBden;
		double AHB_Nbar = -AHB * fHBdenPrime/fHBden;
		//VW part:
		double Ginv = 1 - lambda*b*Nbar[i]; if(Ginv<=0.0) { grad_Nbar[i] = NaN; Aex[i]=NaN; return; }
		double AVW = prefacVW1*log(Ginv) -  Nbar[i]*prefacVW2;
		double AVW_Nbar = T*alpha/Ginv - prefacVW2;
		//FMT part:
		double n3 = zi3*Nbar[i], den = 1./(1-n3);
		double AFMT_Nbar = T*zi3 * (den*den*den)*2*(2-n3);
		double AFMT = T * (den*den)*n3*(4-3*n3); //corresponds to Carnahan-Starling EOS
		//Total
		grad_Nbar[i] = AHB_Nbar + AVW_Nbar - AFMT_Nbar;
		Aex[i] = AHB + AVW - AFMT;
	}
};

#endif // JDFTX_FLUID_FEX_H2O_SCALAREOS_INTERNAL_H
