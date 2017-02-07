/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

This file is part of Fluid1D.

Fluid1D is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Fluid1D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Fluid1D.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#include <fluid/Fex_H2O_BondedVoids.h>
#include <fluid/Fex_LJ.h>
#include <fluid/FluidMixture.h>
#include <core/Units.h>
#include <core/Operators.h>

//EOS functional fit parameters:
const double RV0 = 1.290*Angstrom;
const double TV = 258.7*Kelvin;
const double kappa = 1.805e5*Kelvin*pow(Angstrom,3);
const double RO = 1.419*Angstrom;
const double sigmaU = 2.62*Angstrom;

//SPC/E O-H length
const double rOH = 1.0*Angstrom;

Fex_H2O_BondedVoids::Fex_H2O_BondedVoids(FluidMixture& fluidMixture)
: Fex(fluidMixture),
RV(RV0*exp(-T/TV)),
propO(gInfo,  RO,0.0, +0.8476,&siteChargeKernel),
propH(gInfo, 0.0,0.0, -0.4238,&siteChargeKernel),
propV(gInfo,  RV,0.0, 0,0),
molecule("H2O",
		&propO, vector3<>(0,0,0),
		&propH, vector3<>(-1,-1,1)*(rOH/sqrt(3)), vector3<>(+1,+1,+1)*(rOH/sqrt(3)),
		&propV, vector3<>(+1,+1,-1)*((RO+RV)/sqrt(3)), vector3<>(-1,-1,-1)*((RO+RV)/sqrt(3)) )
{
	//Initialize the kernels:
	setQkernel(siteChargeKernel, gInfo, 1.385*Angstrom);
	setLJatt(Ua, gInfo, -9.0/(32*sqrt(2)*M_PI*pow(sigmaU,3)), sigmaU);
}

double Fex_H2O_BondedVoids::get_aDiel() const
{	return 1 - T/(7.35e3*Kelvin);
}

double Fex_H2O_BondedVoids::compute(const ScalarFieldTilde* Ntilde, ScalarFieldTilde* grad_Ntilde) const
{	ScalarFieldTilde V = O((-kappa) * (Ua * Ntilde[0]));
	grad_Ntilde[0] += V;
	return 0.5*dot(V,Ntilde[0]);
}
double Fex_H2O_BondedVoids::computeUniform(const double* N, double* grad_N) const
{	grad_N[0] += (-kappa)*Ua.data()[0]*N[0];
	return 0.5*(-kappa)*N[0]*Ua.data()[0]*N[0];
}
void Fex_H2O_BondedVoids::directCorrelations(const double* N, ScalarFieldTildeCollection& C) const
{	ScalarFieldTilde UaTilde(Ua, gInfo);
	C[fluidMixture.corrFuncIndex(0,0,this)] += (-kappa) * UaTilde;
}

