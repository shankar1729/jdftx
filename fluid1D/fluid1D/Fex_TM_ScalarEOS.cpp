/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman

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

#include <fluid/Fex_H2O_ScalarEOS_internal.h>
#include <fluid1D/Fex_TM_ScalarEOS.h>
#include <fluid1D/Fex_LJ.h>
#include <fluid1D/FluidMixture.h>

Fex_TM_ScalarEOS::Fex_TM_ScalarEOS(FluidMixture& fluidMixture, double Tc, double Pc, double omega, double sphereRadius)
: Fex(fluidMixture),
eval(new TaoMasonEOS_eval(T, Tc, Pc, omega, sphereRadius))
{
	//Initialize the kernels:
	setQuarticQkernel(siteChargeKernel, gInfo, 0.33);
	setLJatt(fex_LJatt, gInfo, -9.0/(32*sqrt(2)*M_PI*pow(2*sphereRadius,3)), 2*sphereRadius);
}

Fex_TM_ScalarEOS::~Fex_TM_ScalarEOS()
{	delete eval;
}

double Fex_TM_ScalarEOS::compute(const ScalarFieldTilde* Ntilde, ScalarFieldTilde* grad_Ntilde) const
{	//Compute LJatt weighted density:
	ScalarField Nbar = I(fex_LJatt * Ntilde[0]);
	//Evaluated weighted denisty functional:
	ScalarField Aex(&gInfo), AexPrime(&gInfo);
	serialLoop(eval, gInfo.S, Nbar.data(), Aex.data(), AexPrime.data());
	//Convert gradients:
	ScalarFieldTilde OJAex = O(J(Aex));
	grad_Ntilde[0] += fex_LJatt * Idag(Diag(AexPrime) * Jdag(O(Ntilde[0]))) + OJAex;
	return dot(Ntilde[0], OJAex);
}

double Fex_TM_ScalarEOS::computeUniform(const double* N, double* grad_N) const
{	double AexPrime, Aex;
	(*eval)(0, &N[0], &Aex, &AexPrime);
	grad_N[0] += Aex + N[0]*AexPrime;
	return N[0]*Aex;
}

void Fex_TM_ScalarEOS::directCorrelations(const double* N, ScalarFieldTildeCollection& C) const
{	//Compute upto second derivative of per-particle free energy:
	double AexPrime, Aex; (*eval)(0, &N[0], &Aex, &AexPrime);
	const double dN = N[0]*1e-7;
	double AexPrimePlus,  AexPlus,  Nplus  = N[0]+dN; (*eval)(0, &Nplus,  &AexPlus,  &AexPrimePlus);
	double AexPrimeMinus, AexMinus, Nminus = N[0]-dN; (*eval)(0, &Nminus, &AexMinus, &AexPrimeMinus);
	double AexDblPrime = (AexPrimePlus - AexPrimeMinus) / (2*dN);
	//Accumulate correlations:
	ScalarFieldTilde fex_LJattTilde(fex_LJatt, gInfo); //A scalar field version of kernel to ease arithmetic
	C[fluidMixture.corrFuncIndex(0,0,this)] += (2*AexPrime*fex_LJattTilde + N[0]*AexDblPrime*(fex_LJatt*fex_LJattTilde));
}

double Fex_TM_ScalarEOS::vdwRadius() const
{	return eval->vdwRadius();
}

//--------- Chloroform -----------

const double rCCl_CHCl3 = 1.762*Angstrom;
const double rCH_CHCl3 = 1.073*Angstrom;
const double thetaHCCl_CHCl3 = 107.98 * M_PI/180;

Fex_CHCl3_ScalarEOS::Fex_CHCl3_ScalarEOS(FluidMixture& fluidMixture)
: Fex_TM_ScalarEOS(fluidMixture, 536.6*Kelvin, 5328.68*KPascal, 0.216, 2.06*Angstrom),
propC(gInfo, eval->sphereRadius,0., -0.5609,&siteChargeKernel, true, 0.878*pow(Angstrom,3),&siteChargeKernel),
propH(gInfo, 0.,0.,0.0551,&siteChargeKernel, true, 0.135*pow(Angstrom,3),&siteChargeKernel),
propCl(gInfo, 0.,0.,0.1686,&siteChargeKernel, true, 1.910*pow(Angstrom,3),&siteChargeKernel),
molecule("CHCl3",
	&propC,
		 vector3<>(0,0,0),
	&propH,
		vector3<>(0,0,rCH_CHCl3),
	&propCl,
		vector3<>(0, rCCl_CHCl3*sin(thetaHCCl_CHCl3), rCCl_CHCl3*cos(thetaHCCl_CHCl3)),
		vector3<>(+sqrt(0.75)*rCCl_CHCl3*sin(thetaHCCl_CHCl3), -0.5*rCCl_CHCl3*sin(thetaHCCl_CHCl3), rCCl_CHCl3*cos(thetaHCCl_CHCl3)),
		vector3<>(-sqrt(0.75)*rCCl_CHCl3*sin(thetaHCCl_CHCl3), -0.5*rCCl_CHCl3*sin(thetaHCCl_CHCl3), rCCl_CHCl3*cos(thetaHCCl_CHCl3)) )
{
}

//--------- Carbon tetrachloride -----------

const double rCCl_CCl4 = 1.7829*Angstrom;

Fex_CCl4_ScalarEOS::Fex_CCl4_ScalarEOS(FluidMixture& fluidMixture)
: Fex_TM_ScalarEOS(fluidMixture, 556.4*Kelvin, 4493*KPascal, 0.194, 2.17*Angstrom),
propC(gInfo, eval->sphereRadius,0., -0.6052,&siteChargeKernel, true, 0.878*pow(Angstrom,3),&siteChargeKernel),
propCl(gInfo, 0.,0.,+0.1513,&siteChargeKernel, true, 1.910*pow(Angstrom,3),&siteChargeKernel),
molecule("CCl4",
	&propC,
		 vector3<>(0,0,0),
	&propCl,
		vector3<>(0,0,rCCl_CCl4),
		vector3<>(0, rCCl_CCl4*(sqrt(8.)/3), rCCl_CCl4*(-1./3)),
		vector3<>(+sqrt(0.75)*rCCl_CCl4*(sqrt(8.)/3), -0.5*rCCl_CCl4*(sqrt(8.)/3), rCCl_CCl4*(-1./3)),
		vector3<>(-sqrt(0.75)*rCCl_CCl4*(sqrt(8.)/3), -0.5*rCCl_CCl4*(sqrt(8.)/3), rCCl_CCl4*(-1./3)) )
{
}
