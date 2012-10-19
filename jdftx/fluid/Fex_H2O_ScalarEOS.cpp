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

#include <fluid/Fex_H2O_ScalarEOS_internal.h>
#include <fluid/Fex_H2O_ScalarEOS.h>
#include <fluid/Fex_LJ.h>
#include <core/Units.h>
#include <core/Operators.h>

string rigidMoleculeCDFT_ScalarEOSpaper =  "R. Sundararaman and T.A. Arias, (to be submitted to Comp. Phys. Comm.)";

static const double rOH = 1.0*Angstrom; 
static const double thetaHOH = acos(-1.0/3); 
//static const double rOH = 0.96719*Angstrom; //HACK WARNING
//static const double thetaHOH = 1.8069; //103.53 degrees HACK WARNING

Fex_H2O_ScalarEOS::Fex_H2O_ScalarEOS(FluidMixture& fluidMixture)
: Fex(fluidMixture),
fex_LJatt(gInfo), siteChargeKernel(gInfo),
eval(new ScalarEOS_eval(T)),
propO(gInfo, eval->sphereRadius,0.0, 0.8476,&siteChargeKernel),
propH(gInfo, 0.0,0.0,-0.4238,&siteChargeKernel),
//propO(gInfo, eval->sphereRadius,0.0, 0.81036,&siteChargeKernel),
//propH(gInfo, 0.0,0.0,-0.40518,&siteChargeKernel),
//propO(gInfo, eval->sphereRadius,0.0, 0.66,&siteChargeKernel),  //HACK WARNING
//propH(gInfo, 0.0,0.0,-0.33,&siteChargeKernel), //HACK WARNING
molecule("H2O",
	&propO,
		 vector3<>(0,0,0),
	&propH,
		 vector3<>(0, -rOH*sin(0.5*thetaHOH), rOH*cos(0.5*thetaHOH)),
		 vector3<>(0, +rOH*sin(0.5*thetaHOH), rOH*cos(0.5*thetaHOH)) )
{
	//Initialize the kernels:
	applyFuncGsq(gInfo, setCoulombCutoffKernel, siteChargeKernel.data); siteChargeKernel.set();
	setLJatt(fex_LJatt, -9.0/(32*sqrt(2)*M_PI*pow(2*eval->sphereRadius,3)), 2*eval->sphereRadius);
	Citations::add("Scalar-EOS water functional", rigidMoleculeCDFT_ScalarEOSpaper);
}

double Fex_H2O_ScalarEOS::get_aDiel() const
{
	return 1 - T/(7.35e3*Kelvin); 
	//return 1 - T/(3.9962e3*Kelvin); //HACK WARNING
}

#ifdef GPU_ENABLED
void Fex_H20_ScalarEOS_gpu(int nr, const double* Nbar, double* Fex, double* grad_Nbar, ScalarEOS_eval eval);
#endif
double Fex_H2O_ScalarEOS::compute(const DataGptr* Ntilde, DataGptr* grad_Ntilde) const
{	//Compute LJatt weighted density:
	DataRptr Nbar = I(Ntilde[0]*fex_LJatt), grad_Nbar; nullToZero(grad_Nbar, gInfo);
	//Evaluated weighted density functional:
	DataRptr Aex(DataR::alloc(gInfo,isGpuEnabled()));
	#ifdef GPU_ENABLED
	Fex_H20_ScalarEOS_gpu(gInfo.nr, Nbar->dataGpu(), Aex->dataGpu(), grad_Nbar->dataGpu(), *eval);
	#else
	threadedLoop(eval, gInfo.nr, Nbar->data(), Aex->data(), grad_Nbar->data());
	#endif
	//Convert gradients:
	DataRptr NO = I(Ntilde[0]);
	grad_Ntilde[0] += fex_LJatt*Idag(NO*grad_Nbar) + Idag(Aex);
	return gInfo.dV*dot(NO,Aex);
}

double Fex_H2O_ScalarEOS::computeUniform(const double* N, double* grad_N) const
{	double AexPrime, Aex;
	(*eval)(0, &N[0], &Aex, &AexPrime);
	grad_N[0] += Aex + N[0]*AexPrime;
	return N[0]*Aex;
}
