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

#include <fluid1D/Fex_LJ.h>
#include <fluid1D/FluidMixture.h>
#include <core1D/Operators.h>
#include <gsl/gsl_sf.h>

void setLJatt(SphericalKernel& kernel, const GridInfo& gInfo, double eps, double sigma)
{	kernel.resize(gInfo.S);
	for(int i=0; i<gInfo.S; i++)
	{	//Scale k,r so that well minimum is located at r=1, with depth 1:
		double kScale = sigma*pow(2.0,1.0/6);
		double k = gInfo.G[i]*kScale;
		double kSq = k*k;
		double result;
		if(k<35.0)
		{	result =
				(M_PI/1814400)*(
					- 2*cos(k)*(-564480 + kSq*(301680 + kSq*(24 + kSq*(-2 + kSq))))
					- 2*k*sin(k)*(297360 + kSq*(120 + kSq*(-6 + kSq)))
					+ pow(k,3)*(302400 + pow(k,6))*(M_PI - 2*gsl_sf_Si(k)) )
				+ (4*M_PI/3)*((2.2*gsl_sf_bessel_jl(0,k)) + gsl_sf_bessel_jl(2,k));
		}
		else //Use asymptotic form for large k
		{	result = (288*M_PI)*pow(k,-10) * (
				cos(k)*(-23950080 + kSq*(75880 +kSq*(-287 + kSq)))
				+ k*sin(k)*(1315160 + kSq*(-4585 + kSq*18)) );
		}
		//Scale back to actual values:
		kernel[i] = -eps*pow(kScale,3) * result;
	}
}

void setQkernel(SphericalKernel& kernel, const GridInfo& gInfo, double Rvdw)
{	kernel.resize(gInfo.S);
	for(int i=0; i<gInfo.S; i++)
		kernel[i] = gsl_sf_bessel_j0(gInfo.G[i]*Rvdw);
}

Fex_LJ::Fex_LJ(FluidMixture& fluidMixture, double eps, double sigma, string name, double Q)
: Fex(fluidMixture),
	eps(eps), sigma(sigma),
	Qkernel(Q ? new SphericalKernel : 0),
	prop(gInfo,
		 //See ~/DFT/Water1D/LJfluid-coreSoftness.nb for erf model fits
		 //to the Mayer function for the repulsive part of LJ
		 0.5*sigma*(1.1063165169463278 - 0.8595214683635807*exp(-pow(eps/T,0.19185268827189367)/0.44711718709043674)),
		 0.5*sigma*(0.027545888769872514 + 0.07109941355697276*exp(-pow(eps/T,0.5774162302760579)/1.3735105148356759)),
		 Q,Qkernel),
	molecule(name, &prop,  vector3<>(0,0,0))
{	logPrintf("\tInitializing LJ fluid \"%s\" with eps=%lf K and sigma=%lf A:\n", name.c_str(), eps/Kelvin, sigma/Angstrom);
	logPrintf("\t\tCore radius %lf A with softness %lf A\n", prop.sphereRadius/Angstrom, prop.sphereSigma/Angstrom);
	setLJatt(ljatt, gInfo, eps, sigma);
	if(Qkernel) setQkernel(*Qkernel, gInfo, 0.5*sigma);
}
Fex_LJ::~Fex_LJ()
{	if(Qkernel) delete Qkernel;
}

double Fex_LJ::compute(const ScalarFieldTilde* Ntilde, ScalarFieldTilde* grad_Ntilde) const
{	ScalarFieldTilde V = O(ljatt * Ntilde[0]);
	grad_Ntilde[0] += V;
	return 0.5*dot(V, Ntilde[0]);
}
double Fex_LJ::computeUniform(const double* N, double* grad_N) const
{	grad_N[0] += ljatt.data()[0]*N[0];
	return 0.5*N[0]*ljatt.data()[0]*N[0];
}
void Fex_LJ::directCorrelations(const double* N, ScalarFieldTildeCollection& C) const
{	ScalarFieldTilde ljattTilde(ljatt, gInfo);
	C[fluidMixture.corrFuncIndex(0,0,this)] += ljattTilde;
}


Fmix_LJ::Fmix_LJ(const Fex_LJ& fluid1, const Fex_LJ& fluid2)
: Fmix(fluid1.fluidMixture), fluid1(fluid1), fluid2(fluid2)
{	setLJatt(ljatt, gInfo, sqrt(fluid1.eps*fluid2.eps), 0.5*(fluid1.sigma+fluid2.sigma));
}

string Fmix_LJ::getName() const
{	return fluid1.molecule.name + "<->" + fluid2.molecule.name;
}

double Fmix_LJ::compute(const ScalarFieldTildeCollection& Ntilde, ScalarFieldTildeCollection& grad_Ntilde) const
{	unsigned i1 = fluidMixture.get_offsetDensity(&fluid1);
	unsigned i2 = fluidMixture.get_offsetDensity(&fluid2);
	ScalarFieldTilde V1 = O(ljatt * Ntilde[i1]);
	ScalarFieldTilde V2 = O(ljatt * Ntilde[i2]);
	grad_Ntilde[i1] += V2;
	grad_Ntilde[i2] += V1;
	return dot(V1,Ntilde[i2]);
}
double Fmix_LJ::computeUniform(const std::vector<double>& N, std::vector<double>& grad_N) const
{	unsigned i1 = fluidMixture.get_offsetDensity(&fluid1);
	unsigned i2 = fluidMixture.get_offsetDensity(&fluid2);
	grad_N[i1] += ljatt.data()[0]*N[i2];
	grad_N[i2] += ljatt.data()[0]*N[i1];
	return N[i1]*ljatt.data()[0]*N[i2];
}
void Fmix_LJ::directCorrelations(const std::vector<double>& N, ScalarFieldTildeCollection& C) const
{	ScalarFieldTilde ljattTilde(ljatt, gInfo);
	unsigned i1 = fluidMixture.get_offsetDensity(&fluid1);
	unsigned i2 = fluidMixture.get_offsetDensity(&fluid2);
	C[fluidMixture.corrFuncIndex(i1,i2)] += ljattTilde;
}

