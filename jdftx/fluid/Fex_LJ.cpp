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

#include <fluid/Fex_LJ.h>
#include <fluid/FluidMixture.h>
#include <electronic/operators.h>
#include <core/Operators.h>
#include <gsl/gsl_sf.h>

double setLJatt_calc(double G, double eps, double sigma)
{	//Scale so that well minimum is located at r=1, with depth 1:
	double kScale = sigma*pow(2.0,1.0/6);
	double k = G*kScale;
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
	return -eps*pow(kScale,3) * result;
}
void setLJatt(RadialFunctionG& kernel, const GridInfo& gInfo, double eps, double sigma)
{	kernel.init(0, gInfo.dGradial, gInfo.GmaxGrid, setLJatt_calc, eps, sigma);
}


Fex_LJ::Fex_LJ(const FluidMixture* fluidMixture, const FluidComponent* comp, double eps, double sigmaOverride)
: Fex(fluidMixture, comp), eps(eps), sigma(2.*molecule.sites[0]->Rhs)
{
	if(sigmaOverride) sigma = sigmaOverride;
	logPrintf("     Initializing LJ excess functional with eps=%lf Eh and sigma=%lf bohrs\n", eps, sigma);
	setLJatt(ljatt, gInfo, eps, sigma); 
}
Fex_LJ::~Fex_LJ()
{	ljatt.free();
}

double Fex_LJ::compute(const ScalarFieldTilde* Ntilde, ScalarFieldTilde* Phi_Ntilde) const
{	ScalarFieldTilde V = gInfo.nr * (ljatt * Ntilde[0]);
	Phi_Ntilde[0] += V;
	return 0.5*gInfo.dV*dot(V,Ntilde[0]);
}
double Fex_LJ::computeUniform(const double* N, double* Phi_N) const
{	Phi_N[0] += ljatt(0)*N[0];
	return 0.5*N[0]*ljatt(0)*N[0];
}


Fmix_LJ::Fmix_LJ(FluidMixture* fluidMixture, std::shared_ptr<FluidComponent> fluid1, std::shared_ptr<FluidComponent> fluid2, double eps, double sigma)
: Fmix(fluidMixture), fluid1(fluid1), fluid2(fluid2)
{
          string name1 = fluid1->molecule.name;
          string name2 = fluid2->molecule.name;
	  logPrintf("\n     Initializing attractive LJ mixing functional between %s and %s\n		sigma: %lg Bohr and eps: %lg H.\n",
		     name1.c_str(),name2.c_str(),sigma,eps);
	  setLJatt(ljatt, gInfo, eps, sigma);
}

Fmix_LJ::~Fmix_LJ()
{	ljatt.free();
}

string Fmix_LJ::getName() const
{	return fluid1->molecule.name + "<->" + fluid2->molecule.name;
}

double Fmix_LJ::compute(const ScalarFieldTildeArray& Ntilde, ScalarFieldTildeArray& Phi_Ntilde) const
{	unsigned i1 = fluid1->offsetDensity;
	unsigned i2 = fluid2->offsetDensity;
	ScalarFieldTilde V1 = gInfo.nr * (ljatt * Ntilde[i1]);
	ScalarFieldTilde V2 = gInfo.nr * (ljatt * Ntilde[i2]);
	Phi_Ntilde[i1] += V2;
	Phi_Ntilde[i2] += V1;
	return gInfo.dV*dot(V1,Ntilde[i2]);
}
double Fmix_LJ::computeUniform(const std::vector<double>& N, std::vector<double>& Phi_N) const
{	unsigned i1 = fluid1->offsetDensity;
	unsigned i2 = fluid2->offsetDensity;
	Phi_N[i1] += ljatt(0)*N[i2];
	Phi_N[i2] += ljatt(0)*N[i1];
	return N[i1]*ljatt(0)*N[i2];
}


Fmix_GaussianKernel::Fmix_GaussianKernel(FluidMixture* fluidMixture, std::shared_ptr<FluidComponent> fluid1, std::shared_ptr<FluidComponent> fluid2, double Esolv, double Rsolv)
: Fmix(fluidMixture), fluid1(fluid1), fluid2(fluid2)
{
  string name1 = fluid1->molecule.name;
  string name2 = fluid2->molecule.name;
  logPrintf("     Initializing gaussian kernel mixing functional between %s and %s\n		Rsolv: %lg and Esolv: %lg.\n",name1.c_str(),name2.c_str(),Rsolv,Esolv);
  Kmul = -Esolv*(4*M_PI*pow(Rsolv,3))/3;
  Ksolv.init(0, gInfo.dGradial, gInfo.GmaxGrid, RadialFunctionG::gaussTilde, 1.0, Rsolv/sqrt(2));
  
}

Fmix_GaussianKernel::~Fmix_GaussianKernel()
{

}

string Fmix_GaussianKernel::getName() const
{
    return fluid1->molecule.name + "<->" + fluid2->molecule.name;
}

double Fmix_GaussianKernel::computeUniform(const std::vector< double >& N, std::vector< double >& Phi_N) const
{
		unsigned i1 = fluid1->offsetDensity;
		unsigned i2 = fluid2->offsetDensity;
		Phi_N[i1] += Kmul*Ksolv(0)*N[i2];
		Phi_N[i2] += Kmul*Ksolv(0)*N[i1];
		return N[i1]*Kmul*Ksolv(0)*N[i2];
}

double Fmix_GaussianKernel::compute(const ScalarFieldTildeArray& Ntilde, ScalarFieldTildeArray& Phi_Ntilde) const
{
		unsigned i1 = fluid1->offsetDensity;
		unsigned i2 = fluid2->offsetDensity;
		ScalarFieldTilde V1 = gInfo.nr * Kmul*(Ksolv * Ntilde[i1]);
		ScalarFieldTilde V2 = gInfo.nr * Kmul*(Ksolv * Ntilde[i2]);
		Phi_Ntilde[i1] += V2;
		Phi_Ntilde[i2] += V1;
		return gInfo.dV*dot(V1,Ntilde[i2]);
}
