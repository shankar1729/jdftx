/*-------------------------------------------------------------------
Copyright 2014 Ravishankar Sundararaman, Kendra Letchworth-Weaver

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

#ifndef JDFTX_FLUID_FEX_LJ_H
#define JDFTX_FLUID_FEX_LJ_H

//! @file Fex_LJ.h
//! @brief Lennard-Jones fluids

#include <fluid/Fex.h>
#include <fluid/Fmix.h>

//! Initialize kernel to the attarctive part of a Lennard-Jones potential
void setLJatt(RadialFunctionG& kernel, const GridInfo& gInfo, double eps, double sigma);


//! Lennard Jones fluid treated as a mean field perturbation about a soft FMT core
class Fex_LJ : public Fex
{
public:

	//! Create a fluid of Lennard-Jones particles with well-depth eps
	//! The range parameter sigma is set to the hard sphere diameter of the first site
  Fex_LJ(const FluidMixture*, const FluidComponent*, double eps, double sigmaOverride = 0.0);
    virtual ~Fex_LJ();
	
	double compute(const ScalarFieldTilde* Ntilde, ScalarFieldTilde* Phi_Ntilde) const;
	double computeUniform(const double* N, double* Phi_N) const;
private:
	friend class Fmix_LJ; //allow Fmix_LJ to peek at eps and sigma to select coupling parameters
	double eps, sigma;
	RadialFunctionG ljatt;
};

//! Lennard-Jones interaction functional
class Fmix_LJ : public Fmix
{
public:
	//! Add a lennard-jones coupling between two LJ fluid
	//! Could set parameters using the Lorentz-Berthelot Mixing Rules if desired
	Fmix_LJ(FluidMixture*, std::shared_ptr<FluidComponent> fluid1, std::shared_ptr<FluidComponent> fluid2, double eps, double sigma);
	virtual ~Fmix_LJ();
	string getName() const;

	double compute(const ScalarFieldTildeArray& Ntilde, ScalarFieldTildeArray& Phi_Ntilde) const;
	double computeUniform(const std::vector<double>& N, std::vector<double>& Phi_N) const;
private:
	std::shared_ptr<FluidComponent> fluid1, fluid2;
	RadialFunctionG ljatt;
};

//! Gaussian Kernel interaction functional
class Fmix_GaussianKernel : public Fmix
{
public:
	//! Add a gaussian coupling potential between two fluids
	//!< mixing parameter Esolv: depth of interaction potential in hartree
	//!< mixing parameter Rsolv; width of gaussian kernel interaction potential in bohr
	Fmix_GaussianKernel(FluidMixture*, std::shared_ptr<FluidComponent> fluid1, std::shared_ptr<FluidComponent> fluid2, double Esolv, double Rsolv);
	virtual ~Fmix_GaussianKernel();
	string getName() const;

	double compute(const ScalarFieldTildeArray& Ntilde, ScalarFieldTildeArray& Phi_Ntilde) const;
	double computeUniform(const std::vector<double>& N, std::vector<double>& Phi_N) const;
private:
	std::shared_ptr<FluidComponent> fluid1, fluid2;
	RadialFunctionG Ksolv; //shape function for interaction
	double Kmul; //prefactor to interaction
	
};

#endif // JDFTX_FLUID_FEX_LJ_H
