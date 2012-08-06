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

#ifndef JDFTX_FLUID_FEX_LJ_H
#define JDFTX_FLUID_FEX_LJ_H

//! @file Fex_LJ.h
//! @brief Lennard-Jones fluids

#include <fluid/Fex.h>
#include <fluid/Fmix.h>

//! Initialize kernel to the attarctive part of a Lennard-Jones potential
void setLJatt(RealKernel& kernel, double eps, double sigma);


//! Lennard Jones fluid treated as a mean field perturbation about a soft FMT core
class Fex_LJ : public Fex
{
public:

	//! Create a fluid of Lennard-Jones particles with well-depth eps and diameter sigma
	//! name sets the label of the single site of the molecule
	//! The particles can optionally be charged, the charge kernel is set to a gaussian of width equal to the core radius
	Fex_LJ(FluidMixture& fluidMixture, double eps, double sigma, string name, double Q=0.0);
	~Fex_LJ();

	const Molecule* getMolecule() const { return &molecule; }
	double get_aDiel() const { return 1.0; };
	double compute(const DataGptr* Ntilde, DataGptr* grad_Ntilde) const;
	double computeUniform(const double* N, double* grad_N) const;
private:
	friend class Fmix_LJ; //allow Fmix_LJ to peek at eps and sigma to select coupling parameters
	double eps, sigma;
	RealKernel* Qkernel;
	SiteProperties prop;
	Molecule molecule;
	RealKernel ljatt;
};

//! Interaction functional for two Lennard-Jones fluids
class Fmix_LJ : public Fmix
{
public:
	//! Add a lennard-jones coupling between two LJ fluid
	//! with parameters set using the Lorentz-Berthelot Mixing Rules
	Fmix_LJ(const Fex_LJ& fluid1, const Fex_LJ& fluid2);

	string getName() const;

	double compute(const DataGptrCollection& Ntilde, DataGptrCollection& grad_Ntilde) const;
	double computeUniform(const std::vector<double>& N, std::vector<double>& grad_N) const;
private:
	const Fex_LJ &fluid1, &fluid2;
	RealKernel ljatt;
};

#endif // JDFTX_FLUID_FEX_LJ_H
