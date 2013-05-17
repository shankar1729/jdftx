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

#ifndef JDFTX_FLUID_IDEALGAS_H
#define JDFTX_FLUID_IDEALGAS_H

#include <fluid/Molecule.h>
#include <core/DataCollection.h>
class FluidMixture;
class FluidComponent;

//! Abstract base class for an IdealGas evaluator
class IdealGas
{
public:
	const int nIndep; //!< Number of scalar fields used as independent variables
	const Molecule& molecule; //!< Associated molecule geometry
	const GridInfo& gInfo; //!< grid specifications
	const double T; //!< temperature

	DataRptrCollection V; //!< external site potentials

	//! Initialize and register to be used with excess functional fex in its fluidMixture
	IdealGas(int nIndep, const FluidMixture*, const FluidComponent*);
	virtual ~IdealGas() {}

	//! Create an initial guess for the indep, in presence of V and the extra potential Vex
	//! The initial guess is typically taken to be scale times what would generate the equilibrium
	//! ideal gas density upto caps Elo and Ehi on the molecule energy configurations considered.
	//! This would also be a good place to logPrintf useful statistics about the potential for debugging
	virtual void initState(const DataRptr* Vex, DataRptr* indep, double scale, double Elo=-DBL_MAX, double Ehi=+DBL_MAX) const=0;

	//! Given the independent variables indep, compute the site densities N and G=0 component of polarization density P
	virtual void getDensities(const DataRptr* indep, DataRptr* N, vector3<>& P0) const=0;

	//! Return the ideal gas free energy PhiNI = T Int N - T S + (V-mu).N (where S is implementation dependent)
	//! and accumulate the gradients w.r.t the site densities
	//! Nscale is the factor by which the site densities/moments were scaled after getDensities() in order
	//! to implement fixed N / charge neutrality. Accumulate explicit gradients of PhiNI w.r.t Nscale in Phi_Nscale;
	//! the implicit dependence through N is handled by FluidMixture.
	virtual double compute(const DataRptr* indep, const DataRptr* N, DataRptr* Phi_N, const double Nscale, double& Phi_Nscale) const=0;

	//! Compute Phi_indep, the total gradients w.r.t indep, given th egradients of the entire
	//! functional w.r.t site densities, Phi_N,  and polarization density G=0 Phi_P0.
	//! Nscale will be the same as in compute()
	virtual void convertGradients(const DataRptr* indep, const DataRptr* N, const DataRptr* Phi_N, const vector3<>& Phi_P0, DataRptr* Phi_indep, const double Nscale) const=0;

	double get_Nbulk();

	//! Override the values for bulk density and chemical potential set by fluidMixture::initialize()
	void overrideBulk(double Nbulk, double mu);

protected:
	double Nbulk; //!< equilibirum density of this molecule in the bulk mixture
	double mu; //!< chemical potential for this molecule
	double corrPrefac; //!< prefactor for dipolar rotational correlations
	friend class FluidMixture; //!< FluidMixture::initialize() adjusts Nbulk and mu to get target pressure and mole fractions
};

#endif // JDFTX_FLUID_IDEALGAS_H
