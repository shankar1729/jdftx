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

#include <fluid/Fex.h>
#include <core/DataCollection.h>
#include <float.h>

class FluidMixture;

//! Abstract base class for an IdealGas evaluator
class IdealGas
{
public:
	const int nIndep; //!< Number of scalar fields used as independent variables
	const Molecule *const molecule; //!< Associated molecule geometry
	FluidMixture& fluidMixture; //!< Fluid mixture that this has been registered with
	const GridInfo& gInfo; //!< grid specifications
	const double T; //!< temperature

	DataRptrCollection V; //!< external site potentials

	//! Initialize and register to be used with excess functional fex in its fluidMixture
	//! xBulk is the mole fraction of this component in the bulk
	//! (If all the mole fractions don't add to 1, they will be normalized to do so)
	IdealGas(int nIndep, Fex* fex, double xBulk);

	//! Create an initial guess for the indep, in presence of V and the extra potential Vex
	//! The initial guess is typically taken to be scale times what would generate the equilibrium
	//! ideal gas density upto caps Elo and Ehi on the molecule energy configurations considered.
	//! This would also be a good place to logPrintf useful statistics about the potential for debugging
	virtual void initState(const DataRptr* Vex, DataRptr* indep, double scale, double Elo=-DBL_MAX, double Ehi=+DBL_MAX) const=0;

	//! Given the independent variables indep, compute the site densities N and the cell dipole moment P (not electric, orientation moment)
	virtual void getDensities(const DataRptr* indep, DataRptr* N, vector3<>& P) const=0;

	//! Return the ideal gas free energy PhiNI = T Int N - T S + (V-mu).N (where S is implementation dependent)
	//! and accumulate the gradients w.r.t the site densities
	//! Also accumulate the gradient w.r.t cell dipole moment, grad_P, given P (only required for multipole based ideal gases)
	//! Nscale is the factor by which the site densities/moments were scaled after getDensities() in order
	//! to implement fixed N / charge neutrality. Accumulate explicit gradients of PhiNI w.r.t Nscale in grad_Nscale;
	//! the implicit dependence through N is handled by FluidMixture.
	virtual double compute(const DataRptr* indep, const DataRptr* N, DataRptr* grad_N,
		const vector3<>& P, vector3<>& grad_P, const double Nscale, double& grad_Nscale) const=0;

	//! Compute grad_indep, the total gradients w.r.t indep, given th egradients of the entire
	//! functional w.r.t site densities, grad_N,  and total cell dipole (orientation, not electric) moment grad_P.
	//! Nscale will be the same as in compute()
	virtual void convertGradients(const DataRptr* indep, const DataRptr* N,
		const DataRptr* grad_N, vector3<> grad_P, DataRptr* grad_indep, const double Nscale) const=0;

	//! If Nnorm>=0, this component is switched to the cananoical ensemble (number fixed to Nnorm)
	void set_Nnorm(double Nnorm);

	double get_Nbulk();

	//! Override the values for bulk density and chemical potential set by fluidMixture::setPressure()
	void overrideBulk(double Nbulk, double mu);

private:
	double xBulk; //!< bulk mole fraction (normalized by FluidMixture::setPressure())
	double Nnorm; //!< If non-zero, constrain total number of molecules to be this value (constraint handled by fluidMixture)

protected:
	double Nbulk; //!< equilibirum denisty of this molecule in the bulk mixture
	double mu; //!< chemical potential for this molecule
	friend class FluidMixture; //!< FluidMixture::setPressure() adjusts Nbulk and mu to get target pressure and mole fractions
};

#endif // JDFTX_FLUID_IDEALGAS_H
