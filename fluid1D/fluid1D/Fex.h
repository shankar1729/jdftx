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

#ifndef FLUID1D_FLUID1D_FEX_H
#define FLUID1D_FLUID1D_FEX_H

#include <fluid1D/Molecule.h>
#include <core1D/Data.h>
class FluidMixture;

//! Abstract base class for excess functionals
class Fex
{
public:
	FluidMixture& fluidMixture;
	const GridInfo& gInfo;
	const double T;

	//! Initialize base for use with this fluidMixture
	Fex(FluidMixture& fluidMixture);

	//! Return the underlying molecule speicfication (MUST be the same reference and value each time it's called)
	virtual const Molecule* getMolecule() const=0;

	virtual double get_aDiel() const=0; //! Return the correlation scale-factor for the coulomb term

	//! Return the excess free energy given the basis space site densities
	//! and accumulate the gradient (partial derivative) w.r.t them in grad_Ntilde
	virtual double compute(const ScalarFieldTilde* Ntilde, ScalarFieldTilde* grad_Ntilde) const=0;

	//! Return the uniform fluid excess free energy density given the site densities N
	//! and accumulate the derivative w.r.t them in grad_N. This MUST return the
	//! result corresponding to calling compute() with a uniform scalar field.
	//! This is called several times during FluidMixture::setPressure() to get the desired bulk properties
	virtual double computeUniform(const double* N, double* grad_N) const=0;
};

#endif // FLUID1D_FLUID1D_FEX_H
