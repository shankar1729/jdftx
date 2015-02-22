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

#ifndef JDFTX_FLUID_FEX_H
#define JDFTX_FLUID_FEX_H

#include <fluid/Molecule.h>
#include <core/ScalarField.h>
class FluidMixture;
class FluidComponent;

//! Abstract base class for excess functionals
class Fex
{
public:
	const Molecule& molecule;
	const GridInfo& gInfo;
	const double T;

	Fex(const FluidMixture*, const FluidComponent*);
	virtual ~Fex() {}

	//! Return the excess free energy given the reciprocal space site densities
	//! and accumulate the gradient (functional derivative) w.r.t them in Phi_Ntilde
	virtual double compute(const ScalarFieldTilde* Ntilde, ScalarFieldTilde* Phi_Ntilde) const=0;

	//! Return the uniform fluid excess free energy density given the site densities N
	//! and accumulate the derivative w.r.t them in Phi_N. This MUST return the
	//! result corresponding to calling compute() with a uniform scalar field.
	//! This is called several times during FluidMixture::initialize() to get the desired bulk properties
	virtual double computeUniform(const double* N, double* Phi_N) const=0;
};

#endif // JDFTX_FLUID_FEX_H
