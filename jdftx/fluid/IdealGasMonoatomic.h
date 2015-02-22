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

#ifndef JDFTX_FLUID_IDEALGASMONOATOMIC_H
#define JDFTX_FLUID_IDEALGASMONOATOMIC_H

#include <fluid/IdealGas.h>

//! IdealGas for monoatomic molecules (i.e. no orientation integral)
class IdealGasMonoatomic : public IdealGas
{
public:
	//!Initialize and associate with excess functional fex (and its fluid mixture)
	IdealGasMonoatomic(const FluidMixture*, const FluidComponent*);

	void initState(const ScalarField* Vex, ScalarField* psi, double scale, double Elo, double Ehi) const;
	void getDensities(const ScalarField* psi, ScalarField* N, vector3<>& P0) const;
	double compute(const ScalarField* psi, const ScalarField* N, ScalarField* Phi_N, const double Nscale, double& Phi_Nscale) const;
	void convertGradients(const ScalarField* psi, const ScalarField* N, const ScalarField* Phi_N, const vector3<>&  Phi_P0, ScalarField* Phi_psi, const double Nscale) const;
};

#endif // JDFTX_FLUID_IDEALGASMONOATOMIC_H
