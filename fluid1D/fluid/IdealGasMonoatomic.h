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

#ifndef FLUID1D_FLUID1D_IDEALGASMONOATOMIC_H
#define FLUID1D_FLUID1D_IDEALGASMONOATOMIC_H

#include <fluid/IdealGas.h>

//! IdealGas for monoatomic molecules (i.e. no orientation integral)
class IdealGasMonoatomic : public IdealGas
{
public:
	//!Initialize and associate with excess functional fex (and its fluid mixture)
	IdealGasMonoatomic(Fex* fex, double xBulk);

	void initState(const ScalarField* Vex, ScalarField* psi, double scale, double Elo, double Ehi) const;
	void getDensities(const ScalarField* psi, ScalarField* N, double& P) const;
	double compute(const ScalarField* psi, const ScalarField* N, ScalarField* grad_N,
		const double& P, double& grad_P, const double Nscale, double& grad_Nscale) const;
	void convertGradients(const ScalarField* psi, const ScalarField* N,
		const ScalarField* grad_N, double grad_P, ScalarField* grad_psi, const double Nscale) const;
};

#endif // FLUID1D_FLUID1D_IDEALGASMONOATOMIC_H
