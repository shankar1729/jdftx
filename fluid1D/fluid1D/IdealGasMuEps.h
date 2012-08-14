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

#ifndef FLUID1D_FLUID1D_IDEALGASMUEPS_H
#define FLUID1D_FLUID1D_IDEALGASMUEPS_H

#include <fluid1D/IdealGas.h>
#include <fluid/SO3quad.h>
#include <fluid1D/TranslationOperator.h>

//! IdealGas for polyatomic molecules with the monopole-dipole 'MuEps' independent variables
class IdealGasMuEps : public IdealGas
{
public:
	double Eexternal; //!< External uniform electric field (the Mu-Eps fluid can be uniformly polarized!)

	//!Initialize and associate with excess functional fex (and its fluid mixture)
	//!Also specify the orientation quadrature and translation operator used for the orientation integrals
	IdealGasMuEps(Fex* fex, double xBulk, const SO3quad& quad, const TranslationOperator& trans);

	void initState(const ScalarField* Vex, ScalarField* mueps, double scale, double Elo, double Ehi) const;
	void getDensities(const ScalarField* mueps, ScalarField* N, double& P) const;
	double compute(const ScalarField* mueps, const ScalarField* N, ScalarField* grad_N,
		const double& P, double& grad_P, const double Nscale, double& grad_Nscale) const;
	void convertGradients(const ScalarField* mueps, const ScalarField* N,
		const ScalarField* grad_N, double grad_P, ScalarField* grad_mueps, const double Nscale) const;

private:
	const SO3quad& quad; //!< quadrature for orientation integral
	const TranslationOperator& trans; //!< translation operator for orientation integral
	int site0mult; //!< multiplicity of the first site
	double S; //!< cache the entropy, because it is most efficiently computed during getDensities()
};

#endif // FLUID1D_FLUID1D_IDEALGASMUEPS_H
