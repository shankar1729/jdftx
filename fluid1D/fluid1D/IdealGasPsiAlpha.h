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

#ifndef FLUID1D_FLUID1D_IDEALGASPSIALPHA_H
#define FLUID1D_FLUID1D_IDEALGASPSIALPHA_H

#include <fluid1D/IdealGas.h>
#include <fluid/SO3quad.h>
#include <fluid1D/TranslationOperator.h>

//! IdealGas for polyatomic molecules with the effective potential 'psi_alpha' independent variables
class IdealGasPsiAlpha : public IdealGas
{
public:
	//!Initialize and associate with excess functional fex (and its fluid mixture)
	//!Also specify the orientation quadrature and translation operator used for the orientation integrals
	IdealGasPsiAlpha(Fex* fex, double xBulk, const SO3quad& quad, const TranslationOperator& trans);

	void initState(const ScalarField* Vex, ScalarField* psi, double scale, double Elo, double Ehi) const;
	void getDensities(const ScalarField* psi, ScalarField* N, double& P) const;
	double compute(const ScalarField* psi, const ScalarField* N, ScalarField* grad_N,
		const double& P, double& grad_P, const double Nscale, double& grad_Nscale) const;
	void convertGradients(const ScalarField* psi, const ScalarField* N,
		const ScalarField* grad_N, double grad_P, ScalarField* grad_psi, const double Nscale) const;

private:
	const SO3quad& quad; //!< quadrature for orientation integral
	const TranslationOperator& trans; //!< translation operator for orientation integral
	std::vector<int> siteToIndep; //!< index of independent variable by site number
	std::vector<int> densityToIndep; //!< index of independent variable by site density index
	std::vector<int> indepMult; //!< multiplicity (number of sites) of each independent variable
	int indepMultTot; //!< sum of multiplicities of indep (i.e. number of real sites on molecule)
};

#endif // FLUID1D_FLUID1D_IDEALGASPSIALPHA_H
