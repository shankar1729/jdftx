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

#ifndef JDFTX_FLUID_IDEALGASMUEPS_H
#define JDFTX_FLUID_IDEALGASMUEPS_H

#include <fluid/IdealGas.h>
#include <fluid/SO3quad.h>
#include <fluid/TranslationOperator.h>

//! IdealGas for polyatomic molecules with the monopole-dipole 'MuEps' independent variables
class IdealGasMuEps : public IdealGas
{
public:
	vector3<> Eexternal; //!< External uniform electric field (the Mu-Eps fluid can be uniformly polarized!)

	//!Initialize and associate with excess functional fex (and its fluid mixture)
	//!Also specify the orientation quadrature and translation operator used for the orientation integrals
	IdealGasMuEps(Fex* fex, double xBulk, const SO3quad& quad, const TranslationOperator& trans);

	void initState(const DataRptr* Vex, DataRptr* mueps, double scale, double Elo, double Ehi) const;
	void getDensities(const DataRptr* mueps, DataRptr* N, vector3<>& P) const;
	double compute(const DataRptr* mueps, const DataRptr* N, DataRptr* grad_N,
		const vector3<>& P, vector3<>& grad_P, const double Nscale, double& grad_Nscale) const;
	void convertGradients(const DataRptr* mueps, const DataRptr* N,
		const DataRptr* grad_N, vector3<> grad_P, DataRptr* grad_mueps, const double Nscale) const;

private:
	const SO3quad& quad; //!< quadrature for orientation integral
	const TranslationOperator& trans; //!< translation operator for orientation integral
	int site0mult; //!< multiplicity of the first site
	double S; //!< cache the entropy, because it is most efficiently computed during getDensities()
};

#endif // JDFTX_FLUID_IDEALGASMUEPS_H
