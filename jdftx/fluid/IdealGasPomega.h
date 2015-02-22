/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

#ifndef JDFTX_FLUID_IDEALGASPOMEGA_H
#define JDFTX_FLUID_IDEALGASPOMEGA_H

#include <fluid/IdealGas.h>
#include <fluid/SO3quad.h>
#include <fluid/TranslationOperator.h>
#include <core/VectorField.h>

//! IdealGas for polyatomic molecules with the orientation densities 'P_omega' as independent variables
//! This is also the base class for IdealGas's which used compressed representations of Pomega
class IdealGasPomega : public IdealGas
{
public:
	//!Initialize and associate with excess functional fex (and its fluid mixture)
	//!Also specify the orientation quadrature and translation operator used for the orientation integrals
	IdealGasPomega(const FluidMixture*, const FluidComponent*, const SO3quad& quad, const TranslationOperator& trans, unsigned nIndepOverride=0);

	void initState(const ScalarField* Vex, ScalarField* indep, double scale, double Elo, double Ehi) const;
	void getDensities(const ScalarField* indep, ScalarField* N, vector3<>& P0) const;
	double compute(const ScalarField* indep, const ScalarField* N, ScalarField* Phi_N, const double Nscale, double& Phi_Nscale) const;
	void convertGradients(const ScalarField* indep, const ScalarField* N, const ScalarField* Phi_N, const vector3<>& Phi_P0, ScalarField* Phi_indep, const double Nscale) const;

protected:
	const SO3quad& quad; //!< quadrature for orientation integral
	const TranslationOperator& trans; //!< translation operator for orientation integral
	vector3<> pMol; //!< molecule dipole moment in reference frame
	int oStart, oStop; //!< portion of orientation loop handled by current process
	
	virtual string representationName() const;
	
	//These functions are called once for each orientation:
	virtual void initState_o(int o, const matrix3<>& rot, double scale, const ScalarField& Eo, ScalarField* state) const;
	virtual void getDensities_o(int o, const matrix3<>& rot, const ScalarField* state, ScalarField& logPomega_o) const;
	virtual void convertGradients_o(int o, const matrix3<>& rot, const ScalarField& Phi_logPomega_o, ScalarField* Phi_state) const;
	
private:
	double S; //!< cache the entropy, because it is most efficiently computed during getDensities()
	double Ecorr; VectorField Ecorr_P; //!< cache the correlation correction and its derivatives, since they are most efficiently computed during getDensities()
};

#endif // JDFTX_FLUID_IDEALGASPOMEGA_H
