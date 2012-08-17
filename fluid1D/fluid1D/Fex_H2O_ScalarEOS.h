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

#ifndef FLUID1D_FLUID1D_FEX_H2O_SCALAREOS_H
#define FLUID1D_FLUID1D_FEX_H2O_SCALAREOS_H

#include <fluid1D/Fex.h>

class Fex_H2O_ScalarEOS : public Fex
{
public:
	//! Create water with the ScalarEOS functional (can choose soft or hard sphere version)
	Fex_H2O_ScalarEOS(FluidMixture& fluidMixture);

	const Molecule* getMolecule() const { return &molecule; }
	double get_aDiel() const;
	double compute(const ScalarFieldTilde* Ntilde, ScalarFieldTilde* grad_Ntilde) const;
	double computeUniform(const double* N, double* grad_N) const;
	void directCorrelations(const double* N, ScalarFieldTildeCollection& C) const;
private:
	SphericalKernel fex_LJatt, siteChargeKernel;
	struct ScalarEOS_eval* eval;
	SiteProperties propO;
	SiteProperties propH;
	Molecule molecule;
};


#endif // FLUID1D_FLUID1D_FEX_H2O_SCALAREOS_H
