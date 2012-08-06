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

#ifndef JDFTX_FLUID_FEX_H2O_LISCHNER10_H
#define JDFTX_FLUID_FEX_H2O_LISCHNER10_H
#include <fluid/Fex.h>

class Fex_H2O_Lischner10 : public Fex
{
public:
	//! Create water with the excess functional from:
	//! J. Lischner and T. A. Arias, J. Phys. Chem. B 114, 1946 (2010).
	Fex_H2O_Lischner10(FluidMixture& fluidMixture);

	const Molecule* getMolecule() const { return &molecule; }
	double get_aDiel() const;
	double compute(const DataGptr* Ntilde, DataGptr* grad_Ntilde) const;
	double computeUniform(const double* N, double* grad_N) const;
private:
	RealKernel COO, COH, CHH, fex_gauss, siteChargeKernel;
	SiteProperties propO;
	SiteProperties propH;
	Molecule molecule;
};


#endif // JDFTX_FLUID_FEX_H2O_LISCHNER10_H
