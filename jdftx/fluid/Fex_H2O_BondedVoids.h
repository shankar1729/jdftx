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

#ifndef JDFTX_FLUID_FEX_H2O_BONDEDVOIDS_H
#define JDFTX_FLUID_FEX_H2O_BONDEDVOIDS_H
#include <fluid/Fex.h>

class Fex_H2O_BondedVoids : public Fex
{
public:
	//! Create water with the BondedVoids functional as per [arXiv:1112.1442]
	//! which includes a minor update (in the site charge kernel) to
	//! [R. Sundararaman et al., J. Chem. Phys 137, 044107 (2012)]
	Fex_H2O_BondedVoids(const FluidMixture*, const FluidComponent*);
    virtual ~Fex_H2O_BondedVoids();
	
	double compute(const ScalarFieldTilde* Ntilde, ScalarFieldTilde* Phi_Ntilde) const;
	double computeUniform(const double* N, double* Phi_N) const;

	static const double RV0, TV, kappa, RO, sigmaU; //!< Functional parameters
private:
	RadialFunctionG Ua;
};


#endif // JDFTX_FLUID_FEX_H2O_BONDEDVOIDS_H
