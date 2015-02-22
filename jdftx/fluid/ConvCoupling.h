/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver

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


#ifndef JDFTX_ELECTRONIC_CONVCOUPLING_H
#define JDFTX_ELECTRONIC_CONVCOUPLING_H

#include <fluid/Fmix.h>
#include <fluid/Molecule.h>
#include <fluid/FluidMixture.h>
#include <electronic/common.h>

//! Convolution coupling between electrons and fluids
class ConvCoupling : public Fmix
{
public:
	ConvCoupling(FluidMixture* fluidMixture, const ExCorr& exCorr);
	
	//! Set explicit system properties
	//! @param nCavity "Cavity-effective" density of the explicit system (explicit electrons + chargeball)
	void setExplicit(const ScalarFieldTilde& nCavityTilde);

	//! Main energy and gradients function
	double energyAndGrad(const ScalarFieldTildeArray& Ntilde, ScalarFieldTildeArray* Phi_Ntilde=0, ScalarFieldTilde* Phi_nCavityTilde=0) const;
	
	//Interface to fluid side (virtual functions from Fmix):
	double computeUniform(const std::vector<double>& N, std::vector<double>& Phi_N) const;
	double compute(const ScalarFieldTildeArray& Ntilde, ScalarFieldTildeArray& Phi_Ntilde) const;
	string getName() const;
	double Vxc_bulk;
	const ExCorr& exCorr;

private:
	const std::vector<const FluidComponent*>& component;
	ScalarField nCavity;
};

#endif // JDFTX_ELECTRONIC_CONVCOUPLING_H
