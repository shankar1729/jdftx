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

#include <fluid/ConvCoupling.h>
#include <electronic/ExCorr.h>
#include <fluid/FluidMixture.h>
#include <fluid/Molecule.h>
#include <core/ScalarFieldIO.h>
#include <core/ScalarField.h>
#include <electronic/RadialFunction.h>
#include <electronic/operators.h>

ConvCoupling::ConvCoupling(FluidMixture* fluidMixture, const ExCorr& exCorr)
: Fmix(fluidMixture), exCorr(exCorr), component(fluidMixture->getComponents())
{
	Citations::add("Convolution-coupling for Joint Density Functional Theory",
		"K. Letchworth-Weaver, R. Sundararaman and T.A. Arias, (under preparation)");
}

double ConvCoupling::computeUniform(const std::vector<double>& N, std::vector<double>& Phi_N) const
{	return 0.; //No electronic systen to couple to in the bulk fluid
}

double ConvCoupling::compute(const ScalarFieldTildeArray& Ntilde, ScalarFieldTildeArray& Phi_Ntilde) const 
{	return energyAndGrad(Ntilde, &Phi_Ntilde);
}

string ConvCoupling::getName() const
{	return "ConvCoupling";
}

void ConvCoupling::setExplicit(const ScalarFieldTilde& nCavityTilde)
{	this->nCavity = I(nCavityTilde);
}

double ConvCoupling::energyAndGrad(const ScalarFieldTildeArray& Ntilde, ScalarFieldTildeArray* Phi_Ntilde, ScalarFieldTilde* Phi_nCavityTilde) const
{
	//Compute model electron density of fluid:
	ScalarFieldTilde nFluidTilde;
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const FluidComponent& c = *component[ic];
		for(unsigned i=0; i<c.molecule.sites.size(); i++)
		{	const Molecule::Site& s = *(c.molecule.sites[i]);
			if(s.elecKernel)
				nFluidTilde += s.elecKernel * Ntilde[c.offsetDensity+i];
		}
	}
	ScalarField nFluid = I(nFluidTilde);
	
	//Calculate exchange, correlation, and kinetic energy
	ScalarField nTot = nFluid + nCavity;
	ScalarField Vxc_tot, Vxc_fluid, Vxc_cavity;
	double Phi =
		+ exCorr(nTot, &Vxc_tot, true)
		- exCorr(nFluid, &Vxc_fluid, true)
		- exCorr(nCavity, &Vxc_cavity,  true);
	
	//Accumulate electronic-side gradient if required:
	if(Phi_nCavityTilde)
		*Phi_nCavityTilde += J(Vxc_tot - Vxc_cavity);
	
	//Accumulate fluid-side gradients if required:
	if(Phi_Ntilde)
	{	ScalarFieldTilde Phi_nFluidTilde = Idag(Vxc_tot - Vxc_fluid);
		//Propagate to fluid densities:
		for(unsigned ic=0; ic<component.size(); ic++)
		{	const FluidComponent& c = *component[ic];
			for(unsigned i=0; i<c.molecule.sites.size(); i++)
			{	const Molecule::Site& s = *(c.molecule.sites[i]);
				if(s.elecKernel)
					Phi_Ntilde->at(c.offsetDensity+i) += s.elecKernel * Phi_nFluidTilde;
			}
		}
	}
	
	return Phi;
}
