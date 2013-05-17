/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman, Kendra Letchworth Weaver 

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

#include <electronic/ExCorr.h>
#include <fluid/FluidMixture.h>
#include <fluid/Molecule.h>
#include <core/DataIO.h>
#include <core/Data.h>
#include <electronic/operators.h>
#include <electronic/VanDerWaals.h>
#include <electronic/VDWCoupling.h>

VDWCoupling::VDWCoupling(FluidMixture* fluidMixture, const std::shared_ptr<VanDerWaals>& vdW, double vdwScale)
: Fmix(fluidMixture), vdW(vdW), vdwScale(vdwScale)
{
	//create a list of fluid atomic numbers to calculate vdW interactions
	const std::vector<const FluidComponent*>& component = fluidMixture->getComponents();
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const FluidComponent& c = *component[ic];
		for(unsigned i=0; i<c.molecule.sites.size(); i++)
			atomicNumber.push_back(c.molecule.sites[i]->atomicNumber);
	}
}

double VDWCoupling::computeUniform(const std::vector< double >& N, std::vector< double >& Phi_N) const
{	return 0.; //No electronic system to couple to in the bulk fluid
}

double VDWCoupling::compute(const DataGptrCollection& Ntilde, DataGptrCollection& Phi_Ntilde) const
{	return energyAndGrad(Ntilde, &Phi_Ntilde);
}

double VDWCoupling::energyAndGrad(const DataGptrCollection& Ntilde, DataGptrCollection* Phi_Ntilde, IonicGradient* forces) const
{	return vdW->energyAndGrad(Ntilde, atomicNumber, vdwScale, Phi_Ntilde, forces);
}

string VDWCoupling::getName() const
{	return "VDWCoupling";
}

