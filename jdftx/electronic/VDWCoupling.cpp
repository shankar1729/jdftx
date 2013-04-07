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

VDWCoupling::VDWCoupling(FluidMixture& fluidMixture, std::shared_ptr< VanDerWaals > vdW)
	: Fmix(fluidMixture), vdW(vdW)
{
	//create a list of fluid atomic numbers to calculate vdW interactions
	
	for(unsigned ic=0; ic<fluidMixture.get_nComponents(); ic++)
		{
			const FluidMixture::Component& c = fluidMixture.get_component(ic);
			for(int j=0; j<c.molecule->nIndices; j++)
			{
				const SiteProperties& s = *c.indexedSite[j];
				atomicNumber.push_back(s.atomicNumber); 
			}
		}
}

VDWCoupling::~VDWCoupling()
{
}

double VDWCoupling::computeUniform(const std::vector< double >& N, std::vector< double >& grad_N) const
{
	return 0.0; //No electronic system to couple to in the bulk fluid
}

double VDWCoupling::compute(const DataGptrCollection& Ntilde, DataGptrCollection& grad_Ntilde) const
{	double vdwScale = vdW->getScaleFactor(exCorr->getName(), *scaleFac);
	return vdW->energyAndGrad(Ntilde, atomicNumber, vdwScale, &grad_Ntilde);
}

double VDWCoupling::computeElectronic(const DataRptrCollection* N, IonicGradient* forces) 
{
	DataGptrCollection Ntilde(fluidMixture.get_nDensities());
	for(unsigned i=0; i<fluidMixture.get_nDensities(); i++)
		Ntilde[i] = J((*N)[i]);
	
	double vdwScale = vdW->getScaleFactor(exCorr->getName(), *scaleFac);
	return vdW->energyAndGrad(Ntilde, atomicNumber, vdwScale, 0, forces);
}

void VDWCoupling::dumpDebug(const char* filenamePattern) const
{

}


string VDWCoupling::getName() const
{	
	return "VDWCoupling";
}

