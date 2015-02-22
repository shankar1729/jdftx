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


#ifndef JDFTX_ELECTRONIC_VDWCOUPLING_H
#define JDFTX_ELECTRONIC_VDWCOUPLING_H

#include <fluid/Fmix.h>
#include <fluid/Molecule.h>
#include <fluid/FluidMixture.h>
#include <electronic/common.h>
#include <electronic/VanDerWaals.h>

//! Van der Waals coupling between atoms from electronic DFT and fluid density fields
class VDWCoupling : public Fmix
{	
public:
	VDWCoupling(FluidMixture* fluidMixture, const std::vector< std::vector< vector3<> > >& atpos, const std::shared_ptr<VanDerWaals>& vdW, double vdwScale); 
	
	//! Main energy and gradients function
	double energyAndGrad(const ScalarFieldTildeArray& Ntilde, ScalarFieldTildeArray* Phi_Ntilde=0, IonicGradient* forces=0) const;
	
	//Interface to fluid side (virtual functions from Fmix):
	double computeUniform(const std::vector<double>& N, std::vector<double>& Phi_N) const;
	double compute(const ScalarFieldTildeArray& Ntilde, ScalarFieldTildeArray& Phi_Ntilde) const;
	string getName() const;

private:
	const std::vector< std::vector< vector3<> > >& atpos;
	const std::shared_ptr<VanDerWaals>& vdW;
	double vdwScale;
	std::vector<int> atomicNumber; 
};

#endif // JDFTX_ELECTRONIC_VDWCOUPLING_H