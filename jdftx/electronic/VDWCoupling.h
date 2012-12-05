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
	VDWCoupling(FluidMixture& fluidMixture, std::shared_ptr<VanDerWaals> vdW); 
	
	//! Destructor
	~VDWCoupling();
	
	double computeUniform(const std::vector<double>& N, std::vector<double>& grad_N) const;
	
	double compute(const DataGptrCollection& Ntilde, DataGptrCollection& grad_Ntilde) const;

	string getName() const;
	
	//! Compute ionic gradients (for the electronic side) and return van der Waals coupling energy
	double computeElectronic(const DataRptrCollection* N, IonicGradient* forces=0);

	//! For debugging: Output any coupling-related debug info here
	//! (replace %s in filename pattern with something meaningful and output to that filename)
	//!   It will be called from the electronic code during FluidDebug dump, if enabled
	//! DO NOT hack up any of the other functions!
	void dumpDebug(const char* filenamePattern) const;
	
	const ExCorr* exCorr;
	
	private:
		
	std::shared_ptr<VanDerWaals> vdW;
	std::vector<int> atomicNumber; 
	
	
};

#endif // JDFTX_ELECTRONIC_VDWCOUPLING_H