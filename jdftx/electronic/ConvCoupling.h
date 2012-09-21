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
	ConvCoupling(FluidMixture& fluidMixture);
	
	//! Destructor
	~ConvCoupling();
	
	double computeUniform(const std::vector<double>& N, std::vector<double>& grad_N) const;
	double compute(const DataGptrCollection& Ntilde, DataGptrCollection& grad_Ntilde) const;
	string getName() const;
	
	//! Set the electron density kernel for a particular site to an exponential
	void setExponentialKernel(SiteProperties& s);
	
	//! Set the electron density kernel for a particular site to an exponential
	void setExpCusplessKernel(SiteProperties& s);
	
	//! Set the electron density kernel for a particular site from a binary file (must be created on correct grid)
	void setBinaryKernel(SiteProperties& s);
	
	//! Set the electron density kernel for a particular site from a radial function (interpolates onto grid)
	void setRadialKernel(SiteProperties& s);
	
	//! Set explicit system properties
	//! @param nCavityTilde "Cavity-effective" density of the explicit system (explicit electrons + chargeball)
	void setExplicit(const DataGptr& nCavityTilde);
	
	//! Set explicit system properties
	//! @param nCavity "Cavity-effective" density of the explicit system (explicit electrons + chargeball)
	void setExplicit(const DataRptr& nCavity);

	//! Compute coupling and gradients (for the electronic side)
	double computeElectronic(DataGptr* grad_nCavityTilde=0);

	//! For debugging: Output any coupling-related debug info here
	//! (replace %s in filename pattern with something meaningful and output to that filename)
	//!   It will be called from the electronic code during FluidDebug dump, if enabled
	//! DO NOT hack up any of the other functions!
	void dumpDebug(const char* filenamePattern) const;
	
	//DataRptr nFluid;
	//DataRptr nCavity;
	//DataRptr grad_nCavity;
	const ExCorr* exCorr;
	
	FluidMixture::ConvCouplingData CouplingData;
	
	private:
	DataRptr nCavity;
	
};

#endif // JDFTX_ELECTRONIC_CONVCOUPLING_H
