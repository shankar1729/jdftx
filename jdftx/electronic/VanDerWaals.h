/*-------------------------------------------------------------------
Copyright 2012 Deniz Gunceler

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

#ifndef JDFTX_ELECTRONIC_VANDERWAALS_H
#define JDFTX_ELECTRONIC_VANDERWAALS_H

#include <electronic/RadialFunction.h>
#include <electronic/common.h>
#include <core/DataCollection.h>
#include <core/Coulomb.h>

class VanDerWaals
{
public:
	void setup(const Everything &everything);
	
	//! Van der Waal correction energy for a collection of discrete atoms at fixed locations
	//! Corresponding forces are accumulated to Atom::force for each atom
	//! The specified exchange-correlation functional selects the appropriate scale factor
	double energyAndGrad(std::vector<Atom>& atoms, string exCorrName) const;
	
	//! Van der Waal correction to the interaction energy between the explicit atoms
	//! (from IonInfo) and the continuous fields Ntilde with specified atomic numbers.
	//! The specified exchange-correlation functional selects the appropriate scale factor
	//! The gradient w.r.t site densities is accumulated to grad_Ntilde (if non-null) and
	//! the negative gradient w.r.t discrete atom positions is accumulated to forces (if non-null)
	double energyAndGrad(const DataGptrCollection& Ntilde, const std::vector<int>& atomicNumber, string exCorrName,
		DataGptrCollection* grad_Ntilde=0, IonicGradient* forces=0) const;
	
	~VanDerWaals();
	
private:
	const Everything* e;
	
	//! C6 and R0 parameters for the VDW interactions
	struct AtomParams
	{	double C6;
		double R0;
		//!Construct given C6 in J-nm^6/mol and R0 in Angstrom
		AtomParams(double SI_C6=0., double SI_R0=0.);
	};
	std::vector<AtomParams> atomParams; //!< Returns the C6 coefficient and R0
	std::map<string,double> scalingFactor; //!< ExCorr dependent scale factor
	
	//! Get cached RadialFunctionG for interaction kernel between species of
	//! atomic numnbers Z1 and Z2. The radial function will be created (and cached)
	//! the first time this function is called with a specific pair of Z's.
	//! Note that both Z1 and Z2 should be supported non-zero atomic numbers.
	const RadialFunctionG& getRadialFunction(int Z1, int Z2) const;
	
	std::map<std::pair<int,int>,RadialFunctionG> radialFunctions;
};

#endif // JDFTX_ELECTRONIC_VANDERWAALS_H
