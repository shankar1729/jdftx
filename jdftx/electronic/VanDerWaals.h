
/*-------------------------------------------------------------------
Copyright 2012 Deniz Gunceler, Kendra Letchworth Weaver

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

#include <core/ScalarFieldArray.h>
#include <core/Coulomb.h>

//! @addtogroup LongRange
//! @{

//! Abstract base class of pair-potential dispersion correction methods.
class VanDerWaals
{
public:
	VanDerWaals(const Everything& e) : e(e) {}
	virtual ~VanDerWaals() {}
	
	//! Retrieve the scale factor for a specified exchange-correlation functional.
	//! Use scaleOverride, if supplied and supported, to override the functional's default. 
	//! Quit with an appropriate error message if functional is not parametrized.
	virtual double getScaleFactor(string exCorrName, double scaleOverride=0.) const=0;

	//! Van der Waal correction energy for a collection of discrete atoms at fixed locations.
	//! Corresponding forces should be accumulated to Atom::force for each atom.
	//! If E_RRT is non-null, accumulate contributions to the symmetric lattice derivative (stress * volume).
	virtual double energyAndGrad(std::vector<Atom>& atoms, const double scaleFac, matrix3<>* E_RRT=0) const=0;
	
protected:
	const Everything& e;
};

//! @}
#endif // JDFTX_ELECTRONIC_VANDERWAALS_H
