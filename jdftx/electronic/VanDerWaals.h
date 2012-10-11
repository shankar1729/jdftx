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

#include <electronic/common.h>
#include <core/vector3.h>
#include <core/Util.h>
#include <core/Coulomb.h>

class VanDerWaals
{
public:
	void setup(const Everything &everything);
	
	//! Returns the Van der Waals energy and gradient (in cartesian coordinates)
	//! s6 is the scaling parameter that depends on the Ex-Corr functional used
	double energyAndGrad(std::vector<Atom>& atoms, string EXCorr);
	
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
};

#endif // JDFTX_ELECTRONIC_VANDERWAALS_H
