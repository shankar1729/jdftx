/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

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

#ifndef JDFTX_FLUID_SO3QUAD_H
#define JDFTX_FLUID_SO3QUAD_H

//! @addtogroup so3quad
//! @{

/** @file SO3quad.h
@brief Quadratures for SO(3)/Zn
*/

#include <fluid/S2quad.h>

//! @brief Quadrature for SO(3)/Zn
class SO3quad
{
public:
	//! Initialize SO3/Zn quadrature from an S2 quadrature decsription
	SO3quad(const S2quad&, int Zn);
	
	//! Initialize SO3/Zn quadrature from an S2 quadrature name
	//! (nBeta, nAlpha and nGamma are used only for Euler quadrature)
	SO3quad(S2quadType type, int Zn, unsigned nBeta=0, unsigned nAlpha=0, unsigned nGamma=0);
	
	void print(); //!< Summarize euler angles and weights for the quadrature
	int nOrientations() const; //!< get cardinality of sampling set
	vector3<> euler(int iOrientation) const; //!< get euler angles for the iOrientation'th node
	double weight(int iOrientation) const; //!< get weight for the iOrientation'th node

private:
	int nS1byZn; //!< actual number of S1 samples (reduced by symmetry)
	int nS1; //!< effective number of S1 samples (counting symmetric images)
	std::vector<vector3<> > eulerS2; //!< S2 quadrature points
	std::vector<double> weightS2; //!< S2 quadrature weights
	void setup(const S2quad&, int Zn); //!< Initialize SO3/Zn quadrature from an S2 quadrature decsription
};

//! @}
#endif // JDFTX_FLUID_SO3QUAD_H
