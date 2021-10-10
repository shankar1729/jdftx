/*-------------------------------------------------------------------
Copyright 2021 Ravishankar Sundararaman.

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

#ifndef JDFTX_ELECTRONIC_VANDERWAALSD3_H
#define JDFTX_ELECTRONIC_VANDERWAALSD3_H

#include <electronic/VanDerWaals.h>
#include <core/RadialFunction.h>

//! @addtogroup LongRange
//! @{

//! DFT-D3 pair potential dispersion correction \cite Dispersion-D3
class VanDerWaalsD3 : public VanDerWaals
{
public:
	VanDerWaalsD3(const Everything &e);
	~VanDerWaalsD3();
	
	//Implement virtual functions of VanDerWaals abstract base class
	virtual double getScaleFactor(string exCorrName, double scaleOverride=0.) const;
	virtual double energyAndGrad(std::vector<Atom>& atoms, const double scaleFac, matrix3<>* E_RRT=0) const;
};

//! @}
#endif // JDFTX_ELECTRONIC_VANDERWAALSD3_H
