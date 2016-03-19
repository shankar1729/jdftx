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

#ifndef JDFTX_ELECTRONIC_IONDYNAMICSPARAMS_H
#define JDFTX_ELECTRONIC_IONDYNAMICSPARAMS_H

#include <core/Units.h>

enum ConfiningPotentialType
{	ConfineLinear,
	ConfineQuadratic,
	ConfineCubic,
	ConfineSmoothLinear,
	ConfineNone
};

//! @addtogroup md
//! @{

//! @brief Parameters to control the Verlet algorithm
struct IonDynamicsParams
{	double dt, tMax, kT, alpha;
	ConfiningPotentialType confineType;
	std::vector<double> confineParameters;
	//! Set the default values
	IonDynamicsParams(): dt(1.0*fs), tMax(0.0) ,kT(0.001), alpha(0.0), confineType(ConfineNone){}
};

//! @}

#endif // JDFTX_ELECTRONIC_IONDYNAMICSPARAMS_H
