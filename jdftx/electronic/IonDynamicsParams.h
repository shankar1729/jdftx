/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth-Weaver

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

//! @addtogroup IonicSystem
//! @{
//! @file IonDynamicsParams.h IonDynamicsParams and related definitions

//! Type of confining potential
enum ConfiningPotentialType
{	ConfineLinear,
	ConfineQuadratic,
	ConfineCubic,
	ConfineSmoothLinear,
	ConfineNone
};

//! Strategy for eliminating drift 
enum DriftRemovalType
{      DriftNone, //!< No drift correction
       DriftVelocity, //!< Zero total velocity
       DriftMomentum //!< Zero total momentum
};

//! Parameters to control IonicDynamics
struct IonDynamicsParams
{	double dt; //!< time step
	double tMax; //!< maximum time
	double kT; //!< temperature
	double alpha; //!< velocity scaling parameter
	DriftRemovalType driftType; //!< drift removal strategy
	ConfiningPotentialType confineType; //!< confinement potential type
	std::vector<double> confineParameters; //!< parameters controlling confinement potential
	
	//! Set the default values
	IonDynamicsParams(): dt(1.0*fs), tMax(0.0) ,kT(0.001), alpha(0.0), driftType(DriftMomentum), confineType(ConfineNone){}
};

//! @}
#endif // JDFTX_ELECTRONIC_IONDYNAMICSPARAMS_H
