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

#ifndef JDFTX_ELECTRONIC_IONICDYNAMICSPARAMS_H
#define JDFTX_ELECTRONIC_IONICDYNAMICSPARAMS_H

#include <core/Units.h>

//! @addtogroup IonicSystem
//! @{
//! @file IonicDynamicsParams.h Struct IonicDynamicsParams

//! Parameters to control IonicDynamics
struct IonicDynamicsParams
{	double dt; //!< time step [Eh^-1]
	int nSteps; //!< number of steps
	enum StatMethod { StatNone, Berendsen, NoseHoover } statMethod; //!< Method for thermo- and/or baro-stat
	double T0; //!< initial temperature or set point temperature if StatMethod != StatNone [Eh]
	double P0; //!< pressure set point [Eh/a0^3] (NAN if not barostatting (hydrostatic))
	matrix3<> stress0; //!< stress set point [Eh/a0^3] (NAN if not barostatting (anisotropic))
	double tDampT; //!< thermostat damping time [Eh^-1]
	double tDampP; //!< barostat damping time [Eh^-1]
	int chainLengthT; //!< Nose-Hoover chain length for thermostat
	int chainLengthP; //!< Nose-Hoover chain length for barostat
	double B0; //!< characteristic bulk modulus for Berendsen barostat (default: water bulk modulus)
	
	IonicDynamicsParams() : dt(1.*fs), nSteps(0), statMethod(StatNone),
		T0(298*Kelvin), P0(NAN), stress0(NAN,NAN,NAN),
		tDampT(50.*fs), tDampP(100.*fs),
		chainLengthT(3), chainLengthP(3), B0(2.2E9*Pascal) {}
};

//! @}
#endif // JDFTX_ELECTRONIC_IONICDYNAMICSPARAMS_H
