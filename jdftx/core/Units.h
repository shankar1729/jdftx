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

#ifndef JDFTX_CORE_UNITS_H
#define JDFTX_CORE_UNITS_H

//! @addtogroup Utilities
//! @{

//! @file Units.h Commonly used measurement units in terms of atomic units

#include <cmath>

//Energy, temperature units in Hartrees:
const double eV = 1/27.21138505; //!< eV in Hartrees
const double Joule = 1/4.35974434e-18; //!< Joule in Hartrees
const double KJoule = 1000*Joule; //!< KJoule in Hartrees
const double Kcal = KJoule * 4.184; //!< Kcal in Hartrees
const double Kelvin = 1.3806488e-23*Joule; //!< Kelvin in Hartrees
const double invcm = 1./219474.6313705; //!< Inverse cm in Hartrees

//Length units in bohrs:
const double Angstrom = 1/0.5291772; //!< Angstrom in bohrs
const double meter = 1e10*Angstrom; //!< meter in bohrs
const double liter = 1e-3*pow(meter,3);  //!< liter in cubic bohrs

//Mass units in electron masses:
const double amu = 1822.88839; //!< atomic mass unit in electron masses
const double kg = 1./9.10938291e-31; //!< kilogram in electron masses

//Dimensionless:
const double mol = 6.0221367e23; //!< mole in number (i.e. Avogadro number)

//Commonly used derived units:
const double Newton = Joule/meter;  //!< Newton in Hartree/bohr
const double Pascal = Newton/(meter*meter); //!< Pascal in Hartree/bohr^3
const double KPascal = 1000*Pascal;  //!< KPa in Hartree/bohr^3
const double Bar = 100*KPascal;   //!< bar in Hartree/bohr^3
const double mmHg = 133.322387415*Pascal;  //!< mm Hg in Hartree/bohr^3

//Time
const double sec = sqrt((kg*meter)/Newton); //!< second in inverse Hartrees
const double invSec = 1./sec; //!< inverse second in Hartrees
const double fs = sec*1.0e-15; //!< femtosecond in inverse Hartrees

//Electrical:
const double Coul = Joule/eV; //!< Coulomb in electrons
const double Volt = Joule/Coul; //!< Volt in Hartrees
const double Ampere = Coul/sec; //!< Ampere in electrons/inverse Hartree
const double Ohm = Volt/Ampere; //!< Ohm in inverse conductance quanta

//! @}
#endif //JDFTX_CORE_UNITS_H
