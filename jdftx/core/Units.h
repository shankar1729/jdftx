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

//! @file Units.h
//! @brief Commonly used measurement units in terms of atomic units

#include <cmath>

//Energy, temperature units in Hartrees:
const double eV = 1/27.2116; //!< @f$ eV / E_h @f$
const double Joule = 1/4.3597482e-18; //!< @f$ J / E_h @f$
const double KJoule = 1000*Joule; //!< @f$ KJ / E_h @f$
const double Kcal = KJoule * 4.184; //!< @f$ Kcal / E_h @f$
const double Kelvin = 1.380658e-23*Joule;  //!< @f$ k_B K / E_h @f$

//Length units in bohrs:
const double Angstrom = 1/0.5291772; //!< @f$ \AA / a_0 @f$
const double meter = 1e10*Angstrom; //!< @f$ m / a_0 @f$
const double liter = 1e-3*pow(meter,3);  //!< @f$ l / a_0^3 @f$

//Dimensionless:
const double mol = 6.0221367e23; //!< @f$ mol/\# @f$ (= Avogadro number)

//Commonly used derived units:
const double Newton = Joule/meter;  //!< @f$ N / (E_h/a_0) @f$
const double Pascal = Newton/(meter*meter); //!< @f$ Pa / (E_h/a_0^2) @f$
const double KPascal = 1000*Pascal;  //!< @f$ KPa / (E_h/a_0^2) @f$
const double Bar = 100*KPascal;   //!< @f$ bar / (E_h/a_0^2) @f$
const double mmHg = 133.322387415*Pascal;  //!< @f$ mmHg / (E_h/a_0^2) @f$

#endif //JDFTX_CORE_UNITS_H
