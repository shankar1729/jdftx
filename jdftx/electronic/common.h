/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman
Copyright 1996-2003 Sohrab Ismail-Beigi

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

#ifndef JDFTX_ELECTRONIC_COMMON_H
#define JDFTX_ELECTRONIC_COMMON_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdarg>
#include <core/scalar.h>
#include <core/Util.h>

class matrix;			// general purpose matrix.
class diagMatrix; 		// real diagonal matrix
class ColumnBundle;	// stores wavefunctions
class QuantumNumber;	// a set of quantum numbers, eg k, spin etc
class Control;			// contains control information for the run.
class GridInfo;			// Lattice information (vectors, reciprocal, etc.)
class Basis;			// G vector and realspace grid indices.
class SpeciesInfo;		// information for an ion species
class IonInfo;			// collections of all SpeciesInfos
class Symmetries;		// symmetry info
class ExactExchange;	// Evaluates exact exchange (fock term)
class ExCorr;			// Ionic dynamics information.
class ElecInfo;			// parameters for electronic states
class ElecVars;			// collection of electronic variables.
class Energies;			// collection of all energies.
class Everything;		// A big collection of most of the above structures

struct ElecGradient;
struct IonicGradient;
class InverseKohnSham;
class VanDerWaals;     // pair-potential van der waals corrections
class Vibrations;

#endif // JDFTX_ELECTRONIC_COMMON_H
