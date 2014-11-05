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

#ifndef JDFTX_CORE_RANDOM_H
#define JDFTX_CORE_RANDOM_H

#include <core/scalar.h>

//! @file Random.h Random number generation

namespace Random
{
	void seed(int i); //seed random number generator
	double uniform(double start=0.0, double end=1.0); //!< uniform random numbers between 0 and 1
	double normal(double mean=0.0, double sigma=1.0, double cap=0.0); //!< normal random numbers with mean, sigma and an optional cap if non-zero
	complex normalComplex(double sigma=1.0); //!< normal complex number with mean 0 and deviation sigma
}

#endif //JDFTX_CORE_RANDOM_H
