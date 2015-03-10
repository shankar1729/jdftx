/*-------------------------------------------------------------------
Copyright 2015 Ravishankar Sundararaman

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

#ifndef JDFTX_CORE_PULAYPARAMS_H
#define JDFTX_CORE_PULAYPARAMS_H

#include <cstdio>

//! Parameters to control Pulay mixing
struct PulayParams
{
	FILE* fpLog; //!< Stream to log iterations to
	const char* linePrefix; //!< prefix for each output line of Pulay (default "Pulay: ")
	const char* energyLabel; //!< Label for the minimized quantity (default "E")
	const char* energyFormat; //!< printf format for the minimized quantity (default "%22.15le")

	int nIterations; //!< maximum iterations (single point calculation if 0)
	double energyDiffThreshold; //!< convergence threshold for energy difference between successive iterations
	double residualThreshold; //!< convergence threshold on the residual

	int history; //!< Number of past residuals and vectors that are cached and used for mixing
	double mixFraction;  //!< Mixing fraction for total density / potential
	double qKerker; //!< Wavevector controlling Kerker preconditioning (if negative, auto-set to Gmin)
	double qMetric; //!< Wavevector controlling the metric for overlaps (if negative, auto-set to Gmin)
	
	PulayParams()
	: fpLog(stdout), linePrefix("Pulay: "), energyLabel("E"), energyFormat("%22.15le"),
		nIterations(50), energyDiffThreshold(1e-8), residualThreshold(1e-7),
		history(10), mixFraction(0.5), qKerker(0.8), qMetric(0.8)
	{
	}
};

#endif //JDFTX_CORE_PULAYPARAMS_H