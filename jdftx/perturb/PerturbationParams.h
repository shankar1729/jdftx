/*-------------------------------------------------------------------
Copyright 2023 Brandon Li

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

#ifndef PERTURB_PERTURBATIONPARAMS_H_
#define PERTURB_PERTURBATIONPARAMS_H_

#include <cstdio>

//! @addtogroup Algorithms
//! @{

//! Parameters to control the minimization algorithm
struct PerturbationParams
{
	int nIterations;
	FILE* fpLog; //!< Stream to log iterations to
	const char* linePrefix; //!< prefix for each output line of minimizer, useful for nested minimizations (default "CG\t")
	
	enum Algorithms
	{	CG,
		MINRES
	} algorithm;
	
	double residualTol;
	double residualDiffThreshold;
	bool CGBypass;
	bool recomputeResidual;
	
	//! Set the default values
	PerturbationParams() 
	: nIterations(0), fpLog(stdout), linePrefix("Linear solver:\t"), algorithm(MINRES), residualTol(1e-4), residualDiffThreshold(1e-4), CGBypass(false), recomputeResidual(false) {}
};

//! @}

#endif /* PERTURB_PERTURBATIONPARAMS_H_ */
