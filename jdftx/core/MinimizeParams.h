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

#ifndef JDFTX_CORE_MINIMIZEPARAMS_H
#define JDFTX_CORE_MINIMIZEPARAMS_H

#include <cstdio>

//! @addtogroup optimization
//! @{

//! @brief Parameters to control the Nonlinear Conjugate Gradients algorithm
struct MinimizeParams
{
	//! Search direction update scheme
	enum DirectionUpdateScheme
	{	PolakRibiere, //!< Polak-Ribiere (preconditioned) conjugate gradients (default)
		FletcherReeves, //!< Fletcher-Reeves (preconditioned) conjugate gradients
		HestenesStiefel, //!< Hestenes-Stiefel (preconditioned) conjugate gradients
		SteepestDescent //!< Steepest Descent (always along negative (preconditioned) gradient)
	} dirUpdateScheme;

	//! Line minimization method
	enum LinminMethod
	{	Relax, //!< move by a constant multiple (=alphaTstart) of the search direction (not recommended for CG)
		Quad, //!< use the energy at a test step location to find the minimum along the line (default)
		Cubic //!< use the energy and gradient at a test step location to find the minimum along the line
	} linminMethod;

	int nIterations; //!< Maximum number of iterations (default 100)
	int nDim; //!< Dimension of optimization space; used only for knormThreshold (default 1)
	FILE* fpLog; //!< Stream to logPrintf iterations to
	const char* linePrefix; //!< prefix for each output line of minimizer, useful for nested minimizations (default "CG\t")
	const char* energyLabel; //!< Label for the minimized quantity (default "E")
	double knormThreshold; //!< stop when norm of residual against preconditioner falls below this (default: 0)
	double energyDiffThreshold; //!< stop when energy change is below this for nEnergyDiff successive iterations (default: 0)
	int nEnergyDiff; //!< number of successive iterations for energyDiffThreshold check (default: 2)
	
	int nDirResetNum; //!< Switch to stepeest descents if search direction is reset nDirResetNum out of the last nDirResetDen iterations (default: 2)
	int nDirResetDen; //!< Switch to stepeest descents if search direction is reset nDirResetNum out of the last nDirResetDen iterations (default: 5)
	
	double alphaTstart; //!< initial value for the test-step size (default: 1.0)
	double alphaTmin; //!< minimum value of the test-step size (algorithm gives up when difficulties cause alphaT to fall below this value) (default:1e-10)
	bool updateTestStepSize; //!< set alphaT=alpha after every iteration if true (default: true)

	double alphaTreduceFactor; //!< Factor to reduce alphaT on difficulties (default 0.1)
	double alphaTincreaseFactor; //!< Max ratio of alpha to alphaT, increase alphaT by this factor otherwise (default 3.0)
	int nAlphaAdjustMax; //!< maximum number of times to alpha adjust attempts (default 3)

	bool fdTest; //!< whether to perform a finite difference test before each minimization (default false)
	
	//! Set the default values
	MinimizeParams() 
	: dirUpdateScheme(PolakRibiere), linminMethod(Quad),
		nIterations(100), nDim(1), fpLog(stdout), linePrefix("CG\t"), energyLabel("E"),
		knormThreshold(0), energyDiffThreshold(0), nEnergyDiff(2), nDirResetNum(2), nDirResetDen(5),
		alphaTstart(1.0), alphaTmin(1e-10), updateTestStepSize(true),
		alphaTreduceFactor(0.1), alphaTincreaseFactor(3.0), nAlphaAdjustMax(3),
		fdTest(false) {}
};

//! @}

#endif // JDFTX_CORE_MINIMIZEPARAMS_H
