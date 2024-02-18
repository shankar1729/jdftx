/*-------------------------------------------------------------------
Copyright 2024 Ravishankar Sundararaman

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

#ifndef JDFTX_FLUID_NONLINEAR_COMMON_H
#define JDFTX_FLUID_NONLINEAR_COMMON_H

#include <fluid/FluidSolverParams.h>

namespace NonlinearPCMeval { struct Screening; struct Dielectric; } //Forward declaration of helper classes

//! @addtogroup Solvation
//! @{

//! Shared implementations between nonlinear models (NonlinearPCM and CANON)
class NonlinearCommon
{
protected:
	NonlinearCommon(const FluidSolverParams& fsp, double epsBulk);
	virtual ~NonlinearCommon();

	double pMol, ionNbulk, ionZ;
	NonlinearPCMeval::Screening* screeningEval; //!< Internal helper class for Screening from PCM_internal
	NonlinearPCMeval::Dielectric* dielectricEval; //!< Internal helper class for Dielectric from PCM_internal
	RadialFunctionG gLookup, xLookup; //!< lookup tables for transcendental solutions involved in the dielectric and ionic SCF method
	RadialFunctionG dielEnergyLookup, ionEnergyLookup; //!< lookup tables for energy during nonlinear Poisson-Boltzmann solve
};

//! @}
#endif // JDFTX_FLUID_NONLINEAR_COMMON_H
