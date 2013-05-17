/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver, Deniz Gunceler

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

#ifndef JDFTX_ELECTRONIC_FLUIDSOLVER_H
#define JDFTX_ELECTRONIC_FLUIDSOLVER_H

//! @file FluidSolver.h
//! Common interface for all the fluids to the electronic code

#include <core/Data.h>
#include <electronic/FluidSolverParams.h>
#include <electronic/IonicMinimizer.h>

//! Abstract base class for the fluid solvers
struct FluidSolver
{
	const Everything& e;
	const FluidSolverParams& fsp;
	double epsBulk, epsInf; //!< bulk dielectric constants of fluid
	double k2factor; //!< prefactor to screening term (0 => no ionic screening)
	
	//! Abstract base class constructor - do not use directly - see FluidSolver::createSolver
	FluidSolver(const Everything &e, const FluidSolverParams& fsp);
	virtual ~FluidSolver() {}

	//! Set total explicit charge density and effective electron density to use in cavity formation (i.e. including charge balls)
	//!  and set list of explicit atoms to use in van der Waals corrections
	virtual void set(const DataGptr& rhoExplicitTilde, const DataGptr& nCavityTilde)=0;

	//! Compute gradients with respect to electronic side variables, and return fluid+coupling free energy
	//! Any extra forces on explicit ions due to the fluid should be stored in extraForces
	virtual double get_Adiel_and_grad(DataGptr& Adiel_rhoExplicitTilde, DataGptr& Adiel_nCavityTilde, IonicGradient& extraForces) const =0;

	//! Dump relevant fluid densities (eg. NO and NH) to file(s)
	//! the provided pattern will have a single %s which may be substituted
	//! Fluid solver implementations may override to dump fluid densities, no dumping by default
	virtual void dumpDensities(const char* filenamePattern) const {};

	//! Dump fluid debugging quantities (eg. fluid potential and effective electron density for 3.0)
	//! the provided pattern will have a single %s which may be substituted
	//! Fluid solver implementations may override to dump fluid debug stuff, no dumping by default
	virtual void dumpDebug(const char* filenamePattern) const {};


	//------------Fluid solver implementations must provide these pure virtual functions

	//! Specify whether fluid requires a gummel loop (true) or is minimized each time (false)
	virtual bool needsGummel()=0;

	//! Initialize fluid state from a file
	virtual void loadState(const char* filename)=0;

	//! Save fluid state to a file
	virtual void saveState(const char* filename) const=0;

	//! Minimize fluid side (holding explicit electronic system fixed)
	virtual void minimizeFluid()=0;
};

//! Create and return a JDFTx solver (the solver can be freed using delete)
FluidSolver* createFluidSolver(const Everything& e, const FluidSolverParams& params);

#endif // JDFTX_ELECTRONIC_FLUIDSOLVER_H
