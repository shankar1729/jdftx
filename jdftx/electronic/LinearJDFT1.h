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

#ifndef JDFTX_ELECTRONIC_LINEARJDFT1_H
#define JDFTX_ELECTRONIC_LINEARJDFT1_H

#include <electronic/FluidJDFTx.h>
#include <core/Minimize.h>


//! Compute the shape function (0 to 1) given the cavity-determining electron density
void JDFT1_shapeFunc(const DataRptr& nCavity, DataRptr& shape, const double nc, const double sigma);

//!Compute derivative with respect to cavity-determining electron density, given derivative with respect to shape function
void JDFT1_shapeFunc_grad(const DataRptr& nCavity, const DataRptr& grad_shape, DataRptr& grad_nCavity, const double nc, const double sigma);



//! Some extra quantities specific to LinearJDFT1 added here
struct LinearJDFT1params : public FluidSolverParams
{
	double k2factor; //!< Prefactor to kappaSq (computed from T, ionicConcentration and ionicZelectrolyte)

	//! A copy constructor to set base-class variables
	LinearJDFT1params(const FluidSolverParams& p) : FluidSolverParams(p) {}
};


class LinearJDFT1 : public FluidSolver, public LinearSolvable<DataGptr>
{
public:
	LinearJDFT1(const Everything& e, const FluidSolverParams& fsp); //!< Parameters same as createFluidSolver()
	bool needsGummel() { return false; }

	DataGptr hessian(const DataGptr&); //!< Implements #LinearSolvable::hessian for the dielectric poisson equation
	DataGptr precondition(const DataGptr&); //!< Implements a modified inverse kinetic preconditioner

	//! Set the explicit system charge density and effective cavity-formation electron density:
	void set(const DataGptr& rhoExplicitTilde, const DataGptr& nCavityTilde);

	void minimizeFluid(); //!< Converge using linear conjugate gradients

	//! Get the minimized free energy and the electronic n-gradient
	double get_Adiel_and_grad(DataGptr& grad_rhoExplicitTilde, DataGptr& grad_nCavityTilde);

	void loadState(const char* filename); //!< Load state from file
	void saveState(const char* filename) const; //!< Save state to file

	void dumpDensities(const char* filenamePattern) const;
	
private:
	DataRptr nCavity, epsilon, kappaSq;
	DataGptr rhoExplicitTilde;
	LinearJDFT1params params;
	RealKernel Kkernel; DataRptr epsInv; // for preconditioner
};

#endif // JDFTX_ELECTRONIC_LINEARJDFT1_H
