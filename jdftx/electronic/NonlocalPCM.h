/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_NONLOCALPCM_H
#define JDFTX_ELECTRONIC_NONLOCALPCM_H

#include <electronic/FluidSolver.h>
#include <core/Minimize.h>

class NonlocalPCM : public FluidSolver, public LinearSolvable<DataGptr>
{
public:
	NonlocalPCM(const Everything& e, const FluidSolverParams& fsp); //!< Parameters same as createFluidSolver()
	bool needsGummel() { return false; }

	DataGptr chi(const DataGptr&) const; //!< Apply the non-local chi (i.e. compute induced charge density given a potential)
	DataGptr hessian(const DataGptr&); //!< Implements #LinearSolvable::hessian for the non-local poisson-like equation
	DataGptr precondition(const DataGptr&); //!< Implements a modified inverse kinetic preconditioner

	//! Set the explicit system charge density and effective cavity-formation electron density:
	void set(const DataGptr& rhoExplicitTilde, const DataGptr& nCavityTilde);

	void minimizeFluid(); //!< Converge using linear conjugate gradients

	//! Get the minimized free energy and the electronic n-gradient
	double get_Adiel_and_grad(DataGptr& grad_rhoExplicitTilde, DataGptr& grad_nCavityTilde, IonicGradient& extraForces) const;

	void loadState(const char* filename); //!< Load state from file
	void saveState(const char* filename) const; //!< Save state to file

	void dumpDensities(const char* filenamePattern) const;
	
private:
	const FluidSolverParams& params;
	DataRptr nProduct; //convolution of densities that determines cavity
	DataRptr shape; //cavity shape function (0 to 1)
	DataGptr rhoExplicitTilde; //charge density of explicit system
	std::vector< std::shared_ptr<struct MultipoleResponse> > response; //array of multipolar components in chi
	RealKernel nFluid; //electron density model for the fluid
	RealKernel Kkernel; DataRptr epsInv; double epsBulk; //for preconditioner
};

#endif // JDFTX_ELECTRONIC_NONLOCALPCM_H
