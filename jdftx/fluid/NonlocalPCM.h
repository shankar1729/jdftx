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

#include <fluid/PCM.h>
#include <core/Minimize.h>

class NonlocalPCM : public PCM, public LinearSolvable<DataGptr>
{
public:
	NonlocalPCM(const Everything& e, const FluidSolverParams& fsp); //!< Parameters same as createFluidSolver()
    virtual ~NonlocalPCM();
	bool needsGummel() { return false; }

	DataGptr hessian(const DataGptr&) const; //!< Implements #LinearSolvable::hessian for the non-local poisson-like equation
	DataGptr precondition(const DataGptr&) const; //!< Implements a modified inverse kinetic preconditioner
	double sync(double x) const { mpiUtil->bcast(x); return x; } //!< All processes minimize together; make sure scalars are in sync to round-off error
	
	//! Set the explicit system charge density and effective cavity-formation electron density:
	void set(const DataGptr& rhoExplicitTilde, const DataGptr& nCavityTilde);

	void minimizeFluid(); //!< Converge using linear conjugate gradients

	//! Get the minimized free energy and the electronic n-gradient
	double get_Adiel_and_grad(DataGptr& grad_rhoExplicitTilde, DataGptr& grad_nCavityTilde, IonicGradient& extraForces) const;

	void loadState(const char* filename); //!< Load state from file
	void saveState(const char* filename) const; //!< Save state to file
protected:
	void printDebug(FILE* fp) const;
private:
	double sigmaVdw; //!< gaussian width for weight functions below
	RadialFunctionG& wCavity; //unit-norm weight function (~ electron density / Ztot), points to Sf[0] in PCM
	RadialFunctionG wDiel; //unit-norm weight function for dielectric response
	RadialFunctionG Kkernel; DataRptr epsInv; //for preconditioner
	DataRptr shapeDiel; double A_eta;
};

#endif // JDFTX_ELECTRONIC_NONLOCALPCM_H
