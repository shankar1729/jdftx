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

#ifndef JDFTX_ELECTRONIC_LINEARPCM_H
#define JDFTX_ELECTRONIC_LINEARPCM_H

#include <fluid/PCM.h>
#include <core/Minimize.h>

class LinearPCM : public PCM, public LinearSolvable<DataGptr>
{
public:
	LinearPCM(const Everything& e, const FluidSolverParams& fsp); //!< Parameters same as createFluidSolver()
    virtual ~LinearPCM();
	bool needsGummel() { return false; }

	DataGptr hessian(const DataGptr&) const; //!< Implements #LinearSolvable::hessian for the dielectric poisson equation
	DataGptr precondition(const DataGptr&) const; //!< Implements a modified inverse kinetic preconditioner

	void minimizeFluid(); //!< Converge using linear conjugate gradients
	void loadState(const char* filename); //!< Load state from file
	void saveState(const char* filename) const; //!< Save state to file

protected:
	void set_internal(const DataGptr& rhoExplicitTilde, const DataGptr& nCavityTilde);
	double get_Adiel_and_grad_internal(DataGptr& grad_rhoExplicitTilde, DataGptr& grad_nCavityTilde, IonicGradient& extraForces) const;
	friend class NonlinearPCM;
private:
	RadialFunctionG Kkernel; DataRptr epsInv; // for preconditioner
};

#endif // JDFTX_ELECTRONIC_LINEARPCM_H
