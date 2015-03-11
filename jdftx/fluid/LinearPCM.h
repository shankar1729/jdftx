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

class LinearPCM : public PCM, public LinearSolvable<ScalarFieldTilde>
{
public:
	LinearPCM(const Everything& e, const FluidSolverParams& fsp); //!< Parameters same as createFluidSolver()
    virtual ~LinearPCM();
	bool needsGummel() { return false; }

	ScalarFieldTilde hessian(const ScalarFieldTilde&) const; //!< Implements #LinearSolvable::hessian for the dielectric poisson equation
	ScalarFieldTilde precondition(const ScalarFieldTilde&) const; //!< Implements a modified inverse kinetic preconditioner

	void minimizeFluid(); //!< Converge using linear conjugate gradients
	void loadState(const char* filename); //!< Load state from file
	void saveState(const char* filename) const; //!< Save state to file

protected:
	void set_internal(const ScalarFieldTilde& rhoExplicitTilde, const ScalarFieldTilde& nCavityTilde);
	double get_Adiel_and_grad_internal(ScalarFieldTilde& grad_rhoExplicitTilde, ScalarFieldTilde& grad_nCavityTilde, IonicGradient* extraForces, bool electricOnly) const;
private:
	RadialFunctionG Kkernel; ScalarField epsInv; // for preconditioner
	void updatePreconditioner(const ScalarField& epsilon, const ScalarField& kappaSq);
	
	//Optionally override epsilon and kappaSq (when used as the inner solver in NonlinearPCM's SCF):
	friend class NonlinearPCM;
	ScalarField epsilonOverride, kappaSqOverride;
	void override(const ScalarField& epsilon, const ScalarField& kappaSq);
};

#endif // JDFTX_ELECTRONIC_LINEARPCM_H
