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

#ifndef JDFTX_ELECTRONIC_SALSA_H
#define JDFTX_ELECTRONIC_SALSA_H

#include <fluid/PCM.h>
#include <core/Minimize.h>

class SaLSA : public PCM, public LinearSolvable<ScalarFieldTilde>
{
public:
	SaLSA(const Everything& e, const FluidSolverParams& fsp); //!< Parameters same as createFluidSolver()
    virtual ~SaLSA();
	bool needsGummel() { return false; }

	ScalarFieldTilde chi(const ScalarFieldTilde&) const; //!< Apply the non-local chi (i.e. compute induced charge density given a potential)
	ScalarFieldTilde hessian(const ScalarFieldTilde&) const; //!< Implements #LinearSolvable::hessian for the non-local poisson-like equation
	ScalarFieldTilde precondition(const ScalarFieldTilde&) const; //!< Implements a modified inverse kinetic preconditioner
	double sync(double x) const; //!< All processes minimize together; make sure scalars are in sync to round-off error

	void minimizeFluid(); //!< Converge using linear conjugate gradients
	void loadState(const char* filename); //!< Load state from file
	void saveState(const char* filename) const; //!< Save state to file
	void dumpDensities(const char* filenamePattern) const; //!< dump cavity shape functions

protected:
	void set_internal(const ScalarFieldTilde& rhoExplicitTilde, const ScalarFieldTilde& nCavityTilde);
	double get_Adiel_and_grad_internal(ScalarFieldTilde& grad_rhoExplicitTilde, ScalarFieldTilde& grad_nCavityTilde, IonicGradient* extraForces, bool electricOnly) const;

private:
	std::vector< std::shared_ptr<struct MultipoleResponse> > response; //array of multipolar components in chi
	int rStart, rStop; //MPI division of response array
	RadialFunctionG nFluid; //electron density model for the fluid
	RadialFunctionG Kkernel; ScalarField epsInv; //for preconditioner
	ScalarFieldArray siteShape; //shape functions for sites
};

#endif // JDFTX_ELECTRONIC_SALSA_H
