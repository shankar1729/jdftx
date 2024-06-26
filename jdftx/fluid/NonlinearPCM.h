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


#ifndef JDFTX_ELECTRONIC_NONLINEARPCM_H
#define JDFTX_ELECTRONIC_NONLINEARPCM_H

#include <fluid/PCM.h>
#include <core/Minimize.h>

namespace NonlinearPCMeval { struct Screening; struct Dielectric; } //Forward declaration of helper classes

//! @addtogroup Solvation
//! @{
//! @file NonlinearPCM.h NonlinearPCM and helper classes

//! Nonlinear solvation models: shared electrostatic part implementation
class NonlinearPCM : public PCM, public Minimizable<ScalarFieldTilde>
{
public:
	ScalarFieldTilde phiTot; //!< State of the solver = total electrostatic potential

	//! See createFluidSolver()
	NonlinearPCM(const Everything& e, const FluidSolverParams& params);
    virtual ~NonlinearPCM();
	
	bool prefersGummel() const { return false; }

	void loadState(const char* filename); //!< Load state from file
	void saveState(const char* filename) const; //!< Save state to file
	void dumpDensities(const char* filenamePattern) const;
	void minimizeFluid(); //!< Converge using nonlinear conjugate gradients

	// Interface for Minimizable:
	void step(const ScalarFieldTilde& dir, double alpha);
	double compute(ScalarFieldTilde* grad, ScalarFieldTilde* Kgrad);
	bool report(int iter);

protected:
	void set_internal(const ScalarFieldTilde& rhoExplicitTilde, const ScalarFieldTilde& nCavityTilde);
	double get_Adiel_and_grad_internal(ScalarFieldTilde& Adiel_rhoExplicitTilde, ScalarFieldTilde& Adiel_nCavityTilde, IonicGradient* extraForces, matrix3<>* Adiel_RRT) const;

private:
	int iterLast; //latest iteration number in fluid minimize (used to report iteration count when inner log hidden)
	double A0; //constant energy term added during potential optimization
	double pMol, ionNbulk, ionZ;
	NonlinearPCMeval::Screening* screeningEval; //!< Internal helper class for Screening from PCM_internal
	NonlinearPCMeval::Dielectric* dielectricEval; //!< Internal helper class for Dielectric from PCM_internal
	RadialFunctionG dielEnergyLookup, ionEnergyLookup; //!< lookup tables for energy during nonlinear Poisson-Boltzmann solve
	std::shared_ptr<RealKernel> preconditioner;
	
	//Extra quantities for CANON alone:
	bool isNonlocal; //!< whether nonlocal extensions to NonlinearPCM (CANON) are active
	ScalarFieldTilde nCavityNetTilde; //!< input nCavity + full core before convolution
	ScalarFieldTilde rhoLiquidTilde0; //!< built-in charge density in liquid
	RadialFunctionG Nw0, w1; //!< Solvent l=0 and l=1 weight functions (Nw0 includes Nbulk)
};

//! @}
#endif // JDFTX_ELECTRONIC_NONLINEARPCM_H
