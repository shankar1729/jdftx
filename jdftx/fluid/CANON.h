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

#ifndef JDFTX_FLUID_CANON_H
#define JDFTX_FLUID_CANON_H

#include <fluid/PCM.h>
#include <core/Pulay.h>

//! @addtogroup Solvation
//! @{

//! Nonlocal CANON \cite CANON solvation model implementation (electrostatic part)
class CANON : public PCM, public Pulay<ScalarFieldTilde>
{
public:
	CANON(const Everything& e, const FluidSolverParams& fsp); //!< Parameters same as createFluidSolver()
    virtual ~CANON();
	bool prefersGummel() const { return false; }

	void minimizeFluid(); //!< Converge using Pulay on electrostatic potential
	void loadState(const char* filename); //!< Load state from file
	void saveState(const char* filename) const; //!< Save state to file
	void dumpDensities(const char* filenamePattern) const; //!< dump cavity shape functions

protected:
	//Interface for FluidSolver:
	void set_internal(const ScalarFieldTilde& rhoExplicitTilde, const ScalarFieldTilde& nCavityTilde);
	double get_Adiel_and_grad_internal(ScalarFieldTilde& grad_rhoExplicitTilde, ScalarFieldTilde& grad_nCavityTilde, IonicGradient* extraForces, matrix3<>* Adiel_RRT) const;
	void getSusceptibility_internal(const std::vector<complex>& omega, std::vector<SusceptibilityTerm>& susceptibility, ScalarFieldArray& sArr, bool elecOnly) const;

	//Interface for Pulay<ScalarFieldTilde>:
	double cycle(double dEprev, std::vector<double>& extraValues);
	void axpy(double alpha, const ScalarFieldTilde& X, ScalarFieldTilde& Y) const { ::axpy(alpha, X, Y); }
	double dot(const ScalarFieldTilde& X, const ScalarFieldTilde& Y) const { return ::dot(X, Y); }
	size_t variableSize() const { return gInfo.nG * sizeof(complex); }
	void readVariable(ScalarFieldTilde& X, FILE* fp) const;
	void writeVariable(const ScalarFieldTilde& X, FILE* fp) const;
	ScalarFieldTilde getVariable() const;
	void setVariable(const ScalarFieldTilde&);
	ScalarFieldTilde precondition(const ScalarFieldTilde&) const;
	ScalarFieldTilde applyMetric(const ScalarFieldTilde&) const;

private:
	ScalarFieldTilde phiTot; //!< total electrostatic potential (internal state of fluid)
	VectorField eps; //!< Internal effective electric field
	ScalarField muPlus, muMinus; //!< Effective chemical potentials of cations and anions
	ScalarFieldTilde nCavityNetTilde; //!< input nCavity + full core before convolution
	RadialFunctionG w0, w1; //!< Solvent l=0 and l=1 weight functions
	std::shared_ptr<RealKernel> preconditioner, metric; //!< Pulay mixing and metric kernels
};

//! @}
#endif // JDFTX_FLUID_CANON_H
