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

//! @addtogroup Fluid
//! @{
//! @file FluidSolver.h Common interface for all the fluids to the electronic code

#include <core/ScalarField.h>
#include <fluid/FluidSolverParams.h>
#include <electronic/IonicMinimizer.h>

//! Abstract base class for the fluid solvers
struct FluidSolver
{
	const Everything& e;
	const GridInfo& gInfo; //relevant gInfo for fluid (uses embedded grid when coulomb truncation is enabled)
	const FluidSolverParams& fsp;
	double epsBulk, epsInf; //!< bulk dielectric constants of fluid
	double k2factor; //!< prefactor to screening term (0 => no ionic screening)
	std::vector<std::vector< vector3<> > > atpos; //!atomic positions per species in the relevant coordinate system (depending on embedding option)
	ScalarFieldTilde A_rhoNonES; //!Any non-electrostatic contributions to A_rhoExplicitTilde (removed from dumped d_fluid / d_tot)
	
	//! Abstract base class constructor - do not use directly - see FluidSolver::createSolver
	FluidSolver(const Everything &e, const FluidSolverParams& fsp);
	virtual ~FluidSolver() {}

	double ionWidthMuCorrection() const; //!< correction to electron chemical potential due to finite ion width in fluid interaction

	//Whether a gummel loop is in use
	inline bool useGummel() const
	{	switch(fsp.solveFrequency)
		{	case FluidFreqGummel: return true;
			case FluidFreqInner: return false;
			default: return prefersGummel(); //pick by fluid type
		}
	}
	
	//! Set total explicit charge density and effective electron density to use in cavity formation (i.e. including charge balls)
	//! and set list of explicit atoms to use in van der Waals corrections
	//! This base-class wrapper handles grid embedding (if necessary) and calls set_internal of the derived class
	void set(const ScalarFieldTilde& rhoExplicitTilde, const ScalarFieldTilde& nCavityTilde);

	//! Compute gradients with respect to electronic side variables (if non-null), and return fluid+coupling free energy
	//! Any extra forces on explicit ions due to the fluid should be stored in extraForces (if non-null)
	//! This base-class wrapper handles grid embedding (if necessary) and calls set_internal of the derived class
	double get_Adiel_and_grad(ScalarFieldTilde* Adiel_rhoExplicitTilde=0, ScalarFieldTilde* Adiel_nCavityTilde=0, IonicGradient* extraForces=0) const;

	//! Structure to represent frequency-dependent fluid susceptibility response
	//! [of the form prefactor[iOmega]*s[iSite] * 4*pi/(2*l+1) * Sum_m(Ylm(Ghat)*Ylm(Ghat')) * w(G)*w(G')]
	struct SusceptibilityTerm
	{	int iSite; //!< index to site / shape function (into sTilde)
		int l; //!< angular momentum of response
		const RadialFunctionG* w; //!< weight function, none if null
		std::vector<complex> prefactor; //!< frequency-dependent prefactor
	};
	//! Get susceptibility for a list of specified (complex) frequencies
	void getSusceptibility(const std::vector<complex>& omega, std::vector<SusceptibilityTerm>& susceptibility, ScalarFieldTildeArray& sTilde, bool elecOnly) const;

	virtual double bulkPotential() {return 0.0;}

	//! Dump relevant fluid densities (eg. NO and NH) to file(s)
	//! the provided pattern will have a single %s which may be substituted
	//! Fluid solver implementations may override to dump fluid densities, no dumping by default
	virtual void dumpDensities(const char* filenamePattern) const {};

	//! Dump fluid debugging quantities (eg. fluid potential and effective electron density for 3.0)
	//! the provided pattern will have a single %s which may be substituted
	//! Fluid solver implementations may override to dump fluid debug stuff, no dumping by default
	virtual void dumpDebug(const char* filenamePattern) const {};

	//------------Fluid solver implementations must provide these pure virtual functions

	//! Specify whether fluid prefers a gummel loop (true) or is minimized each time (false)
	virtual bool prefersGummel() const=0;

	//! Initialize fluid state from a file
	virtual void loadState(const char* filename)=0;

	//! Save fluid state to a file
	virtual void saveState(const char* filename) const=0;

	//! Minimize fluid side (holding explicit electronic system fixed)
	virtual void minimizeFluid()=0;
	
protected:
	//! Fluid-dependent implementation of set()
	virtual void set_internal(const ScalarFieldTilde& rhoExplicitTilde, const ScalarFieldTilde& nCavityTilde)=0;

	//! Fluid-dependent implementation of get_Adiel_and_grad()
	virtual double get_Adiel_and_grad_internal(ScalarFieldTilde& Adiel_rhoExplicitTilde, ScalarFieldTilde& Adiel_nCavityTilde, IonicGradient* extraForces) const =0;
	
	//! Fluid-dependent implementation of getSusceptibility()
	virtual void getSusceptibility_internal(const std::vector<complex>& omega, std::vector<SusceptibilityTerm>& susceptibility, ScalarFieldArray& sArr, bool elecOnly) const;
};

//! Create and return a JDFTx solver (the solver can be freed using delete)
FluidSolver* createFluidSolver(const Everything& e, const FluidSolverParams& params);

//! @}
#endif // JDFTX_ELECTRONIC_FLUIDSOLVER_H
