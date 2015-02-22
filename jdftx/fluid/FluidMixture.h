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

#ifndef JDFTX_FLUID_FLUIDMIXTURE_H
#define JDFTX_FLUID_FLUIDMIXTURE_H

#include <fluid/FluidComponent.h>
#include <fluid/Fmix.h>
#include <core/Units.h>
#include <core/Minimize.h>
#include <core/EnergyComponents.h>

//! @brief Mixture of fluids that provides the total free energy functional for minimization
//! Constructing Fex and IdealGas objects require a FluidMixture reference, to which they add themselves.
//! The FluidMixture object is ready to use after initialize() is called.
class FluidMixture : public Minimizable<ScalarFieldArray>
{
public:
	const GridInfo& gInfo;
	const double T; //!< Temperature
	bool verboseLog; //!< print energy components etc. if enabled (off by default)
	vector3<> Eexternal; //!< External uniform electric field

	FluidMixture(const GridInfo& gInfo, const double T=298*Kelvin);
	virtual ~FluidMixture();

	//! Call after initializing and adding all the components of the fluid mixture
	//! This calculates the bulk equilibrium densities for all components and corresponding chemical potentials
	//! so as to achieve pressure P with the specified mole fractions.
	//! The initial guess for the densities is taken from FluidComponent::Nbulk, and this may be altered to select a particular phase
	//! The dielectric constants for mixtures is estimated as a weighted combination of the components, and may be overriden by non-zero overrides
	void initialize(double P=1.01325*Bar, double epsBulkOverride=0., double epsInfOverride=0.);

	//! Compute the boiling pressure (vapor-liquid equilbirum) at current temperature,
	//! given guesses for the total liquid and vapor densities, and optionally retrieve equilibrium vapor densities
	double getBoilingPressure(double NliqGuess, double NvapGuess, std::vector<double>* Nvap=0) const;
	
	unsigned get_nIndep() const { return nIndepIdgas + (polarizable ? 3 : 0); }  //!< get the number of scalar fields used as independent variables
	unsigned get_nDensities() const { return nDensities; } //!< get the total number of site densities

	const std::vector<const FluidComponent*>& getComponents() const; //!< access component list
	ScalarFieldArray state;

	//! External charge density.
	//! The interaction energy between internal charges and rhoExternal is included in getFreeEnergy() etc.
	//! When charged components are present, the unit cell will be neutral including this charge.
	ScalarFieldTilde rhoExternal;
	bool useMFKernel; //!If true, use the mean field kernel for external coulomb interactions as well 
	double Qtol; //!< tolerance for unit cell neutrality (default 1e-12)

	//! Initialize the independent variables
	//! @param scale scale the state that would produce the equilibrium ideal gas densities by this amount to ge tthe guess
	//! @param Elo Lower cap on the individiual molecule energy configurations used in the estimate
	//! @param Ehi Upper cap on the individiual molecule energy configurations used in the estimate
	void initState(double scale = 0.0, double Elo=-DBL_MAX, double Ehi=+DBL_MAX);

	//! Load the state from a single binary file
	void loadState(const char* filename);

	//! Save the state to a single binary file
	void saveState(const char* filename) const;

	//! Optional outputs for operator() and getFreeEnergy(), retrieve results for all non-null pointers
	struct Outputs
	{	
		ScalarFieldArray* N; //!< site densities
		vector3<>* electricP; //!< total electric dipole moment in cell (useful only with multipole-based IdealGas's)
		ScalarFieldTilde* Phi_rhoExternal; //!< derivative of free energy w.r.t rhoExternal
		ScalarFieldArray* psiEff; //!< Estimate ideal gas effective potentials (useful only when no electric field or potential on non-indep sites)
		EnergyComponents* Phi; //!< retrieve fluid energy components
		
		//! initialize all the above to null
		Outputs( ScalarFieldArray* N=0, vector3<>* electricP=0,
			ScalarFieldTilde* Phi_rhoExternal=0, ScalarFieldArray* psiEff=0, EnergyComponents* Phi=0);
	};

	//! @brief Free energy and gradient evaluation
	//! @param[in] indep Current value of the independent variables
	//! @param[out] Phi_indep Gradient w.r.t indep at the current value of indep
	//! @param[out] outputs optional outputs, see Outputs
	//! @return Free energy difference (compared to uniform fluid) at this value of indep
	//! Implemented in FluidMixtureCompute.cpp
	double operator()(const ScalarFieldArray& indep, ScalarFieldArray& Phi_indep, Outputs outputs=Outputs()) const;

	//! @brief Get the free energy, densities and moments for the current state
	//! @param[out] outputs optional outputs, see Outputs
	//! @return Free energy at current state
	double getFreeEnergy(Outputs outputs=Outputs()) const;

	//! Advance step along direction dir by scale alpha (interface for minimize())
	void step(const ScalarFieldArray& dir, double alpha);

	//! Return energy at current state, and optionally the gradient (interface for minimize())
	double compute(ScalarFieldArray* grad);

	ScalarFieldArray precondition(const ScalarFieldArray& grad);
	double sync(double x) const; //!< All processes minimize together; make sure scalars are in sync to round-off error
	
private:
	unsigned nIndepIdgas; //!< number of scalar fields used as independent variables for the component ideal gases
	unsigned nDensities; //!< total number of site densities
	bool polarizable;  //!< whether an additional vector field is required due to polarizable components
	double p;
	double Crot, Cpol; //!< dielectric correlation prefactors (set by initialize())
	
	std::vector<const FluidComponent*> component; //!< array of fluid components
	std::vector<const Fmix*> fmixArr; //!< array of mixing functionals

	std::vector<double> Kindep; //!< preconditioner weight for components
	double Keps; //!< preconditioner weight for polarizability independent variable
	
	void addComponent(FluidComponent*); //!< Called by FluidComponent::addToFluidMixture() to add a component
	friend class FluidComponent;

	void addFmix(const Fmix* fmix); //!< Called by Fmix::Fmix() to add a mixing functional
	friend class Fmix;

	//! Compute the uniform excess free energy and gradient w.r.t component molecular densities (implemented in FluidMixtureCompute.cpp)
	double computeUniformEx(const std::vector< double >& Nmol, std::vector< double >& Phi_Nmol) const;

	//! Compute the pressure of the uniform fluid mixture of total molecular density Ntot
	double compute_p(double Ntot) const;
	
	friend class BoilingPressureSolver;
	static ScalarFieldTilde coulomb(const ScalarFieldTilde& rho) { return (-4*M_PI) * Linv(O(rho)); }
};

#endif // JDFTX_FLUID_FLUIDMIXTURE_H
