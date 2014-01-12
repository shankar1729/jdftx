/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

This file is part of Fluid1D.

Fluid1D is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Fluid1D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Fluid1D.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#ifndef FLUID1D_FLUID1D_FLUIDMIXTURE_H
#define FLUID1D_FLUID1D_FLUIDMIXTURE_H

#include <fluid1D/IdealGas.h>
#include <fluid1D/Fmix.h>
#include <core/Units.h>
#include <core/Minimize.h>

//! @brief Mixture of fluids that provides the total free energy functional for minimization
//! Constructing Fex and IdealGas objects require a FluidMixture reference, to which they add themselves.
//! The FluidMixture object is ready to use after setPressure() is called.
class FluidMixture : public Minimizable<ScalarFieldCollection>
{
public:
	const GridInfo& gInfo;
	const double T; //!< Temperature
	bool verboseLog; //!< print energy components etc. if enabled (off by default)
	
	double Eexternal; //!< External uniform electric field (Note: do not use IdealGasPsiAlpha with polar fluids and a next external field)
	double Kpol; //!< preconditioning factor for polarization (sometimes useful to tune)
	
	FluidMixture(const GridInfo& gInfo, const double T=298*Kelvin);

	//! Call after initializing and adding all the components of the fluid mixture
	//! This calculates the bulk equilibrium densities for all components and corresponding chemical potentials
	//! so as to achieve pressure P with the specified mole fractions.
	//! Nguess is the initial guess for the total molecular density used for this adjustment,
	//! This may be used to select a particular phase (eg. use a low Nguess to pick up the vapor mixture)
	void setPressure(double P=1.01325*Bar, double Nguess=5e-3);

	//! Set the pressure to the boiling point at given temperature
	//! The liquid phase is constrained to have the supplied mole fractions
	//! @param Nvap If non-null, conatins the equilibirum vapor densities on return
	//! @param tol Tolerance on the solve (on Delta(P)/NT and Delta(mu)/T)
	//! @param NliqGuess A guess for the total liquid density at equilibrium
	//! @param NvapGuess A guess for the total vapor density at equilibrium
	//! @return The pressure corresponding to chemical and physical equilibrium
	double setBoilingPressure(std::vector<double>* Nvap=0, double tol=1e-8, double NliqGuess=5e-3, double NvapGuess=1e-6);
	
	unsigned get_nIndep() const { return nIndepIdgas + (polarizable ? 1 : 0); }  //!< get the number of scalar fields used as independent variables
	unsigned get_nDensities() const { return nDensities; } //!< get the total number of site densities

	unsigned get_offsetDensity(const Fex* fex) const; //!< get the offset for the site densities for a particular component in the fluid

	//! One component of the fluid mixture
	struct Component
	{	IdealGas* idealGas; //!< Indep <-> Density converter and entropy calculator
		const Fex* fex; //!< Excess functional (in excess to sphere mixture and long-range)
		const Molecule* molecule; //!< Molecule geometry etc.
		unsigned offsetIndep; //!< Offset to the independent variables of this component
		unsigned offsetDensity; //!< Offset to the site densities that belong to this component
		std::vector<SiteProperties*> indexedSite; //!< site properties indexed by density index
		std::vector<int> indexedSiteMultiplicity; //!< number of sites at each index
	};

	unsigned get_nComponents() const; //!< get number of components
	const Component& get_component(unsigned c) const; //!< get component number c
	ScalarFieldCollection state;

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
		ScalarFieldCollection* N; //!< site densities
		double* electricP; //!< total electric dipole moment in cell (useful only with multipole-based IdealGas's)
		ScalarFieldCollection* psiEff; //!< Estimate ideal gas effective potentials (useful only when no electric field or potential on non-indep sites)
		
		//! initialize all the above to null
		Outputs( ScalarFieldCollection* N=0, double* electricP=0, ScalarFieldCollection* psiEff=0);
	};

	//! @brief Free energy and gradient evaluation
	//! @param[in] indep Current value of the independent variables
	//! @param[out] grad_indep Gradient w.r.t indepo at the current value of indep
	//! @param[out] outputs optional outputs, see Outputs
	//! @return Free energy difference (compared to uniform fluid) at this value of indep
	double operator()(const ScalarFieldCollection& indep, ScalarFieldCollection& grad_indep, Outputs outputs=Outputs()) const;

	//! @brief Get the free energy, densities and moments for the current state
	//! @param[out] outputs optional outputs, see Outputs
	//! @return Free energy at current state
	double getFreeEnergy(Outputs outputs=Outputs()) const;

	//! Advance step along direction dir by scale alpha (interface for minimize())
	void step(const ScalarFieldCollection& dir, double alpha);

	//! Return energy at current state, and optionally the gradient (interface for minimize())
	double compute(ScalarFieldCollection* grad);

	//! Preconditioner: scale by inverse quadrature weights on each channel
    ScalarFieldCollection precondition(const ScalarFieldCollection& grad);


	//! Get the direct correlations c = (-1/T) times the second variational derivative of the excess free energy
	//! evaluated at the uniform fluid mixture with component densities Nmol (implemented in FluidMixtureCorrFunc.cpp)
	//! The output is nDensities*(nDensities+1)/2 in order (11),(21),(22),(31),(32),(33),...
	ScalarFieldTildeCollection getDirectCorrelations(const std::vector<double>& Nmol) const;
	
	//! Get the partial pair correlations g evaluated at the uniform fluid mixture
	//! with component densities Nmol (implemented in FluidMixtureCorrFunc.cpp)
	//! The output is nDensities*(nDensities+1)/2 in order (11),(21),(22),(31),(32),(33),...
	ScalarFieldCollection getPairCorrelations(const std::vector<double>& Nmol) const;
	
	//! Indexing scheme for correlation functions (implemented in FluidMixtureCorrFunc.cpp)
	//! i1 and i2 are absolute density indices of fex is null
	//! and are component relative indices if fex is non-null
	int corrFuncIndex(unsigned i1, unsigned i2, const Fex* fex=0) const;

private:
	unsigned nDensities; //!< total number of site densities
	unsigned nIndepIdgas; //!< number of scalar fields used as independent variables for the component ideal gases
	bool polarizable; //!< whether additional scalar fields are required due to polarizable components
	double p;

	std::vector<Component> component; //!< array of fluid components
	std::vector<const Fmix*> fmixArr; //!< array of mixing functionals

	void addComponent(IdealGas* idealGas, const Fex* fex); //!< Called by IdealGas::IdealGas() to add a component
	friend class IdealGas;

	void addFmix(const Fmix* fmix); //!< Called by Fmix::Fmix() to add a mixing functional
	friend class Fmix;

	//! Compute the uniform excess free energy and gradient w.r.t component molecular densities
	double computeUniformEx(const std::vector< double >& Nmol, std::vector< double >& grad_Nmol) const;

	//! Compute the pressure of the uniform fluid mixture of total molecular density Ntot
	double compute_p(double Ntot) const;
	
	friend class BoilingPressureSolver;
};

#endif // FLUID1D_FLUID1D_FLUIDMIXTURE_H
