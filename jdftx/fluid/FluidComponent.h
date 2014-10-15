/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman, Kendra Letchworth-Weaver

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

#ifndef JDFTX_FLUID_FLUIDCOMPONENT_H
#define JDFTX_FLUID_FLUIDCOMPONENT_H

#include <fluid/S2quad.h>
#include <fluid/Molecule.h>

//!Named fluid components for which bulk properties / geometries / excess functionals are available
struct FluidComponent
{
	//!Component names
	enum Name
	{	//Neutral solvent molecules:
		H2O, //!< Water
		CHCl3, //!< Chloroform
		CCl4, //!< Carbon tetrachloride
		CH3CN, //!< Acetonitrile
		DMC, //!< Dimethyl carbonate
		EC, //!< Ethylene carbonate
		PC, //!< Propylene carbonate
		DMF, //!< Dimethylformamide
		THF, //!< Tetrahydrofuran
		EthylEther, //!< Diethyl ether (http://en.wikipedia.org/wiki/Diethyl_ether_(data_page))
		Chlorobenzene, //! (http://en.wikipedia.org/wiki/Chlorobenzene)
		Isobutanol, //! (http://en.wikipedia.org/wiki/Isobutanol)
		CarbonDisulfide, //! (http://en.wikipedia.org/wiki/Carbon_disulfide),
		DMSO, //! Dimethyl sulfoxide
		CH2Cl2, //! Dichloromethane (methyl chloride)
		Ethanol,
		Methanol,
		Octanol,
		Glyme,
		EthyleneGlycol,
		CustomSolvent, //!< Custom solvent (molecule and other parameters initialized manually)
		//Cations:
		Sodium, //!< Na+
		HydratedSodium, //!< Na+ cluster with 6 octahedrally coordinated H2O molecules
		Potassium, //!< K+
		HydratedPotassium, //!< K+ cluster with 6 octahedrally coordinated H2O molecules
		Hydronium, //!< H3O+
		HydratedHydronium, //!< H3O+ cluster with 4 H2O molecules
		CustomCation, //!< Custom cation (molecule and other parameters initialized manually)
		//Anions:
		Chloride, //!< Cl-
		Fluoride, //!< Fl-
		Perchlorate, //!< ClO4-
		Hydroxide, //!< OH-
		HydratedHydroxide, //!< OH- cluster with 4 H2O molecules
		CustomAnion //!< Custom anion (molecule and other parameters initialized manually)
	};
	const Name name;
	
	//!Type of component - used to determine role of component in simpler PCMs
	enum Type
	{	Solvent, Cation, Anion
	};
	const Type type;
	static Type getType(Name name);
	
	//!Excess functional choices
	enum Functional
	{	ScalarEOS, //!< Generic hard sphere + weighted density functional constrained to equation of state based on [R. Sundararaman and T.A. Arias, arXiv:1302.0026]
		FittedCorrelations, //!< H2O functional from [J. Lischner and T. A. Arias, J. Phys. Chem. B 114, 1946 (2010)] (DEPRECATED)
		BondedVoids, //!< H2O functional from [R. Sundararaman, K. Letchworth-Weaver and T.A. Arias, J. Chem. Phys. 137, 044107 (2012)]
		MeanFieldLJ, //!< Hard sphere + mean field Lennard-Jones perturbation (useful for ions in solution)
		FunctionalNone //!< No excess functional beyond hard spheres / electrostatics (or fex may be initialized manually)
	};
	const Functional functional;
	double epsLJ; //!< Lennard-Jones well depth for mean-field LJ functional
	
	//!Ideal gas representation (does not affect simple fluids which always use IdealGasMonoatomic)
	enum Representation
	{	Pomega, //!< directly work with orientation probability density
		PsiAlpha, //!< site-potential representation
		MuEps //!< multipole density representation truncated at l=1 (default)
	}
	representation;
	
	S2quadType s2quadType; //!< Quadrature on S2 that generates the SO(3) quadrature (default: Octahedron) HACK: Kendra made default 7design24
	unsigned quad_nBeta, quad_nAlpha, quad_nGamma; //!< Subdivisions for euler angle outer-product quadrature
	
	enum TranslationMode
	{	ConstantSpline,
		LinearSpline, //!< default and recommended
		Fourier
	}
	translationMode; //!< type of translation operator used for sampling rigid molecule geometry
	
	//Bulk solvent properties (used by various PCM's):
	double epsBulk; //!< bulk dielectric constant
	double Nbulk; //!< bulk number-density of molecules in bohr^-3 (used as initial guess in mixtures)
	double pMol; //!< dipole moment of each molecule in e-bohr
	double epsInf; //!< optical-frequency dielectric constant
	double Pvap; //!< vapor pressure in Eh/bohr^3
	double sigmaBulk; //!< bulk surface tension in Eh/bohr^2
	double Rvdw; //!< effective van der Waals radius of liquid (derived from equation of state) in bohrs
	double Res; //!< electrostatic radius of solvent (derived from nonlocal response) in bohrs
	
	double Nnorm; //!< If Nnorm>=0, this component is switched to the canonical ensemble (number fixed to Nnorm)
	
	//Molecule geometry and site properties:
	Molecule molecule;

	double pureNbulk(double T) const; //!< get density of solvent component in the pure phase at temperature T (returns 1 mol/liter for ions)

	FluidComponent(Name name, double T, Functional functional); //!< set default properties
	
	//Extra properties when participating in a classical density functional FluidMixture:
	std::shared_ptr<class SO3quad> quad; //!< orientation quadrature
	std::shared_ptr<class TranslationOperator> trans; //!< translation operator
	std::shared_ptr<class IdealGas> idealGas; //!< Indep <-> Density converter and entropy calculator
	std::shared_ptr<class Fex> fex; //!< Excess functional (in excess to sphere mixture and long-range)
	std::shared_ptr<class ScalarEOS> eos; //!< Equation of states for ScalarEOS functionals
	unsigned offsetIndep; //!< Offset to the independent variables of this component
	unsigned offsetDensity; //!< Offset to the site densities that belong to this component
	void addToFluidMixture(class FluidMixture* fluidMixture); //!< Initialize fex, idealGas and register with fluidMixture
};

#endif // JDFTX_FLUID_FLUIDCOMPONENT_H
