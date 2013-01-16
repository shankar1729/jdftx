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

//! @file FluidSolver.h
//! Common interface for all the fluids to the electronic code

#include <cstdio>
#include <core/GridInfo.h>
#include <core/DataIO.h>
#include <core/MinimizeParams.h>
#include <core/Coulomb.h>
#include <electronic/ExCorr.h>
#include <electronic/IonicMinimizer.h>
#include <fluid/FluidMixture.h>
#include <fluid/S2quad.h>
#include <fluid/Molecule.h>
#include <fluid/Fex_HardSphereIon.h>
#include <fluid/Fex_H2O_ScalarEOS.h>

typedef enum
{
	FluidNone, //!< No fluid
	FluidLinear, //!< Linear local-dielectric fluid [K. Letchworth-Weaver and T. A. Arias, Phys. Rev. B 86, 075140 (2012)]
	FluidLinearPCM, //!< Reparametrized linear local-dielectric fluid including non-electrostatic terms (EXPERIMENTAL)
	FluidNonlinearPCM, //!< Nonlinear local-dielectric fluid including non-electrostatic terms (EXPERIMENTAL)
	FluidNonlocalPCM, //!< Nonlocal Polarizable Continuum Model (EXPERIMENTAL)
	FluidFittedCorrelations, //!< Functional from [J. Lischner and T. A. Arias, J. Phys. Chem. B 114, 1946 (2010)]
	FluidScalarEOS, //!< Scalar EOS functional
	FluidScalarEOSCustom, //!< Scalar EOS functional with custom geometry or sites 
	FluidBondedVoids, //!< Functional from [R. Sundararaman, K. Letchworth-Weaver and T.A. Arias, J. Chem. Phys. 137, 044107 (2012)]
	FluidHSIonic //!< Functional of optionally charged hard spheres (EXPERIMENTAL)
}
FluidType; //!< Fluid type


//! Parameters controlling Non-local PCM 
struct NonlocalPCMparams
{	int lMax; //!< Maximum angular momentum to include
	//TODO: completely control NonlocalPCM from commands via this struct
	
	NonlocalPCMparams() : lMax(2) {}
};



//! Extra parameters for fluids:
struct FluidSolverParams
{
	double T; //!< temperature:
	bool verboseLog; //!< for Linear control iteration progress, for rest control sub-iteration progress output
	vector3<> zeroRef; //!<coordinates of zero reference of potential;
	
	//For the 1.0 solvers (Linear/Nonlinear) alone:
	double epsilonBulk; //!< Bulk dielectric constant
	double ionicConcentration; //!< Concentration of ions
	int ionicZelectrolyte; //!< Magnitude of charge on each ion
	double ionicRadiusPlus; //!< Hard-sphere radius of anion (Plus is wrt electron-positive convention)
	double ionicRadiusMinus; //!< Hard-sphere radius of cation (Minus is wrt electron-positive convention)
	double nc; //!< critical density for the 1.0 fluid shape function
	double sigma; //!< smoothing factor for the  1.0 fluid shape function
	double cavityTension; //! Surface tension (hartree per bohr^2) of the cavity
	double cavityPressure;  //! Volume dependent energy/tension of the cavity (simulates vdw)
	
	//For Nonlocal PCM (SaLSA) alone:
	NonlocalPCMparams npcmParams;
	
	//For Nonlinear1.0 alone:
	bool linearDielectric; //!< If true, work in the linear dielectric response limit
	bool linearScreening; //!< If true, work in the linearized Poisson-Boltzman limit for the ions
	double Nbulk; //!< Bulk number-density of molecules in bohr^-3
	double pMol; //!< Dipole moment of each molecule in e-bohr
	double epsInf; //! Optical-frequency dielectric constant
	
	//For Explicit Fluid JDFT alone:
	ConvolutionCouplingSiteModel convCouplingH2OModel; //!< selects parameter set for convolution coupling water
	double convCouplingScale; //!< scales von Weisacker correction in kinetic energy functional for fluid coupling 
	double VDWCouplingScale; //!< scales van der Waals correction for fluid coupling
	S2quadType s2quadType; //!< Quadrature on S2 that generates the SO(3) quadrature
	unsigned quad_nBeta, quad_nAlpha, quad_nGamma; //!< Subdivisions for euler angle outer-product quadrature
	
	double nuclearWidth; //!< Gaussian width of fluid site nuclear charge densities for coupling
	
//Kendra: probably get rid of block below and replace with H2OSites	
	double oxygenWidth; //!< Exponential width of electron distribution of oxygen atom in fluid
	double hydrogenWidth; //!< Exponential width of electron distribution of hydrogen atom in fluid
	double oxygenSiteCharge; //!< Oxygen site charge in water convolution coupling
	double hydrogenSiteCharge; //!< Hydrogen site charge in water convolution coupling
	double oxygenZnuc; //!< Oxygen nuclear charge in water convolution coupling
	double hydrogenZnuc; //!< Hydrogen nuclear charge in water convolution coupling
	string oxygenFilename; //!<  Filename which contains electron distribution of oxygen atom in fluid
	string hydrogenFilename; //!< Filename which contains electron distribution of hydrogen atom in fluid
	
	ExCorr exCorr; //!< Fluid exchange, correlation. and kinetic energy functional

	//For water in Explicit Fluid JDFT
	std::vector<H2OSite> H2OSites; //!< vector of H2OSite objects to be initialized in fluidMixture
	
	//For Ionic fluid in Explicit Fluid JDFT
	std::vector<HardSphereIon> hSIons; //!< vector of hard sphere ion objects to be initialized in fluidMixture
	
	FluidSolverParams() : verboseLog(false) {}
};


//! Abstract base class for the fluid solvers
struct FluidSolver
{
	const Everything& e;

	//! Abstract base class constructor - do not use directly - see FluidSolver::createSolver
	FluidSolver(const Everything &everything);
		
	virtual ~FluidSolver() {}

	//! Set total explicit charge density and effective electron density to use in cavity formation (i.e. including charge balls)
	//!  and set list of explicit atoms to use in van der Waals corrections
	virtual void set(const DataGptr& rhoExplicitTilde, const DataGptr& nCavityTilde)=0;

	//! Compute gradients with respect to electronic side variables, and return fluid+coupling free energy
	//! Any extra forces on explicit ions due to the fluid should be stored in extraForces
	virtual double get_Adiel_and_grad(DataGptr& grad_rhoExplicitTilde, DataGptr& grad_nCavityTilde, IonicGradient& extraForces)=0;

	//grad_rhoExplicitTilde is d_fluid and grad_nCavityTilde is V_cavity

	//! Dump relevant fluid densities (eg. NO and NH) to file(s)
	//! the provided pattern will have a single %s which may be substituted
	//! Fluid solver implementations may override to dump fluid densities, no dumping by default
	virtual void dumpDensities(const char* filenamePattern) const {};

	//! Dump fluid debugging quantities (eg. fluid potential and effective electron density for 3.0)
	//! the provided pattern will have a single %s which may be substituted
	//! Fluid solver implementations may override to dump fluid debug stuff, no dumping by default
	virtual void dumpDebug(const char* filenamePattern) const {};


	//------------Fluid solver implementations must provide these pure virtual functions

	//! Specify whether fluid requires a gummel loop (true) or is minimized each time (false)
	virtual bool needsGummel()=0;

	//! Initialize fluid state from a file
	virtual void loadState(const char* filename)=0;

	//! Save fluid state to a file
	virtual void saveState(const char* filename) const=0;

	//! Minimize fluid side (holding explicit electronic system fixed)
	virtual void minimizeFluid()=0;
};

//! Create and return a JDFTx solver (the solver can be freed using delete)
//! @param type Type of fluid
//! @param Nx Samples along first axis (outer index, slowest)
//! @param Ny Samples along second axis
//! @param Nz Samples along third axis (inner index, fastest)
//! @param R Lattice vectors
//! @param params Extra parameters, some functional specific
FluidSolver* createFluidSolver(FluidType type, const Everything& e, FluidSolverParams& params);

#endif // JDFTX_ELECTRONIC_FLUIDSOLVER_H
