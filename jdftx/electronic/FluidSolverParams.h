/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman, Kendra Letchworth Weaver, Deniz Gunceler

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

#ifndef JDFTX_ELECTRONIC_FLUIDSOLVERPARAMS_H
#define JDFTX_ELECTRONIC_FLUIDSOLVERPARAMS_H

//! @file FluidSolverParams.h
//! Parameters describing the fluids in the electronic code

#include <fluid/S2quad.h>
#include <fluid/Fex_HardSphereIon.h>
#include <electronic/ExCorr.h>

enum FluidType
{
	FluidNone, //!< No fluid

	FluidLinearPCM, //!< Linear local-dielectric fluid, optionally including non-electrostatic terms
	FluidNonlinearPCM, //!< Nonlinear local-dielectric fluid including non-electrostatic terms
	FluidNonlocalPCM, //!< Nonlocal Polarizable Continuum Model (EXPERIMENTAL)
	FluidFittedCorrelations, //!< Functional from [J. Lischner and T. A. Arias, J. Phys. Chem. B 114, 1946 (2010)]
	FluidScalarEOS, //!< Scalar EOS functional
	FluidScalarEOSCustom, //!< Scalar EOS functional with custom geometry or sites 
	FluidBondedVoids, //!< Functional from [R. Sundararaman, K. Letchworth-Weaver and T.A. Arias, J. Chem. Phys. 137, 044107 (2012)]
	FluidHSIonic //!< Functional of optionally charged hard spheres (EXPERIMENTAL)
};

enum PCMVariant
{	PCM_SGA13, //!< Local-response dielectric fluid or electrolyte with weighted-density cavotation and dispersion [R. Sundararaman, D. Gunceler and T.A. Arias, (under preparation)]
	PCM_GLSSA13, //!< Local-response dielectric fluid or electrolyte with empirical cavity tension [D. Gunceler, K. Letchworth-Weaver, R. Sundararaman, K.A. Schwarz and T.A. Arias, arXiv:1301.6189]
	PCM_LA12, //!< Linear local-response electrolyte [K. Letchworth-Weaver and T.A. Arias, Phys. Rev. B 86, 075140 (2012)]
	PCM_PRA05 //!< Linear local-response dielectric fluid [S.A. Petrosyan SA, A.A. Rigos and T.A. Arias, J Phys Chem B. 109, 15436 (2005)]
};

//! Parameters controlling Non-local PCM 
struct NonlocalPCMparams
{	int lMax; //!< Maximum angular momentum to include
	//TODO: completely control NonlocalPCM from commands via this struct
	
	NonlocalPCMparams() : lMax(2) {}
};



//! Extra parameters for fluids:
struct FluidSolverParams
{
	FluidType fluidType;
	PCMVariant pcmVariant;
	
	double T; //!< temperature
	bool verboseLog; //!< whether iteration progress is printed for Linear PCM's, and whether sub-iteration progress is printed for others
	vector3<> zeroRef; //!<coordinates of zero reference of potential;
	
	//Bulk solvent properties (used by various PCM's):
	double epsBulk; //!< bulk dielectric constant
	double Nbulk; //!< bulk number-density of molecules in bohr^-3
	double pMol; //!< dipole moment of each molecule in e-bohr
	double epsInf; //!< optical-frequency dielectric constant
	double Pvap; //!< vapor pressure in Eh/bohr^3
	double sigmaBulk; //!< bulk surface tension in Eh/bohr^2
	double Rvdw; //!< effective van der Waals radius of liquid (derived from equation of state) in bohrs
	
	//PCM fit parameters:
	double nc; //!< critical density for the PCM cavity shape function
	double sigma; //!< smoothing factor for the PCM cavity shape function
	double cavityTension; //!< effective surface tension (including dispersion etc.) of the cavity (hartree per bohr^2)
	
	//For the PCM's alone:
	double ionicConcentration; //!< Concentration of ions
	int ionicZelectrolyte; //!< Magnitude of charge on each ion
	double ionicRadiusPlus; //!< Hard-sphere radius of anion (Plus is wrt electron-positive convention)
	double ionicRadiusMinus; //!< Hard-sphere radius of cation (Minus is wrt electron-positive convention)
	
	//For Nonlocal PCM (SaLSA) alone:
	NonlocalPCMparams npcmParams;
	
	//Debug parameters for Nonlinear PCM's:
	bool linearDielectric; //!< If true, work in the linear dielectric response limit
	bool linearScreening; //!< If true, work in the linearized Poisson-Boltzman limit for the ions
	
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
	
	FluidSolverParams();
	
	//! Names of solvents with predefined PCM parameters
	enum SolventName
	{	H2O, //!< Water
		CHCl3, //!< Chloroform
		CCl4, //!< Carbon tetrachloride
		DMC, //!< Dimethyl carbonate
		EC, //!< Ethylene carbonate
		PC, //!< Propylene carbonate
		DMF, //!< Dimethylformamide
		THF //!< Tetrahydrofuran
	} solventName; //!< solvent to load default parameters for
	
	//! Set predefined parameters for solventName (for a specific model)
	void setPCMparams();
	
	string initWarnings; //!< warnings emitted during parameter initialization, if any
	
	bool needsVDW() const; //!< whether pair-potential vdW corrections are required
};

#endif // JDFTX_ELECTRONIC_FLUIDSOLVERPARAMS_H