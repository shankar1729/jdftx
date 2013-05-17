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

#include <electronic/ExCorr.h>
#include <fluid/FluidComponent.h>

enum FluidType
{
	FluidNone, //!< No fluid
	FluidLinearPCM, //!< Linear local-dielectric fluid, optionally including non-electrostatic terms
	FluidNonlinearPCM, //!< Nonlinear local-dielectric fluid including non-electrostatic terms
	FluidNonlocalPCM, //!< Nonlocal polarizable continuum model (EXPERIMENTAL)
	FluidClassicalDFT //!< Classical density functional description of fluid (EXPERIMENTAL)
};

enum PCMVariant
{	PCM_SLSA13, //!< Non-local PCM, use only with fluid type NonlocalPCM [R. Sundararaman, K. Letchworth-Weaver, K. Schwarz and T.A. Arias (under preparation)]
	PCM_SGA13, //!< Local-response dielectric fluid or electrolyte with weighted-density cavotation and dispersion [R. Sundararaman, D. Gunceler and T.A. Arias, (under preparation)]
	PCM_GLSSA13, //!< Local-response dielectric fluid or electrolyte with empirical cavity tension [D. Gunceler, K. Letchworth-Weaver, R. Sundararaman, K.A. Schwarz and T.A. Arias, arXiv:1301.6189]
	PCM_LA12, //!< Linear local-response electrolyte [K. Letchworth-Weaver and T.A. Arias, Phys. Rev. B 86, 075140 (2012)]
	PCM_PRA05 //!< Linear local-response dielectric fluid [S.A. Petrosyan SA, A.A. Rigos and T.A. Arias, J Phys Chem B. 109, 15436 (2005)]
};

//! Extra parameters for fluids:
struct FluidSolverParams
{
	FluidType fluidType;
	PCMVariant pcmVariant;
	
	double T; //!< temperature
	double P; //!< pressure
	double epsBulkOverride, epsInfOverride; //!< Override default dielectric constants if non-zero
	bool verboseLog; //!< whether iteration progress is printed for Linear PCM's, and whether sub-iteration progress is printed for others
	
	const std::vector< std::shared_ptr<FluidComponent> >& components; //!< list of all fluid components
	const std::vector< std::shared_ptr<FluidComponent> >& solvents; //!< list of solvent components
	const std::vector< std::shared_ptr<FluidComponent> >& cations; //!< list of cationic components
	const std::vector< std::shared_ptr<FluidComponent> >& anions; //!< list of anionic components
	
	void addComponent(const std::shared_ptr<FluidComponent>& component); //!< Add component to the component list as well as one of solvents, anions or cations as appropriate
	
	//Fit parameters:
	double nc; //!< critical density for the PCM cavity shape function
	double sigma; //!< smoothing factor for the PCM cavity shape function
	double cavityTension; //!< effective surface tension (including dispersion etc.) of the cavity (hartree per bohr^2)
	double vdwScale; //!< overall scale factor for Grimme pair potentials (or damping range scale factor for vdw-TS when implemented)
	
	//For Nonlocal PCM (SaLSA) alone:
	int lMax;
	
	//Debug parameters for Nonlinear PCM's:
	bool linearDielectric; //!< If true, work in the linear dielectric response limit
	bool linearScreening; //!< If true, work in the linearized Poisson-Boltzman limit for the ions

	//For Explicit Fluid JDFT alone:
	ExCorr exCorr; //!< Fluid exchange-correlation and kinetic energy functional

	string initWarnings; //!< warnings emitted during parameter initialization, if any
	
	FluidSolverParams();
	void setPCMparams(); //!< Set predefined parameters for solventName (for a specific model)
	bool needsVDW() const; //!< whether pair-potential vdW corrections are required
	bool ionicScreening() const; //!< whether list of fluid components includes ionic species for Debye screening
private:
	std::vector< std::shared_ptr<FluidComponent> > components_, solvents_, cations_, anions_; //internal mutable versions of the public component lists
};

#endif // JDFTX_ELECTRONIC_FLUIDSOLVERPARAMS_H
