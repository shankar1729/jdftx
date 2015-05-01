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
#include <core/PulayParams.h>

enum FluidType
{
	FluidNone, //!< No fluid
	FluidLinearPCM, //!< Linear local-dielectric fluid, optionally including non-electrostatic terms
	FluidNonlinearPCM, //!< Nonlinear local-dielectric fluid including non-electrostatic terms
	FluidSaLSA, //!< Spherically-averaged liquid susceptibility ansatz (nonlocal PCM)
	FluidClassicalDFT //!< Classical density functional description of fluid (EXPERIMENTAL)
};

	//!Mixing functional choices
enum FMixFunctional
{
	FMixNone, //!< No Interaction
	LJPotential, //!<Interaction with attractive part of Lennard Jones with sigma/eps potential only
	GaussianKernel, //!< Interaction with gaussian kernel with Rsolv/Esolv
};

//! Parameters needed to mix fluids
struct FmixParams
{
	std::shared_ptr<FluidComponent> fluid1,fluid2;
	FMixFunctional FmixType; //!<Type of Fmix to be used (GaussianKernel or LJPotential)
	double energyScale,lengthScale; //!<Energy scale (eps for LJ potential) and range parameter (sigma for LJ potential)
};

enum PCMVariant
{	PCM_SaLSA, //!< Use only with fluid type SaLSA [R. Sundararaman, K. Schwarz, K. Letchworth-Weaver, and T.A. Arias, JCP 142, 054102 (2015)]
	PCM_CANDLE, //!< Charge-asymmetric nonlocally-determined local-electric (CANDLE) solvation model [R. Sundararaman and W.A. Goddard III, JCP 142, 064107 (2015)]
	PCM_SGA13, //!< Local-response dielectric fluid or electrolyte with weighted-density cavitation and dispersion [R. Sundararaman, D. Gunceler and T.A. Arias, JCP 141, 134105 (2014)]
	PCM_GLSSA13, //!< Local-response dielectric fluid or electrolyte with empirical cavity tension [D. Gunceler, K. Letchworth-Weaver, R. Sundararaman, K.A. Schwarz and T.A. Arias, MSMSE 21, 074005 (2013)]
	PCM_LA12, //!< Linear local-response electrolyte [K. Letchworth-Weaver and T.A. Arias, Phys. Rev. B 86, 075140 (2012)]
	PCM_PRA05, //!< Linear local-response dielectric fluid [S.A. Petrosyan SA, A.A. Rigos and T.A. Arias, J Phys Chem B. 109, 15436 (2005)]
	PCM_SCCS_g09,      //!< g09 parametrization of SCCS local linear model for water [Andreussi et al. J. Chem. Phys. 136, 064102 (2012)]
	PCM_SCCS_g03,      //!< g03 parametrization of SCCS local linear model for water [Andreussi et al. J. Chem. Phys. 136, 064102 (2012)]
	PCM_SCCS_g03p,     //!< g03' parametrization of SCCS local linear model for water [Andreussi et al. J. Chem. Phys. 136, 064102 (2012)]
	PCM_SCCS_g09beta,  //!< g09+beta parametrization of SCCS local linear model for water [Andreussi et al. J. Chem. Phys. 136, 064102 (2012)]
	PCM_SCCS_g03beta,  //!< g03+beta parametrization of SCCS local linear model for water [Andreussi et al. J. Chem. Phys. 136, 064102 (2012)]
	PCM_SCCS_g03pbeta, //!< g03'+beta parametrization of SCCS local linear model for water [Andreussi et al. J. Chem. Phys. 136, 064102 (2012)]
	PCM_SCCS_cation,   //!< cations-only parametrization of SCCS local linear model for water [Dupont et al., J. Chem. Phys. 139, 214110 (2013)]
	PCM_SCCS_anion     //!< anions-only parametrization of SCCS local linear model for water [Dupont et al., J. Chem. Phys. 139, 214110 (2013)]
};

//! Check for any of the SCCS cases:
#define case_PCM_SCCS_any \
	case PCM_SCCS_g09: \
	case PCM_SCCS_g03: \
	case PCM_SCCS_g03p: \
	case PCM_SCCS_g09beta: \
	case PCM_SCCS_g03beta: \
	case PCM_SCCS_g03pbeta: \
	case PCM_SCCS_cation: \
	case PCM_SCCS_anion
inline bool isPCM_SCCS(PCMVariant pcmVariant) { switch(pcmVariant) { case_PCM_SCCS_any: return true; default: return false; } }


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
	
	//For CANDLE alone:
	double Ztot; //!< number of valence electrons
	double eta_wDiel; //!< control electrostatic weight function (gaussian convolved by delta(r-eta) at l=1) (fit parameter)
	double sqrtC6eff; //!< (effective C6 parameter in J-nm^6/mol)^(1/2) for the entire molecule (fit parameter) (vdwScale unnecessary and not used due to this)
	double pCavity; //!< sensitivity of cavity to surface electric field to emulate charge asymmetry [e-bohr/Eh]  (fit parameter)
	
	//For SCCS alone:
	double rhoMin, rhoMax; //!< start and end of transition
	double rhoDelta; //!< Delta used for "quantum surface"
	double cavityPressure; //!< volume term (used in some parametrizations)
	
	//For SaLSA alone:
	int lMax;
	
	//Debug parameters for Nonlinear PCM's:
	bool linearDielectric; //!< If true, work in the linear dielectric response limit
	bool linearScreening; //!< If true, work in the linearized Poisson-Boltzman limit for the ions
	bool nonlinearSCF; //!< whether to use an SCF method for nonlinear PCMs
	PulayParams scfParams; //!< parameters controlling Pulay mixing for SCF version of nonlinear PCM
	
	//For Explicit Fluid JDFT alone:
	ExCorr exCorr; //!< Fluid exchange-correlation and kinetic energy functional
        std::vector<FmixParams> FmixList; //!< Tabulates which components interact through an additional Fmix

	string initWarnings; //!< warnings emitted during parameter initialization, if any
	
	FluidSolverParams();
	void setPCMparams(); //!< Set predefined parameters for solventName (for a specific model)
	void setCDFTparams(); //!< Set predefined parameters for solventName (for a classical DFT model)
	bool needsVDW() const; //!< whether pair-potential vdW corrections are required
	bool ionicScreening() const; //!< whether list of fluid components includes ionic species for Debye screening
private:
	std::vector< std::shared_ptr<FluidComponent> > components_, solvents_, cations_, anions_; //internal mutable versions of the public component lists
};

#endif // JDFTX_ELECTRONIC_FLUIDSOLVERPARAMS_H
