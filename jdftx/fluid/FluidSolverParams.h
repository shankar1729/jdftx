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

//! @addtogroup Fluid
//! @{
//! @file FluidSolverParams.h
//! Parameters describing the fluids in the electronic code

#include <electronic/ExCorr.h>
#include <fluid/FluidComponent.h>
#include <core/PulayParams.h>

enum FluidType
{
	FluidNone, //!< No fluid
	FluidLinearPCM, //!< Linear local-dielectric fluid, optionally including non-electrostatic terms \cite PCM-Kendra
	FluidNonlinearPCM, //!< Nonlinear local-dielectric fluid including non-electrostatic terms \cite NonlinearPCM
	FluidSaLSA, //!< Spherically-averaged liquid susceptibility ansatz \cite SaLSA
	FluidClassicalDFT //!< Classical density functional description of fluid \cite PolarizableCDFT \cite RigidCDFT
};

enum FluidSolveFrequency
{
	FluidFreqInner, //!< Solve fluid every electronic step
	FluidFreqGummel, //!< Use a Gummel iteration
	FluidFreqDefault //!< Decide based on fluid type (Inner for linear fluids, Gummel for rest)
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
{	PCM_SaLSA, //!< Use only with fluid type FluidSaLSA \cite SaLSA
	PCM_CANON, //!< Use only with fluid type FluidNonlinearPCM \cite CANON
	PCM_CANDLE, //!< Charge-asymmetric nonlocally-determined local-electric (CANDLE) solvation model \cite CANDLE
	PCM_SGA13, //!< Local-response dielectric fluid or electrolyte with weighted-density cavitation and dispersion \cite CavityWDA
	PCM_GLSSA13, //!< Local-response dielectric fluid or electrolyte with empirical cavity tension \cite NonlinearPCM
	PCM_LA12, //!< Linear local-response electrolyte \cite PCM-Kendra
	PCM_SoftSphere, //!< Soft-sphere continuum solvation model \cite PCM-SoftSphere
	PCM_FixedCavity, //!< Fixed-cavity electrostatic-only continuum solvation model
	PCM_SCCS_g09,      //!< g09 parametrization of SCCS local linear model for water \cite PCM-SCCS
	PCM_SCCS_g03,      //!< g03 parametrization of SCCS local linear model for water \cite PCM-SCCS
	PCM_SCCS_g03p,     //!< g03' parametrization of SCCS local linear model for water  \cite PCM-SCCS
	PCM_SCCS_g09beta,  //!< g09+beta parametrization of SCCS local linear model for water \cite PCM-SCCS
	PCM_SCCS_g03beta,  //!< g03+beta parametrization of SCCS local linear model for water \cite PCM-SCCS
	PCM_SCCS_g03pbeta, //!< g03'+beta parametrization of SCCS local linear model for water \cite PCM-SCCS
	PCM_SCCS_cation,   //!< cations-only parametrization of SCCS local linear model for water \cite PCM-SCCS-charged
	PCM_SCCS_anion     //!< anions-only parametrization of SCCS local linear model for water \cite PCM-SCCS-charged
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
	vector3<> epsBulkTensor; //!< Override default dielectric constants with a tensor if non-zero (assuming Cartesian coords are principal axes, LinearPCM only)
	bool verboseLog; //!< whether iteration progress is printed for Linear PCM's, and whether sub-iteration progress is printed for others
	FluidSolveFrequency solveFrequency;
	
	const std::vector< std::shared_ptr<FluidComponent> >& components; //!< list of all fluid components
	const std::vector< std::shared_ptr<FluidComponent> >& solvents; //!< list of solvent components
	const std::vector< std::shared_ptr<FluidComponent> >& cations; //!< list of cationic components
	const std::vector< std::shared_ptr<FluidComponent> >& anions; //!< list of anionic components
	
	void addComponent(const std::shared_ptr<FluidComponent>& component); //!< Add component to the component list as well as one of solvents, anions or cations as appropriate
	
	//Fit parameters:
	double nc; //!< critical density for the PCM cavity shape function
	double sigma; //!< smoothing factor for the PCM cavity shape function (dimensionless for most, but in bohrs for SoftSphere)
	double cavityTension; //!< effective surface tension (including dispersion etc.) of the cavity (hartree per bohr^2)
	double vdwScale; //!< overall scale factor for Grimme pair potentials (or damping range scale factor for vdw-TS when implemented)
	
	//For CANDLE alone:
	double Ztot; //!< number of valence electrons (also used by CANON)
	double eta_wDiel; //!< electrostatic cavity expansion widthin bohrs (fit parameter)
	double sqrtC6eff; //!< effective C6 parameter in J-nm^6/mol)^(1/2) for the entire molecule (fit parameter) (vdwScale unnecessary and not used due to this); also used by CANON
	double pCavity; //!< sensitivity of cavity to surface electric field to emulate charge asymmetry [e-a0/Eh]  (fit parameter)
	
	//For SCCS alone:
	double rhoMin, rhoMax; //!< start and end of transition
	double rhoDelta; //!< Delta used for "quantum surface"
	double cavityPressure; //!< volume term (used in some parametrizations)
	
	//For SaLSA alone:
	int lMax;
	
	//For CANON alone:
	double Res; //! Electrostatic radius used for dielectric nonlocality in CANON
	double Zcenter; //!< Charge at center used to determine asymmetry in CANON
	
	//For soft sphere model alone:
	double getAtomicRadius(const class SpeciesInfo& sp) const; //!< get the solute atom radius for the soft-sphere solvation model given species
	double cavityScale; //!< radius scale factor
	double ionSpacing; //!< extra spacing from dielectric to ionic cavity (in bohrs); also used by CANON
	
	//For fixed-cavity model alone:
	string cavityFile; //!< filename of cavity to read in
	
	//Cavity masking parameters (erf slab centered at zMask0 with half-width zMaskH where fluid is excluded):
	double zMask0; //z center in lattice coordinates for cavity mask
	double zMaskH; //half-width in z lattice coordinates for cavity mask
	double zMaskIonH; //half-width in z lattice coordinates for ionic cavity mask
	double zMaskSigma; //smoothness of z-mask in bohrs
	
	//Debug parameters for Nonlinear PCM's:
	bool linearDielectric; //!< If true, work in the linear dielectric response limit
	bool linearScreening; //!< If true, work in the linearized Poisson-Boltzman limit for the ions
	bool nonlinearSCF; //!< whether to use an SCF method for nonlinear PCMs
	double screenOverride; //! overrides screening factor with this value
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

//! @}
#endif // JDFTX_ELECTRONIC_FLUIDSOLVERPARAMS_H
