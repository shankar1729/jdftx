/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman, Deniz Gunceler

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

#ifndef JDFTX_ELECTRONIC_PCM_H
#define JDFTX_ELECTRONIC_PCM_H

#include <fluid/FluidSolver.h>
#include <core/RadialFunction.h>
#include <core/EnergyComponents.h>
#include <core/Coulomb.h>

//! @addtogroup Solvation
//! @{

//! Base class for all PCMs
class PCM : public FluidSolver
{
public:
	PCM(const Everything& e, const FluidSolverParams& fsp);
	virtual ~PCM();
	
	void dumpDensities(const char* filenamePattern) const; //!< dump cavity shape functions
	void dumpDebug(const char* filenamePattern) const; //!< generate fluidDebug text file with common info to all PCMs

protected:
	EnergyComponents Adiel; //!< PCM energy components
	ScalarFieldTilde rhoExplicitTilde; //!< Charge density of explicit (electronic) system
	ScalarField nCavity, tauCavity, nCavityEx[2]; //!< Cavity determining electron density (or product for SaLSA, or KE density for SG14tauVW, and expanded electron densities for the SGA13 variant)
	ScalarFieldArray shape; //!< Electrostatic cavity shape function. Second component, if any, is separate ionic cavity
	ScalarField shapeVdw; //!< Separate cavitation/dispersion shape function for the SGA13 variant
	
	virtual void printDebug(FILE* fp) const {} //!< over-ride to get extra PCM-specific output in fluidDebug text file
	
	void updateCavity(); //!< update shape function(s) from nCavity, and energies dependent upon shape alone
	
	//! Propagate A_shape (+ cached Acavity_shape) and accumulate to gradients w.r.t nCavity and rhoExplicitTilde.
	//! Set fluid force contributions (for atom-sphere cavities) if non-null.
	void propagateCavityGradients(const ScalarFieldArray& A_shape, ScalarField& A_nCavity, ScalarFieldTilde& A_rhoExplicitTilde, IonicGradient* forces) const;
	
	//! Accumulate extra fluid forces (vdw and full-core forces, when applicable)
	void accumExtraForces(IonicGradient* forces, const ScalarFieldTilde& A_nCavityTilde) const;
	
	ScalarFieldTilde getFullCore() const; //!< get full core correction for PCM variants that need them
private:
	ScalarField Acavity_shape, Acavity_shapeVdw; //!< Cached gradients of cavitation (and dispersion) energies w.r.t shape functions (assumed Acavity does not depend on ionic cavity)
	double A_nc, A_tension, A_vdwScale, A_eta_wDiel, A_pCavity, A_cavityScale; //!< Cached derivatives w.r.t fit parameters (accessed via dumpDebug() for PCM fits)
	double Rex[2]; //!< radii for cavity expansion (SGA13 only)
	RadialFunctionG wExpand[2]; //!< weight function for cavity expansion (SGA13 only)
	RadialFunctionG wCavity; //!< weight function for nonlocal cavitation energy
	std::vector<vector3<>> atposAll; //!< all solute atomic positions (SoftSphere only)
	std::vector<vector3<int>> latticeReps; //!< lattice repetitions needed to treat periodic boundaries correctly for SoftSphere
	std::vector<double> Rall; //!< all solute atomic radii, with overall scale factor pre-multiplied (SoftSphere only)
	std::vector<double> RallIonic; //!< Rall + ion spacing adde, used for separate ionic cavity, if any (SoftSphere only)
	ScalarFieldArray zMask; //optional cavity mask function
	int nShape; //natural number of shape functions of the solvation model (2 if ionspacing is used to make ionic cavity, else 1)
	bool fixedCavityMasked; //!< whether mask has already been applied to fixed cavity
protected:
	std::vector<RadialFunctionG> Sf; //!< spherically-averaged structure factors for each solvent site
	std::vector<int> atomicNumbers; //!< atomic number for each solvent site (for dispersion interactions)
	static ScalarFieldTilde coulomb(const ScalarFieldTilde& rho) { return (-4*M_PI) * Linv(O(rho)); }
	friend struct ChargedDefect;
};

//! @}
#endif // JDFTX_ELECTRONIC_PCM_H
