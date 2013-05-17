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

#include <electronic/FluidSolver.h>
#include <electronic/RadialFunction.h>
#include <core/EnergyComponents.h>

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
	DataGptr rhoExplicitTilde; //!< Charge density of explicit (electronic) system
	DataRptr nCavity, nCavityEx[2]; //!< Cavity determining electron density (or product for nonlocalPCM, and expanded electron densities for the SGA13 variant)
	DataRptr shape, shapeVdw; //!< Electrostatic cavity shape function (and separate cavitation/dispersion shape function for the SGA13 variant)
	std::shared_ptr<IonicGradient> vdwForces; //!< cached forces due to dispersion terms (if any)
	
	virtual void printDebug(FILE* fp) const {} //!< over-ride to get extra PCM-specific output in fluidDebug text file
	
	void updateCavity(); //!< update shape function(s) from nCavity, and energies dependent upon shape alone
	void propagateCavityGradients(const DataRptr& A_shape, DataRptr& A_nCavity) const; //!< propagate A_shape (+ cached Acavity_shape) to those A_nCavity (set result, not accumulate)
private:
	DataRptr Acavity_shape, Acavity_shapeVdw; //!< Cached gradients of cavitation (and dispersion) energies w.r.t shape functions
	double A_nc, A_tension, A_vdwScale; //!< Cached derivatives w.r.t fit parameters (accessed via dumpDebug() for PCM fits)
	double Rex[2]; //!< radii for cavity expansion (SGA13 only)
	RadialFunctionG wExpand[2]; //!< weight function for cavity expansion (SGA13 only)
	RadialFunctionG wCavity; //!< weight function for nonlocal cavitation energy
protected:
	std::vector<RadialFunctionG> Sf; //!< spherically-averaged structure factors for each solvent site
};

#endif // JDFTX_ELECTRONIC_PCM_H
