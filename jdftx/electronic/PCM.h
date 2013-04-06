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
#include <core/EnergyComponents.h>

//! Base class for all PCMs
class PCM : public FluidSolver
{
public:
	PCM(const Everything& e, const FluidSolverParams& fsp);

	void dumpDebug(const char* filenamePattern) const; //!< generate fluidDebug text file with common info to all PCMs

protected:
	FluidSolverParams fsp;
	double k2factor; //!< Prefactor to kappaSq
	
	EnergyComponents Adiel; //!< PCM energy components
	DataGptr rhoExplicitTilde; //!< Charge density of explicit (electronic) system
	DataRptr nCavity; //!< Cavity determining electron density (includes core and chargeball contributions)
	DataRptr shape, shapeEx; //!< Cavity shape function (and expanded shape function for the SGA13 variant)

	virtual void printDebug(FILE* fp) const {} //!< over-ride to get extra PCM-specific output in fluidDebug text file
	
	void updateCavity(); //!< update shape function(s) from nCavity, and energies dependent upon shape alone
	void propagateCavityGradients(const DataRptr& A_shape, DataRptr& A_nCavity) const; //!< propagate gradients w.r.t shape function (including cached cavity gradients) to those w.r.t nCavity
private:
	DataRptr Acavity_shape, Acavity_shapeEx; //!< Cached gradients of cavitation (and dispersion) energies w.r.t shape functions
	double A_nc, A_tension, A_vdwScale; //!< Cached deerivatives w.r.t fit parameters (accessed via dumpDebug() for PCM fits)
};

#endif // JDFTX_ELECTRONIC_PCM_H