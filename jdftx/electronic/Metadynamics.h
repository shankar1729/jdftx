/*-------------------------------------------------------------------
Copyright 2025 Andrew Diggs, Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_METADYNAMICS_H
#define JDFTX_ELECTRONIC_METADYNAMICS_H

#include <core/Coulomb.h>
#include <core/GridInfo.h>

//! Metadynamics on a bond between two atoms
struct MetadynamicsBond
{
	int iSp1, at1; //!< species-name and 0-based index of first atom
	int iSp2, at2; //!< species-name and 0-based index of second atom
	double resolution; //!< spatial resolution of bias potential in bohrs
	double energy_per_step; //!< integral of energy in Eh-a0 added to bias potential profile each time step
	string state_filename; //!< file name to load/save bias potential from
	diagMatrix potential; //!< bias potential accumulated during metadynamics (in Hartrees)
	
	MetadynamicsBond() : iSp1(-1), at1(-1), iSp2(-1), at2(-1), resolution(0.), energy_per_step(0.) {}
	
	void initialize(); //!< initialize potential, including load from file if specified
	void save(); //!< save current bias potential to file if specified
	
	//Update bias potential, accumulate corresponding forces and return corresponding potential energy
	double energyAndGrad(const GridInfo& gInfo, std::vector<Atom>& atoms);
};

#endif // JDFTX_ELECTRONIC_METADYNAMICS_H

