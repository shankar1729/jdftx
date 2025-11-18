/*-------------------------------------------------------------------
Copyright 2022 Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_IONIC_GAUSSIAN_POTENTIAL_H
#define JDFTX_ELECTRONIC_IONIC_GAUSSIAN_POTENTIAL_H

#include <core/Coulomb.h>
#include <core/GridInfo.h>

//! Gaussian potential and forces on atoms, centered at origin
struct IonicGaussianPotential
{
	int iSpecies; //!< which species it applies to
	double U0; //!< peak amplitude in Hartrees
	double sigma; //!< width (standard deviation) in bohrs
	enum Geometry {Spherical, Cylindrical, Planar} geometry;
	vector3<> center; //!< center of the Gaussian potential in reduced coordinates (0.5,0.5,0.5) is center of cell

	IonicGaussianPotential() : iSpecies(-1), U0(0.), sigma(0.), geometry(Spherical), center(0.0, 0.0, 0.0) {}

	//Compute energy and accumulate forces, and optionally stresses, due to external potential:
	double energyAndGrad(const GridInfo& gInfo, std::vector<Atom>& atoms, matrix3<>* E_RRTptr) const;
};

#endif // JDFTX_ELECTRONIC_IONIC_GAUSSIAN_POTENTIAL_H
