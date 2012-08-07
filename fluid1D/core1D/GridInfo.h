/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

This file is part of Fluid1D.

Fluid1D is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Fluid1D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Fluid1D.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#ifndef FLUID1D_CORE1D_GRIDINFO_H
#define FLUID1D_CORE1D_GRIDINFO_H

/** @file GridInfo.h
@brief Geometry of the simulation grid
*/

#include <fftw3.h>
#include <vector>

/** @brief Simulation grid descriptor
//! @ingroup griddata
*/
class GridInfo
{
public:
	
	const enum CoordinateSystem
	{	Spherical, //!< 1D spherical, `even extension' about r = rMax
		Cylindrical, //!< 1D cylindrical (all quantities per bohr length in z), `even extension' about rho = rMax
		Planar //!< 1D planar (all quantities per bohr^2 area in xy), even extension about both ends in z (0 and rMax)
	} coord; //!< Coordinate system
	const int S; //!< Sample count
	const double rMax; //!< Length or maximum radius of simulation grid
	
	//! Setup simulation grid
	//! @param coord Coordinate system
	//! @param S Sample count (equal to basis function count for all implemented bases)
	//! @param hMean Mean grid spacing, defined by rMax/S
	GridInfo(CoordinateSystem coord, int S, double hMean);
	~GridInfo();
	
	std::vector<double> r; //!< Nodes of quadrature grid
	std::vector<double> G; //!< Momenta of basis functions
	std::vector<double> w; //!< Weights on quadrature grid
	std::vector<double> wTilde; //!< Normalization of basis functions

	fftw_plan planPlanarI, planPlanarIdag, planPlanarID, planPlanarIDdag; //!< FFTW plans for planar transforms
	std::vector<double> matI, matID, matIDD; //!< Dense SxS row-major matrices for spherical/cylindrical transforms
};

#endif // FLUID1D_CORE1D_DATA_H
