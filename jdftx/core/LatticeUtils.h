/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

#ifndef JDFTX_CORE_LATTICEUTILS_H
#define JDFTX_CORE_LATTICEUTILS_H

//! @file GridInfo.h Miscellaneous utilities relating to properties of Bravais Lattices

#include <core/GridInfo.h>
#include <vector>

//! Relative threshold for symmetry detection
static const double symmThreshold = 1e-4;
static const double symmThresholdSq = symmThreshold * symmThreshold;

//! Given a set of lattice vectors in the columns of R,
//! return the minimal lattice vectors (shortest linear combinations)
//! and optionally retrieve the integer transmission matrix
//! (Rreduced = R * transmission) and its inverse.
matrix3<> reduceLatticeVectors(const matrix3<>& R,
	matrix3<int>* transmission=0, matrix3<int>* invTransmission=0);

//! Find the symmetries of a Bravais lattice
//! Optionally retrieve the reduced lattice vectors
//! and transmission matrices which are computed as
//! an intermediate step in determining symmetries.
std::vector<matrix3<int>> getSymmetries(const matrix3<>& R, matrix3<>* Rreduced=0,
	matrix3<int>* transmission=0, matrix3<int>* invTransmission=0);

//! Supercell corresponding to a given k-point mesh
struct Supercell
{
	const GridInfo& gInfo; //!< unit cell and corresponding grids
	std::vector<vector3<>> kmesh; //!< closure of kmeshReduced under symmetry group sym
	matrix3<> Rsuper; //!< super-cell lattice vectors
	matrix3<int> super; //!< linear combinations to get Rsuper (Rsuper = R * super)
	
	//! Construct supercell given the unit cell and grid definition,
	//! symmetry-reduced k-point mesh and list of symmetries
	Supercell(const GridInfo& gInfo,
		const std::vector<vector3<>>& kmeshReduced,
		const std::vector<matrix3<int>>& sym);
};

#endif // JDFTX_CORE_LATTICEUTILS_H
