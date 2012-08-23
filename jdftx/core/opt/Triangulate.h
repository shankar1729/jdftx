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

#ifndef JDFTX_CORE_OPT_TRIANGULATE_H
#define JDFTX_CORE_OPT_TRIANGULATE_H

//! @file Triangulate.h Wrapper for the periodic 3D triangulation routines in the CGAL library

#include <core/vector3.h>
#include <vector>

namespace Triangulation
{
	//! Vertex in triangulation
	struct Vertex
	{	vector3<> pos; //!< position of the vertex in lattice coordinates (including offset)
		size_t index; //!< index of corresponding point in the fundamental domain in the input array
	};

	//! Tetrahedral cell in triangulation
	struct Cell
	{	Vertex vertex[4];
	};

	/**
	Compute a 3D periodic Delaunay triangulation
	@param points A list of points in [-0.5,0.5)^3 (lattice coordinates)
	@param cells List of tetrahedra in the computed triangulation
	*/
	void triangulate(const std::vector< vector3<> >& points, std::vector<Cell>& cells);
}

#endif // JDFTX_CORE_OPT_TRIANGULATE_H
