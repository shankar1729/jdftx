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

#include <core/opt/Triangulate.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <map>
#include <core/string.h>

namespace Triangulation
{
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef CGAL::Periodic_3_triangulation_traits_3<K> GT;
	typedef CGAL::Periodic_3_Delaunay_triangulation_3<GT> PDT;

	void triangulate(const std::vector< vector3<> >& points, std::vector<Cell>& cells)
	{	//Domain is a unit cube centered at origin: 
		PDT::Iso_cuboid domain(-0.5,-0.5,-0.5, 0.5,0.5,0.5);
		//Convert the points to the CGAL format and maintain an index map to the original set:
		std::vector<PDT::Point> cgalPoints(points.size());
		std::map< vector3<>, size_t > indexMap;
		for(size_t i=0; i<points.size(); i++)
		{	cgalPoints[i] = PDT::Point(points[i][0],points[i][1], points[i][2]);
			indexMap[points[i]] = i;
		}
		//Create the triangulation:
		PDT T(cgalPoints.begin(), cgalPoints.end(), domain);

		//Loop over the tetrahedra:
		for(PDT::Periodic_tetrahedron_iterator t=T.periodic_tetrahedra_begin(PDT::UNIQUE); t!=T.periodic_tetrahedra_end(PDT::UNIQUE); t++)
		{	Cell cell;
			for(int j=0; j<4; j++)
			{	PDT::Point origPoint = (*t)[j].first; //in orig domain
				PDT::Point netPoint = T.point((*t)[j]); //with offset
				vector3<> orig(origPoint[0], origPoint[1], origPoint[2]);
				vector3<> net(netPoint[0], netPoint[1], netPoint[2]);
				//Locate point in original set:
				std::map< vector3<>, size_t >::iterator iter = indexMap.find(orig);
				if(iter == indexMap.end())
				{	ostringstream oss;
					oss << "Error in triangulation: point (" << orig[0] << "," << orig[1] << "," <<  orig[2] << ") not in original set.\n";
					throw oss.str();
				}
				cell.vertex[j].pos = net;
				cell.vertex[j].index = iter->second;
			}
			cells.push_back(cell);
		}
	}
}
