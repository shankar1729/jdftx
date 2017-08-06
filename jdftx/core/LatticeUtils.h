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

//! @addtogroup Geometry
//! @{

//! @file LatticeUtils.h Miscellaneous utilities relating to properties of Bravais Lattices

#include <core/GridInfo.h>
#include <core/Util.h>

//! Relative threshold for symmetry detection
extern double symmThreshold, symmThresholdSq;

//! Given a set of lattice vectors in the columns of R,
//! return the minimal lattice vectors (shortest linear combinations)
//! and optionally retrieve the integer transmission matrix
//! (Rreduced = R * transmission) and its inverse.
matrix3<> reduceLatticeVectors(const matrix3<>& R,
	matrix3<int>* transmission=0, matrix3<int>* invTransmission=0);

//! Find the symmetries of a Bravais lattice, where some
//! of the lattice directions may optionally be truncated.
//! Optionally retrieve the reduced lattice vectors
//! and transmission matrices which are computed as
//! an intermediate step in determining symmetries.
std::vector<matrix3<int>> getSymmetries(const matrix3<>& R,
	vector3<bool> isTruncated=vector3<bool>(false,false,false),
	matrix3<>* Rreduced=0, matrix3<int>* transmission=0, matrix3<int>* invTransmission=0);

//! Supercell corresponding to a given k-point mesh
struct Supercell
{
	const GridInfo& gInfo; //!< unit cell and corresponding grids
	std::vector<vector3<>> kmesh; //!< closure of kmeshReduced under symmetry group sym
	matrix3<> Rsuper; //!< super-cell lattice vectors
	matrix3<int> super; //!< linear combinations to get Rsuper (Rsuper = R * super)
	
	//! Transformation from reduced k-point mesh to full one. The resulting k-point
	//! in the full mesh is (~sym[iSym]) * kmeshReduced[iReduced] * invert + offset
	struct KmeshTransform
	{	unsigned iReduced; //!< corresponding reduced index
		unsigned iSym; //!< symmetry matrix for transforming reduced to current value
		int invert; //!< sign of transformation to include inversion symmetry in k-space
		vector3<int> offset; //!< additional translation to get to kmesh from reduced one
	};
	std::vector<KmeshTransform> kmeshTransform;
	
	//! Construct supercell given the unit cell and grid definition,
	//! symmetry-reduced k-point mesh and list of symmetries
	Supercell(const GridInfo& gInfo,
		const std::vector<vector3<>>& kmeshReduced,
		const std::vector<SpaceGroupOp>& sym, const std::vector<int>& invertList);
};

//! Get a list of unit cells in a supercell, with padding at the boundaries to maintain a Wigner-Seitz
//! supercell range (smoothed by rSmooth) for all pairs of lattice coordinates between arrays x1 and x2.
//! Returns a list of cell positions in lattice coordinates, along with the weights for matrix elements
//! connecting each pair from x1 and x2 (which will add to 1 over each set of equivalent cells)
//! Optionally write the cell map to a file, if fname is non-null
std::map<vector3<int>, class matrix> getCellMap(const matrix3<>& R, const matrix3<>& Rsup, const vector3<bool>& isTruncated,
	const std::vector<vector3<>>& x1, const std::vector<vector3<>>& x2, double rSmooth, string fname=string());

//Helper function for PeriodicLookup< vector3<> > used in Supercell::Supercell
inline vector3<> getCoord(const vector3<>& pos) { return pos; }

//! O(1) look-up table for finding periodic image within symmThreshold
//! Needs a companion function vector3<> getCoord(T) that returns the lattice coordinates corresponding to type T
template<typename T> class PeriodicLookup
{	const std::vector<T>& points;
	vector3<int> S; //lookup mesh sample count
	std::vector< std::vector<size_t> > indices; //list of indices into points, for each  lookup mesh cell
	
	inline size_t meshIndex(vector3<int> iv) const
	{	for(int k=0; k<3; k++) //wrap to [0,S)
		{	if(iv[k] < 0) iv[k] += S[k];
			if(iv[k] >= S[k]) iv[k] -= S[k];
		}
		return iv[0]+S[0]*size_t(iv[1]+S[1]*iv[2]);
	}

public:
	//!Initialize given an array of points and a metric corresponding to the lattice coordinates
	//!(optionally override the number of points, if the points array is to be dynamically grown)
	PeriodicLookup(const std::vector<T>& points, matrix3<> metric, size_t nPointsTarget=0) : points(points)
	{	//Set S such that prod(S) ~ nPoints and the sample counts are proportional to dimension length
		vector3<> Stmp; for(int k=0; k<3; k++) Stmp[k] = sqrt(metric(k,k));
		Stmp *= pow(std::max(points.size(), nPointsTarget)/(Stmp[0]*Stmp[1]*Stmp[2]), 1./3); //normalize so that product is nPoints
		for(int k=0; k<3; k++)
		{	S[k] = std::max(1, int(round(Stmp[k])));
			assert(symmThreshold*S[k] < 0.5);
		}
		//Initialize indices:
		indices.resize(S[0]*S[1]*S[2]);
		for(size_t iPoint=0; iPoint<points.size(); iPoint++)
			addPoint(iPoint, points[iPoint]);
	}
	
	//! Call every time a point is added to the end of the points vector externally.
	//! Note that points may only be added to the end of the vector (otherwise some of the existing indices are invalidated!)
	void addPoint(size_t iPoint, const T& point)
	{	vector3<> v = getCoord(point);
		vector3<int> iv;
		for(int k=0; k<3; k++)
		{	v[k] -= floor(v[k]); //in [0,1)
			iv[k] = int(floor(v[k]*S[k] + 0.5)); //in [0, S]
		}
		indices[meshIndex(iv)].push_back(iPoint);
	}
	
	//!Return index of a point within symmThresold of x, and string::npos if none found
	//!(Optionally also only return points whose tag matches the specified value)
	template<typename Tag = double, typename TagEquality = std::equal_to<Tag> >
	size_t find(vector3<> v, Tag tag = Tag(), const std::vector<Tag>* tagArr=0, TagEquality tagEquality = std::equal_to<Tag>()) const
	{	vector3<int> ivMin, ivMax;
		for(int k=0; k<3; k++)
		{	v[k] -= floor(v[k]); //in [0,1)
			ivMin[k] = int(floor((v[k]-symmThreshold)*S[k] + 0.5)); //in [-1, S]
			ivMax[k] = int(floor((v[k]+symmThreshold)*S[k] + 0.5)); //in [0, S]
		}
		vector3<int> iv;
		for(iv[0]=ivMin[0]; iv[0]<=ivMax[0]; iv[0]++)
		for(iv[1]=ivMin[1]; iv[1]<=ivMax[1]; iv[1]++)
		for(iv[2]=ivMin[2]; iv[2]<=ivMax[2]; iv[2]++)
			for(size_t index: indices[meshIndex(iv)])
				if(circDistanceSquared(v, getCoord(points[index])) < symmThresholdSq
					&& (!tagArr || tagEquality(tag, tagArr->at(index))) )
					return index;
		return string::npos;
	}
};

//! @}
#endif // JDFTX_CORE_LATTICEUTILS_H
