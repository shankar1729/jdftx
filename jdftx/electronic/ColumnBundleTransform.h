/*-------------------------------------------------------------------
Copyright 2014 Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_COLUMNBUNDLETRANSFORM_H
#define JDFTX_ELECTRONIC_COLUMNBUNDLETRANSFORM_H

#include <electronic/Basis.h>
#include <electronic/matrix.h>

//! Handle transformation of ColumnBundles upon symmetry operations
class ColumnBundleTransform
{
public:
	struct BasisWrapper
	{	BasisWrapper(const Basis& basis);
		const Basis& basis; //!< reference to underlying basis
		vector3<int> iGbox; //!< half-size of index look-up table
		vector3<int> pitch; //!< pitch of each dimension in index look-up table
		std::vector<int> table; //!< index look-up table
	};
	
	/**
	Define transformations between basis/k 'C' indexed continuously and basis/k 'D' indexed
	discontinuously i.e. scatter will occur from C to D, whereas gather will occur from D to C.
	The D basis is wrapped so that a lookup table generated during construction
	can be reused if many transformations to a common D are required.
	The symmetry operation sym is in covariant lattice coordinates,
	and explicit k-point inversion is specified by invert = +/- 1.
	Optionally, super specifies a transform to a supercell (D corresponds to a supercell of C);
	columnbundles will however not be automatically re-normalized for the supercell
	Note that the transformation is kD = kC * sym * invert * super + offset (where offset is determined automatically).
	If the ColumnBundles are spinorial (nSpinor=2), corresponding spin-space transformations will also be applied
	*/
	ColumnBundleTransform(const vector3<>& kC, const Basis& basisC, const vector3<>& kD, const BasisWrapper& basisDwrapper,
		int nSpinor, const matrix3<int>& sym, int invert, const matrix3<int>& super = matrix3<int>(1,1,1));
	~ColumnBundleTransform();
	
	//Non-copyable:
	ColumnBundleTransform(const ColumnBundleTransform&)=delete;
	ColumnBundleTransform& operator=(const ColumnBundleTransform&)=delete;
	
	void scatterAxpy(complex alpha, const ColumnBundle& C_C, int bC, ColumnBundle& C_D, int bD) const; //!< scatter-accumulate a single column
	void gatherAxpy(complex alpha, const ColumnBundle& C_D, int bD, ColumnBundle& C_C, int bC) const; //!< gather-accumulate a single column
	
	void scatterAxpy(complex alpha, const ColumnBundle& C_C, ColumnBundle& C_D, int bDstart, int bDstep) const; //!< scatter-accumulate all columns of C_C
	void gatherAxpy(complex alpha, const ColumnBundle& C_D, int bDstart, int bDstep, ColumnBundle& C_C) const; //!< gather-accumulate all columns of C_C

private:
	const vector3<>& kC; const Basis& basisC;
	const vector3<>& kD; const Basis& basisD;
	int nSpinor, invert;
	
	//Index array:
	std::vector<int> index; int* indexPref;
	#ifdef GPU_ENABLED
	int* indexGpu;
	#endif

	matrix spinorRot; //spinor space rotation
};

#endif //JDFTX_ELECTRONIC_COLUMNBUNDLETRANSFORM_H