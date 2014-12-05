/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman
Copyright 1996-2003 Sohrab Ismail-Beigi

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

#ifndef JDFTX_ELECTRONIC_BASIS_H
#define JDFTX_ELECTRONIC_BASIS_H

#include <electronic/common.h>
#include <core/vector3.h>
#include <core/GridInfo.h>
#include <vector>

class Basis
{
public:
	const GridInfo* gInfo; //!< pointer to the grid specs
	const IonInfo* iInfo; //!< pointer to the ion information (basis is conceptually ultrasoft-pseudopotential dependent)
	
	size_t nbasis; //!< number of basis elements (i.e. G-vectors)
	vector3<int> *iGarr; //!< the (integer) G-vectors for the basis in recip. lattice coords
	int *index; //!< indices of the basis functions in the FFT boxes used
	#ifdef GPU_ENABLED
	vector3<int> *iGarrGpu; //!< GPU copy of G-vector index coefficients
	int *indexGpu; //!< copy of index array on the GPU
	#endif
	int* indexPref; //!< points to indexGpu in GPU mode and index otherwise
	vector3<int> *iGarrPref; //!< points to iGarrGpu in GPU mode and iGarr otherwise
	
	std::vector<int> head; //!< short list of low G basis locations (used for phase fixing)
	
	Basis();
	~Basis();

	Basis(const Basis&); //!< copy by reference
	Basis& operator=(const Basis&); //!< copy by reference

	//! Setup the indices and integer G-vectors within Ecut for kpoint k
	void setup(const GridInfo& gInfo, const IonInfo& iInfo, double Ecut, const vector3<> k);

	//! Create a custom basis with an arbitrary indexing scheme
	void setup(const GridInfo& gInfo, const IonInfo& iInfo, const std::vector<int>& indexVec);
	
private:
	bool ownsData; //whether the pointers were allocated by this object and need to be freed
	void setup(const GridInfo& gInfo, const IonInfo& iInfo,
		const std::vector<int>& indexVec,
		const std::vector< vector3<int> >& iGvec); //set the data arrays from vectors
};

#endif // JDFTX_ELECTRONIC_BASIS_H
