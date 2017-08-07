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

#include <core/ManagedMemory.h>

class GridInfo;
class IonInfo;

//! @addtogroup ElecSystem
//! @{

//! Wavefunction basis
class Basis
{
public:
	const GridInfo* gInfo; //!< pointer to the grid specs
	const IonInfo* iInfo; //!< pointer to the ion information (basis is conceptually ultrasoft-pseudopotential dependent)
	
	size_t nbasis; //!< number of basis elements (i.e. G-vectors)
	IndexVecArray iGarr;
	IndexArray index;
	std::vector<int> head; //!< short list of low G basis locations (used for phase fixing)
	
	Basis();
	Basis(const Basis&); //!< copy by reference
	Basis& operator=(const Basis&); //!< copy by reference

	//! Setup the indices and integer G-vectors within Ecut for kpoint k
	void setup(const GridInfo& gInfo, const IonInfo& iInfo, double Ecut, const vector3<> k);

	//! Create a custom basis with an arbitrary indexing scheme
	void setup(const GridInfo& gInfo, const IonInfo& iInfo, const std::vector<int>& indexVec);
	
private:
	void setup(const GridInfo& gInfo, const IonInfo& iInfo,
		const std::vector<int>& indexVec,
		const std::vector< vector3<int> >& iGvec); //set the data arrays from vectors
};

//! @}
#endif // JDFTX_ELECTRONIC_BASIS_H
