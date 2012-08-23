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

#ifndef JDFTX_ELECTRONIC_DOS_H
#define JDFTX_ELECTRONIC_DOS_H

#include <electronic/common.h>
#include <core/vector3.h>

//! (Weighted-) density of states calculator
class DOS
{
public:
	
	//! Weight-function definition
	struct Weight
	{	//!Weight-function type
		enum Type
		{	Total, //!< total density of states in unit cell
			Slice, //!< density of states within an arbitrary planar slice
			Sphere, //!< density of states inside an arbitrary sphere
			AtomSlice, //!< density of states in a planar slice centered on an atom
			AtomSphere, //!< density of states in a sphere centered on an atom
			File, //!< density of states with an arbitrary weight function read from a file
			Delim //!< delimiter used by command parser (not a valid weight type)
		};
		
		Type type; //!< weight function type
		vector3<int> direction; //!< lattice plane specification for slice mode
		
		vector3<> center; //!< center of slice/sphere in lattice coordinates
		double radius; //!< radius for sphere modes or half-width for slice modes (bohrs)
		
		size_t specieIndex; //!< Specie index for atom-centered modes
		size_t atomIndex; //!< Atom index for atom-centered modes
		
		string filename; //!< Weight-function filename
	};
	
	std::vector<Weight> weights; //!< list of weight functions (default: total DOS only)
	bool occupied; //!< whether to weight the density of states by fillings (default: no)
	double Etol; //!< tolerance for identifying eigenvalues (energy resolution) (default: 1e-6)
	
	DOS();
	void setup(const Everything&); //!< initialize
	void dump(); //!< dump density of states to file (filename obtained from Dump)
private:
	const Everything* e;
};

#endif // JDFTX_ELECTRONIC_DOS_H