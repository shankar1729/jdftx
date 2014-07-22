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
#include <electronic/ElecInfo.h>
#include <core/vector3.h>

//! (Weighted-) density of states calculator
class DOS
{
public:
	
	
	//! Weight-function definition
	struct Weight
	{	
		//!Treatment of band/state fillings
		enum FillingMode
		{	Complete, //!< Do not weight density of states by fillings
			Occupied, //!< Weight density of states by fillings
		};
		FillingMode fillingMode; //!< treatment of band/state fillings
		
		//!Treatment of spin in non-collinear calculations:
		vector3<> Mhat; //!< spin polarization direction (unit vector if active, zero if inactive)
		
		//!Weight-function type
		enum Type
		{	Total, //!< total density of states in unit cell
			Slice, //!< density of states within an arbitrary planar slice
			Sphere, //!< density of states inside an arbitrary sphere
			AtomSlice, //!< density of states in a planar slice centered on an atom
			AtomSphere, //!< density of states in a sphere centered on an atom
			File, //!< density of states with an arbitrary weight function read from a file
			Orbital, //!< atomic-orbital projected density of states
			OrthoOrbital //!< Lowdin orthonormalized-atomic-orbital projected density of states
		};
		Type type; //!< weight function type
		
		vector3<int> direction; //!< lattice plane specification for slice mode
		
		vector3<> center; //!< center of slice/sphere in lattice coordinates
		double radius; //!< radius for sphere modes or half-width for slice modes (bohrs)
		
		size_t specieIndex; //!< Specie index for atom-centered modes
		size_t atomIndex; //!< Atom index for atom-centered modes
		
		string filename; //!< Weight-function input filename for File mode
		
		struct OrbitalDesc
		{	int l, m, s; unsigned n; //!< pseudo-atom quantum numbers (m = l+1 is used to signify total l-contribution)
			SpinType spinType; //!< SpinOrbit => relativistic spec, SpinZ => Up/Dn, SpinNone => no spin (collinear) (Note different meanings from eInfo.spinType)
			OrbitalDesc(int l=0, int m=0, int s=0, SpinType spinType=SpinNone) : l(l), m(m), s(s), n(0), spinType(spinType) {}
			void parse(string desc); //!< set values from a string (throws a string exception on invalid input)
			operator string() const; //!< convert to a string description
			bool operator<(const OrbitalDesc&) const; //!< so that OrbitalDesc can be used as a key
		};
		OrbitalDesc orbitalDesc;
		
		string getDescription(const Everything&) const; //!< return a descriptive string
	};
	
	std::vector<Weight> weights; //!< list of weight functions (default: total DOS only)
	double Etol; //!< tolerance for identifying eigenvalues (energy resolution) (default: 1e-6)
	double Esigma; //!< optional gaussian width in spectrum
	
	DOS();
	void setup(const Everything&); //!< initialize
	void dump(); //!< dump density of states to file (filename obtained from Dump)
private:
	const Everything* e;
};

#endif // JDFTX_ELECTRONIC_DOS_H