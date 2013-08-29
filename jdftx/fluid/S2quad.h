/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

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

#ifndef JDFTX_FLUID_SO3PLATONIC_H
#define JDFTX_FLUID_SO3PLATONIC_H

//! @addtogroup so3quad
//! @{

/** @file S2quad.h
@brief Quadratures on S2 used to generate quadratures on SO(3)
*/

#include <vector>
#include <core/vector3.h>
#include <core/string.h>
#include <core/Util.h>

//! Abstract base class for a S2 quadrature definition (used to generate the SO3 qudarature)
struct S2quad
{
	std::vector<vector3<> > euler; //!< list of S2 quadrature nodes in first two components with suggested S1 phase in third component
	std::vector<double> weight; //!< suggested S2 quadrature weights (same length as euler, will be normalized to add up to 1)
	virtual int jMax() const=0; //!< max angular momentum that the S2 quadrature is exact to (will be checked)
	virtual int nS1() const=0; //!< suggested number of samples for the S1 sector (will be rounded up to nearest multiple of Zn)
	virtual string name() const=0; //!< A short descriptive name for the S2 quadrature
	virtual ~S2quad() {}
protected:
	//! Utility to add nodes, optionally modifying weights/S1 phases
	//! @param n direction vector for node (need not be normalized)
	//! @param relWeight relative weight (weights will be normalized to add up to 1 automatically), uniform weights by default
	//! @param s1phase Phase for the S1 part in going to SO(3) quadrature, 0 by default.
	void add(vector3<> n, double relWeight=1., double s1phase=0.);
};

//! List of available quadratures
enum S2quadType
{	QuadEuler,
	QuadTetrahedron,
	QuadOctahedron,
	QuadIcosahedron,
	Quad7design_24,
	Quad8design_36,
	Quad9design_48,
	Quad10design_60,
	Quad11design_70,
	Quad12design_84,
	Quad13design_94,
	Quad14design_108,
	Quad15design_120,
	Quad16design_144,
	Quad17design_156,
	Quad18design_180,
	Quad19design_204,
	Quad20design_216,
	Quad21design_240,
};

//! Names for available quadratures (used in input files etc.)
static EnumStringMap<S2quadType> S2quadTypeMap
(	QuadEuler,       "Euler",
	QuadTetrahedron, "Tetrahedron", 
	QuadOctahedron,  "Octahedron",
	QuadIcosahedron, "Icosahedron",
	Quad7design_24,   "7design24",
	Quad8design_36,   "8design36",
	Quad9design_48,   "9design48",
	Quad10design_60,  "10design60",
	Quad11design_70,  "11design70",
	Quad12design_84,  "12design84",
	Quad13design_94,  "13design94",
	Quad14design_108, "14design108",
	Quad15design_120, "15design120",
	Quad16design_144, "16design144",
	Quad17design_156, "17design156",
	Quad18design_180, "18design180",
	Quad19design_204, "19design204",
	Quad20design_216, "20design216",
	Quad21design_240, "21design240"
);


//! @brief Outer-product quadrature on ZYZ euler angles
class EulerProduct : public S2quad
{
public:
	//! Specify sampling along each auler angle
	//! @param nBeta Number of Gauss-legendre samples in cos(beta) on [-1,1] (must be non-zero)
	//! @param nAlpha Number of uniform samples in alpha on [0,2pi) (taken to be 2*nBeta if 0)
	//! @param nGamma Number of uniform samples in gamma on [0,2pi) (taken to be 2*nBeta if 0)
	//! Note nGamma will be reduced by Zn symmetry in SO(3) quadrature generation
	EulerProduct(unsigned nBeta, unsigned nAlpha=0, unsigned nGamma=0);
	int jMax() const;
	int nS1() const;
	string name() const;
	
private:
	unsigned nBeta, nAlpha, nGamma;
};


//--------------- Platonic Solid rotation groups ------------------

//! @brief Tetrahedron rotation group
class Tetrahedron : public S2quad
{
public:
	Tetrahedron();
	int jMax() const { return 2; }
	int nS1() const { return 3; }
	string name() const { return "Tetrahedron"; }
};

//! @brief Octahedron rotation group
class Octahedron : public S2quad
{
public:
	Octahedron();
	int jMax() const { return 3; }
	int nS1() const { return 4; }
	string name() const { return "Octahedron"; }
};

//! @brief Icosahedron rotation group
class Icosahedron : public S2quad
{
public:
	Icosahedron();
	int jMax() const { return 5; }
	int nS1() const { return 5; } //!< This is a special case that doesn't have nS1=jMax+1
	string name() const { return "Icosahedron"; }
};


//--------------- Hardin and Sloane Spherical designs -----------------
// Auto-converted from http://www2.research.att.com/~njas/sphdesigns

//! @brief Spherical 7-design with 24 nodes
class S2_7design_24 : public S2quad
{
public:
	S2_7design_24();
	int jMax() const { return 7; }
	int nS1() const { return 8; }
	string name() const { return "7-design"; }
};

//! @brief Spherical 8-design with 36 nodes
class S2_8design_36 : public S2quad
{
public:
	S2_8design_36();
	int jMax() const { return 8; }
	int nS1() const { return 9; }
	string name() const { return "8-design"; }
};

//! @brief Spherical 9-design with 48 nodes
class S2_9design_48 : public S2quad
{
public:
	S2_9design_48();
	int jMax() const { return 9; }
	int nS1() const { return 10; }
	string name() const { return "9-design"; }
};

//! @brief Spherical 10-design with 60 nodes
class S2_10design_60 : public S2quad
{
public:
	S2_10design_60();
	int jMax() const { return 10; }
	int nS1() const { return 11; }
	string name() const { return "10-design"; }
};

//! @brief Spherical 11-design with 70 nodes
class S2_11design_70 : public S2quad
{
public:
	S2_11design_70();
	int jMax() const { return 11; }
	int nS1() const { return 12; }
	string name() const { return "11-design"; }
};

//! @brief Spherical 12-design with 84 nodes
class S2_12design_84 : public S2quad
{
public:
	S2_12design_84();
	int jMax() const { return 12; }
	int nS1() const { return 13; }
	string name() const { return "12-design"; }
};

//! @brief Spherical 13-design with 94 nodes
class S2_13design_94 : public S2quad
{
public:
	S2_13design_94();
	int jMax() const { return 13; }
	int nS1() const { return 14; }
	string name() const { return "13-design"; }
};

//! @brief Spherical 14-design with 108 nodes
class S2_14design_108 : public S2quad
{
public:
	S2_14design_108();
	int jMax() const { return 14; }
	int nS1() const { return 15; }
	string name() const { return "14-design"; }
};

//! @brief Spherical 15-design with 120 nodes
class S2_15design_120 : public S2quad
{
public:
	S2_15design_120();
	int jMax() const { return 15; }
	int nS1() const { return 16; }
	string name() const { return "15-design"; }
};

//! @brief Spherical 16-design with 144 nodes
class S2_16design_144 : public S2quad
{
public:
	S2_16design_144();
	int jMax() const { return 16; }
	int nS1() const { return 17; }
	string name() const { return "16-design"; }
};

//! @brief Spherical 17-design with 156 nodes
class S2_17design_156 : public S2quad
{
public:
	S2_17design_156();
	int jMax() const { return 17; }
	int nS1() const { return 18; }
	string name() const { return "17-design"; }
};

//! @brief Spherical 18-design with 180 nodes
class S2_18design_180 : public S2quad
{
public:
	S2_18design_180();
	int jMax() const { return 18; }
	int nS1() const { return 19; }
	string name() const { return "18-design"; }
};

//! @brief Spherical 19-design with 204 nodes
class S2_19design_204 : public S2quad
{
public:
	S2_19design_204();
	int jMax() const { return 19; }
	int nS1() const { return 20; }
	string name() const { return "19-design"; }
};

//! @brief Spherical 20-design with 216 nodes
class S2_20design_216 : public S2quad
{
public:
	S2_20design_216();
	int jMax() const { return 20; }
	int nS1() const { return 21; }
	string name() const { return "20-design"; }
};

//! @brief Spherical 21-design with 240 nodes
class S2_21design_240 : public S2quad
{
public:
	S2_21design_240();
	int jMax() const { return 21; }
	int nS1() const { return 22; }
	string name() const { return "21-design"; }
};

//! @}

#endif // JDFTX_FLUID_SO3PLATONIC_H
