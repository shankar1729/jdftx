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

#ifndef JDFTX_CORE_COULOMB_H
#define JDFTX_CORE_COULOMB_H

#include <core/GridInfo.h>
#include <core/string.h>
#include <memory>
#include <set>


struct CoulombParams
{	//! Truncation geometry
	enum Geometry
	{	Periodic, //!< Fully periodic calculation (default)
		Slab, //!< Truncated along one lattice direction, periodic in two
		Wire, //!< Truncated along two lattice directions, periodic in one
		Cylindrical, //!< Cylindrical truncation, with 1D periodicity along axis
		Isolated, //!< Isolated system (all directions truncated)
		Spherical //!< Spherical isolation in all directions
	};
	Geometry geometry; //!< Truncation geometry
	int iDir; //!< Truncated lattice direction for Slab or periodic direction for Wire
	double Rc; //!< Truncation radius for cylindrical / spherical modes (0 => in-radius of Wigner-Seitz cell)
	
	double ionMargin; //!< margin around ions when checking localization constraints
	
	bool embed; //!< whether to embed in double-sized box (along truncated directions) to compute Coulomb interactions
	vector3<> embedCenter; //!< 'center' of the system, when it is embedded into the larger box (in lattice coordinates)
	bool embedFluidMode; //!< if true, don't truncate, just evaluate coulomb interactions in the larger box (fluid screening does the image separation instead)
	
	vector3<> Efield; //!< electric field (in Cartesian coordinates, atomic units [Eh/e/a0])
	
	//Parameters for computing exchange integrals:
	//! Regularization method for G=0 singularities in exchange
	enum ExchangeRegularization
	{	None, //!< No regularization (3D periodic or non-periodic systems only)
		AuxiliaryFunction, //!< Auxiliary function method (3D periodic systems only) [P. Carrier et al, PRB 75, 205126 (2007)]
		ProbeChargeEwald, //!< Ewald sum on a probe charge per unit cell (3D/2D/1D periodic systems) [J. Paier et al, JCP 122, 234102 (2005)]
		SphericalTruncated, //!< Wigner-Seitz volume spherical truncation  [J. Spencer et al, PRB 77, 193110 (2008)]
		WignerSeitzTruncated //!< Wigner-Seitz cell truncation [R. Sundararaman et al, arXiv:1302.6204]
	};
	ExchangeRegularization exchangeRegularization; //!< exchange regularization method
	std::set<double> omegaSet; //!< set of exchange erf-screening parameters
	std::shared_ptr<class Supercell> supercell; //!< Description of k-point supercell for exchange
	
	CoulombParams();
	
	//! Create a Coulomb object corresponding to the parameters of this class
	std::shared_ptr<class Coulomb> createCoulomb(const GridInfo& gInfo) const;
	
	//! Get a list of which directions are truncated:
	vector3<bool> isTruncated() const;
};


//! Information required for pair-potential evaluations
struct Atom
{	double Z; //!< charge
	vector3<> pos; //!< position in lattice coordinates (covariant)
	vector3<> force; //!< force in lattice coordinates (contravariant)
	int atomicNumber; //!< atomic number
	int sp; //!< species index for atom; used for VDW
	
	Atom(double Z, vector3<> pos, vector3<> force=vector3<>(0,0,0), int atomicNumber=0, int sp=0)
	: Z(Z), pos(pos), force(force), atomicNumber(atomicNumber), sp(sp)
	{
	}
};


//! Abstract base class for Ewald summation in arbitrary dimension
class Ewald
{
public:
	//!Get the energy of a point charge configurtaion, and accumulate corresponding forces
	//!The implementation will shift each Atom::pos by lattice vectors to bring it to
	//!the fundamental zone (or Wigner-Seitz cell as appropriate)
	virtual double energyAndGrad(std::vector<Atom>& atoms) const=0;
};


//! Abstract base class for the (optionally truncated) Coulomb interaction
class Coulomb
{
public:
	//! Special point-charge handling mode when using embedded truncation
	enum PointChargeMode
	{	PointChargeNone, //!< smooth charge distributions, no handling required
		PointChargeLeft, //!< point charge distribution on left end of operator (use for pointCharge^ K smooth)
		PointChargeRight //!< point charge distribution on right end of operator (use for smooth^ K pointCharge, hermitian conjugate of PointChargeLeft)
	};
	
	//! Apply Coulomb kernel (destructible input).
	//! Pass appropriate pointChargeMode when applying to nucelar densities for special handling of high-frequency
	//! components necessary in the non-translationally invariant scheme i.e. when params.embed==true
	ScalarFieldTilde operator()(ScalarFieldTilde&&, PointChargeMode pointChargeMode=PointChargeNone) const;
	
	//! Apply Coulomb kernel (implemented in base class using virtual destructible input version)
	ScalarFieldTilde operator()(const ScalarFieldTilde&, PointChargeMode pointChargeMode=PointChargeNone) const;
	
	//! Create the appropriate Ewald class, if required, and call Ewald::energyAndGrad
	//! Includes interaction with Efield, if present (Requires embedded truncation)
	double energyAndGrad(std::vector<Atom>& atoms) const; 

	//! Generate the potential due to the Efield (if any) (Requires embedded truncation)
	ScalarField getEfieldPotential() const;
	
	//! Apply regularized coulomb kernel for exchange integral with k-point difference kDiff
	//! and optionally screened with range parameter omega (destructible input)
	complexScalarFieldTilde operator()(complexScalarFieldTilde&&, vector3<> kDiff, double omega) const;

	//! Apply regularized coulomb kernel for exchange integral with k-point difference kDiff
	//! and optionally screened with range parameter omega (destructible input)
	complexScalarFieldTilde operator()(const complexScalarFieldTilde&, vector3<> kDiff, double omega) const;

private:
	const GridInfo& gInfoOrig; //!< original grid
protected:
	const CoulombParams& params;
	const GridInfo& gInfo; //!< embedding grid, which is 2x larger in truncated directions if params.embed == true
	std::shared_ptr<Ewald> ewald;
	std::map<double, std::shared_ptr<struct ExchangeEval>> exchangeEval;
	friend class ExchangeEval;
	
	Coulomb(const GridInfo& gInfoOrig, const CoulombParams& params);
	virtual ~Coulomb();
	
	//! Call to initialize exchangeEval if exact exchange is required
	//! NOTE: this must be called from the end of each derived class constructor
	void initExchangeEval();
	
	//!Apply the Coulomb operator (on optionally embedded grid) with appropriate truncation
	//!Embedding is handled in base class wrapper functions above
	virtual ScalarFieldTilde apply(ScalarFieldTilde&&) const=0;
	
	//!Each implementation must create and return the corresponding Ewald evaluator
	//!for the supplied lattice vectors R which may correspond to a supercell of
	//!gInfo.R along the periodic directions (the truncated directions will be identical)
	//!The number of atoms may be used for choosing the optimum gaussian width sigma
	virtual std::shared_ptr<Ewald> createEwald(matrix3<> R, size_t nAtoms) const=0;

private:
	//Data for mesh-embedded truncation:
	GridInfo gInfoEmbed; //!< embedding grid - internal object initialized only if params.embed == true
	vector3<int> ivCenter; //!< params.embedCenter in original grid mesh coordinates, rounded to closest grid point
	vector3<> xCenter; //!< meshCenter in original grid lattice coordinates (differs from params.embedCenter only due to grid-point rounding)
	vector3<> embedScale; //!< lattiec coordinates scale factor for moving from original to embedding mesh
	int* embedIndex; //!< list of indices from original mesh to larger mesh
	std::vector< std::pair<int,int*> > symmIndex; //!< list of number of equivalence classes and corresponding indices per cardinality for boundary symmetrization
	struct WignerSeitz* wsOrig; //!< Wigner-seitz cell of original mesh
	double ionWidth; //!< Range separation parameter for dealing with point charges in the embedded method
	RealKernel* ionKernel;
	ScalarFieldTilde embedExpand(const ScalarFieldTilde& in) const; //!< expand to embedding grid and symmetrize boundaries
	complexScalarFieldTilde embedExpand(complexScalarFieldTilde&& in) const; //!< expand to embedding grid and symmetrize boundaries
	ScalarFieldTilde embedShrink(const ScalarFieldTilde& in) const; //!< symmetrize boundaries and shrink to original grid (dagger of embedExpand)
	complexScalarFieldTilde embedShrink(complexScalarFieldTilde&& in) const; //!< symmetrize boundaries and shrink to original grid (dagger of embedExpand)
	friend class FluidSolver;
	friend class SlabEpsilon;
	friend class ChargedDefect;
};

#endif // JDFTX_CORE_COULOMB_H
