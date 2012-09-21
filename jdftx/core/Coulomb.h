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


struct ExchangeRegularization
{	//! Regularization method (for G=0 singularities in exchange)
	enum Method
	{	None, //!< No regularization (3D periodic or non-periodic systems only)
		AuxiliaryFunction, //!< Auxiliary function method (3D periodic systems only) [P. Carrier et al, PRB 75, 205126 (2007)]
		SphericalTruncated, //!< Wigner-Seitz volume spherical truncation  [J. Spencer et al, PRB 77, 193110 (2008)]
		WignerSeitzTruncated //!< Wigner-Seitz cell truncation [R. Sundararaman et al, (under preparation)]
	};
	Method method; //!< Regularization method
	
	double sigma; //!< Gaussian width for Wigner-Seitz truncated method
	string filename; //!< File to cache computed exchange kernel (used only by method WignerSeitzTruncated)
};


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
	int iDir; //!< Truncated lattice direction for Slab; Periodic direction for Wire
	double borderWidth; //!< Border width in smoothed Wigner-Seitz truncation (used by Wire, Isolated)
	double ionMargin; //!< Extra margin for the separation between ions from the truncation boundary (all truncated modes)
	double Rc; //!< Truncation radius for cylindrical / spherical modes (0 => in-radius of Wigner-Seitz cell)
	string filename; //!< File to cache computed kernel in (used only by Wire, Isolated)
	
	//Parameters for computing exchange integrals:
	ExchangeRegularization exchangeRegularization; //!< exchange regularization method
	std::set<double> omegaSet; //!< set of exchange erf-screening parameters
	std::shared_ptr<class Supercell> supercell; //!< Description of k-point supercell for exchange
	
	//! Create a Coulomb object corresponding to the parameters of this class
	std::shared_ptr<class Coulomb> createCoulomb(const GridInfo& gInfo) const;
};


//! Information required for pair-potential evaluations
struct Atom
{	double Z; //!< charge
	int atomicNumber; //!< atomic number
	vector3<> pos; //!< position in lattice coordinates (covariant)
	vector3<> force; //!< force in lattice coordinates (contravariant)
};


//! Abstract base class for the (optionally truncated) Coulomb interaction
class Coulomb
{
public:
	//! Apply Coulomb kernel (destructible input): implemented in derived classes
	virtual DataGptr operator()(DataGptr&&) const=0;
	
	//! Apply Coulomb kernel (implemented in base class using virtual destructible input version)
	DataGptr operator()(const DataGptr&) const;
	
	//!Get the energy of a point charge configurtaion, and accumulate corresponding forces
	//!The implementation will shift each Atom::pos by lattice vectors to bring it to
	//!the fundamental zone (or Wigner-Seitz cell as appropriate)
	virtual double energyAndGrad(std::vector<Atom>& atoms) const=0; 

	//! Apply regularized coulomb kernel for exchange integral with k-point difference kDiff
	//! and optionally screened with range parameter omega (destructible input)
	complexDataGptr operator()(complexDataGptr&&, vector3<> kDiff, double omega) const;

	//! Apply regularized coulomb kernel for exchange integral with k-point difference kDiff
	//! and optionally screened with range parameter omega (destructible input)
	complexDataGptr operator()(const complexDataGptr&, vector3<> kDiff, double omega) const;

protected:
	const GridInfo& gInfo;
	const CoulombParams& params;
	std::map<double, std::shared_ptr<struct ExchangeEval>> exchangeEval;
	
	Coulomb(const GridInfo& gInfo, const CoulombParams& params);
	
	//! Call to initialize exchangeEval if exact exchange is required
	//! NOTE: this must be called from the end of each derived class constructor
	void initExchangeEval();
};

#endif // JDFTX_CORE_COULOMB_H
