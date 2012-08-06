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

#ifndef JDFTX_ELECTRONIC_RADIALSCHRODINGER_H
#define JDFTX_ELECTRONIC_RADIALSCHRODINGER_H

#include <electronic/RadialFunction.h>
#include <electronic/matrix.h>
#include <core/Minimize.h>
#include <core/scalar.h>
#include <cstring>
#include <vector>
#include <map>

//! Radial schrodinger equation solver (non-interacting eigen-problem for an atom)
class RadialSchrodinger
{
public:
	//! Create a radial atom on grid specified by rArr and drArr, with potential V-Z/r
	//! iMatch sets the grid point at which matching occurs, and is set to the middle if 0 or unspecified
	RadialSchrodinger(const std::vector<double>& rArr, const std::vector<double>& drArr,
		const std::vector<double>& V, double Z, size_t iMatch=0.);
	
	//! Determine the optimum eigenfunction fillings (indexed by l and then nNodes)
	//! Optionally retrieve the corresponding eigenvalues (if Eptr is non-null)
	std::vector< std::vector<double> > getFillings(double nElectrons,
		std::vector< std::vector<double> >* Eptr=0);
	
	//! Optional outputs from compute: retrieve all non-null quantities
	struct Outputs
	{
		std::vector<double>* n; //!< electron density
		std::vector<double>* Dn; //!< density gradient (radial direction, others 0 by symmetry)
		std::vector<double>* tau; //!< kinetic energy density
		std::vector< std::vector< std::vector<complex> > >* z; //!< eigenfunctions indexed by l and nNodes. Each z = (u, u') (See schEqn)
		std::vector< std::vector<double> >* E; //!< eigenvalues  indexed by l and nNodes
		
		Outputs(std::vector<double>* n=0, std::vector<double>* Dn=0, std::vector<double>* tau=0,
			std::vector< std::vector< std::vector<complex> > >* z=0, std::vector< std::vector<double> >* E=0)
		: n(n), Dn(Dn), tau(tau), z(z), E(E) {}
	};
	
	//! Compute total non-interacting energy of atom, given its fillings, and collect any optional outputs
	double compute(const std::vector< std::vector<double> >& F, Outputs outputs);
	
private:
	const std::vector<double>& rArr; //!< radial grid
	const std::vector<double>& drArr; //!< intgeration weights for radial grid, divided by 4*pi*r^2
	std::vector<double> Vsub; //!< sampled part of V, including midway points
	double Z; //!< nuclear charge which contributes -Z/r to the potential in addition to V
	size_t iMatch; //!< index of point at which inward and outward solves are matched
	
	std::vector< std::map<int,double> > nodesEmin, nodesEmax; //!< min and max E so far recorded at a given node count, for each l.
	std::vector< std::map<double,double> > cachedEerr; //!< Cached matching error for attempted E's for each l
	std::vector< std::map<int,double> > cachedE; //!< Cached eigenvalues for each l and node count

	//! Use complex arithmetic to represent the second order schrodinger equation
	//!  Real(z) stores u(r) = r^(1-l) R(r) and Imag(z) stores u'(r)
	//!  iSub is the index into the sub-sampled potential array
	inline complex schEqn(double r, int iSub, complex z, int l, double E) const
	{	double rInv = r ? 1./r : 0.;
		double V = Vsub[iSub] - Z*rInv;
		return complex(z.imag(), 2*(z.real()*(V-E) - z.imag()*(l+1)*rInv));
	}
	
	//! Inward and outward integration of the Schrodinger equation
	//! @param l Angular momentum
	//! @param E Eigenvalue
	//! @param zSave Retrieve z = (u, u') with u = R/r^l at all points (if non-null)
	//! @return Number of nodes in radial wavefunction (this value is also cached in the nodesEmin and nodesEmax maps)
	//! The returned wavefunction is normalized according to Integral 4 pi r^2 R(r) = 1
	//! The inward and outward portions are matched at iMatch'th grid point (can be set using construtcor)
	//! The matching error is stored in cachedEerr and can be retrieved using getEerr()
	int solveSchEqn(int l, double E, complex* zSave=0);
	
	//! Binary-search (interval halving or doubling as appropriate!)
	//! till given nNodes is encountered at particular l,
	//! and update nodesEmin[l] and nodesEmax[l]
	//! NOTE: The maps must be seeded with two calls to solveSchEqn(l, E)
	//! with distinct E's at this l, before calling this routine.
	void locateNodeCount(int l, int nNodes);
	
	//! Compute the error function whose roots are the eigenvalues
	double getEerr(int l, double E);
	
	//! Zero-in on a bracketed eigenvalue using Ridder's method
	double findE(int l, double Elo, double Ehi, double tol);

	//! Find eigenvalue (and optionally eigenfunction) for a specific node count and angular momentum
	double getEig(int l, int nNodes, complex* zEvec=0);
};


//! Invert Kohn-Sham equations to get effective potential for a given electron density
class InvertKS : public Minimizable<diagMatrix>
{
public:
	//! Initialize for a given target electron density nTarget
	InvertKS(const RadialFunctionR& nTarget);
	RadialFunctionR getTau(); //!< Get the orbital KE density

	//Functions for Minimizable interface:
	void step(const diagMatrix& dV, double alpha);
	double compute(diagMatrix* E_V);
	diagMatrix precondition(const diagMatrix& grad);
	
private:
	diagMatrix V; //!< Effective Kohn-Sham potential
	double nElectrons; //!< number of electrons
	const std::vector<double> &r, &dr, &n0; //!< Radial grid, integration weights and target electron density
	std::vector<double> n; //!< non-interacting electron density
	std::vector< std::vector<double> > F; //!< non-interacting fillings
};

#endif // JDFTX_ELECTRONIC_RADIALSCHRODINGER_H
