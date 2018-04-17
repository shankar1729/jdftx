/*-------------------------------------------------------------------
Copyright 2017 Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_TETRAHEDRALDOS_H
#define JDFTX_ELECTRONIC_TETRAHEDRALDOS_H

#include <core/matrix.h>
#include <core/string.h>
#include <vector>
#include <array>

//! Evaluate DOS using the tetrahedron method
class TetrahedralDOS
{
public:
	
	const int nSpins; //!< number of separate spin channels
	const int nBands; //!< number of bands
	const int nWeights; //!< number of weighted-DOS's beig calculated (also counting total DOS i.e. weight = 1)
	const int nReduced; //!< number of reduced k-points
	const int nStates; //!< nReduced * nSpins
	
	//! Initialize calculator for a given uniform k-point mesh kmesh and indices iReduced to reduced mesh.
	//! If iReduced is empty, then nReduced = kmesh.size() and eigenvalues / weights must be provided on the full mesh
	//! R is the unit cell lattice vectors, while super are its linear combinations in the k-point sampled supercell.
	//! nSpins, nBands and nWeights initialize corresponding members of class.
	//! weightSum sets the net contribution from the BZ integral (can include spin-degeneracy factors here, if any).
	TetrahedralDOS(std::vector<vector3<>> kmesh, std::vector<int> iReduced,
		const matrix3<>& R, const matrix3<int>& super,
		int nSpins, int nBands, int nWeights, double weightSum=1.);
	
	inline double& e(int iState, int iBand) { return eigs[iState + nStates*iBand]; } //!< access eigenvalue
	inline const double& e(int iState, int iBand) const { return eigs[iState + nStates*iBand]; } //!< access eigenvalue (const version)
	inline double& w(int iWeight, int iState, int iBand) { return weights[iWeight + nWeights*(iState + nStates*iBand)]; } //!< access weight
	inline const double& w(int iWeight, int iState, int iBand) const { return weights[iWeight + nWeights*(iState + nStates*iBand)]; } //!< access weight (const version)
	
	void setEigs(const std::vector<diagMatrix>& E); //!< set all eigenvalues together (instead of using e())
	void setWeights(int iWeight, const std::vector<diagMatrix>& weights); //!< set all weights for given iWeight together (instead of using e()); all weights are initially 1
	
	//! Replace clusters of eigenvalues that differ by less than Etol by a single value equal to their mean
	void weldEigenvalues(double Etol);

	typedef std::pair<double, std::vector<double> > LsplineElem; //!< Single rnergy and DOS values with all weights at that energy
	typedef std::vector<LsplineElem> Lspline; //!< Set of all energy and DOS values as a linear spline

	//! Generate the density of states for a given spin channel
	//! Etol sets the width of the delta-function DOS of bands that are completely flat (potentially welded within Etol)
	Lspline getDOS(int iSpin, double Etol) const;

	//! Apply gaussian smoothing of width Esigma
	Lspline gaussSmooth(const Lspline& in, double Esigma) const;

	//! Write the given DOS Lspline to a file (optionally with a header)
	void printDOS(const Lspline& dos, string filename, string header=string());

private:
	//Tetrahedron in Brillouin zone triangulation:
	struct Tetrahedron
	{	std::array<int,4> q; //state indices for the vertices
		double V; //volume of tetrahedron
	};
	std::vector<Tetrahedron> tetrahedra;
	std::vector<double> eigs; //flat array of eigenvalues (inner index state, outer index bands)
	std::vector<double> weights; //flat array of DOS weights (inner index weight function, middle index state, and outer index bands)
	
	
	//! Accumulate contribution from one tetrahedron (exactly a cubic spline for linear interpolation)
	//! to the weighted DOS for all weight functions (from a single band)
	void accumTetrahedron(const Tetrahedron& t, int iBand, int iSpin, struct Cspline& wdos) const;

	//! Coalesce overlapping splines: convert an arbitrary set of spline pieces into a regular ordered piecewise spline
	void coalesceIntervals(struct Cspline& cspline) const;

	//! Convert cubic splines to integrated linear splines which handle discontinuities and singularities better.
	//! The Cspline object must be coalesced before passing to this function
	Lspline convertLspline(const struct Cspline& cspline) const;

	//! Collect contributions from multiple linear splines (one for each band)
	Lspline mergeLsplines(const std::vector<Lspline>& lsplines) const;
};

#endif //JDFTX_ELECTRONIC_TETRAHEDRALDOS_H
