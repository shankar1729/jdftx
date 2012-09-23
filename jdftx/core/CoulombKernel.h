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

#ifndef JDFTX_CORE_COULOMBKERNEL_H
#define JDFTX_CORE_COULOMBKERNEL_H

#include <core/WignerSeitz.h>
#include <core/matrix3.h>
#include <core/string.h>

//! @file CoulombKernel.h Shared routines for creating, saving and loading Wigner-Seitz truncated Coulomb kernels

//! Description of a Wigner-Seitz truncated coulomb kernel and corresponding I/O routines
//! The file format contains a header with geometry and truncation details
//! and will be compressed using lattice symmetries, if possible.
struct CoulombKernelDesc
{	const matrix3<> R; //!< lattice vectors
	const vector3<int> S; //!< sample count
	const vector3<bool> isTruncated; //!< whether corresponding lattice direction is truncated
	const vector3<> sigmaBorder; //!< border smearing width for each direction (non-orthogonal directions must have same sigma)
	double omega; //!< erf-screening parameter (used for screened exchange kernels)
	
	CoulombKernelDesc(const matrix3<> R, const vector3<int> S,
		const vector3<bool> isTruncated, const vector3<> sigmaBorder, double omega=0.);
	
	//! Save kernel to file.
	//! data must contain kernel in Fourier space in fftw c2r order (S[0]*S[1]*(1+S[2]/2) entries)
	//! data will be symmetrized appropriately, irrespective of whether filename is null.
	void saveKernel(double* data, string filename) const;
	
	//! Load kernel from file if available. Returns true if kernel loaded successfully.
	//! data must be allocated for S[0]*S[1]*(1+S[2]/2) entries (fftw c2r order)
	bool loadKernel(double* data, string filename) const;
	
	//! Compute the kernel in the fully isolated case (all entries of isTruncated = true)
	//! data must be allocated for S[0]*S[1]*(1+S[2]/2) entries (fftw c2r order)
	//! ws is the Wigner-Seitz cell corresponding to lattice vectors R
	void computeIsolatedKernel(double* data, const WignerSeitz& ws) const;
	
	static const double nSigmasPerWidth; //!< number of gaussian widths within border width
	
	//! Compute maximum sigma such that gaussians of sigma and sigmaOther
	//! with centers separated by L have negligible overlap at working precision
	static double getMaxSigma(double L, double sigmaOther);
	
private:
	std::vector<matrix3<int>> sym; //!< symmetry matrices
	//Various indiviudally optimized cases of computeIsolatedKernel
	void computeNonOrtho(double* data, const WignerSeitz& ws) const; //!< No lattice direction orthogonal to other two
	void computeRightPrism(double* data, const WignerSeitz& ws) const; //!< 1 lattice direction orthogonal to other 2
	void computeCylinder(double* data, const WignerSeitz& ws) const; //!< Capped cylinder (for exchange in Cylindrical truncation)
};

#endif // JDFTX_CORE_COULOMBKERNEL_H
