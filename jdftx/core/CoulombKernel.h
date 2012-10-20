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

//! Wigner-Seitz truncated coulomb kernel generator
struct CoulombKernel
{	const matrix3<> R; //!< lattice vectors
	const vector3<int> S; //!< sample count
	const vector3<bool> isTruncated; //!< whether corresponding lattice direction is truncated
	double omega; //!< erf-screening parameter (used for screened exchange kernels)
	
	CoulombKernel(const matrix3<> R, const vector3<int> S, const vector3<bool> isTruncated, double omega=0.);
	
	//! Initialize the truncated kernel in data
	//! data must be allocated for S[0]*S[1]*(1+S[2]/2) entries (fftw c2r order).
	//! ws is the Wigner-Seitz cell corresponding to lattice vectors R.
	//!      Supported modes include fully truncated (Isolated or Wigner-Seitz
	//! truncated exchange kernel) and one direction periodic (Wire geometry).
	void compute(double* data, const WignerSeitz& ws) const;
	
	static const double nSigmasPerWidth; //!< number of sigmas at which gaussian is negligible at working precision
	
private:
	//Various indiviudally optimized cases of computeKernel:
	void computeIsolated(double* data, const WignerSeitz& ws) const; //!< Fully truncated
	void computeWire(double* data, const WignerSeitz& ws) const; //!< 1 periodic direction
};

#endif // JDFTX_CORE_COULOMBKERNEL_H
