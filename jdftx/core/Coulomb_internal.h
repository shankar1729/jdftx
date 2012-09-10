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

#ifndef JDFTX_CORE_COULOMB_INTERNAL_H
#define JDFTX_CORE_COULOMB_INTERNAL_H

//! @file Coulomb_internal.h Shared inline functions for anlaytical truncated Coulomb kernels

#include <core/matrix3.h>

//! Periodic coulomb interaction (4 pi/G^2)
struct CoulombPeriodic_calc
{	__hostanddev__ double operator()(const vector3<int>& iG, const matrix3<>& GGT) const
	{	double Gsq = GGT.metric_length_squared(iG);
		return Gsq ? (4*M_PI)/Gsq : 0.;
	}
};

//! Slab-truncated coulomb interaction
struct CoulombSlab_calc
{	int iDir; double hlfL;
	CoulombSlab_calc(int iDir, double hlfL) : iDir(iDir), hlfL(hlfL) {}
	__hostanddev__ double operator()(const vector3<int>& iG, const matrix3<>& GGT) const
	{	double Gsq = GGT.metric_length_squared(iG);
		double Gplane = Gsq - GGT(iDir,iDir) * iG[iDir]*iG[iDir]; //G along the non-truncated directions
		Gplane = Gplane>0. ? sqrt(Gplane) : 0.; //safe sqrt to prevent NaN from roundoff errors
		return (4*M_PI) * (Gsq ? (1. - exp(-Gplane*hlfL) * cos(M_PI*iG[iDir]))/Gsq : -0.5*hlfL*hlfL);
	}
};

//! Sphere-truncated coulomb interaction
struct CoulombSpherical_calc
{	double Rc;
	CoulombSpherical_calc(double Rc) : Rc(Rc) {}
	__hostanddev__ double operator()(const vector3<int>& iG, const matrix3<>& GGT) const
	{	double Gsq = GGT.metric_length_squared(iG);
		return Gsq ? (4*M_PI) * (1. - cos(Rc*sqrt(Gsq)))/Gsq : (2*M_PI)*Rc*Rc;
	}
};

#ifdef GPU_ENABLED
void coulombAnalytic_gpu(vector3<int> S, const matrix3<>& GGT, const CoulombPeriodic_calc& calc, complex* data);
void coulombAnalytic_gpu(vector3<int> S, const matrix3<>& GGT, const CoulombSlab_calc& calc, complex* data);
void coulombAnalytic_gpu(vector3<int> S, const matrix3<>& GGT, const CoulombSpherical_calc& calc, complex* data);
#endif
void coulombAnalytic(vector3<int> S, const matrix3<>& GGT, const CoulombPeriodic_calc& calc, complex* data);
void coulombAnalytic(vector3<int> S, const matrix3<>& GGT, const CoulombSlab_calc& calc, complex* data);
void coulombAnalytic(vector3<int> S, const matrix3<>& GGT, const CoulombSpherical_calc& calc, complex* data);

//! Compute erf(x)/x (with x~0 handled properly)
__hostanddev__ double erf_by_x(double x)
{	double xSq = x*x;
	if(xSq<1e-6) return (1./sqrt(M_PI))*(2. - xSq*(2./3 + 0.2*xSq));
	else return erf(x)/x;
}

#endif // JDFTX_CORE_COULOMB_INTERNAL_H
