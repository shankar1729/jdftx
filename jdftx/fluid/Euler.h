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

#ifndef JDFTX_FLUID_EULER_H
#define JDFTX_FLUID_EULER_H

//! @addtogroup so3quad
//! @{

/** @file Euler.h
@brief Various Euler angle related utilities

The Euler angle convention used here is R(alpha,beta,gamma) = R_Z(alpha) * R_Y(beta) * R_Z(gamma)
All routines are inlined for performance.
*/

#include <cmath>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_pow_int.h>
#include <core/matrix3.h>

//! @brief Get the position of the new Z-axis, given alpha and beta (does not depend on gamma)
inline vector3<> polarUnitVector(double alpha, double beta)
{	register double sBeta = sin(beta);
	return vector3<>(sBeta*cos(alpha), sBeta*sin(alpha), cos(beta));
}

//! @brief Calculates euler[0]=alpha and euler[1]=beta given newZ, the direction of the Z-axis in the rotated frame
inline void getEulerAxis(const vector3<>& newZ, vector3<>& euler)
{	euler[1] = acos(newZ[2]/newZ.length());
	if(euler[1]*(M_PI-euler[1]) < 1e-6) //beta = 0 or pi
	{	euler[0] = 0.0; //alpha is meaningless
	}
	else
	{	euler[0] = atan2(newZ[1], newZ[0]);
	}
}

//! @brief Calculate rotation matrices from euler angles
inline matrix3<> matrixFromEuler(const vector3<>& euler)
{	return rotation(euler[0], 2) //R_Z(alpha)
		* rotation(euler[1], 1) //R_Y(beta)
		* rotation(euler[2], 2); //R_Z(gamma)
}

//! @brief Calculate euler angles from rotation matrices
inline vector3<> eulerFromMatrix(const matrix3<>& mat)
{	if(fabs(mat(2,2))>1.-1e-7) //beta is 0 or pi and only gamma+/-alpha is determined (set alpha=0 w.l.og)
		return vector3<>(0., //alpha (set 0 w.l.o.g)
			(mat(2,2)>0 ? 0. : M_PI), //beta
			atan2(-mat(1,0),mat(1,1)) ); //gamma
	else
		return vector3<>(	atan2(mat(1,2),-mat(0,2)), //alpha
					acos(mat(2,2)), //beta
					atan2(mat(2,1),mat(2,0)) ); //gamma
}


//! @brief Wigner d-function d^j_{m1 m2}(beta) for ZYZ Euler angles
inline double wigner_d(const int j, const int m1, const int m2, const double beta)
{	int t1 = j+m1;
	int t2 = j-m2;
	int t3 = m2-m1;
	int t4 = 2*j + m1 - m2;

	int sMin = (-t3 < 0) ? 0 : -t3;
	int sMax = (t2 < t1) ? t2 : t1;

	double cosbeta = cos(beta * 0.5);
	double sinbeta = sin(beta * 0.5);
	double a = 0.5*(gsl_sf_lnfact(t1) + gsl_sf_lnfact(t2) + gsl_sf_lnfact(j-m1) + gsl_sf_lnfact(j+m2));

	register double result = 0.0;
	register double sign = (sMin % 2 ? -1.0 : 1.0);
	for(int s = sMin; s <= sMax; s++)
	{	result += copysign(gsl_sf_pow_int(cosbeta, t4-2*s) * gsl_sf_pow_int(sinbeta, t3+2*s)
				* exp(a - (gsl_sf_lnfact(t1 - s) + gsl_sf_lnfact(t2 - s) + gsl_sf_lnfact(t3 + s) + gsl_sf_lnfact(s))), sign);
		sign = -sign;
	}
	return result;
}


//! @}

#endif // JDFTX_FLUID_EULER_H
