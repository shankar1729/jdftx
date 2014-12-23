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

#ifndef JDFTX_CORE_SCALAR_H
#define JDFTX_CORE_SCALAR_H

#include <cmath>

#ifndef __device__ //in .cpp files
	#define __hostanddev__ inline
#else //in .cu files
	#define __hostanddev__ inline __device__ __host__
	#define __in_a_cu_file__
	#include <cuda_runtime.h>
#endif

//! Ceiling of a positive integer division, templated over int types
template<class T> T ceildiv(T num, T den) { return (num + den - 1)/den; }

//! Return largest multiple of den smaller than num, templated over int types
template<class T> T floorMultiple(T num, T den) { return (num/den)*den; }

#ifdef __APPLE__
#ifndef __in_a_cu_file__
inline void sincos(double x, double* s, double* c)
{	*s = sin(x);
	*c = cos(x);
}
#endif
#endif

//! Complex number (need to define our own because we need operators for gpu code as well)
struct complex
{	double x, y;

	//Accessors
	__hostanddev__ double& real() { return x; }
	__hostanddev__ double& imag() { return y; }
	__hostanddev__ const double& real() const { return x; }
	__hostanddev__ const double& imag() const { return y; }

	//Constructors
	__hostanddev__ complex(double x=0, double y=0) : x(x), y(y) {}
	#ifdef __in_a_cu_file__
	__hostanddev__ complex(const double2& c) : x(c.x), y(c.y) {} //!< convert from cuda complex
	__hostanddev__ operator double2() const { double2 ret; ret.x=x; ret.y=y; return ret;} //!< convert to cuda complex
	#endif

	//Arithmetic
	__hostanddev__ complex& operator+=(const complex& c) { x+=c.x; y+=c.y; return *this; }
	__hostanddev__ complex& operator+=(double r) { x+=r; return *this; }
	__hostanddev__ complex operator+(const complex& c) const { return complex(x+c.x, y+c.y); }
	__hostanddev__ complex operator+(double r) const { return complex(x+r, y); }
	__hostanddev__ complex& operator-=(const complex& c) { x-=c.x; y-=c.y; return *this; }
	__hostanddev__ complex& operator-=(double r) { x-=r; return *this; }
	__hostanddev__ complex operator-(const complex& c) const { return complex(x-c.x, y-c.y); }
	__hostanddev__ complex operator-(double r) const { return complex(x-r, y); }
	__hostanddev__ complex operator-() const { return complex(-x, -y); }
	__hostanddev__ complex& operator*=(const complex& c) { return (*this = *this * c); }
	__hostanddev__ complex& operator*=(double r) { x*=r; y*=r; return *this; }
	__hostanddev__ complex operator*(const complex& c) const { return complex(x*c.x-y*c.y, y*c.x+x*c.y); }
	__hostanddev__ complex operator*(double r) const { return complex(x*r, y*r); }
	__hostanddev__ complex& operator/=(const complex& c) { return (*this = *this / c); }
	__hostanddev__ complex& operator/=(double r) { return (*this *= 1.0/r); }
	__hostanddev__ complex operator/(const complex& c) const { return complex(x*c.x+y*c.y, y*c.x-x*c.y) / c.norm(); }
	__hostanddev__ complex operator/(double r) const { return *this * (1.0/r); }

	__hostanddev__ double norm() const { return x*x + y*y; }
	__hostanddev__ double abs() const { return sqrt(norm()); }
	__hostanddev__ double arg() const { return atan2(y,x); }
	__hostanddev__ complex conj() const { return complex(x,-y); }
};

__hostanddev__ double real(const complex& c) { return c.real(); }
__hostanddev__ double imag(const complex& c) { return c.imag(); }
__hostanddev__ double norm(const complex& c) { return c.norm(); }
__hostanddev__ double abs(const complex& c) { return c.abs(); }
__hostanddev__ double arg(const complex& c) { return c.arg(); }
__hostanddev__ complex conj(const complex& c) { return c.conj(); }
__hostanddev__ double conj(const double& c) { return c; } //provided to ease templating over complex and double

__hostanddev__ complex operator+(double r, const complex& c) { return c+r; }
__hostanddev__ complex operator-(double r, const complex& c) { return -c+r; }
__hostanddev__ complex operator*(double r, const complex& c) { return c*r; }


//! Compute cis(x) = exp(iota x) = cos(x) + iota sin(x)
__hostanddev__ complex cis(double x)
{	double s, c; sincos(x, &s, &c);
	return complex(c, s);
}

#endif // JDFTX_CORE_SCALAR_H
