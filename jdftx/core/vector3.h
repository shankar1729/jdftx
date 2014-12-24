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

#ifndef JDFTX_CORE_VECTOR3_H
#define JDFTX_CORE_VECTOR3_H

//! @file vector3.h
//! 3-vector and fields there of

#include <core/scalar.h>
#include <vector>
#include <cstdio>

#define LOOP3(code) { for(int k=0; k<3; k++) { code } }

//! Generic 3-vector
template<typename scalar=double> class vector3
{
	scalar v[3];
public:
	//Accessors
	__hostanddev__ scalar& operator[](int k) { return v[k]; }
	__hostanddev__ const scalar& operator[](int k) const { return v[k]; }
	__hostanddev__ scalar& x() { return v[0]; }
	__hostanddev__ scalar& y() { return v[1]; }
	__hostanddev__ scalar& z() { return v[2]; }
	__hostanddev__ const scalar& x() const { return v[0]; }
	__hostanddev__ const scalar& y() const { return v[1]; }
	__hostanddev__ const scalar& z() const { return v[2]; }

	//Constructor
	__hostanddev__ explicit vector3(scalar a=0, scalar b=0, scalar c=0) { v[0]=a; v[1]=b; v[2]=c; }
	vector3(std::vector<scalar> a) { LOOP3(v[k]=a[k];) }
	template<typename scalar2> __hostanddev__ explicit vector3(const vector3<scalar2>& a) { LOOP3(v[k]=a[k];) } 
	
	//Arithmetic:
	__hostanddev__ vector3 operator+(const vector3 &a) const { return vector3(v[0]+a[0], v[1]+a[1], v[2]+a[2]); }
	__hostanddev__ vector3 operator+=(const vector3 &a) { LOOP3( v[k]+=a[k]; ) return *this; }
	__hostanddev__ vector3 operator+(const scalar a) const { return vector3(v[0]+a, v[1]+a, v[2]+a); }
	__hostanddev__ vector3 operator+=(const scalar a) { LOOP3( v[k]+=a; ) return *this; }

	__hostanddev__ vector3 operator-() const { return vector3(-v[0],-v[1],-v[2]); }
	__hostanddev__ vector3 operator-(const vector3 &a) const { return vector3(v[0]-a[0], v[1]-a[1], v[2]-a[2]); }
	__hostanddev__ vector3 operator-=(const vector3 &a) { LOOP3( v[k]-=a[k]; ) return *this; }

	__hostanddev__ vector3 operator/(scalar s) const { return (*this)*(1.0/s); }
	__hostanddev__ vector3& operator/=(scalar s) { return (*this)*=(1.0/s); }

	__hostanddev__ scalar length_squared() const { return v[0]*v[0] + v[1]*v[1] + v[2]*v[2]; }
	__hostanddev__ scalar length() const { return sqrt(length_squared()); }

	void print(FILE* fp, const char *format) const { std::fprintf(fp, "[ "); LOOP3( fprintf(fp, format, v[k]); ) std::fprintf(fp, " ]\n"); }
	__hostanddev__ bool operator==(const vector3& w) const { LOOP3( if(v[k] != w[k]) return false; ) return true; }
	__hostanddev__ bool operator<(const vector3& w) const { LOOP3( if(v[k]!=w[k]) return v[k]<w[k]; ) return false; }
};

template<typename scalar> __hostanddev__ vector3<scalar> operator+(scalar s, const vector3<scalar>& a) { return vector3<scalar>(a[0]+s, a[1]+s, a[2]+s); }
__hostanddev__ vector3<int> operator+(int s, const vector3<int>& a) { return vector3<int>(a[0]+s, a[1]+s, a[2]+s); }
__hostanddev__ vector3<int> operator+(const vector3<int>& a, int s) { return vector3<int>(a[0]+s, a[1]+s, a[2]+s); }
//Mixed additions (special cases)
template<typename scalar> __hostanddev__ vector3<scalar> operator+(const vector3<scalar>& a, int s) { return vector3<scalar>(a[0]+s, a[1]+s, a[2]+s); }
template<typename scalar> __hostanddev__ vector3<scalar> operator+(int s, const vector3<scalar>& a) { return vector3<scalar>(a[0]+s, a[1]+s, a[2]+s); }
__hostanddev__ vector3<> operator+(vector3<> a, vector3<int> b) { return vector3<>(a[0]+b[0], a[1]+b[1], a[2]+b[2]); }
__hostanddev__ vector3<> operator+(vector3<int> a, vector3<> b) { return vector3<>(a[0]+b[0], a[1]+b[1], a[2]+b[2]); }

template<typename scalar> __hostanddev__ vector3<scalar>& operator*=(vector3<scalar>& a, scalar s) { LOOP3(a[k]*=s;) return a; }
template<typename scalar> __hostanddev__ vector3<scalar>& operator*=(vector3<scalar>& a, double s) { LOOP3(a[k]*=s;) return a; }
template<typename scalar> __hostanddev__ vector3<scalar>& operator*=(vector3<scalar>& a, int s) { LOOP3(a[k]*=s;) return a; }
__hostanddev__ vector3<>& operator*=(vector3<>& a, double s) { LOOP3(a[k]*=s;) return a; }
__hostanddev__ vector3<>& operator*=(vector3<>& a, int s) { LOOP3(a[k]*=s;) return a; }
__hostanddev__ vector3<int>& operator*=(vector3<int>& a, int s) { LOOP3(a[k]*=s;) return a; }

template<typename scalar> __hostanddev__  vector3<scalar> operator*(scalar s, const vector3<scalar> &a) { vector3<scalar> v; LOOP3(v[k]=a[k]*s;) return v; }
template<typename scalar> __hostanddev__  vector3<scalar> operator*(double s, const vector3<scalar> &a) { vector3<scalar> v; LOOP3(v[k]=a[k]*s;) return v; }
template<typename scalar> __hostanddev__  vector3<scalar> operator*(const vector3<scalar> &a, double s) { vector3<scalar> v; LOOP3(v[k]=a[k]*s;) return v; }
template<typename scalar> __hostanddev__  vector3<scalar> operator*(int s, const vector3<scalar> &a) { vector3<scalar> v; LOOP3(v[k]=a[k]*s;) return v; }
template<typename scalar> __hostanddev__  vector3<scalar> operator*(const vector3<scalar> &a, int s) { vector3<scalar> v; LOOP3(v[k]=a[k]*s;) return v; }
__hostanddev__  vector3<> operator*(double s, const vector3<> &a) { vector3<> v; LOOP3(v[k]=a[k]*s;) return v; }
__hostanddev__  vector3<> operator*(const vector3<> &a, double s) { vector3<> v; LOOP3(v[k]=a[k]*s;) return v; }
__hostanddev__  vector3<> operator*(double s, const vector3<int> &a) { vector3<> v; LOOP3(v[k]=a[k]*s;) return v; }
__hostanddev__  vector3<> operator*(const vector3<int> &a, double s) { vector3<> v; LOOP3(v[k]=a[k]*s;) return v; }
__hostanddev__  vector3<complex> operator*(complex s, const vector3<> &a) { vector3<complex> v; LOOP3(v[k]=a[k]*s;) return v; }
__hostanddev__  vector3<complex> operator*(const vector3<> &a, complex s) { vector3<complex> v; LOOP3(v[k]=a[k]*s;) return v; }
__hostanddev__  vector3<complex> operator*(complex s, const vector3<int> &a) { vector3<complex> v; LOOP3(v[k]=a[k]*s;) return v; }
__hostanddev__  vector3<complex> operator*(const vector3<int> &a, complex s) { vector3<complex> v; LOOP3(v[k]=a[k]*s;) return v; }

//dot product of similar generic vectors, and all combinations with real vectors
template<typename scalar> __hostanddev__ scalar dot(const vector3<scalar>& a, const vector3<scalar>& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
template<typename scalar> __hostanddev__ scalar dot(const vector3<double>& a, const vector3<scalar>& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
template<typename scalar> __hostanddev__ scalar dot(const vector3<int>& a, const vector3<scalar>& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
template<typename scalar> __hostanddev__ scalar dot(const vector3<scalar>& a, const vector3<double>& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
template<typename scalar> __hostanddev__ scalar dot(const vector3<scalar>& a, const vector3<int>& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
__hostanddev__ double dot(const vector3<double>& a, const vector3<double>& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
__hostanddev__ double dot(const vector3<int>& a, const vector3<double>& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
__hostanddev__ double dot(const vector3<double>& a, const vector3<int>& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
__hostanddev__ int dot(const vector3<int>& a, const vector3<int>& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

//! cross product
template<typename scalar> __hostanddev__ vector3<scalar> cross(const vector3<scalar> &a, const vector3<scalar> &b)
{	return vector3<scalar>(
		a[1]*b[2] - a[2]*b[1],
		a[2]*b[0] - a[0]*b[2],
		a[0]*b[1] - a[1]*b[0] );
}

//! box product / triple product
template<typename scalar> __hostanddev__ scalar box(const vector3<scalar>& a, const vector3<scalar>& b, const vector3<scalar>& c)
{	return dot(a,cross(b,c));
}

//! Round vector3<> to vector3<int> (and optionally retrieve error)
__hostanddev__ vector3<int> round(const vector3<>& v, double* err=0)
{	vector3<int> out;
	LOOP3( out[k] = round(v[k]); )
	if(err) *err = (v+(-out)).length();
	return out;
} 

//! Return squared distance between a and b, interpreted as coordinates on a unit S1xS1xS1 embedded in R^6
//! i.e. sum of squares of distances between each coordinate put on a circle. This is useful for checking
//! distances in a periodic cell safely.
__hostanddev__ double circDistanceSquared(const vector3<>& a, const vector3<>& b)
{	return (cis(2*M_PI*a[0]) - cis(2*M_PI*b[0])).norm()
		+ (cis(2*M_PI*a[1]) - cis(2*M_PI*b[1])).norm()
		+ (cis(2*M_PI*a[2]) - cis(2*M_PI*b[2])).norm();
}


//! GCD of integers, templated over integer types
template<typename T> T gcd(T x, T y)
{	while(y != 0)
	{	T yPrev = y;
		y = x % y;
		x = yPrev;
	}
	return x;
}

//! Reduce an integer vector by its gcd
template<typename T> vector3<T> gcdReduce(const vector3<T>& d)
{	T g = gcd(gcd(d[0], d[1]), d[2]);
	return vector3<T>(d[0]/g, d[1]/g, d[2]/g);
}


//! Load vector from a constant vector field
template<typename scalar> __hostanddev__ vector3<scalar> loadVector(const vector3<const scalar*>& vArr, int i)
{	return vector3<scalar>( vArr[0][i], vArr[1][i], vArr[2][i] );
}
//! Load vector from a vector field
template<typename scalar> __hostanddev__ vector3<scalar> loadVector(const vector3<scalar*>& vArr, int i)
{	return vector3<scalar>( vArr[0][i], vArr[1][i], vArr[2][i] );
}
//! Store vector to a vector field
template<typename scalar> __hostanddev__ void storeVector(const vector3<scalar>& v, vector3<scalar*>& vArr, int i)
{	LOOP3( vArr[k][i] = v[k]; )
}
//! Accumulate vector onto a vector field
template<typename scalar> __hostanddev__ void accumVector(const vector3<scalar>& v, vector3<scalar*>& vArr, int i)
{	LOOP3( vArr[k][i] += v[k]; )
}

#undef LOOP3
#endif // JDFTX_CORE_VECTOR3_H
