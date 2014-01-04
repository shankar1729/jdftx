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

#ifndef JDFTX_CORE_MATRIX3_H
#define JDFTX_CORE_MATRIX3_H

#include <core/vector3.h>

template<typename scalar=double> class matrix3
{
	scalar m[3][3];

public:
	//accessors:
	__hostanddev__ scalar& operator()(int i, int j) { return m[i][j]; }
	__hostanddev__ const scalar& operator()(int i, int j) const { return m[i][j]; }
	__hostanddev__ vector3<scalar> row(int i) const { return vector3<scalar>(m[i][0], m[i][1], m[i][2]); }
	__hostanddev__ vector3<scalar> column(int i) const { return vector3<scalar>(m[0][i], m[1][i], m[2][i]); }

	__hostanddev__ void set_row(int i, const vector3<scalar>& v) { for(int j=0; j<3; j++) m[i][j] = v[j]; }
	__hostanddev__ void set_rows(const vector3<scalar>& v0, const vector3<scalar>& v1, const vector3<scalar>& v2)
	{	for(int j=0; j<3; j++) { m[0][j] = v0[j]; m[1][j] = v1[j]; m[2][j] = v2[j]; }
	}
	__hostanddev__ void set_col(int j, const vector3<scalar>& v) { for(int i=0; i<3; i++) m[i][j] = v[i]; }
	__hostanddev__ void set_cols(const vector3<scalar>& v0, const vector3<scalar>& v1, const vector3<scalar>& v2)
	{	for(int i=0; i<3; i++) { m[i][0] = v0[i]; m[i][1] = v1[i]; m[i][2] = v2[i]; }
	}

	//constructors:
	explicit __hostanddev__ matrix3(scalar d0=0, scalar d1=0, scalar d2=0)
	{	m[0][0] = d0; m[1][1] = d1, m[2][2] = d2;
		m[0][1] = m[0][2] = m[1][0] = m[1][2] = m[2][0] = m[2][1] = 0.0;
	}
	__hostanddev__ matrix3(
		scalar m00, scalar m01, scalar m02,
		scalar m10, scalar m11, scalar m12,
		scalar m20, scalar m21, scalar m22 )
	{	m[0][0] = m00; m[0][1] = m01; m[0][2] = m02;
		m[1][0] = m10; m[1][1] = m11; m[1][2] = m12;
		m[2][0] = m20; m[2][1] = m21; m[2][2] = m22;
	}
	template<typename scalar2> explicit __hostanddev__ matrix3(const matrix3<scalar2>& n)
	{	for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				m[i][j] = scalar(n(i,j));
	}
	
	//arithmetic operators
	__hostanddev__ matrix3<scalar> operator-() const
	{	matrix3<scalar> n;
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				n(i,j) = -m[i][j];
		return n;
	}
	__hostanddev__ matrix3<scalar> operator+(const matrix3<scalar> &n) const
	{	matrix3<scalar> ret;
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				ret(i,j) = m[i][j] + n(i,j);
		return ret;
	}
	__hostanddev__ matrix3<scalar>& operator+=(const matrix3<scalar> &n)
	{	for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				m[i][j] += n(i,j);
		return *this;
	}
	__hostanddev__ matrix3<scalar> operator-(const matrix3<scalar> &n) const
	{	matrix3<scalar> ret;
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				ret(i,j) = m[i][j] - n(i,j);
		return ret;
	}
	__hostanddev__ matrix3<scalar>& operator-=(const matrix3<scalar> &n)
	{	for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				m[i][j] -= n(i,j);
		return *this;
	}
	__hostanddev__ matrix3<scalar>& operator*=(scalar s)
	{	for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				m[i][j] *= s;
		return *this;
	}
	__hostanddev__ matrix3<scalar> operator*(scalar s) const
	{	matrix3<scalar> ret;
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				ret(i,j) = m[i][j] * s;
		return ret;
	}
	
	//!Compute the length of a vector, using this matrix as the metric
	#define METRIC_LENGTH_SQUARED \
		return v[0]*v[0]*m[0][0] + v[1]*v[1]*m[1][1] + v[2]*v[2]*m[2][2] \
			+ 2*(v[0]*v[1]*m[0][1] + v[0]*v[2]*m[0][2] + v[1]*v[2]*m[1][2]);
	__hostanddev__ double metric_length_squared(const vector3<double> &v) const { METRIC_LENGTH_SQUARED }
	__hostanddev__ scalar metric_length_squared(const vector3<int> &v) const { METRIC_LENGTH_SQUARED }
	#undef METRIC_LENGTH_SQUARED
	
	__hostanddev__ matrix3<scalar> operator/(scalar s) const { return (*this) * (1.0/s); }
	__hostanddev__ matrix3<scalar>& operator/(scalar s) { return (*this) *= (1.0/s); }

	//! transpose
	__hostanddev__ matrix3<scalar> operator~() const
	{	matrix3<scalar> ret;
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				ret(i,j) = m[j][i];
		return ret;
	}

	void print(FILE* fp, const char *format, bool brackets=true) const
	{	for(int i=0; i<3; i++)
		{	if(brackets) fprintf(fp, "[ ");
			for(int j=0; j<3; j++) fprintf(fp, format, m[i][j]);
			if(brackets) fprintf(fp, " ]\n"); else fprintf(fp, "\n");
		}
	}

	//Comparison operators
	__hostanddev__ bool operator==(const matrix3<scalar>& n) const
	{	for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				if(m[i][j] != n(i,j))
					return false;
		return true;
	}
	__hostanddev__ bool operator!=(const matrix3<scalar>& n) const
	{	return ! (*this == n);
	}
};

//Multiplies:
template<typename scalar> __hostanddev__ matrix3<scalar> operator*(scalar s, const matrix3<scalar> &m) { return m*s; }

template<typename scalar> __hostanddev__ matrix3<scalar> outer(const vector3<scalar> &a, const vector3<scalar> &b)
{	matrix3<scalar> m;
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			m(i,j) = a[i] * b[j];
	return m;
}

#define MUL_MAT_VEC(retType) \
	vector3<retType> ret; \
	for(int i=0; i<3; i++) \
		for(int j=0; j<3; j++) \
			ret[i] += m(i,j) * v[j]; \
	return ret;
template<typename scalar> __hostanddev__ vector3<scalar> operator*(const matrix3<scalar>& m, const vector3<scalar> &v)
{	MUL_MAT_VEC(scalar)
}
template<typename scalar> __hostanddev__ vector3<scalar> operator*(const matrix3<scalar>& m, const vector3<int> &v)
{	MUL_MAT_VEC(scalar)
}
template<typename scalar> __hostanddev__ vector3<scalar> operator*(const matrix3<int>& m, const vector3<scalar> &v)
{	MUL_MAT_VEC(scalar)
}
__hostanddev__ vector3<int> operator*(const matrix3<int>& m, const vector3<int> &v)
{	MUL_MAT_VEC(int)
}
#undef MUL_MAT_VEC

#define MUL_VEC_MAT(retType) \
	vector3<retType> ret; \
	for(int i=0; i<3; i++) \
		for(int j=0; j<3; j++) \
			ret[j] += v[i] * m(i,j); \
	return ret;
template<typename scalar> __hostanddev__ vector3<scalar> operator*(const vector3<scalar> &v, const matrix3<scalar>& m)
{	MUL_VEC_MAT(scalar)
}
template<typename scalar> __hostanddev__ vector3<scalar> operator*(const vector3<int> &v, const matrix3<scalar>& m)
{	MUL_VEC_MAT(scalar)
}
template<typename scalar> __hostanddev__ vector3<scalar> operator*(const vector3<scalar> &v, const matrix3<int>& m)
{	MUL_VEC_MAT(scalar)
}
__hostanddev__ vector3<int> operator*(const vector3<int> &v, const matrix3<int>& m)
{	MUL_VEC_MAT(int)
}
#undef MUL_VEC_MAT

#define MUL_MAT_MAT(retType) \
	matrix3<retType> ret; \
	for(int i=0; i<3; i++) \
		for(int j=0; j<3; j++) \
			for(int k=0; k<3; k++) \
				ret(i,j) += m(i,k) * n(k,j); \
	return ret;
template<typename scalar> __hostanddev__ matrix3<scalar> operator*(const matrix3<scalar> &m, const matrix3<scalar>& n)
{	MUL_MAT_MAT(scalar)
}
template<typename scalar> __hostanddev__ matrix3<scalar> operator*(const matrix3<scalar> &m, const matrix3<int>& n)
{	MUL_MAT_MAT(scalar)
}
template<typename scalar> __hostanddev__ matrix3<scalar> operator*(const matrix3<int> &m, const matrix3<scalar>& n)
{	MUL_MAT_MAT(scalar)
}
__hostanddev__ matrix3<int> operator*(const matrix3<int> &m, const matrix3<int>& n)
{	MUL_MAT_MAT(int)
}
#undef MUL_MAT_MAT
template<typename scalar> __hostanddev__ matrix3<scalar>& operator*=(matrix3<scalar> &m, const matrix3<scalar>& n)
{ return (m = m * n); }


template<typename scalar> __hostanddev__ matrix3<scalar> Diag(vector3<scalar> v)
{ return matrix3<scalar>(v[0],v[1],v[2]); }

template<typename scalar> __hostanddev__ scalar trace(const matrix3<scalar> &m) { return m(0,0)+m(1,1)+m(2,2); }
template<typename scalar> __hostanddev__ scalar det(const matrix3<scalar> &m) { return box(m.row(0),m.row(1),m.row(2)); }
template<typename scalar> __hostanddev__ matrix3<scalar> adjugate(const matrix3<scalar> &m)
{	matrix3<scalar> adj;
	adj.set_cols(cross(m.row(1),m.row(2)), cross(m.row(2),m.row(0)), cross(m.row(0),m.row(1)));
	return adj;
}
__hostanddev__ matrix3<> inv(const matrix3<> &m)
{	return (1./det(m)) * adjugate(m);
}
__hostanddev__ double nrm2(const matrix3<>& m) { return sqrt(trace((~m)*m)); } //!< 2-norm of matrix

//! Create a rotation matrix
__hostanddev__ matrix3<> rotation(double theta, int axis)
{	double s, c; sincos(theta, &s, &c);
	switch(axis)
	{	case 0:
			return matrix3<>(
			 1,  0,  0,
			 0,  c,  s,
			 0, -s,  c );
		case 1:
			return matrix3<>(
			 c,  0, -s,
			 0,  1,  0,
			 s,  0,  c );
		default:
			return matrix3<>(
			 c,  s,  0,
			-s,  c,  0,
			 0,  0,  1 );
	}
}

#endif //JDFTX_CORE_MATRIX3_H
