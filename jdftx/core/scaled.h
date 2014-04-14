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

#ifndef JDFTX_CORE_SCALED_H
#define JDFTX_CORE_SCALED_H

//! @file scaled.h Template to avoid (delay) scaling operations on linear objects

template<typename T> struct scaled
{	const T& data;
	double scale;
	scaled(const T& data, double scale=1.0) : data(data), scale(scale) {}
	operator T() const { T ret(data); return ret *= scale; }
	scaled<T>& operator*=(double s) { scale *= s; return *this; }
};

template<typename T> T& operator+=(T& y, const scaled<T>& x) { if(y) axpy(x.scale, x.data, y); else y = x.scale * x.data; return y; }
template<typename T> T& operator-=(T& y, const scaled<T>& x) { if(y) axpy(-x.scale, x.data, y); else y = (-x.scale) * x.data; return y; }
template<typename T> T operator+(const scaled<T>& x, const scaled<T>& y) { T ret(x); ret += y; return ret; }
template<typename T> T operator-(const scaled<T>& x, const scaled<T>& y) { T ret(x); ret -= y; return ret; }

template<typename T> scaled<T> operator-(const scaled<T>& x)  { return scaled<T>(x.data, -x.scale); }
template<typename T> scaled<T> operator*(double s, const scaled<T>& x) { return scaled<T>(x.data, x.scale * s); }
template<typename T> scaled<T> operator*(const scaled<T>& x, double s) { return scaled<T>(x.data, x.scale * s); }

#endif
