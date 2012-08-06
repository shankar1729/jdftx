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

#ifndef JDFTX_CORE_BSPLINE_H
#define JDFTX_CORE_BSPLINE_H

//! @addtogroup excessfunctionals
//! @{

//! @file Bspline.h
//! @brief Fast Bezier spline lookup table

#include <cstdio>

//! Wraps a cubic Bezier spline lookup table around function and derivative evaluations
class Bspline
{
protected:
	virtual double getValue(double)=0; //!< Function that needs to be converted to a lookup table (must implement in derived class)
	virtual double getDeriv(double)=0; //!< Derivative of the function that needs to be converted to a lookup table (must implement in derived class)
public:
	//! Initialize the lookup table to cover [xStart:xEnd)
	//! (getValue and getDeriv will get called outside this range)
	Bspline(double xStart, double xEnd, int nIntervals);
	~Bspline();

	void init(); //!< Initialize the cubic spline coefficients based on getValue() and getDeriv()

	//! Print a report on estimated max errors in function and deriv due to interpolation to stream fp.
	//! If plotFileName is provided, a 3 column data [x getValue(x) value(x)] for x in [xStart xEnd] will be printed to it
	void maxErrorEstimate(FILE* fp, const char* plotFileName=0);

	inline double value(double x); //!< Calculate function using lookup table
	inline double deriv(double x); //!< Calculate derivative using lookup table
	inline void value_deriv(double x, double& value, double& deriv); //!< Calculate function and derivative using lookup table
private:
	int nCoeff;
	double xStart, h, invh;
	double* coeff;
};

//! @}

//###################################################################################################
//####  Implementation  ####
//##########################
//!@cond


#include <cmath>

inline double Bspline::value(double x)
{	register double t = invh*(x-xStart);
	register int it = floor(t);
	if(it<0 || it>=nCoeff) return getValue(x);
	register double* a = coeff + 4*it;
	t -= it;
	return a[0] + t*(a[1] + t*(a[2] + t*a[3]));
}

inline double Bspline::deriv(double x)
{	register double t = invh*(x-xStart);
	register int it = floor(t);
	if(it<0 || it>=nCoeff) return getDeriv(x);
	register double* a = coeff + 4*it;
	t -= it;
	return invh*(a[1] + t*(2.0*a[2] + t*3.0*a[3]));
}

inline void Bspline::value_deriv(double x, double& value, double& deriv)
{	register double t = invh*(x-xStart);
	register int it = floor(t);
	if(it<0 || it>=nCoeff) { value=getValue(x); deriv=getDeriv(x); return; }
	register double* a = coeff + 4*it;
	t -= it;
	value = a[0] + t*(a[1] + t*(a[2] + t*a[3]));
	deriv = invh*(a[1] + t*(2.0*a[2] + t*3.0*a[3]));
}

//! @endcond

#endif // JDFTX_CORE_BSPLINE_H
