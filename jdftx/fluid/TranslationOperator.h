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

#ifndef JDFTX_FLUID_TRANSLATIONOPERATOR_H
#define JDFTX_FLUID_TRANSLATIONOPERATOR_H

//! @file TranslationOperator.h
//! Options for translating site densities for rigid molecule ideal gas implementations

#include <core/GridInfo.h>
#include <core/ScalarField.h>

//! Abstract base class for translation operators
class TranslationOperator
{
public:
	const GridInfo& gInfo;

	TranslationOperator(const GridInfo& gInfo);
	virtual ~TranslationOperator() {}

	//! Compute @f$ y += alpha T_t(x) @f$ ,
	//! where @f$ T_t @f$  is the translation operator @f$ T_t(x(r)) = x(r+t) @f$ modulo the lattice vectors
	//! T must conserve integral(x) and satisfy @f$ T^{\dagger}_t = T_{-t} @f$ exactly for gradient correctness
	//! Note that @f$ T^{-1}_t = T_{-t} @f$ may only be approximately true for some implementations.
	virtual void taxpy(const vector3<>& t, double alpha, const ScalarField& x, ScalarField& y) const=0;
};

//! Translation operator which works in real space using interpolating splines
class TranslationOperatorSpline : public TranslationOperator
{
public:
	//! Types of interpolating spline available for translation
	const enum SplineType
	{	Constant, //!< 0th order spline i.e. nearest neighbour interpolation
		Linear //!< 1st order spline i.e. linear interpolation
	} splineType;

	TranslationOperatorSpline(const GridInfo& gInfo, SplineType splineType);
	void taxpy(const vector3<>& t, double alpha, const ScalarField& x, ScalarField& y) const;
};

//! The exact translation operator in PW basis, although much slower and with potential ringing issues
class TranslationOperatorFourier : public TranslationOperator
{
public:
	TranslationOperatorFourier(const GridInfo& gInfo);
	void taxpy(const vector3<>& t, double alpha, const ScalarField& x, ScalarField& y) const;
};

#endif // JDFTX_FLUID_TRANSLATIONOPERATOR_H
