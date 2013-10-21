/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

This file is part of Fluid1D.

Fluid1D is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Fluid1D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Fluid1D.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#ifndef FLUID1D_FLUID1D_TRANSLATIONOPERATOR_H
#define FLUID1D_FLUID1D_TRANSLATIONOPERATOR_H

//! @file TranslationOperator.h
//! Options for translating site densities for rigid molecule ideal gas implementations

#include <core1D/GridInfo.h>
#include <core1D/Data.h>
#include <core/vector3.h>
#include <memory>
#include <map>
#include <thread>

//! Abstract base-class for reduced-dimensionality generalized translation operator
//! Implements application of operator @f$ S_{\vec{a}(\vec{r})} @f$ defined by
//! @f$ S_{\vec{a}} f(\vec{r}) = f(\vec{r}+\vec{a}(\vec{r})) @f$,
//! and its hermitian conjugate
class TranslationOperator
{
public:
	//! Perform @f$ y += \alpha S_{\vec{a}}(x) @f$
	virtual void S_axpy(const vector3<>& a, double alpha, const ScalarField& x, ScalarField& y) const=0;

	//! Perform @f$ y += \alpha S^{\dag}_{\vec{a}}(x) @f$
	virtual void Sdag_axpy(const vector3<>& a, double alpha, const ScalarField& x, ScalarField& y) const=0;
};


//! Implementation of TranslationOperator based on linear-spline sampling
class TranslationOperatorLspline : public TranslationOperator
{
public:
	const GridInfo& gInfo;
	TranslationOperatorLspline(const GridInfo& gInfo);
	void S_axpy(const vector3<>& a, double alpha, const ScalarField& x, ScalarField& y) const;
	void Sdag_axpy(const vector3<>& a, double alpha, const ScalarField& x, ScalarField& y) const;

private:
	std::map< vector3<>, std::shared_ptr<struct LsplineMatrix> > lsplineMatrix; //!< precomputed matrix list
	std::mutex lsplineMutex; //!< for thread safety
	std::shared_ptr<struct LsplineMatrix> getMatrix(const vector3<>& a) const;
};

#endif // FLUID1D_FLUID1D_TRANSLATIONOPERATOR_H
