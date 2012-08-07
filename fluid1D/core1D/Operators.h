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

#ifndef FLUID1D_CORE1D_OPERATORS_H
#define FLUID1D_CORE1D_OPERATORS_H

#include <core1D/Data.h>
#include <core/scaled.h>


//------------------------------ Linear Unary operators ------------------------------

ScalarFieldTilde O(const ScalarFieldTilde&); //!<Inner product operator
ScalarFieldTilde O(ScalarFieldTilde&&); //!<Inner product operator

ScalarField I(const ScalarFieldTilde&); //!< Forward transform: basis -> real space
ScalarField ID(const ScalarFieldTilde&); //!< Forward transform of vector derivative: basis -> real space
ScalarField IDD(const ScalarFieldTilde&); //!< Forward transform of tensor second derivative: basis -> real space
ScalarField Jdag(const ScalarFieldTilde&); //!< Transpose of inverse transform transpose: basis -> real space

ScalarFieldTilde J(const ScalarField&); //!< Inverse transform: Real space -> basis
ScalarFieldTilde Idag(const ScalarField&); //!< Transpose of forward transform: Real space -> basis
ScalarFieldTilde IDdag(const ScalarField&); //!< Transpose of forward transform of vector derivative: Real space -> basis
ScalarFieldTilde IDDdag(const ScalarField&); //!< Transpose of forward transform of tensor second derivative: Real space -> basis

ScalarField JdagOJ(const ScalarField&); //!< Evaluate Jdag(O(J())), which avoids 2 fourier transforms in PW-like bases
ScalarField JdagOJ(ScalarField&&); //!< Evaluate Jdag(O(J())), which avoids 2 fourier transforms in PW-like bases

ScalarFieldTilde L(const ScalarFieldTilde&); //!< Laplacian
ScalarFieldTilde L(ScalarFieldTilde&&); //!< Laplacian
ScalarFieldTilde Linv(const ScalarFieldTilde&); //!< Inverse Laplacian
ScalarFieldTilde Linv(ScalarFieldTilde&&); //!< Inverse Laplacian

//------------------------------ Nonlinear Unary operators ------------------------------

ScalarField exp(const ScalarField&); //!< Elementwise exponential
ScalarField log(const ScalarField&); //!< Elementwise logarithm
ScalarField sqrt(const ScalarField&); //!< Elementwise square root
ScalarField inv(const ScalarField&); //!< Elementwise reciprocal
ScalarField pow(const ScalarField&, double alpha); //!< Elementwise power

//------------------------------ Multiplication operators ------------------------------

//scale and unary-
ScalarField& operator*=(ScalarField& X, double s);
scaled<ScalarField> operator*(double s, const ScalarField &Y);
scaled<ScalarField> operator*(const ScalarField &Y, double s);
scaled<ScalarField> operator-(const ScalarField &Y);
ScalarFieldTilde& operator*=(ScalarFieldTilde& X, double s);
scaled<ScalarFieldTilde> operator*(double s, const ScalarFieldTilde &Y);
scaled<ScalarFieldTilde> operator*(const ScalarFieldTilde &Y, double s);
scaled<ScalarFieldTilde> operator-(const ScalarFieldTilde &Y);

//Convolution:
ScalarFieldTilde operator*(const SphericalKernel&, const ScalarFieldTilde&); //!< Convolution (preserve input)
ScalarFieldTilde operator*(const SphericalKernel&, ScalarFieldTilde&&); //!< Convolution (destructible input)

//------------------------------ Linear combine operators ------------------------------

ScalarField& operator+=(ScalarField& Y, const ScalarField &X);
ScalarField& operator-=(ScalarField& Y, const ScalarField &X);
ScalarField operator+(const ScalarField &Y1, const ScalarField &Y2);
ScalarField operator-(const ScalarField &Y1,const ScalarField &Y2);
ScalarFieldTilde& operator+=(ScalarFieldTilde& Y, const ScalarFieldTilde &X);
ScalarFieldTilde& operator-=(ScalarFieldTilde& Y, const ScalarFieldTilde &X);
ScalarFieldTilde operator+(const ScalarFieldTilde &Y1, const ScalarFieldTilde &Y2);
ScalarFieldTilde operator-(const ScalarFieldTilde &Y1,const ScalarFieldTilde &Y2);

//Real-space scalar additions:
ScalarField& operator+=(ScalarField&, double); //!< Increment by scalar
ScalarField operator+(double, const ScalarField&); //!< Add scalar (preserve inputs)
ScalarField operator+(const ScalarField&, double); //!< Add scalar (preserve inputs)
ScalarField& operator-=(ScalarField&, double); //!< Decrement by scalar
ScalarField operator-(double, const ScalarField&); //!< Subtract from scalar (preserve inputs)
ScalarField operator-(const ScalarField&, double); //!< Subtract scalar (preserve inputs)

//------------------------- Initialization utilties ---------------------------------

void initRandom(ScalarField&, double cap=0.0); //!< initialize element-wise with a unit-normal random number (with a cap if cap>0)
void initRandomFlat(ScalarField&); //!< initialize element-wise with a unit-flat [0:1) random number


#endif // FLUID1D_CORE1D_OPERATORS_H
