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

#ifndef JDFTX_CORE_OPERATORS_H
#define JDFTX_CORE_OPERATORS_H

//! @addtogroup griddata
//! @{

/** @file Operators.h
@brief Operators on #ScalarField's and #ScalarFieldTilde's
*/

#include <core/ScalarField.h>
#include <core/BlasExtra.h>
#include <core/GridInfo.h>
#include <core/LoopMacros.h>

//----------------- Real / complex conversion ------------------
ScalarField Real(const complexScalarField&); //!< real part of a complex scalar field (real-space)
ScalarFieldTilde Real(const complexScalarFieldTilde&); //!< real part of a complex scalar field (reciprocal space)
ScalarField Imag(const complexScalarField&); //!< imaginary part of a complex scalar field (real-space)
ScalarFieldTilde Imag(const complexScalarFieldTilde&); //!< imaginary part of a complex scalar field (reciprocal space)
complexScalarField Complex(const ScalarField&); //!< convert real to complex scalar field with zero imaginary part (real-space)
complexScalarField Complex(const ScalarField& re, const ScalarField& im); //!< construct complex scalar field fromr eal and imaginary parts (real-space)
complexScalarFieldTilde Complex(const ScalarFieldTilde&); //!< convert real to complex scalar field with zero imaginary part (reciprocal-space)

//------------------------------ Linear Unary operators ------------------------------

ScalarFieldTilde O(const ScalarFieldTilde&); //!<Inner product operator (diagonal in PW basis)
ScalarFieldTilde O(ScalarFieldTilde&&); //!<Inner product operator (diagonal in PW basis)
complexScalarFieldTilde O(const complexScalarFieldTilde&); //!<Inner product operator (diagonal in PW basis)
complexScalarFieldTilde O(complexScalarFieldTilde&&); //!<Inner product operator (diagonal in PW basis)

//Transform operators:
//  compat: force GPU transform output to be FFTW compatible (affects only how redundant Nyquist frequency components in the C->R transforms are handled)
//  nThreads: maximum number of threads to use; use to prevent thread nesting when performing transforms of vector fields etc.

ScalarField I(const ScalarFieldTilde&, bool compat=false, int nThreads=0); //!< Forward transform: PW basis -> real space (preserve input)
ScalarField I(ScalarFieldTilde&&, bool compat=false, int nThreads=0); //!< Forward transform: PW basis -> real space (destructible input)
complexScalarField I(const complexScalarFieldTilde&, int nThreads=0); //!< Forward transform: PW basis -> real space (preserve input)
complexScalarField I(complexScalarFieldTilde&&, int nThreads=0); //!< Forward transform: PW basis -> real space (destructible input)

ScalarFieldTilde J(const ScalarField&, int nThreads=0); //!< Inverse transform: Real space -> PW basis
complexScalarFieldTilde J(const complexScalarField&, int nThreads=0); //!< Inverse transform: Real space -> PW basis (preserve input)
complexScalarFieldTilde J(complexScalarField&&, int nThreads=0); //!< Inverse transform: Real space -> PW basis (destructible input)

ScalarFieldTilde Idag(const ScalarField&, int nThreads=0); //!< Forward transform transpose: Real space -> PW basis
complexScalarFieldTilde Idag(const complexScalarField&, int nThreads=0); //!< Forward transform transpose: Real space -> PW basis (preserve input)
complexScalarFieldTilde Idag(complexScalarField&&, int nThreads=0); //!< Forward transform transpose: Real space -> PW basis (destructible input)

ScalarField Jdag(const ScalarFieldTilde&, bool compat=false, int nThreads=0); //!< Inverse transform transpose: PW basis -> real space (preserve input)
ScalarField Jdag(ScalarFieldTilde&&, bool compat=false, int nThreads=0); //!< Inverse transform transpose: PW basis -> real space (destructible input)
complexScalarField Jdag(const complexScalarFieldTilde&, int nThreads=0); //!< Inverse transform transpose: PW basis -> real space (preserve input)
complexScalarField Jdag(complexScalarFieldTilde&&, int nThreads=0); //!< Inverse transform transpose: PW basis -> real space (destructible input)

ScalarField JdagOJ(const ScalarField&); //!< Evaluate Jdag(O(J())), which avoids 2 fourier transforms in PW basis (preserve input)
ScalarField JdagOJ(ScalarField&&); //!< Evaluate Jdag(O(J())), which avoids 2 fourier transforms in PW basis (destructible input)
complexScalarField JdagOJ(const complexScalarField&); //!< Evaluate Jdag(O(J())), which avoids 2 fourier transforms in PW basis (preserve input)
complexScalarField JdagOJ(complexScalarField&&); //!< Evaluate Jdag(O(J())), which avoids 2 fourier transforms in PW basis (destructible input)

ScalarFieldTilde L(const ScalarFieldTilde&); //!< Laplacian
ScalarFieldTilde L(ScalarFieldTilde&&); //!< Laplacian
complexScalarFieldTilde L(const complexScalarFieldTilde&); //!< Laplacian
complexScalarFieldTilde L(complexScalarFieldTilde&&); //!< Laplacian

ScalarFieldTilde Linv(const ScalarFieldTilde&); //!< Inverse Laplacian
ScalarFieldTilde Linv(ScalarFieldTilde&&); //!< Inverse Laplacian
complexScalarFieldTilde Linv(const complexScalarFieldTilde&); //!< Inverse Laplacian
complexScalarFieldTilde Linv(complexScalarFieldTilde&&); //!< Inverse Laplacian

void zeroNyquist(RealKernel& Gdata); //!< zeros out all the nyquist components of a real G-kernel
void zeroNyquist(ScalarFieldTilde& Gptr); //!< zeros out all the nyquist components of a G-space data array
void zeroNyquist(ScalarField& Rptr); //!< zeros out all the nyquist components of an R-space data array

//------------------------------ Nonlinear Unary operators ------------------------------

ScalarField exp(const ScalarField&); //!< Elementwise exponential (preserve input)
ScalarField exp(ScalarField&&); //!< Elementwise exponential (destructible input)

ScalarField log(const ScalarField&); //!< Elementwise logarithm (preserve input)
ScalarField log(ScalarField&&); //!< Elementwise logarithm (destructible input)

ScalarField sqrt(const ScalarField&); //!< Elementwise square root (preserve input)
ScalarField sqrt(ScalarField&&); //!< Elementwise square root (destructible input)

ScalarField inv(const ScalarField&); //!< Elementwise reciprocal (preserve input)
ScalarField inv(ScalarField&&); //!< Elementwise reciprocal (destructible input)

ScalarField pow(const ScalarField&, double alpha); //!< Elementwise power (preserve input)
ScalarField pow(ScalarField&&, double alpha); //!< Elementwise power (destructible input)

#define Tptr std::shared_ptr<T> //!< shorthand for writing the template operators (undef'd at end of header)

template<class T> Tptr clone(const Tptr& X) { if(X) return X->clone(); else return 0; } //!< Clone (NOTE: operator= is by reference for the ScalarField classes)

//------------------------------ Multiplication operators ------------------------------

template<class T> Tptr& operator*=(Tptr& in, double scaleFac) { if(in) in->scale *= scaleFac; return in; } //!< Scale
template<class T> Tptr operator*(const Tptr& in, double scaleFac) { Tptr out(in->clone()); return out *= scaleFac; } //!< Scalar multiply (preserve input)
template<class T> Tptr operator*(double scaleFac, const Tptr& in) { Tptr out(in->clone()); return out *= scaleFac; } //!< Scalar multiply (preserve input)
template<class T> Tptr operator*(Tptr&& in, double scaleFac) { return in *= scaleFac; } //!< Scalar multiply (destructible input)
template<class T> Tptr operator*(double scaleFac, Tptr&& in) { return in *= scaleFac; } //!< Scalar multiply (destructible input)

//! Generic elementwise conjugate for complex data:
template<class T> Tptr conj(Tptr&& in)
{	callPref(eblas_dscal)(in->nElem, -1., ((double*)in->dataPref(false))+1, 2); //negate the imaginary parts
	return in;
}
template<class T> Tptr conj(const Tptr& in) { return conj(clone(in)); }

//! Generic elementwise multiply for complex data:
template<class T> Tptr& operator*=(Tptr& in, const Tptr& other)
{	in->scale *= other->scale;
	callPref(eblas_zmul)(in->nElem, other->dataPref(false), 1, in->dataPref(false), 1);
	return in;
}
ScalarField& operator*=(ScalarField& in, const ScalarField& other); //!< Elementwise multiply for real data
template<class T> Tptr operator*(const Tptr& in1, const Tptr& in2) { Tptr out(in1->clone()); return out *= in2; } //!< Elementwise multiply (preserve inputs)
template<class T> Tptr operator*(const Tptr& in1, Tptr&& in2) { return in2 *= in1; } //!< Elementwise multiply (destructible input)
template<class T> Tptr operator*(Tptr&& in1, const Tptr& in2) { return in1 *= in2; } //!< Elementwise multiply (destructible input)
template<class T> Tptr operator*(Tptr&& in1, Tptr&& in2) { return in1 *= in2; } //!< Elementwise multiply (destructible inputs)

//Extra operators in R-space alone for mixed complex-real elementwise multiplications:
complexScalarField& operator*=(complexScalarField&, const ScalarField&); //!< elementwise multiply 
complexScalarField operator*(const complexScalarField&, const ScalarField&);//!< elementwise multiply (preserve inputs)
complexScalarField operator*(const ScalarField&, const complexScalarField&);//!< elementwise multiply (preserve inputs)
complexScalarField operator*(complexScalarField&&, const ScalarField&);//!< elementwise multiply (destructible inputs)
complexScalarField operator*(const ScalarField&, complexScalarField&&);//!< elementwise multiply (destructible inputs)

//Extra operators in G-space alone for real kernel multiplication:
ScalarFieldTilde& operator*=(ScalarFieldTilde&, const RealKernel&); //!< Elementwise multiply
ScalarFieldTilde operator*(const RealKernel&, const ScalarFieldTilde&); //!< Elementwise multiply (preserve inputs)
ScalarFieldTilde operator*(const ScalarFieldTilde&, const RealKernel&); //!< Elementwise multiply (preserve inputs)
ScalarFieldTilde operator*(const RealKernel&, ScalarFieldTilde&&); //!< Elementwise multiply (destructible input)
ScalarFieldTilde operator*(ScalarFieldTilde&&, const RealKernel&); //!< Elementwise multiply (destructible input)


//------------------------------ Linear combine operators ------------------------------

//!Generic axpy for complex data types (Note: null pointers are treated as zero)
template<typename T> void axpy(double alpha, const Tptr& X, Tptr& Y)
{	if(X)
	{	if(Y)
		{	if(Y->scale == 0.0) { Y = X * alpha; }
			else callPref(eblas_zaxpy)(X->nElem, alpha*X->scale/Y->scale, X->dataPref(false), 1, Y->dataPref(false), 1);
		}
		else Y = X * alpha;
	}
	//if X is null, nothing needs to be done, Y remains unchanged
}
void axpy(double alpha, const ScalarField& X, ScalarField& Y); //!< Real data Linear combine: Y += alpha * X (Note: null pointers are treated as zero)
template<class T> Tptr& operator+=(Tptr& in, const Tptr& other) { axpy(+1.0, other, in); return in; } //!< Increment
template<class T> Tptr& operator-=(Tptr& in, const Tptr& other) { axpy(-1.0, other, in); return in; } //!< Decrement
template<class T> Tptr operator+(const Tptr& in1, const Tptr& in2) { Tptr out(in1->clone()); return out += in2; } //!< Add (preserve inputs)
template<class T> Tptr operator+(const Tptr& in1, Tptr&& in2) { return in2 += in1; } //!< Add (destructible input)
template<class T> Tptr operator+(Tptr&& in1, const Tptr& in2) { return in1 += in2; } //!< Add (destructible input)
template<class T> Tptr operator+(Tptr&& in1, Tptr&& in2) { return in1 += in2; } //!< Add (destructible inputs)
template<class T> Tptr operator-(const Tptr& in1, const Tptr& in2) { Tptr out(in1->clone()); return out -= in2; } //!< Subtract (preserve inputs)
template<class T> Tptr operator-(const Tptr& in1, Tptr&& in2) { return (in2 -= in1) *= -1.0; } //!< Subtract (destructible input)
template<class T> Tptr operator-(Tptr&& in1, const Tptr& in2) { return in1 -= in2; } //!< Subtract (destructible input)
template<class T> Tptr operator-(Tptr&& in1, Tptr&& in2) { return in1 -= in2; } //!< Subtract (destructible inputs)
template<class T> Tptr operator-(const Tptr& in) { return (-1.0)*in; } //!< Negate
template<class T> Tptr operator-(Tptr&& in) { return in*=(-1.0); } //!< Negate
//Extra operators in R-space alone for scalar additions:
ScalarField& operator+=(ScalarField&, double); //!< Increment by scalar
ScalarField operator+(double, const ScalarField&); //!< Add scalar (preserve inputs)
ScalarField operator+(const ScalarField&, double); //!< Add scalar (preserve inputs)
ScalarField operator+(double, ScalarField&&); //!< Add scalar (destructible input)
ScalarField operator+(ScalarField&&, double); //!< Add scalar (destructible input)
ScalarField& operator-=(ScalarField&, double); //!< Decrement by scalar
ScalarField operator-(double, const ScalarField&); //!< Subtract from scalar (preserve inputs)
ScalarField operator-(const ScalarField&, double); //!< Subtract scalar (preserve inputs)
ScalarField operator-(double, ScalarField&&); //!< Subtract from scalar (destructible input)
ScalarField operator-(ScalarField&&, double); //!< Subtract scalar (destructible input)


//------------------------------ Norms and dot products ------------------------------

template<typename T> complex dot(const Tptr& X, const Tptr& Y) //!< Generic inner product for complex types
{	return callPref(eblas_zdotc)(X->nElem, X->dataPref(), 1, Y->dataPref(), 1); 
}
template<typename T> double nrm2(const Tptr& X) //!< Generic 2-norm for complex types
{	return callPref(eblas_dznrm2)(X->nElem, X->dataPref(), 1); 
}
template<typename T> complex sum(const Tptr& X) //!< Generic sum for complex types
{	FieldData dataOne(X->gInfo, "complexScalarField", 1, 2, false); *((complex*)dataOne.data()) = 1.0;
	return callPref(eblas_zdotc)(X->nElem, (complex*)dataOne.dataPref(), 0, X->dataPref(), 1); 
}
//Special handling for real scalar fields:
double dot(const ScalarField&, const ScalarField&); //!< Inner product
double dot(const ScalarFieldTilde&, const ScalarFieldTilde&); //!< Inner product
double nrm2(const ScalarField&); //!< 2-norm
double nrm2(const ScalarFieldTilde&); //!< 2-norm
double sum(const ScalarField&); //!< Sum of elements
double sum(const ScalarFieldTilde&); //!< Equivalent to dot() with a ScalarFieldTilde of all 1s (NOTE: sum(X) != sum(I(X)))

double integral(const ScalarField&); //!< Integral in the unit cell (dV times sum())
double integral(const ScalarFieldTilde&); //!< Integral in the unit cell (just fetches the G=0 component with correct prefactor)
complex integral(const complexScalarField&); //!< Integral in the unit cell (dV times sum())
complex integral(const complexScalarFieldTilde&); //!< Integral in the unit cell (just fetches the G=0 component with correct prefactor)

//------------------------------ Grid conversion utilities ------------------------------
ScalarFieldTilde changeGrid(const ScalarFieldTilde&, const GridInfo& gInfoNew); //Fourier up/down-sample to get to the new grid
ScalarField changeGrid(const ScalarField&, const GridInfo& gInfoNew); //Fourier up/down-sample to get to the new grid (real-space version)
complexScalarFieldTilde changeGrid(const complexScalarFieldTilde&, const GridInfo& gInfoNew); //Fourier up/down-sample to get to the new grid
complexScalarField changeGrid(const complexScalarField&, const GridInfo& gInfoNew); //Fourier up/down-sample to get to the new grid (real-space version)

//------------------------------ Initialization utilities ------------------------------

#include <string.h>

template<typename T> void initZero(Tptr& X) { X->zero(); }
template<typename T> void initZero(Tptr& X, const GridInfo& gInfo) { if(X) X->zero(); else nullToZero(X, gInfo); }
template<typename T> void nullToZero(Tptr& X, const GridInfo& gInfo) { if(!X) { X=T::alloc(gInfo,isGpuEnabled()); initZero(X); } } //!< If X is null, allocate and initialize to 0
void initRandom(ScalarField&, double cap=0.0); //!< initialize element-wise with a unit-normal random number (with a cap if cap>0)
void initRandomFlat(ScalarField&); //!< initialize element-wise with a unit-flat [0:1) random number
void initGaussianKernel(RealKernel&, double x0); //!< initialize to gaussian kernel exp(-(G x0/2)^2)
void initTranslation(ScalarFieldTilde&, const vector3<>& r);  //!< initialize to translation operator exp(-i G.r)
ScalarFieldTilde gaussConvolve(const ScalarFieldTilde&, double sigma); //!< convolve with a gaussian
ScalarFieldTilde gaussConvolve(ScalarFieldTilde&&, double sigma); //!< convolve with a gaussian (destructible input)

//! Evaluate a function f(i, Gsq, args...) at each point in reciprocal space indexed by i
template<typename Func, typename... Args> void applyFuncGsq(const GridInfo& gInfo, const Func& f, Args... args);

//! Evaluate a function f(i, r, args...) at each point in real space index by i
template<typename Func, typename... Args> void applyFunc_r(const GridInfo& gInfo, const Func& f, Args... args);

//------------------------------ Debug utilities ------------------------------

void printStats(const ScalarField& X, const char* name, FILE* fp=stdout); //!< Print mean and standard deviation of array with specified name (debug utility)


//! Check the symmetry of a linear operator
//! @tparam Callable An operator with function call signature Vec Callable(const Vec&)
//! @tparam Vec Any operand type representing an element of a vector space
template<typename Callable, typename Vec> void checkSymmetry(Callable* func, const Vec& v1, const Vec& v2, const char* funcName)
{	double dot1 = dot(v1, (*func)(v2));
	double dot2 = dot(v2, (*func)(v1));
	double dotDiff = fabs(dot1-dot2);
	logPrintf("Relative error in symmetry of %s: %le\n", funcName, dotDiff/sqrt(fabs(dot1)*fabs(dot2)));
}

//! @}

#undef Tptr

//!@cond


template<typename Func, typename... Args>
void applyFuncGsq_sub(size_t iStart, size_t iStop, const vector3<int> S, const matrix3<> GGT, const Func* f, Args... args)
{	THREAD_halfGspaceLoop( (*f)(i, GGT.metric_length_squared(iG), args...); )
}
template<typename Func, typename... Args> void applyFuncGsq(const GridInfo& gInfo, const Func& f, Args... args)
{	threadLaunch(applyFuncGsq_sub<Func,Args...>, gInfo.nG, gInfo.S, gInfo.GGT, &f, args...);
}

template<typename Func, typename... Args>
void applyFunc_r_sub(size_t iStart, size_t iStop, const vector3<int> S, const vector3<> h[3], const Func* f, Args... args)
{	THREAD_rLoop
	(	vector3<> ri = iv[0]*h[0] + iv[1]*h[1] + iv[2]*h[2];
		(*f)(i, ri, args...);
	)
}
template<typename Func, typename... Args> void applyFunc_r(const GridInfo& gInfo, const Func& f, Args... args)
{	threadLaunch(applyFunc_r_sub<Func,Args...>, gInfo.nr, gInfo.S, gInfo.h, f, args...);
}
//!@endcond
#endif //JDFTX_CORE_OPERATORS_H
