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
@brief Operators on #DataRptr's and #DataGptr's
*/

#include <core/Data.h>
#include <core/BlasExtra.h>
#include <core/GridInfo.h>
#include <core/LoopMacros.h>

//----------------- Real / complex conversion ------------------
DataRptr Real(const complexDataRptr&); //!< real part of a complex scalar field (real-space)
DataGptr Real(const complexDataGptr&); //!< real part of a complex scalar field (reciprocal space)
DataRptr Imag(const complexDataRptr&); //!< imaginary part of a complex scalar field (real-space)
DataGptr Imag(const complexDataGptr&); //!< imaginary part of a complex scalar field (reciprocal space)
complexDataRptr Complex(const DataRptr&); //!< convert real to complex scalar field with zero imaginary part (real-space)
complexDataGptr Complex(const DataGptr&); //!< convert real to complex scalar field with zero imaginary part (reciprocal-space)

//------------------------------ Linear Unary operators ------------------------------

DataGptr O(const DataGptr&); //!<Inner product operator (diagonal in PW basis)
DataGptr O(DataGptr&&); //!<Inner product operator (diagonal in PW basis)
complexDataGptr O(const complexDataGptr&); //!<Inner product operator (diagonal in PW basis)
complexDataGptr O(complexDataGptr&&); //!<Inner product operator (diagonal in PW basis)

DataRptr I(const DataGptr&, bool compat=false); //!< Forward transform: PW basis -> real space (preserve input)
DataRptr I(DataGptr&&, bool compat=false); //!< Forward transform: PW basis -> real space (destructible input)
complexDataRptr I(const complexDataGptr&); //!< Forward transform: PW basis -> real space (preserve input)
complexDataRptr I(complexDataGptr&&); //!< Forward transform: PW basis -> real space (destructible input)

DataGptr J(const DataRptr&); //!< Inverse transform: Real space -> PW basis
complexDataGptr J(const complexDataRptr&); //!< Inverse transform: Real space -> PW basis (preserve input)
complexDataGptr J(complexDataRptr&&); //!< Inverse transform: Real space -> PW basis (destructible input)

DataGptr Idag(const DataRptr&); //!< Forward transform transpose: Real space -> PW basis
complexDataGptr Idag(const complexDataRptr&); //!< Forward transform transpose: Real space -> PW basis (preserve input)
complexDataGptr Idag(complexDataRptr&&); //!< Forward transform transpose: Real space -> PW basis (destructible input)

DataRptr Jdag(const DataGptr&, bool compat=false); //!< Inverse transform transpose: PW basis -> real space (preserve input)
DataRptr Jdag(DataGptr&&, bool compat=false); //!< Inverse transform transpose: PW basis -> real space (destructible input)
complexDataRptr Jdag(const complexDataGptr&); //!< Inverse transform transpose: PW basis -> real space (preserve input)
complexDataRptr Jdag(complexDataGptr&&); //!< Inverse transform transpose: PW basis -> real space (destructible input)

DataRptr JdagOJ(const DataRptr&); //!< Evaluate Jdag(O(J())), which avoids 2 fourier transforms in PW basis (preserve input)
DataRptr JdagOJ(DataRptr&&); //!< Evaluate Jdag(O(J())), which avoids 2 fourier transforms in PW basis (destructible input)
complexDataRptr JdagOJ(const complexDataRptr&); //!< Evaluate Jdag(O(J())), which avoids 2 fourier transforms in PW basis (preserve input)
complexDataRptr JdagOJ(complexDataRptr&&); //!< Evaluate Jdag(O(J())), which avoids 2 fourier transforms in PW basis (destructible input)

DataGptr L(const DataGptr&); //!< Laplacian
DataGptr L(DataGptr&&); //!< Laplacian
complexDataGptr L(const complexDataGptr&); //!< Laplacian
complexDataGptr L(complexDataGptr&&); //!< Laplacian

DataGptr Linv(const DataGptr&); //!< Inverse Laplacian
DataGptr Linv(DataGptr&&); //!< Inverse Laplacian
complexDataGptr Linv(const complexDataGptr&); //!< Inverse Laplacian
complexDataGptr Linv(complexDataGptr&&); //!< Inverse Laplacian

void zeroNyquist(RealKernel& Gdata); //!< zeros out all the nyquist components of a real G-kernel
void zeroNyquist(DataGptr& Gptr); //!< zeros out all the nyquist components of a G-space data array
void zeroNyquist(DataRptr& Rptr); //!< zeros out all the nyquist components of an R-space data array

//------------------------------ Nonlinear Unary operators ------------------------------

DataRptr exp(const DataRptr&); //!< Elementwise exponential (preserve input)
DataRptr exp(DataRptr&&); //!< Elementwise exponential (destructible input)

DataRptr log(const DataRptr&); //!< Elementwise logarithm (preserve input)
DataRptr log(DataRptr&&); //!< Elementwise logarithm (destructible input)

DataRptr sqrt(const DataRptr&); //!< Elementwise square root (preserve input)
DataRptr sqrt(DataRptr&&); //!< Elementwise square root (destructible input)

DataRptr inv(const DataRptr&); //!< Elementwise reciprocal (preserve input)
DataRptr inv(DataRptr&&); //!< Elementwise reciprocal (destructible input)

DataRptr pow(const DataRptr&, double alpha); //!< Elementwise power (preserve input)
DataRptr pow(DataRptr&&, double alpha); //!< Elementwise power (destructible input)

#define Tptr std::shared_ptr<T> //!< shorthand for writing the template operators (undef'd at end of header)

template<class T> Tptr clone(const Tptr& X) { return X->clone(); } //!< Clone (NOTE: operator= is by reference for Data*ptr)

//------------------------------ Multiplication operators ------------------------------

template<class T> Tptr& operator*=(Tptr& in, double scaleFac) { in->scale *= scaleFac; return in; } //!< Scale
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
DataRptr& operator*=(DataRptr& in, const DataRptr& other); //!< Elementwise multiply for real data
template<class T> Tptr operator*(const Tptr& in1, const Tptr& in2) { Tptr out(in1->clone()); return out *= in2; } //!< Elementwise multiply (preserve inputs)
template<class T> Tptr operator*(const Tptr& in1, Tptr&& in2) { return in2 *= in1; } //!< Elementwise multiply (destructible input)
template<class T> Tptr operator*(Tptr&& in1, const Tptr& in2) { return in1 *= in2; } //!< Elementwise multiply (destructible input)
template<class T> Tptr operator*(Tptr&& in1, Tptr&& in2) { return in1 *= in2; } //!< Elementwise multiply (destructible inputs)

//Extra operators in R-space alone for mixed complex-real elementwise multiplications:
complexDataRptr& operator*=(complexDataRptr&, const DataRptr&); //!< elementwise multiply 
complexDataRptr operator*(const complexDataRptr&, const DataRptr&);//!< elementwise multiply (preserve inputs)
complexDataRptr operator*(const DataRptr&, const complexDataRptr&);//!< elementwise multiply (preserve inputs)
complexDataRptr operator*(complexDataRptr&&, const DataRptr&);//!< elementwise multiply (destructible inputs)
complexDataRptr operator*(const DataRptr&, complexDataRptr&&);//!< elementwise multiply (destructible inputs)

//Extra operators in G-space alone for real kernel multiplication:
DataGptr& operator*=(DataGptr&, const RealKernel&); //!< Elementwise multiply
DataGptr operator*(const RealKernel&, const DataGptr&); //!< Elementwise multiply (preserve inputs)
DataGptr operator*(const DataGptr&, const RealKernel&); //!< Elementwise multiply (preserve inputs)
DataGptr operator*(const RealKernel&, DataGptr&&); //!< Elementwise multiply (destructible input)
DataGptr operator*(DataGptr&&, const RealKernel&); //!< Elementwise multiply (destructible input)


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
void axpy(double alpha, const DataRptr& X, DataRptr& Y); //!< Real data Linear combine: Y += alpha * X (Note: null pointers are treated as zero)
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
DataRptr& operator+=(DataRptr&, double); //!< Increment by scalar
DataRptr operator+(double, const DataRptr&); //!< Add scalar (preserve inputs)
DataRptr operator+(const DataRptr&, double); //!< Add scalar (preserve inputs)
DataRptr operator+(double, DataRptr&&); //!< Add scalar (destructible input)
DataRptr operator+(DataRptr&&, double); //!< Add scalar (destructible input)
DataRptr& operator-=(DataRptr&, double); //!< Decrement by scalar
DataRptr operator-(double, const DataRptr&); //!< Subtract from scalar (preserve inputs)
DataRptr operator-(const DataRptr&, double); //!< Subtract scalar (preserve inputs)
DataRptr operator-(double, DataRptr&&); //!< Subtract from scalar (destructible input)
DataRptr operator-(DataRptr&&, double); //!< Subtract scalar (destructible input)


//------------------------------ Norms and dot products ------------------------------

template<typename T> complex dot(const Tptr& X, const Tptr& Y) //!< Generic inner product for complex types
{	return callPref(eblas_zdotc)(X->nElem, X->dataPref(), 1, Y->dataPref(), 1); 
}
template<typename T> double nrm2(const Tptr& X) //!< Generic 2-norm for complex types
{	return callPref(eblas_dznrm2)(X->nElem, X->dataPref(), 1); 
}
template<typename T> complex sum(const Tptr& X) //!< Generic sum for complex types
{	Data dataOne(X->gInfo, 1, 2, false); *((complex*)dataOne.data()) = 1.0;
	return callPref(eblas_zdotc)(X->nElem, (complex*)dataOne.dataPref(), 0, X->dataPref(), 1); 
}
//Special handling for real scalar fields:
double dot(const DataRptr&, const DataRptr&); //!< Inner product
double dot(const DataGptr&, const DataGptr&); //!< Inner product
double nrm2(const DataRptr&); //!< 2-norm
double nrm2(const DataGptr&); //!< 2-norm
double sum(const DataRptr&); //!< Sum of elements
double sum(const DataGptr&); //!< Equivalent to dot() with a DataGptr of all 1s (NOTE: sum(X) != sum(I(X)))

double integral(const DataRptr&); //!< Integral in the unit cell (dV times sum())
double integral(const DataGptr&); //!< Integral in the unit cell (just fetches the G=0 component with correct prefactor)
complex integral(const complexDataRptr&); //!< Integral in the unit cell (dV times sum())
complex integral(const complexDataGptr&); //!< Integral in the unit cell (just fetches the G=0 component with correct prefactor)

//------------------------------ Grid conversion utilities ------------------------------
DataGptr changeGrid(const DataGptr&, const GridInfo& gInfoNew); //Fourier up/down-sample to get to the new grid
DataRptr changeGrid(const DataRptr&, const GridInfo& gInfoNew); //Fourier up/down-sample to get to the new grid (real-space version)
complexDataGptr changeGrid(const complexDataGptr&, const GridInfo& gInfoNew); //Fourier up/down-sample to get to the new grid
complexDataRptr changeGrid(const complexDataRptr&, const GridInfo& gInfoNew); //Fourier up/down-sample to get to the new grid (real-space version)

//------------------------------ Initialization utilities ------------------------------

#include <string.h>

template<typename T> void initZero(Tptr& X) { X->zero(); }
template<typename T> void initZero(Tptr& X, const GridInfo& gInfo) { if(X) X->zero(); else nullToZero(X, gInfo); }
template<typename T> void nullToZero(Tptr& X, const GridInfo& gInfo) { if(!X) { X=T::alloc(gInfo,isGpuEnabled()); initZero(X); } } //!< If X is null, allocate and initialize to 0
void initRandom(DataRptr&, double cap=0.0); //!< initialize element-wise with a unit-normal random number (with a cap if cap>0)
void initRandomFlat(DataRptr&); //!< initialize element-wise with a unit-flat [0:1) random number
void initGaussianKernel(RealKernel&, double x0); //!< initialize to gaussian kernel exp(-(G x0/2)^2)
void initTranslation(DataGptr&, const vector3<>& r);  //!< initialize to translation operator exp(-i G.r)
DataGptr gaussConvolve(const DataGptr&, double sigma); //!< convolve with a gaussian
DataGptr gaussConvolve(DataGptr&&, double sigma); //!< convolve with a gaussian (destructible input)

//! Evaluate a function f(i, Gsq, args...) at each point in reciprocal space indexed by i
template<typename Func, typename... Args> void applyFuncGsq(const GridInfo& gInfo, const Func& f, Args... args);

//! Evaluate a function f(i, r, args...) at each point in real space index by i
template<typename Func, typename... Args> void applyFunc_r(const GridInfo& gInfo, const Func& f, Args... args);

//------------------------------ Debug utilities ------------------------------

void printStats(const DataRptr& X, const char* name, FILE* fp=stdout); //!< Print mean and standard deviation of array with specified name (debug utility)


//! Check the symmetry of a linear operator
//! @tparam Callable An operator with function call signature Vec Callable(const Vec&)
//! @tparam Vec Any operand type representing an element of a vector space
template<typename Callable, typename Vec> void checkSymmetry(Callable* func, const Vec& v1, const Vec& v2, const char* funcName)
{	double dot1 = dot(v1, (*func)(v2));
	double dot2 = dot(v2, (*func)(v1));
	double dotDiff = fabs(dot1-dot2);
	printf("Relative error in symmetry of %s: %le\n", funcName, dotDiff/sqrt(fabs(dot1)*fabs(dot2)));
}

//! @}

#undef Tptr

//!@cond


template<typename Func, typename... Args>
void applyFuncGsq_sub(size_t iStart, size_t iStop, const vector3<int> S, const matrix3<> GGT, const Func* f, Args... args)
{	THREAD_halfGspaceLoop( (*f)(i, GGT.metric_length_squared(iG), args...); )
}
template<typename Func, typename... Args> void applyFuncGsq(const GridInfo& gInfo, const Func& f, Args... args)
{	threadLaunch(shouldThreadOperators() ? 0 : 1, //0 => max threads
		applyFuncGsq_sub<Func,Args...>, gInfo.nG, gInfo.S, gInfo.GGT, &f, args...);
}

template<typename Func, typename... Args>
void applyFunc_r_sub(size_t iStart, size_t iStop, const vector3<int> S, const vector3<> h[3], const Func* f, Args... args)
{	THREAD_rLoop
	(	vector3<> ri = iv[0]*h[0] + iv[1]*h[1] + iv[2]*h[2];
		(*f)(i, ri, args...);
	)
}
template<typename Func, typename... Args> void applyFunc_r(const GridInfo& gInfo, const Func& f, Args... args)
{	threadLaunch(shouldThreadOperators() ? 0 : 1, //0 => max threads
		applyFunc_r_sub<Func,Args...>, gInfo.nr, gInfo.S, gInfo.h, f, args...);
}
//!@endcond
#endif //JDFTX_CORE_OPERATORS_H
