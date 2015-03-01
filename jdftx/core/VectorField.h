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

#ifndef JDFTX_CORE_VECTORFIELD_H
#define JDFTX_CORE_VECTORFIELD_H

//! @addtogroup griddata
//! @{

/** @file VectorField.h
@brief Generic multiplet of data arrays (and specialized to triplets for vector fields in real/reciprocal space)
*/

#include <core/ScalarField.h>
#include <core/GridInfo.h>
#include <core/Operators.h>
#include <core/vector3.h>
#include <vector>

#define Tptr std::shared_ptr<T> //!< shorthand for writing the template operators (undef'd at end of header)
#define TptrMul ScalarFieldMultiplet<T,N> //!< shorthand for the template operators/functions (undef'd at end of file)
#define RptrMul ScalarFieldMultiplet<ScalarFieldData,N> //!< shorthand for real-space-only template operators/functions (undef'd at end of file)
#define GptrMul ScalarFieldMultiplet<ScalarFieldTildeData,N> //!< shorthand for reciprocal-space-only template operators/functions (undef'd at end of file)

//!@cond
#define Nloop(code) for(int i=0; i<N; i++) {code} //loop over components (used repeatedly below)
//!@endcond


//! @brief Generic multiplet object with overloaded arithmetic
//! @tparam T ScalarFieldData or ScalarFieldTildeData
//! @tparam N Number of elements in multiplet
template<class T, int N> struct ScalarFieldMultiplet
{	std::vector<Tptr> component; //!< the array of components (also accessible via #operator[])

	//! @brief Construct multiplet from an array of data sets (or default: initialize to null)
	//! @param in Pointer to array, or null to initialize each component to null
	ScalarFieldMultiplet(const Tptr* in=0) : component(N) { Nloop( component[i] = (in ? in[i] : 0); ) }

	//! @brief Construct a multiplet with allocated data
	//! @param gInfo Simulation grid info / memory manager to use to allocate the data
	ScalarFieldMultiplet(const GridInfo& gInfo, bool onGpu=false) : component(N) { Nloop( component[i] = Tptr(T::alloc(gInfo,onGpu)); ) }

	Tptr& operator[](int i) { return component[i]; } //!< Retrieve a reference to the i'th component (no bound checks)
	const Tptr& operator[](int i) const { return component[i]; } //!< Retrieve a const reference to the i'th component (no bound checks)
	ScalarFieldMultiplet clone() const { TptrMul out; Nloop( out[i] = component[i]->clone(); ) return out; } //!< Clone data (note assignment will be reference for the actual data)

	std::vector<typename T::DataType*> data(); //!< Get the component data pointers in an std::vector
	std::vector<const typename T::DataType*> data() const; //!< Get the component data pointers in an std::vector (const version)
	std::vector<const typename T::DataType*> const_data() const { return data(); } //!< Get the component data pointers in an std::vector (const version)
	#ifdef GPU_ENABLED
	std::vector<typename T::DataType*> dataGpu(); //!< Get the component GPU data pointers in an std::vector
	std::vector<const typename T::DataType*> dataGpu() const; //!< Get the component GPU data pointers in an std::vector (const version)
	std::vector<const typename T::DataType*> const_dataGpu() const { return dataGpu(); } //!< Get the component GPU data pointers in an std::vector (const version)
	#endif

	//Preferred data access (on GPU if available, on CPU otherwise)
	#ifdef GPU_ENABLED
	std::vector<typename T::DataType*> dataPref() { return dataGpu(); }
	std::vector<const typename T::DataType*> dataPref() const { return dataGpu(); }
	std::vector<const typename T::DataType*> const_dataPref() const { return dataGpu(); }
	#else
	std::vector<typename T::DataType*> dataPref()  { return data(); }
	std::vector<const typename T::DataType*> dataPref() const { return data(); }
	std::vector<const typename T::DataType*> const_dataPref() const { return data(); }
	#endif

	
	explicit operator bool() const { bool ret=true; Nloop(ret = ret && component[i];) return ret; } //!< Cast to bool: true if all components are non-null
	void loadFromFile(const char* fileName); //!< Load all components from a single binary file
	void saveToFile(const char* fileName) const; //!< Save all components from a single binary file
};
typedef ScalarFieldMultiplet<ScalarFieldData,3> VectorField; //!< Real space vector field
typedef ScalarFieldMultiplet<ScalarFieldTildeData,3> VectorFieldTilde; //!< Reciprocal space vector field
typedef ScalarFieldMultiplet<ScalarFieldData,5> TensorField; //!< Symmetric traceless tensor: real space field
typedef ScalarFieldMultiplet<ScalarFieldTildeData,5> TensorFieldTilde; //!< Symmetric traceless tensor: reciprocal space field

//Allocation/init/copy etc:
template<class T,int N> TptrMul clone(const TptrMul& X) { return X ? X.clone() : (TptrMul)0; } //!< Clone (NOTE: operator= is by reference for ScalarField multiplets)
template<class T,int N> void initZero(TptrMul& X) { Nloop( initZero(X[i]); ) } //!< Initialize data to 0 and scale factors to 1
template<class T,int N> void nullToZero(TptrMul& X, const GridInfo& gInfo) { Nloop(nullToZero(X[i],gInfo);) } //!< Allocate and initialize each component of X to 0 if null
template<int N> void initRandom(RptrMul& X, double cap=0.0) {Nloop(initRandom(X[i],cap);)} //!< initialize element-wise with a unit-normal random number (with a cap if cap>0)
template<int N> void initRandomFlat(RptrMul& X) {Nloop(initRandomFlat(X[i]);)} //!< initialize element-wise with a unit-flat [0:1) random number
template<int N> void randomize(RptrMul& X) { initRandom(X, 3.); } //!< alternate interface required by Minimizable

//Multiplication of commensurate multiplets
template<class T,int N> TptrMul& operator*=(TptrMul& in, const TptrMul& other) { Nloop(in[i]*=other[i];) return in; } //!< Elementwise multiply each component
template<class T,int N> TptrMul operator*(const TptrMul& in1, const TptrMul& in2) { TptrMul out(in1.clone()); return out *= in2; } //!< Elementwise multiply each component (preserve inputs)
template<class T,int N> TptrMul operator*(TptrMul&& in1, const TptrMul& in2) { return in1 *= in2; } //!< Elementwise multiply each component (destructible input)
template<class T,int N> TptrMul operator*(const TptrMul& in2, TptrMul&& in1) { return in2 *= in1; } //!< Elementwise multiply each component (destructible input)
template<class T,int N> TptrMul operator*(TptrMul&& in1, TptrMul&& in2) { return in1 *= in2; } //!< Elementwise multiply each component (destructible inputs)

//Multiplication of multiplets with singlets
template<class T,int N> TptrMul& operator*=(TptrMul& inM, const Tptr& inS) { Nloop(inM[i]*=inS;) return inM; } //!< Elementwise multiply each component
template<class T,int N> TptrMul operator*(const TptrMul& inM, const Tptr& inS) { TptrMul out(inM.clone()); return out *= inS; } //!< Elementwise multiply each component (preserve inputs)
template<class T,int N> TptrMul operator*(const Tptr& inS, const TptrMul& inM) { TptrMul out(inM.clone()); return out *= inS; } //!< Elementwise multiply each component (preserve inputs)
template<class T,int N> TptrMul operator*(TptrMul&& inM, const Tptr& inS) { return inM *= inS; } //!< Elementwise multiply each component (destructible input)
template<class T,int N> TptrMul operator*(const Tptr& inS, TptrMul&& inM) { return inM *= inS; } //!< Elementwise multiply each component (destructible input)

//Multiplication by scalars:
template<class T,int N> TptrMul& operator*=(TptrMul& in, double scaleFac) { Nloop(in[i]*=scaleFac;) return in; } //!< Scale
template<class T,int N> TptrMul operator*(const TptrMul& in, double scaleFac) { TptrMul out(in.clone()); return out *= scaleFac; } //!< Scalar multiply (preserve input)
template<class T,int N> TptrMul operator*(double scaleFac, const TptrMul& in) { TptrMul out(in.clone()); return out *= scaleFac; } //!< Scalar multiply (preserve input)
template<class T,int N> TptrMul operator*(TptrMul&& in, double scaleFac) { return in *= scaleFac; } //!< Scalar multiply (destructible input)
template<class T,int N> TptrMul operator*(double scaleFac, TptrMul&& in) { return in *= scaleFac; } //!< Scalar multiply (destructible input)
template<class T> ScalarFieldMultiplet<T,3> operator*(vector3<> v, const Tptr& in) { ScalarFieldMultiplet<T,3> out; for(int k=0; k<3; k++) out[k] = v[k] * in; return out; } //!< 3-vector multiply
template<class T> Tptr dot(vector3<> v, const ScalarFieldMultiplet<T,3>& in) { Tptr out; for(int k=0; k<3; k++) out += v[k] * in[k]; return out; } //!< 3-vector multiply

//Linear combine operators:
template<class T,int N> void axpy(double alpha, const TptrMul& X, TptrMul& Y) { Nloop(axpy(alpha, X[i], Y[i]);) } //!< Linear combine Y += alpha * X
template<class T,int N> TptrMul& operator+=(TptrMul& in, const TptrMul& other) { axpy(+1.0, other, in); return in; } //!<Increment
template<class T,int N> TptrMul& operator-=(TptrMul& in, const TptrMul& other) { axpy(-1.0, other, in); return in; } //!<Decrement
template<class T,int N> TptrMul operator+(const TptrMul& in1, const TptrMul& in2) { TptrMul out(in1.clone()); return out += in2; } //!<Add (preserve inputs)
template<class T,int N> TptrMul operator+(const TptrMul& in1, TptrMul&& in2) { return in2 += in1; } //!<Add (destructible input)
template<class T,int N> TptrMul operator+(TptrMul&& in1, const TptrMul& in2) { return in1 += in2; } //!<Add (destructible input)
template<class T,int N> TptrMul operator+(TptrMul&& in1, TptrMul&& in2) { return in1 += in2; } //!<Add (destructible inputs)
template<class T,int N> TptrMul operator-(const TptrMul& in1, const TptrMul& in2) { TptrMul out(in1.clone()); return out -= in2; } //!<Subtract (preserve input)
template<class T,int N> TptrMul operator-(const TptrMul& in1, TptrMul&& in2) { return (in2 -= in1) *= -1.0; } //!<Subtract (destructible input)
template<class T,int N> TptrMul operator-(TptrMul&& in1, const TptrMul& in2) { return in1 -= in2; } //!<Subtract (destructible input)
template<class T,int N> TptrMul operator-(TptrMul&& in1, TptrMul&& in2) { return in1 -= in2; } //!<Subtract (destructible inputs)
template<class T,int N> TptrMul operator-(const TptrMul& in) { return (-1.0)*in; } //!< Negate
template<class T,int N> TptrMul operator-(TptrMul&& in) { return in*=(-1.0); } //!< Negate

//Linear combine with singlets:
template<class T,int N> void axpy(double alpha, const Tptr& X, TptrMul& Y) { Nloop(axpy(alpha, X, Y[i]);) } //!< Linear combine Y += alpha * X
template<class T,int N> TptrMul& operator+=(TptrMul& in, const Tptr& other) { axpy(+1.0, other, in); return in; } //!<Increment
template<class T,int N> TptrMul& operator-=(TptrMul& in, const Tptr& other) { axpy(-1.0, other, in); return in; } //!<Decrement
template<class T,int N> TptrMul operator+(const TptrMul& in1, const Tptr& in2) { TptrMul out(in1.clone()); return out += in2; } //!<Add (preserve inputs)
template<class T,int N> TptrMul operator+(const Tptr& in1, const TptrMul& in2) { TptrMul out(in2.clone()); return out += in1; } //!<Add (preserve inputs)
template<class T,int N> TptrMul operator+(const Tptr& in1, TptrMul&& in2) { return in2 += in1; } //!<Add (destructible input)
template<class T,int N> TptrMul operator+(TptrMul&& in1, const Tptr& in2) { return in1 += in2; } //!<Add (destructible input)
template<class T,int N> TptrMul operator-(const TptrMul& in1, const Tptr& in2) { TptrMul out(in1.clone()); return out -= in2; } //!<Subtract (preserve input)
template<class T,int N> TptrMul operator-(const Tptr& in1, const TptrMul& in2) { TptrMul out(in2.clone()); return (out -= in1) *= -1.0; } //!<Subtract (preserve input)
template<class T,int N> TptrMul operator-(TptrMul&& in1, const Tptr& in2) { return in1 -= in2; } //!<Subtract (destructible input)
template<class T,int N> TptrMul operator-(const Tptr& in1, TptrMul&& in2) { return (in2 -= in1) *= -1.0; } //!<Subtract (destructible input)

//Norms and dot products
template<class T,int N> double dot(const TptrMul& X, const TptrMul& Y) { double ret=0.0; Nloop(ret+=dot(X[i],Y[i]);) return ret; } //!< Inner product
template<class T,int N> double nrm2(const TptrMul& X) { return sqrt(dot(X,X)); } //!< 2-norm
template<class T,int N> double sum(const TptrMul& X) { double ret=0.0; Nloop(ret+=sum(X[i]);) return ret; } //!< Sum of elements
inline vector3<> getGzero(const VectorFieldTilde& X) { vector3<> ret; for(int k=0; k<3; k++) if(X[k]) ret[k]=X[k]->getGzero(); return ret; } //!< return G=0 components
inline void setGzero(const VectorFieldTilde& X, vector3<> v) { for(int k=0; k<3; k++) if(X[k]) X[k]->setGzero(v[k]); } //!< set G=0 components
inline vector3<> sumComponents(const VectorField& X) { return vector3<>(sum(X[0]), sum(X[1]), sum(X[2])); } //!< Sum of elements (component-wise)
inline ScalarField lengthSquared(const VectorField& X) { return X[0]*X[0] + X[1]*X[1] + X[2]*X[2]; } //!< Elementwise length squared
inline ScalarField dotElemwise(const VectorField& X, const VectorField& Y) { return X[0]*Y[0] + X[1]*Y[1] + X[2]*Y[2]; } //!< Elementwise dot

//Extra operators in R-space alone for scalar additions:
template<int N> RptrMul& operator+=(RptrMul& in, double scalar) { Nloop(in[i]+=scalar;) return in; } //!<Increment by scalar
template<int N> RptrMul operator+(double scalar, const RptrMul& in) { RptrMul out(in.clone()); return out += scalar; } //!<Add scalar (preserve input)
template<int N> RptrMul operator+(const RptrMul& in, double scalar) { RptrMul out(in.clone()); return out += scalar; } //!<Add scalar (preserve input)
template<int N> RptrMul operator+(double scalar, RptrMul&& in) { return in += scalar; } //!<Add scalar (destructible input)
template<int N> RptrMul operator+(RptrMul&& in, double scalar) { return in += scalar; } //!<Add scalar (destructible input)

//Extra operators in G-space alone for real kernel multiplications:
template<int N> GptrMul& operator*=(GptrMul& in, const RealKernel& kernel) { Nloop(in[i]*=kernel;) return in; } //!< Multiply by kernel
template<int N> GptrMul operator*(const RealKernel& kernel, const GptrMul& in) { GptrMul out(in.clone()); return out *= kernel; } //!< Multiply by kernel (preserve input)
template<int N> GptrMul operator*(const GptrMul& in, const RealKernel& kernel) { GptrMul out(in.clone()); return out *= kernel; } //!< Multiply by kernel (preserve input)
template<int N> GptrMul operator*(const RealKernel& kernel, GptrMul&& in) { return in *= kernel; } //!< Multiply by kernel (destructible input)
template<int N> GptrMul operator*(GptrMul&& in, const RealKernel& kernel) { return in *= kernel; } //!< Multiply by kernel (destructible input)

//Transform operators:
template<int N> GptrMul O(GptrMul&& X) { Nloop( O((ScalarFieldTilde&&)X[i]); ) return X; } //!< Inner product operator (diagonal in PW basis)
template<int N> GptrMul O(const GptrMul& X) { return O(X.clone()); } //!< Inner product operator (diagonal in PW basis)
template<int N> RptrMul I(GptrMul&& X, bool compat=false); //!< Forward transform: PW basis -> real space (destructible input)
template<int N> GptrMul J(const RptrMul& X); //!< Inverse transform: Real space -> PW basis
template<int N> GptrMul Idag(const RptrMul& X); //!< Forward transform transpose: Real space -> PW basis
template<int N> RptrMul Jdag(GptrMul&& X, bool compat=false); //!< Inverse transform transpose: PW basis -> real space (destructible input)
template<int N> RptrMul Jdag(const GptrMul& X, bool compat=false) { return Jdag(X.clone(), compat); } //!< Inverse transform transpose: PW basis -> real space (preserve input)
template<int N> RptrMul I(const GptrMul& X, bool compat=false) { return I(X.clone(), compat); } //!< Forward transform: PW basis -> real space (preserve input)

//Special operators for triplets (implemented in operators.cpp):
VectorFieldTilde gradient(const ScalarFieldTilde&); //!< compute the gradient of a complex field, returns cartesian components
VectorField gradient(const ScalarField&); //!< compute the gradient of a complex field, returns cartesian components
ScalarFieldTilde divergence(const VectorFieldTilde&); //!< compute the divergence of a vector field specified in cartesian components
ScalarField divergence(const VectorField&); //!< compute the divergence of a vector field specified in cartesian components

//Special operators for symmetric traceless tensors (implemented in operators.cpp)
TensorFieldTilde tensorGradient(const ScalarFieldTilde&); //!< symmetric traceless tensor second derivative of a scalar field
ScalarFieldTilde tensorDivergence(const TensorFieldTilde&); //!<  second derivative contraction of a symmetric traceless tensor field

//Debug:
template<int N> void printStats(const RptrMul&, const char* name, FILE* fpLog=stdout); //!< Print mean and standard deviation of each component array with specified name (debug utility)

//! @}

//###################################################################################################
//####  Implementation  ####
//##########################
//!@cond

template<class T, int N>
std::vector<typename T::DataType*> TptrMul::data()
{	std::vector<typename T::DataType*> ret(N, (typename T::DataType*)0);
	Nloop( if(component[i]) ret[i] = component[i]->data(); )
	return ret;
}
template<class T, int N>
std::vector<const typename T::DataType*> TptrMul::data() const
{	std::vector<const typename T::DataType*> ret(N, (const typename T::DataType*)0);
	Nloop( if(component[i]) ret[i] = component[i]->data(); )
	return ret;
}
#ifdef GPU_ENABLED
template<class T, int N>
std::vector<typename T::DataType*> TptrMul::dataGpu()
{	std::vector<typename T::DataType*> ret(N, (typename T::DataType*)0);
	Nloop( if(component[i]) ret[i] = component[i]->dataGpu(); )
	return ret;
}
template<class T, int N>
std::vector<const typename T::DataType*> TptrMul::dataGpu() const
{	std::vector<const typename T::DataType*> ret(N, (typename T::DataType*)0);
	Nloop( if(component[i]) ret[i] = component[i]->dataGpu(); )
	return ret;
}
#endif

template<class T, int N>
void TptrMul::loadFromFile(const char* filename)
{	//Checks for the correct filesize
	off_t expectedLen = 0;
	Nloop(expectedLen += sizeof(typename T::DataType) * component[i]->nElem;)
	off_t fLen = fileSize(filename);
	if(fLen != expectedLen)
	{	die("\nLength of '%s' was %ld instead of the expected %ld bytes.\n"
				"Hint: Are you really reading the correct file?\n\n",
				filename, (unsigned long)fLen, (unsigned long)expectedLen);
	}
	
	FILE* fp = fopen(filename, "rb");
	if(!fp) die("Could not open %s for reading.\n", filename)
	Nloop(
		if(!component[i]) die("Component %d was null in loadFromFile(\"%s\").\n", i, filename)	
		if(fread(component[i]->data(), sizeof(typename T::DataType), component[i]->nElem, fp) < unsigned(component[i]->nElem))
			die("File ended too soon while reading component %d in loadFromFile(\"%s\").\n", i, filename)
	)
	fclose(fp);
}
template<class T, int N>
void TptrMul::saveToFile(const char* filename) const
{	FILE* fp = fopen(filename, "wb");
	if(!fp) die("Could not open %s for writing.\n", filename)
	Nloop(
		if(!component[i]) die("Component %d was null in saveToFile(\"%s\").\n", i, filename)
		fwrite(component[i]->data(), sizeof(typename T::DataType), component[i]->nElem, fp);
	)
	fclose(fp);
}

namespace ScalarFieldMultipletPrivate
{
	inline ScalarField IcompatTrue(ScalarFieldTilde&& in, int nThreads) { return I((ScalarFieldTilde&&)in, true, nThreads); }
	inline ScalarField IcompatFalse(ScalarFieldTilde&& in, int nThreads) { return I((ScalarFieldTilde&&)in, false, nThreads); }
	inline ScalarField JdagCompatTrue(ScalarFieldTilde&& in, int nThreads) { return Jdag((ScalarFieldTilde&&)in, true, nThreads); }
	inline ScalarField JdagCompatFalse(ScalarFieldTilde&& in, int nThreads) { return Jdag((ScalarFieldTilde&&)in, false, nThreads); }
	
	template<typename FuncOut, typename FuncIn, typename Out, typename In>
	void threadUnary_sub(int iOpThread, int nOpThreads, int nThreadsTot, int N, FuncOut (*func)(FuncIn,int), Out* out, In in)
	{	//Divide jobs amongst operator threads:
		int iStart = (iOpThread*N)/nOpThreads;
		int iStop = ((iOpThread+1)*N)/nOpThreads;
		//Divide total threads amongst operator threads:
		int nThreads = ((iOpThread+1)*nThreadsTot)/nOpThreads - (iOpThread*nThreadsTot)/nOpThreads;
		for(int i=iStart; i<iStop; i++)
			(*out)[i] = func((FuncIn)in[i], nThreads);
	}

	template<typename FuncOut, typename FuncIn, typename Out, typename In>
	void threadUnary(FuncOut (*func)(FuncIn,int), int N, Out* out, In in)
	{	int nThreadsTot = (isGpuEnabled() || !shouldThreadOperators()) ? 1 : nProcsAvailable;
		int nOperatorThreads = std::min(nThreadsTot, N);
		threadLaunch(nOperatorThreads, threadUnary_sub<FuncOut,FuncIn,Out,In>, 0, nThreadsTot, N, func, out, in);
	}
};

template<int N>
RptrMul I(GptrMul&& X, bool compat)
{	using namespace ScalarFieldMultipletPrivate;
	RptrMul out;
	ScalarField (*func)(ScalarFieldTilde&&,int) = compat ? IcompatTrue : IcompatFalse;
	threadUnary<ScalarField,ScalarFieldTilde&&>(func, N, &out, X);
	return out;
}

template<int N>
GptrMul J(const RptrMul& X)
{	using namespace ScalarFieldMultipletPrivate;
	GptrMul out;
	ScalarFieldTilde (*func)(const ScalarField&,int) = J;
	threadUnary(func, N, &out, X);
	return out;
}

template<int N>
GptrMul Idag(const RptrMul& X)
{	using namespace ScalarFieldMultipletPrivate;
	GptrMul out;
	ScalarFieldTilde (*func)(const ScalarField&,int) = Idag;
	threadUnary(func, N, &out, X);
	return out;
}

template<int N>
RptrMul Jdag(GptrMul&& X, bool compat)
{	using namespace ScalarFieldMultipletPrivate;
	RptrMul out;
	ScalarField (*func)(ScalarFieldTilde&&,int) = compat ? JdagCompatTrue : JdagCompatFalse;
	threadUnary<ScalarField,ScalarFieldTilde&&>(func, N, &out, X);
	return out;
}

template<int N> void printStats(const RptrMul& X, const char* namePrefix, FILE* fpLog)
{	char name[256];
	Nloop( sprintf(name, "%s[%d]", namePrefix, i); printStats(X[i], name, fpLog); )
}

#undef Tptr
#undef TptrMul
#undef RptrMul
#undef GptrMul
//!@endcond

#endif // JDFTX_CORE_VECTORFIELD_H
