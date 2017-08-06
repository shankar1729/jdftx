/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Deniz Gunceler

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

#ifndef JDFTX_CORE_SCALARFIELDARRAY_H
#define JDFTX_CORE_SCALARFIELDARRAY_H

//! @addtogroup DataStructures
//! @{

//! @file ScalarFieldArray.h Variable length arrays of ScalarField and ScalarFieldTildeArray, and their operators

#include <core/VectorField.h>
#include <cstdio>

typedef std::vector<ScalarField> ScalarFieldArray; //!< dynamic size collection of real space scalar fields
typedef std::vector<ScalarFieldTilde> ScalarFieldTildeArray; //!< dynamic size collection of reciprocal space scalar fields
#define TptrCollection std::vector<std::shared_ptr<T> > //!< shorthand for templates below (undef'd at end of file)

//! Extract a std::vector of data pointers from a ScalarFieldArray
template<typename T> std::vector<typename T::DataType*> dataPref(TptrCollection& x)
{	std::vector<typename T::DataType*> xData(x.size());
	for(unsigned s=0; s<x.size(); s++)
		xData[s] = x[s] ? x[s]->dataPref() : 0;
	return xData;
}

//! Extract a std::vector of const data pointers from a const ScalarFieldArray
template<typename T> std::vector<const typename T::DataType*> constDataPref(const TptrCollection& x)
{	std::vector<const typename T::DataType*> xData(x.size());
	for(unsigned s=0; s<x.size(); s++)
		xData[s] = x[s] ? x[s]->dataPref() : 0;
	return xData;
}

//! Create a copy of the data (note operator= references same data since Tptr's are pointers!)
template<typename T> TptrCollection clone(const TptrCollection& x)
{	TptrCollection ret(x.size());
	for(unsigned i=0; i<x.size(); i++) if(x[i]) ret[i] = clone(x[i]);
	return ret;
}

//! Scale
template<typename T> TptrCollection& operator*=(TptrCollection& x, double alpha)
{ 	for(unsigned i=0; i<x.size(); i++) if(x[i]) x[i] *= alpha;
	return x;
}
template<class T> TptrCollection operator*(const TptrCollection& in, double scaleFac) { 
	 TptrCollection out(clone(in)); return out *= scaleFac;
}
template<class T> TptrCollection operator*(double scaleFac, const TptrCollection& in) { 
	 TptrCollection out(clone(in)); return out *= scaleFac;
}
template<class T> TptrCollection operator*(TptrCollection&& in, double scaleFac) { return in *= scaleFac; } //!< Add (destructible input)
template<class T> TptrCollection operator*(double scaleFac, TptrCollection&& in) { return in *= scaleFac; } //!< Add (destructible input)

//! y += alpha x
template<typename T> void axpy(double alpha, const TptrCollection& x, TptrCollection& y)
{	assert(x.size()==y.size());
	for(unsigned i=0; i<x.size(); i++) axpy(alpha, x[i], y[i]);
}

//! elementise multiply
inline ScalarFieldArray operator*(const ScalarFieldArray& x, const ScalarFieldArray& y)
{	assert(x.size()==y.size());
	ScalarFieldArray z(x.size());
	for(unsigned i=0; i<x.size(); i++) z[i] = x[i]*y[i];
	return z;
}
inline ScalarFieldArray operator*(ScalarFieldArray&& x, const ScalarFieldArray& y)
{	assert(x.size()==y.size());
	for(unsigned i=0; i<x.size(); i++) x[i] *= y[i];
	return x;
}
inline ScalarFieldArray operator*(const ScalarFieldArray& x, ScalarFieldArray&& y)
{	assert(x.size()==y.size());
	for(unsigned i=0; i<x.size(); i++) y[i] *= x[i];
	return y;
}

inline ScalarFieldArray& operator*=(ScalarFieldArray& y, const ScalarField& x)
{	for(unsigned i=0; i<y.size(); i++) if(y[i]) y[i] *= x;
	return y;
}
inline ScalarFieldArray operator*(const ScalarField& x, ScalarFieldArray&& y) { return y *= x; }
inline ScalarFieldArray operator*(ScalarFieldArray&& y, const ScalarField& x) { return y *= x; }
inline ScalarFieldArray operator*(const ScalarField& x, const ScalarFieldArray& y) { ScalarFieldArray z = clone(y); return z *= x; }
inline ScalarFieldArray operator*(const ScalarFieldArray& y, const ScalarField& x) { ScalarFieldArray z = clone(y); return z *= x; }

//! Increment
template<class T> TptrCollection& operator+=(TptrCollection& in, const TptrCollection& other) { axpy(+1.0, other, in); return in; }
//! Decrement
template<class T> TptrCollection& operator-=(TptrCollection& in, const TptrCollection& other) { axpy(-1.0, other, in); return in; }

// Addition
template<class T> TptrCollection operator+(const TptrCollection& in1, const TptrCollection& in2) { //!< Add (destructible inputs)
	TptrCollection out(clone(in1)); 
	return out += in2; 
} 
template<class T> TptrCollection operator+(const TptrCollection& in1, TptrCollection&& in2) { return in2 += in1; } //!< Add (destructible input)
template<class T> TptrCollection operator+(TptrCollection&& in1, const TptrCollection& in2) { return in1 += in2; } //!< Add (destructible input)
template<class T> TptrCollection operator+(TptrCollection&& in1, TptrCollection&& in2) { return in1 += in2; } //!< Add (destructible inputs)

// Subtraction
template<class T> TptrCollection operator-(const TptrCollection& in1, const TptrCollection& in2) { //!< Add (destructible inputs)
	TptrCollection out(clone(in1)); 
	return out -= in2; 
} 
template<class T> TptrCollection operator-(const TptrCollection& in1, TptrCollection&& in2) { return in2 -= in1; } //!< Add (destructible inputs)
template<class T> TptrCollection operator-(TptrCollection&& in1, const TptrCollection& in2) { return in1 -= in2; } //!< Add (destructible inputs)
template<class T> TptrCollection operator-(TptrCollection&& in1, TptrCollection&& in2) { return in1 -= in2; } //!< Add (destructible inputs)

//! Inner product
template<typename T> double dot(const TptrCollection& x, const TptrCollection& y)
{	assert(x.size()==y.size());
	double ret = 0.0;
	for(unsigned i=0; i<x.size(); i++) if(x[i] && y[i]) ret += dot(x[i], y[i]);
	return ret;
}

//Spherical tensor derivatives
ScalarFieldTildeArray lGradient(const ScalarFieldTilde&, int l); //!< spherical tensor gradient of order l (2l+1 outputs, multiplied by Ylm(Ghat) (iG)^l)
ScalarFieldTilde lDivergence(const ScalarFieldTildeArray&, int l); //!< spherical tensor divergence of order l (2l+1 inputs, multiplied by Ylm(Ghat) (iG)^l, and summed)

//! Initialize (non-null) data to zero
template<typename T> void initZero(TptrCollection& x)
{	for(unsigned i=0; i<x.size(); i++) initZero(x[i]);
}

//! Allocate all null components and initialize them to zero
//! If N is non-zero, resize collection to N scalar fields
template<typename T> void nullToZero(TptrCollection& x, const GridInfo& gInfo, int N=0)
{	if(N) x.resize(N);
	for(unsigned i=0; i<x.size(); i++) nullToZero(x[i], gInfo);
}

//! Initialize to random numbers (uniform on 0 to 1)
template<typename T> void initRandomFlat(TptrCollection& x)
{	for(unsigned i=0; i<x.size(); i++) initRandomFlat(x[i]);
}
//! Initialize to normal random numbers:
template<typename T> void randomize(TptrCollection& x)
{	for(unsigned i=0; i<x.size(); i++) initRandom(x[i], 3.0);
}

//! Load an array of scalar fields from a file
template<typename T> void loadFromFile(TptrCollection& x, const char* filename)
{	//Checks for the correct filesize
	off_t expectedLen = 0;
	for(unsigned i=0; i<x.size(); i++){expectedLen += sizeof(typename T::DataType) * x[i]->nElem;}
	off_t fLen = fileSize(filename);
	if(fLen != expectedLen)
	{	die("\nLength of '%s' was %ld instead of the expected %ld bytes.\n"
				"Hint: Are you really reading the correct file?\n\n",
				filename, (unsigned long)fLen, (unsigned long)expectedLen);
	}
	
	FILE* fp = fopen(filename, "rb");
	if(!fp) die("Could not open %s for reading.\n", filename)
	for(unsigned i=0; i<x.size(); i++)
	{	if(!x[i]) die("x[%d] was null in loadFromFile(x,\"%s\").\n", i, filename)
		if(freadLE(x[i]->data(), sizeof(typename T::DataType), x[i]->nElem, fp) < unsigned(x[i]->nElem))
			die("File ended too soon while reading x[%d] in loadFromFile(x,\"%s\").\n", i, filename)
	}
	fclose(fp);
}

//! Save an array of scalar fields to a file
template<typename T> void saveToFile(const TptrCollection& x, const char* filename)
{	FILE* fp = fopen(filename, "wb");
	if(!fp) die("Could not open %s for writing.\n", filename)
	for(unsigned i=0; i<x.size(); i++)
	{	if(!x[i]) die("x[%d] was null in saveToFile(x,\"%s\").\n", i, filename)
		fwriteLE(x[i]->data(), sizeof(typename T::DataType), x[i]->nElem, fp);
	}
	fclose(fp);
}

//----------------- Transform operators -------------------------

inline ScalarFieldArray I(ScalarFieldTildeArray&& X) //!< Reciprocal to real space
{	using namespace ScalarFieldMultipletPrivate;
	ScalarFieldArray out(X.size());
	ScalarField (*func)(ScalarFieldTilde&&,int) = I;
	threadUnary<ScalarField,ScalarFieldTilde&&>(func, int(X.size()), &out, X);
	return out;
}
inline ScalarFieldArray I(const ScalarFieldTildeArray& X) { return I(clone(X)); } //!< Reciprocal to real space

inline ScalarFieldTildeArray J(const ScalarFieldArray& X) //!< Real to reciprocal space
{	using namespace ScalarFieldMultipletPrivate;
	ScalarFieldTildeArray out(X.size());
	ScalarFieldTilde (*func)(const ScalarField&,int) = J;
	threadUnary(func, int(X.size()), &out, X);
	return out;
}

inline ScalarFieldTildeArray Idag(const ScalarFieldArray& X) //!< Hermitian conjugate of I
{	using namespace ScalarFieldMultipletPrivate;
	ScalarFieldTildeArray out(X.size());
	ScalarFieldTilde (*func)(const ScalarField&,int) = Idag;
	threadUnary(func, int(X.size()), &out, X);
	return out;
}

inline ScalarFieldArray Jdag(ScalarFieldTildeArray&& X) //!< Hermitian conjugate of J
{	using namespace ScalarFieldMultipletPrivate;
	ScalarFieldArray out(X.size());
	ScalarField (*func)(ScalarFieldTilde&&,int) = Jdag;
	threadUnary<ScalarField,ScalarFieldTilde&&>(func, int(X.size()), &out, X);
	return out;
}
inline ScalarFieldArray Jdag(const ScalarFieldTildeArray& X) { return Jdag(clone(X)); } //!< Hermitian conjugate of J

#undef TptrCollection

//! @}
#endif // JDFTX_CORE_SCALARFIELDARRAY_H
