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

#ifndef JDFTX_CORE_SCALARFIELD_H
#define JDFTX_CORE_SCALARFIELD_H


//! @addtogroup DataStructures
//! @{

/** @file ScalarField.h
@brief  Real and complex scalar fields in real and reciprocal space
*/

#include <core/scalar.h>
#include <core/Util.h>
#include <core/ManagedMemory.h>
#include <core/matrix.h>
#include <core/GridInfo.h>
#include <memory>

class GridInfo; //Grid description and memory manager
struct ScalarFieldData; //Real space data storage container for real scalar fields
struct ScalarFieldTildeData; //Reciprocal space data storage container for real scalar fields
struct complexScalarFieldData; //Real space data storage container for complex scalar fields
struct complexScalarFieldTildeData; //Reciprocal space data storage container for complex scalar fields

//DO NOT work with above data types directly, instead use the following smart pointer types:
typedef std::shared_ptr<ScalarFieldData> ScalarField; //!< A smart reference-counting pointer to #ScalarFieldData
typedef std::shared_ptr<ScalarFieldTildeData> ScalarFieldTilde; //!< A smart reference-counting pointer to #ScalarFieldTildeData
typedef std::shared_ptr<complexScalarFieldData> complexScalarField; //!< A smart reference-counting pointer to #complexScalarFieldData
typedef std::shared_ptr<complexScalarFieldTildeData> complexScalarFieldTilde; //!< A smart reference-counting pointer to #complexScalarFieldTildeData

//! ManagedMemory wrapper with gridInfo and pending scale factor for ScalarField* classes
template<typename T> struct FieldData : private ManagedMemory<T>
{
	const int nElem; //!< number of elements
	double scale; //!< overall scale factor of the data array
	const GridInfo& gInfo; //!< simulation grid info

	FieldData(const GridInfo& gInfo, string category, int nElem, bool onGpu=false) : nElem(nElem), scale(1.), gInfo(gInfo)
	{	ManagedMemory<T>::memInit(category, nElem, onGpu);
	}
	
	//! Copy data and scale (used by clone())
	void copyData(const FieldData<T>& other)
	{	scale = other.scale;
		memcpy((ManagedMemory<T>&)(*this), other);
	}
	
	//! Absorb scale factor into data
	void absorbScale() const
	{	if(scale != 1.)
		{	FieldData& X = (FieldData&)(*this); //cast to non-const (this function modifies data, but is logically constant)
			::scale(scale, (ManagedMemory<T>&)X);
			X.scale = 1.;
		}
	}
	
	#define getDataCode(dataLoc) \
		if(shouldAbsorbScale) absorbScale(); \
		return ManagedMemory<T>::dataLoc();
	T* data(bool shouldAbsorbScale=true) { getDataCode(data) } //!< get a pointer to the actual data (after absorbing the scale factor, unless otherwise specified)
	const T* data(bool shouldAbsorbScale=true) const { getDataCode(data) } //!< get a const pointer to the actual data (after absorbing the scale factor, unless otherwise specified)
	#ifdef GPU_ENABLED
	T* dataGpu(bool shouldAbsorbScale=true) { getDataCode(dataGpu) } //!< get a pointer to the actual data (after absorbing the scale factor, unless otherwise specified)
	const T* dataGpu(bool shouldAbsorbScale=true) const { getDataCode(dataGpu) } //!< get a const pointer to the actual data (after absorbing the scale factor, unless otherwise specified)
	#endif
	#undef getDataCode
	
	//Shorthands for fetching gpu data in gpu mode and cpu data in cpu mode:
	#ifdef GPU_ENABLED
	T* dataPref(bool shouldAbsorbScale=true) { return dataGpu(shouldAbsorbScale); }
	const T* dataPref(bool shouldAbsorbScale=true) const { return dataGpu(shouldAbsorbScale); }
	#else
	T* dataPref(bool shouldAbsorbScale=true) { return data(shouldAbsorbScale); }
	const T* dataPref(bool shouldAbsorbScale=true) const { return data(shouldAbsorbScale); }
	#endif
	

	inline void zero() { ManagedMemory<T>::zero(); } //!< initialize to zero
	inline bool isOnGpu() const { return ManagedMemory<T>::isOnGpu(); } //!< Check where the data is (for #ifdef simplicity exposed even when no GPU_ENABLED)

	//Inter-process communication (expose corresponding ManagedMemory functions):
	inline void send(int dest, int tag=0) const { absorbScale(); ManagedMemory<T>::send(dest,tag); } //!< send to another process
	inline void recv(int src, int tag=0) { absorbScale(); ManagedMemory<T>::recv(src,tag); } //!< receive from another process
	inline void bcast(int root=0) { absorbScale(); ManagedMemory<T>::bcast(root); } //!< synchronize across processes (using value on specified root process)
	inline void allReduce(MPIUtil::ReduceOp op, bool safeMode=false) { absorbScale(); ManagedMemory<T>::allReduce(op, safeMode); } //!< apply all-to-all reduction

protected:
	struct PrivateTag {}; //!< Used to prevent direct use of ScalarField constructors, and force the shared_ptr usage
};

/**
@brief Real space real scalar field data
Do not use this data structure directly or from a simple pointer
ScalarFieldData*; work only with #ScalarField's. The public functions of ScalarFieldData
can be accessed with -> from the ScalarField.
*/
struct ScalarFieldData : public FieldData<double>
{	typedef double DataType; //!< Type of data in container (useful for templating)
	ScalarField clone() const; //!< clone the data (NOTE: assigning ScalarField's makes a new reference to the same data)
	static ScalarField alloc(const GridInfo& gInfo, bool onGpu=false); //!< Create real space data
	ScalarFieldData(const GridInfo& gInfo, bool onGpu, PrivateTag); //!< called only by ScalarFieldData::alloc()
	matrix toMatrix() const; //!<convert to (complex) matrix
};


/**
@brief Reciprocal space real scalar field data
Do not use this data structure directly or from a simple pointer
ScalarFieldTildeData*; work only with #ScalarFieldTilde's. The public functions of ScalarFieldTildeData
can be accessed with -> from the ScalarFieldTilde.
*/
struct ScalarFieldTildeData : public FieldData<complex>
{	typedef complex DataType; //!< Type of data in container (useful for templating)
	ScalarFieldTilde clone() const; //!< clone the data (NOTE: assigning ScalarFieldTilde's makes a new reference to the same data)
	static ScalarFieldTilde alloc(const GridInfo& gInfo, bool onGpu=false); //!< Create reciprocal space data
	double getGzero() const; //!< get the G=0 component
	void setGzero(double Gzero); //!< set the G=0 component
	ScalarFieldTildeData(const GridInfo& gInfo, bool onGpu, PrivateTag); //!< called only by ScalarFieldTildeData::alloc()
};


/**
@brief Real space complex scalar field data
Do not use this data structure directly or from a simple pointer
complexScalarFieldData*; work only with #complexScalarField's. The public functions of complexScalarFieldData
can be accessed with -> from the complexScalarField.
*/
struct complexScalarFieldData : public FieldData<complex>
{	typedef complex DataType; //!< Type of data in container (useful for templating)
	complexScalarField clone() const; //!< clone the data (NOTE: assigning complexScalarField's makes a new reference to the same data)
	static complexScalarField alloc(const GridInfo& gInfo, bool onGpu=false); //!< Create real space data
	complexScalarFieldData(const GridInfo& gInfo, bool onGpu, PrivateTag); //!< called only by complexScalarFieldData::alloc()
	matrix toMatrix() const; //!<convert to matrix
};

/**
@brief Reciprocal space complex scalar field data
Do not use this data structure directly or from a simple pointer
complexScalarFieldTildeData*; work only with #complexScalarFieldTilde's. The public functions of complexScalarFieldTildeData
can be accessed with -> from the complexScalarFieldTilde.
*/
struct complexScalarFieldTildeData : public FieldData<complex>
{	typedef complex DataType; //!< Type of data in container (useful for templating)
	complexScalarFieldTilde clone() const; //!< clone the data (NOTE: assigning complexScalarFieldTilde's makes a new reference to the same data)
	static complexScalarFieldTilde alloc(const GridInfo& gInfo, bool onGpu=false); //!< Create reciprocal space data
	complexScalarFieldTildeData(const GridInfo& gInfo, bool onGpu, PrivateTag); //!< called only by complexScalarFieldTildeData::alloc()
};



//! Special class for storing real reciprocal-space kernels encountered ever so often for convolutions
struct RealKernel : public FieldData<double>
{	RealKernel(const GridInfo& gInfo); 
};

//! @}
#endif //JDFTX_CORE_SCALARFIELD_H
