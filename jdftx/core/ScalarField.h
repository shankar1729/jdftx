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


//! @addtogroup griddata
//! @{

/** @file ScalarField.h
@brief  Real and complex scalar fields in real and reciprocal space
*/

#include <memory>
#include <core/scalar.h>
#include <core/Util.h>
#include <core/ManagedMemory.h>

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

//Define shorthands for fetching gpu data in gpu mode and cpu data in cpu mode:
#ifdef GPU_ENABLED
	#define DECLARE_DATA_PREF_ACCESS \
		DataType* dataPref(bool shouldAbsorbScale=true) { return dataGpu(shouldAbsorbScale); } \
		const DataType* dataPref(bool shouldAbsorbScale=true) const { return dataGpu(shouldAbsorbScale); }
#else
	#define DECLARE_DATA_PREF_ACCESS \
		DataType* dataPref(bool shouldAbsorbScale=true) { return data(shouldAbsorbScale); } \
		const DataType* dataPref(bool shouldAbsorbScale=true) const { return data(shouldAbsorbScale); }
#endif

//! Base class for #ScalarFieldData and #ScalarFieldTildeData
struct FieldData : private ManagedMemory
{
	int nElem; //!< number of elements = #gInfo.nr
	double scale; //!< overall scale factor of the data array
	const GridInfo& gInfo; //!< simulation grid info

	void absorbScale() const; //!< absorb scale factor into data
	inline void zero() { ManagedMemory::zero(); } //!< initialize to zero

	typedef void DataType; //!< this base class has no specific data type
	
	void* data(bool shouldAbsorbScale=true); //!< get a pointer to the actual data (after absorbing the scale factor, unless otherwise specified)
	const void* data(bool shouldAbsorbScale=true) const; //!< get a const pointer to the actual data (after absorbing the scale factor, unless otherwise specified)
	#ifdef GPU_ENABLED
	void* dataGpu(bool shouldAbsorbScale=true); //!< get a pointer to the actual data (after absorbing the scale factor, unless otherwise specified)
	const void* dataGpu(bool shouldAbsorbScale=true) const; //!< get a const pointer to the actual data (after absorbing the scale factor, unless otherwise specified)
	#endif

	DECLARE_DATA_PREF_ACCESS
	
	inline bool isOnGpu() const { return ManagedMemory::isOnGpu(); } //!< Check where the data is (for #ifdef simplicity exposed even when no GPU_ENABLED)

	FieldData(const GridInfo& gInfo, string category, int nElem, int nDoublesPerElem=1, bool onGpu=false);
	void copyData(const FieldData& other); //!< copy data and scale (used by clone())

	//Inter-process communication (expose corresponding ManagedMemory functions):
	inline void send(int dest, int tag=0) const { absorbScale(); ManagedMemory::send(dest,tag); } //!< send to another process
	inline void recv(int src, int tag=0) { absorbScale(); ManagedMemory::recv(src,tag); } //!< receive from another process
	inline void bcast(int root=0) { absorbScale(); ManagedMemory::bcast(root); } //!< synchronize across processes (using value on specified root process)
	inline void allReduce(MPIUtil::ReduceOp op, bool safeMode=false) { absorbScale(); ManagedMemory::allReduce(op, safeMode, isReal); } //!< apply all-to-all reduction

protected:
	struct PrivateTag {}; //!< Used to prevent direct use of ScalarField constructors, and force the shared_ptr usage
private:
	bool isReal; //!< whether underlying data type is real (nDoublesPerElem==1 at construction)
};

//Shorthand for defining the data() and dataGpu() functions in derived classes of FieldData
#ifdef GPU_ENABLED
	#define DECLARE_DATA_ACCESS \
		DataType* data(bool shouldAbsorbScale=true) { return (DataType*)FieldData::data(shouldAbsorbScale); } \
		const DataType* data(bool shouldAbsorbScale=true) const { return (const DataType*)FieldData::data(shouldAbsorbScale); } \
		DataType* dataGpu(bool shouldAbsorbScale=true) { return (DataType*)FieldData::dataGpu(shouldAbsorbScale); } \
		const DataType* dataGpu(bool shouldAbsorbScale=true) const { return (const DataType*)FieldData::dataGpu(shouldAbsorbScale); }
#else
	#define DECLARE_DATA_ACCESS \
		DataType* data(bool shouldAbsorbScale=true) { return (DataType*)FieldData::data(shouldAbsorbScale); } \
		const DataType* data(bool shouldAbsorbScale=true) const { return (const DataType*)FieldData::data(shouldAbsorbScale); }
#endif

/**
@brief Real space real scalar field data
Do not use this data structure directly or from a simple pointer
ScalarFieldData*; work only with #ScalarField's. The public functions of ScalarFieldData
can be accessed with -> from the ScalarField.
*/
struct ScalarFieldData : public FieldData
{	typedef double DataType; //!< Type of data in container (useful for templating)
	DECLARE_DATA_ACCESS
	DECLARE_DATA_PREF_ACCESS
	ScalarField clone() const; //!< clone the data (NOTE: assigning ScalarField's makes a new reference to the same data)
	static ScalarField alloc(const GridInfo& gInfo, bool onGpu=false); //!< Create real space data
	ScalarFieldData(const GridInfo& gInfo, bool onGpu, PrivateTag); //!< called only by ScalarFieldData::alloc()
};


/**
@brief Reciprocal space real scalar field data
Do not use this data structure directly or from a simple pointer
ScalarFieldTildeData*; work only with #ScalarFieldTilde's. The public functions of ScalarFieldTildeData
can be accessed with -> from the ScalarFieldTilde.
*/
struct ScalarFieldTildeData : public FieldData
{	typedef complex DataType; //!< Type of data in container (useful for templating)
	DECLARE_DATA_ACCESS
	DECLARE_DATA_PREF_ACCESS
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
struct complexScalarFieldData : public FieldData
{	typedef complex DataType; //!< Type of data in container (useful for templating)
	DECLARE_DATA_ACCESS
	DECLARE_DATA_PREF_ACCESS
	complexScalarField clone() const; //!< clone the data (NOTE: assigning complexScalarField's makes a new reference to the same data)
	static complexScalarField alloc(const GridInfo& gInfo, bool onGpu=false); //!< Create real space data
	complexScalarFieldData(const GridInfo& gInfo, bool onGpu, PrivateTag); //!< called only by complexScalarFieldData::alloc()
};

/**
@brief Reciprocal space complex scalar field data
Do not use this data structure directly or from a simple pointer
complexScalarFieldTildeData*; work only with #complexScalarFieldTilde's. The public functions of complexScalarFieldTildeData
can be accessed with -> from the complexScalarFieldTilde.
*/
struct complexScalarFieldTildeData : public FieldData
{	typedef complex DataType; //!< Type of data in container (useful for templating)
	DECLARE_DATA_ACCESS
	DECLARE_DATA_PREF_ACCESS
	complexScalarFieldTilde clone() const; //!< clone the data (NOTE: assigning complexScalarFieldTilde's makes a new reference to the same data)
	static complexScalarFieldTilde alloc(const GridInfo& gInfo, bool onGpu=false); //!< Create reciprocal space data
	complexScalarFieldTildeData(const GridInfo& gInfo, bool onGpu, PrivateTag); //!< called only by complexScalarFieldTildeData::alloc()
};



//! Special class for storing real reciprocal-space kernels encountered ever so often for convolutions
struct RealKernel
{	const GridInfo& gInfo;
	int nElem; //!< number of elements = gInfo.nG
	double* data; //!< cpu data pointer
	#ifdef GPU_ENABLED
	double* dataGpu; //!< gpu data pointer (unlike above classes, both cpu and gpu stored simultaneously)
	#endif
	double* dataPref; //!< points to data or dataGpu depending on GPU_ENABLED

	RealKernel(const GridInfo&); //!< allocates data (and dataGpu if GPU_ENABLED)
	~RealKernel();
	void set(); //!< call after initializing data on cpu (will copy from data to dataGpu if GPU_ENABLED)
};

#undef DECLARE_DATA_PREF_ACCESS
#undef DECLARE_DATA_ACCESS
//! @}

#endif //JDFTX_CORE_SCALARFIELD_H

