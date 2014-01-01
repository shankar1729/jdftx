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

#ifndef JDFTX_CORE_DATA_H
#define JDFTX_CORE_DATA_H


//! @addtogroup griddata
//! @{

/** @file Data.h
@brief  Data storage containers for real and reciprocal space arrays over the simulation grid.
*/

#include <memory>
#include <core/scalar.h>
#include <core/Util.h>

class GridInfo; //Grid description and memory manager
struct DataR; //Real space data storage container for real scalar fields
struct DataG; //Reciprocal space data storage container for real scalar fields
struct complexDataR; //Real space data storage container for complex scalar fields
struct complexDataG; //Reciprocal space data storage container for complex scalar fields

//DO NOT work with above data types directly, instead use the following smart pointer types:
typedef std::shared_ptr<DataR> DataRptr; //!< A smart reference-counting pointer to #DataR
typedef std::shared_ptr<DataG> DataGptr; //!< A smart reference-counting pointer to #DataG
typedef std::shared_ptr<complexDataR> complexDataRptr; //!< A smart reference-counting pointer to #complexDataR
typedef std::shared_ptr<complexDataG> complexDataGptr; //!< A smart reference-counting pointer to #complexDataG

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

//! Base class for #DataR and #DataG
struct Data
{
	int nElem; //!< number of elements = #gInfo.nr
	double scale; //!< overall scale factor of the data array
	const GridInfo& gInfo; //!< simulation grid info

	void absorbScale() const; //!< absorb scale factor into data
	void zero(); //!< initialize to zero

	typedef void DataType; //!< this base class has no specific data type
	
	void* data(bool shouldAbsorbScale=true); //!< get a pointer to the actual data (after absorbing the scale factor, unless otherwise specified)
	const void* data(bool shouldAbsorbScale=true) const; //!< get a const pointer to the actual data (after absorbing the scale factor, unless otherwise specified)
	#ifdef GPU_ENABLED
	void* dataGpu(bool shouldAbsorbScale=true); //!< get a pointer to the actual data (after absorbing the scale factor, unless otherwise specified)
	const void* dataGpu(bool shouldAbsorbScale=true) const; //!< get a const pointer to the actual data (after absorbing the scale factor, unless otherwise specified)
	#endif

	DECLARE_DATA_PREF_ACCESS
	
	bool isOnGpu() const { return onGpu; } //!< Check where the data is (for #ifdef simplicity exposed even when no GPU_ENABLED)

	Data(const GridInfo& gInfo, int nElem, int nDoublesPerElem, bool onGpu);
	~Data();
	void copyData(const Data& other); //!< copy data and scale (used by clone())

	//Inter-process communication:
	void send(int dest, int tag=0) const; //send to another process
	void recv(int src, int tag=0); //receive from another process
	void bcast(int root=0); //synchronize across processes (using value on specified root process)
	void allReduce(MPIUtil::ReduceOp op, bool safeMode=false); //apply all-to-all reduction (see MPIUtil::allReduce)

private:
	int nDoubles; //!< number of doubles stored in pData (typically either gInfo.nr or ginfo.nG*2)

	void* pData; //!< managed data pointer
	bool onGpu;
	#ifdef GPU_ENABLED
	void toCpu(); //!< move data to the CPU
	void toGpu(); //!< move data to the GPU
	#endif
};

//Shorthand for defining the data() and dataGpu() functions in derived classes of Data
#ifdef GPU_ENABLED
	#define DECLARE_DATA_ACCESS \
		DataType* data(bool shouldAbsorbScale=true) { return (DataType*)Data::data(shouldAbsorbScale); } \
		const DataType* data(bool shouldAbsorbScale=true) const { return (const DataType*)Data::data(shouldAbsorbScale); } \
		DataType* dataGpu(bool shouldAbsorbScale=true) { return (DataType*)Data::dataGpu(shouldAbsorbScale); } \
		const DataType* dataGpu(bool shouldAbsorbScale=true) const { return (const DataType*)Data::dataGpu(shouldAbsorbScale); }
#else
	#define DECLARE_DATA_ACCESS \
		DataType* data(bool shouldAbsorbScale=true) { return (DataType*)Data::data(shouldAbsorbScale); } \
		const DataType* data(bool shouldAbsorbScale=true) const { return (const DataType*)Data::data(shouldAbsorbScale); }
#endif

/**
@brief Real space real scalar field data
Do not use this data structure directly or from a simple pointer
DataR*; work only with #DataRptr's. The public functions of DataR
can be accessed with -> from the DataRptr.
*/
struct DataR : public Data
{	typedef double DataType; //!< Type of data in container (useful for templating)
	DECLARE_DATA_ACCESS
	DECLARE_DATA_PREF_ACCESS
	DataRptr clone() const; //!< clone the data (NOTE: assigning DataRptr's makes a new reference to the same data)
	static DataRptr alloc(const GridInfo& gInfo, bool onGpu=false); //!< Create real space data
private:
	DataR(const GridInfo& gInfo, bool onGpu); //!< called only by DataR::alloc()
};

//Override allReduce for DataType=complex to prevent unsupported operations
#define OVERRIDE_allReduce \
	void allReduce(MPIUtil::ReduceOp op, bool safeMode=false) \
	{	assert(op!=MPIUtil::ReduceProd && op!=MPIUtil::ReduceMax && op!=MPIUtil::ReduceMin); \
		Data::allReduce(op, safeMode); \
	}

/**
@brief Reciprocal space real scalar field data
Do not use this data structure directly or from a simple pointer
DataG*; work only with #DataGptr's. The public functions of DataG
can be accessed with -> from the DataGptr.
*/
struct DataG : public Data
{	typedef complex DataType; //!< Type of data in container (useful for templating)
	DECLARE_DATA_ACCESS
	DECLARE_DATA_PREF_ACCESS
	DataGptr clone() const; //!< clone the data (NOTE: assigning DataGptr's makes a new reference to the same data)
	static DataGptr alloc(const GridInfo& gInfo, bool onGpu=false); //!< Create reciprocal space data
	double getGzero() const; //!< get the G=0 component
	void setGzero(double Gzero); //!< set the G=0 component
	OVERRIDE_allReduce
private:
	DataG(const GridInfo& gInfo, bool onGpu); //!< called only by DataG::alloc()
};


/**
@brief Real space complex scalar field data
Do not use this data structure directly or from a simple pointer
complexDataR*; work only with #complexDataRptr's. The public functions of complexDataR
can be accessed with -> from the complexDataRptr.
*/
struct complexDataR : public Data
{	typedef complex DataType; //!< Type of data in container (useful for templating)
	DECLARE_DATA_ACCESS
	DECLARE_DATA_PREF_ACCESS
	complexDataRptr clone() const; //!< clone the data (NOTE: assigning complexDataRptr's makes a new reference to the same data)
	static complexDataRptr alloc(const GridInfo& gInfo, bool onGpu=false); //!< Create real space data
	OVERRIDE_allReduce
private:
	complexDataR(const GridInfo& gInfo, bool onGpu); //!< called only by complexDataR::alloc()
};

/**
@brief Reciprocal space complex scalar field data
Do not use this data structure directly or from a simple pointer
complexDataG*; work only with #complexDataGptr's. The public functions of complexDataG
can be accessed with -> from the complexDataGptr.
*/
struct complexDataG : public Data
{	typedef complex DataType; //!< Type of data in container (useful for templating)
	DECLARE_DATA_ACCESS
	DECLARE_DATA_PREF_ACCESS
	complexDataGptr clone() const; //!< clone the data (NOTE: assigning complexDataGptr's makes a new reference to the same data)
	static complexDataGptr alloc(const GridInfo& gInfo, bool onGpu=false); //!< Create reciprocal space data
	OVERRIDE_allReduce
private:
	complexDataG(const GridInfo& gInfo, bool onGpu); //!< called only by complexDataG::alloc()
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

#undef OVERRIDE_allReduce
#undef DECLARE_DATA_PREF_ACCESS
#undef DECLARE_DATA_ACCESS
//! @}

#endif //JDFTX_CORE_DATA_H

