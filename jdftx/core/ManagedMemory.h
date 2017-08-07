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

#ifndef JDFTX_CORE_MANAGEDMEMORY_H
#define JDFTX_CORE_MANAGEDMEMORY_H

#include <core/Util.h>
#include <core/vector3.h>
#include <core/BlasExtra.h>

//! @addtogroup DataStructures
//! @{

//! @file ManagedMemory.h Base class and operators for managed-memory objects

//! Base class for managed-memory objects (that could potentially live on GPUs as well) with unspecified data type
class ManagedMemoryBase
{
public:
	static void reportUsage(); //!< print memory usage report

protected:
	ManagedMemoryBase(): nBytes(0),c(0),onGpu(false) {} //!< Initialize a valid state, but don't allocate anything
	~ManagedMemoryBase() { memFree(); }

	void memFree(); //!< Free memory
	void memInit(string category, size_t nBytes, bool onGpu=false); //!< Allocate memory
	void memMove(ManagedMemoryBase&&); //!< Steal the other object's data (used for move constructors/assignment)

	string category; //!< category of managed memory objects to report memory usage under
	size_t nBytes; //!< Size of stored data
	void* c; //!< Actual data storage
	bool onGpu; //!< For reduced #ifdef's, this flag is retained even in the absence of gpu support
	void toCpu() const; //!< move data to the CPU (does nothing without GPU_ENABLED); logically const, but data location may change
	void toGpu() const; //!< move data to the GPU (does nothing without GPU_ENABLED); logically const, but data location may change
};

//! Base class for managed memory of a specified data type
template<typename T> class ManagedMemory : private ManagedMemoryBase
{
protected:
	//Extend expose core functions from base:
    ManagedMemory() : nElem(0) {} //!< Initialize a valid state, but don't allocate anything
    ~ManagedMemory() { memFree(); }
	void memFree(); //!< Free memory
	void memInit(string category, size_t nElem, bool onGpu=false); //!< Allocate memory
	void memMove(ManagedMemory<T>&&); //!< Steal the other object's data (used for move constructors/assignment)

private:
	size_t nElem;

public:
	//! @brief Return a pointer to the actual data
	//! Return a CPU pointer to the actual data, will move data from GPU to CPU if necessary
	//! In GPU mode, care must be taken when calling this from multiple cpu threads
	//! Only the "GPU Owner" thread may call this when the data is actually on the GPU.
	//! Ideally call once from main thread to get data onto the cpu before launching other cpu threads
	T* data() { toCpu(); return (T*)c; }

	//! @brief Return a const pointer to the actual data
	//! Return a CPU pointer to the actual data, will move data from GPU to CPU if necessary
	//! In GPU mode, care must be taken when calling this from multiple cpu threads
	//! Only the "GPU Owner" thread may call this when the data is actually on the GPU.
	//! Ideally call once from main thread to get data onto the cpu before launching other cpu threads
	const T* data() const  { toCpu(); return (const T*)c; }

	#ifdef GPU_ENABLED
	T* dataGpu() { toGpu(); return (T*)c; } //!< Get a GPU data pointer (must be called from GPU owner thread)
	const T* dataGpu() const { toGpu(); return (const T*)c; } //!< Get a const GPU data pointer (must be called from GPU owner thread)
	#endif

	size_t nData() const { return nElem; } //!< number of data points
	bool isOnGpu() const { return onGpu; } //!< Check where the data is (for #ifdef simplicity exposed even when no GPU_ENABLED)

	//Iterator access on CPU:
	T* begin() { return data(); } //!< pointer to start of array
	const T* begin() const { return data(); } //!< const pointer to start of array
	const T* cbegin() const { return data(); } //!< const pointer to start of array
	T* end() { return data()+nElem; } //!< pointer just past end of array
	const T* end() const { return data()+nElem; } //!< const pointer just past end of array
	const T* cend() const { return data()+nElem; } //!< const pointer just past end of array

	//Utilities to automatically select "preferred" data i.e. GPU when it is enabled, CPU otherwise
	//This eases the overload of CPU/GPU functions based on complex vs complex
	#ifdef GPU_ENABLED
	inline T* dataPref() { return dataGpu(); }
	inline const T* dataPref() const { return dataGpu(); }
	#else
	inline T* dataPref() { return data(); }
	inline const T* dataPref() const { return data(); }
	#endif

	//Inter-process communication:
	void send(int dest, int tag=0) const; //!< send to another process
	void recv(int src, int tag=0); //!< receive from another process
	void bcast(int root=0); //!< synchronize across processes (using value on specified root process)
	void allReduce(MPIUtil::ReduceOp op, bool safeMode=false); //!< apply all-to-all reduction (see MPIUtil::allReduce)

	void read(const char *fname); //!< binary read from a file
	void read(FILE *filep); //!< binary read from a stream
	void write(const char *fname) const; //!< binary-write to a file
	void write(FILE *filep) const; //!< binary-write toa stream
	void read_real(const char *fname); //!< binary read real-part from file, setting imaginary parts to 0
	void read_real(FILE *filep); //!< binary read real-part from stream, setting imaginary parts to 0
	void write_real(const char *fname) const; //!< binary write real-parts to file
	void write_real(FILE *filep) const; //!< binary write real-parts to stream
	void dump(const char* fname, bool realPartOnly) const; //!< write as complex or real-part alone and report discarded imaginary part, if any
	void zero(); //!< set all elements to zero
};

//! General ManagedMemory class for one-off use as a CPU <-> GPU array.
//! For each commonly-used physical object, derive  directly from
//! ManagedMemory and implement operators; do not use this wrapper.
template<typename T> struct ManagedArray : public ManagedMemory<T>
{	void init(size_t size, bool onGpu=false); //!< calls memInit with category "misc"
	ManagedArray(const T* ptr=0, size_t N=0); //!< optionally initialize N elements from a pointer
	ManagedArray(const std::vector<T>&); //!< initialize from an std::vector
	ManagedArray& operator=(const ManagedArray&); //!< copy-assignment
	ManagedArray& operator=(ManagedArray&&); //!< move-assignment
};

//! Managed array of integers (indices)
struct IndexArray : public ManagedMemory<int>
{	void init(size_t size, bool onGpu=false) { memInit("IndexArrays", size, onGpu); }
};

//! Managed array of integer vectors
struct IndexVecArray : public ManagedMemory<vector3<int>>
{	void init(size_t size, bool onGpu=false) { memInit("IndexArrays", size, onGpu); }
};


//Some common elementwise / vector-like operations
template<typename T> void memcpy(ManagedMemory<T>&, const ManagedMemory<T>&); //!< copy entire object over
void scale(double alpha, ManagedMemory<double>& y); //! scale y *= alpha
void scale(double alpha, ManagedMemory<complex>& y); //! scale y *= alpha
void scale(complex alpha, ManagedMemory<complex>& y); //! scale y *= alpha
void scale(const ManagedMemory<complex>& x, ManagedMemory<complex>& y); //! scale y *= x elementwise
void axpy(complex alpha, const ManagedMemory<complex>& x, ManagedMemory<complex>& y); //!< standard blas y += alpha*x
double nrm2(const ManagedMemory<complex>&); //!< 2-norm, pretending it is a vector
complex dotc(const ManagedMemory<complex>& a, const ManagedMemory<complex>& b); //!< return a^H b

//! @}

//--------- Template implementations ---------
//! @cond
#include <type_traits>


template<typename T> void ManagedMemory<T>::memFree()
{	ManagedMemoryBase::memFree(); //first invoke base class version
	nElem = 0;
}

template<typename T> void ManagedMemory<T>::memInit(string category, size_t nElem, bool onGpu)
{	this->nElem = nElem;
	ManagedMemoryBase::memInit(category, nElem*sizeof(T), onGpu);
}

template<typename T> void ManagedMemory<T>::memMove(ManagedMemory<T>&& mOther)
{	ManagedMemoryBase::memMove((ManagedMemoryBase&&)mOther); //first invoke base class version
	std::swap(nElem, mOther.nElem);
}

template<typename T> void ManagedMemory<T>::read(const char *fname)
{	intptr_t fsizeExpected = nData() * sizeof(T);
	intptr_t fsize = fileSize(fname);
	if(fsize != fsizeExpected)
		die("Length of '%s' was %" PRIdPTR " instead of the expected %" PRIdPTR " bytes.\n", fname, fsize, fsizeExpected);
	FILE *fp = fopen(fname, "rb");
	if(!fp) die("Error opening %s for reading.\n", fname);
	read(fp);
	fclose(fp);
}
template<typename T> void ManagedMemory<T>::read(FILE *fp)
{	size_t nDone = freadLE(data(), sizeof(T), nData(), fp);
	if(nDone<nData()) die("Error after processing %lu of %lu records.\n", nDone, nData());
}

template<typename T> void ManagedMemory<T>::write(const char *fname) const
{	FILE *fp = fopen(fname,"wb");
	if(!fp) die("Error opening %s for writing.\n", fname);
	write(fp);
	fclose(fp);
}
template<typename T> void ManagedMemory<T>::write(FILE *fp) const
{	size_t nDone = fwriteLE(data(), sizeof(T), nData(), fp);
	if(nDone<nData()) die("Error after processing %lu of %lu records.\n", nDone, nData());
}

template<typename T> void ManagedMemory<T>::read_real(const char *fname)
{	assert((std::is_same<T,complex>::value));
	intptr_t fsizeExpected = nData() * sizeof(double);
	intptr_t fsize = fileSize(fname);
	if(fsize != fsizeExpected)
		die("Length of '%s' was %" PRIdPTR " instead of the expected %" PRIdPTR " bytes.\n", fname, fsize, fsizeExpected);
	FILE *fp = fopen(fname,"rb");
	read_real(fp);
	fclose(fp);
}
template<typename T> void ManagedMemory<T>::read_real(FILE *fp)
{	assert((std::is_same<T,complex>::value));
	double *dataReal = new double[nData()];
	freadLE(dataReal, sizeof(double), nData(), fp);
	complex* thisData = (complex*)this->data();
	for (size_t i=0; i<nData(); i++) thisData[i] = dataReal[i];
	delete[] dataReal;
}

template<typename T> void ManagedMemory<T>::write_real(const char *fname) const
{	assert((std::is_same<T,complex>::value));
	FILE *fp = fopen(fname,"wb");
	if(!fp) die("Error opening %s for writing.\n", fname);
	write_real(fp);
	fclose(fp);
}
template<typename T> void ManagedMemory<T>::write_real(FILE *fp) const
{	assert((std::is_same<T,complex>::value));
	const complex* thisData = (const complex*)this->data();
	double *dataReal = new double[nData()];
	for(size_t i=0; i<nData(); i++) dataReal[i] = thisData[i].real();
	fwriteLE(dataReal, sizeof(double), nData(), fp);
	delete[] dataReal;
}
template<typename T> void ManagedMemory<T>::dump(const char* fname, bool realPartOnly) const
{	logPrintf("Dumping '%s' ... ", fname); logFlush();
	if(realPartOnly)
	{	assert((std::is_same<T,complex>::value));
		write_real(fname);
		//Collect imaginary part:
		double nrm2tot = nrm2(*this); 
		double nrm2im = callPref(eblas_dnrm2)(nData(), ((double*)dataPref())+1, 2); //look only at imaginary parts with a stride of 2
		logPrintf("done. Relative discarded imaginary part: %le\n", nrm2im / nrm2tot);
	}
	else
	{	write(fname);
		logPrintf("done.\n");
	}
}

template<typename T> void ManagedMemory<T>::zero()
{	callPref(eblas_zero)(nData(), dataPref());
}

//Which data to use for MPI operations:
#if defined(GPU_ENABLED) && defined(CUDA_AWARE_MPI)
#define dataMPI dataGpu
#else
#define dataMPI data
#endif
template<typename T> void ManagedMemory<T>::send(int dest, int tag) const
{	assert(mpiUtil->nProcesses()>1);
	mpiUtil->send(dataMPI(), nData(), dest, tag);
}
template<typename T> void ManagedMemory<T>::recv(int src, int tag)
{	assert(mpiUtil->nProcesses()>1);
	mpiUtil->recv(dataMPI(), nData(), src, tag);
}
template<typename T> void ManagedMemory<T>::bcast(int root)
{	if(mpiUtil->nProcesses()>1)
		mpiUtil->bcast(dataMPI(), nData(), root);
}
template<typename T> void ManagedMemory<T>::allReduce(MPIUtil::ReduceOp op, bool safeMode)
{	if(mpiUtil->nProcesses()>1)
		mpiUtil->allReduce(dataMPI(), nData(), op, safeMode);
}
#undef dataMPI

template<typename T> void memcpy(ManagedMemory<T>& a, const ManagedMemory<T>& b)
{	assert(a.nData() == b.nData());
	if(!a.nData()) return; //no data to copy
	callPref(eblas_copy)(a.dataPref(), b.dataPref(), a.nData());
}

//---- Implementations of ManagedArray ----
template<typename T> void ManagedArray<T>::init(size_t size, bool onGpu)
{	ManagedMemory<T>::memInit("misc", size, onGpu);
}

template<typename T> ManagedArray<T>::ManagedArray(const T* ptr, size_t N)
{	if(ptr && N)
	{	init(N);
		eblas_copy(this->data(), ptr, N);
	}
}

template<typename T> ManagedArray<T>::ManagedArray(const std::vector<T>& v)
{	if(v.size())
	{	init(v.size());
		eblas_copy(this->data(), v.data(), v.size());
	}
}

template<typename T> ManagedArray<T>& ManagedArray<T>::operator=(const ManagedArray<T>& other)
{	if(other.nData())
	{	init(other.nData());
		memcpy(*this, other);
	}
	else ManagedMemory<T>::memFree();
	return *this;
}

template<typename T> ManagedArray<T>& ManagedArray<T>::operator=(ManagedArray<T>&& other)
{	ManagedMemory<T>::memMove((ManagedMemory<T>&&)other);
	return *this;
}

//! @endcond
#endif // JDFTX_CORE_MANAGEDMEMORY_H
