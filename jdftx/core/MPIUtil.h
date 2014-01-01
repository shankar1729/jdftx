/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman

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

#ifndef JDFTX_CORE_MPIUTIL_H
#define JDFTX_CORE_MPIUTIL_H

#include <core/string.h>
#include <cstdlib>
#include <cstdio>

#ifdef MPI_ENABLED
#include <mpi.h>
#endif

//! MPI wrapper class
class MPIUtil
{
	int nProcs, iProc;
public:
	int iProcess() const { return iProc; } //!< rank of current process
	int nProcesses() const { return nProcs; }  //!< number of processes
	bool isHead() const { return iProc==0; } //!< whether this is the root process (makes code more readable)

	MPIUtil(int argc, char** argv);
	~MPIUtil();
	void exit(int errCode) const; //!< global exit (kill other MPI processes as well)

	//Point-to-point functions:
	template<typename T> void send(const T* data, size_t nData, int dest, int tag) const; //!< generic array send
	template<typename T> void recv(T* data, size_t nData, int src, int tag) const; //!< generic array receive
	template<typename T> void send(const T& data, int dest, int tag) const; //!< generic scalar send
	template<typename T> void recv(T& data, int src, int tag) const; //!< generic scalar receive
	void send(const bool* data, size_t nData, int dest, int tag) const; //!< send specialization for bool which is not natively supported by MPI
	void recv(bool* data, size_t nData, int src, int tag) const; //!< receive specialization for bool which is not natively supported by MPI
	void send(const string& s, int dest, int tag) const; //!< send string
	void recv(string& s, int src, int tag) const; //!< send string
	
	//Broadcast functions:
	template<typename T> void bcast(T* data, size_t nData, int root=0) const; //!< generic array broadcast
	template<typename T> void bcast(T& data, int root=0) const; //!< generic scalar broadcast
	void bcast(bool* data, size_t nData, int root=0) const; //!< specialization for bool which is not natively supported by MPI
	void bcast(string& s, int root=0) const; //!< broadcast string

	//Reduce functions (safe mode gaurantees identical results irrespective of round-off (but could be slower)):
	enum ReduceOp { ReduceMin, ReduceMax, ReduceSum, ReduceProd, ReduceLAnd, ReduceBAnd, ReduceLOr, ReduceBOr, ReduceLXor, ReduceBXor };
	template<typename T> void allReduce(T* data, size_t nData, ReduceOp op, bool safeMode=false) const; //!< generic array reduction
	template<typename T> void allReduce(T& data, ReduceOp op, bool safeMode=false) const; //!< generic scalar reduction
	void allReduce(bool* data, size_t nData, ReduceOp op, bool safeMode=false) const;  //!< specialization for bool which is not natively supported by MPI

	//File access (tiny subset of MPI-IO, using byte offsets alone, and made to closely resemble stdio):
	#ifdef MPI_ENABLED
	typedef MPI_File File;
	#else
	typedef FILE* File;
	#endif
	void fopenRead(File& fp, const char* fname, size_t fsizeExpected=0, const char* fsizeErrMsg=0) const; //!< open file for reading and optionally check file size
	void fopenWrite(File& fp, const char* fname) const; //!< open file for writing
	void fclose(File& fp) const;
	void fseek(File fp, long offset, int whence) const; //!< syntax consistent with fseek from stdio
	void fread(void *ptr, size_t size, size_t nmemb, File fp) const;
	void fwrite(const void *ptr, size_t size, size_t nmemb, File fp) const;
};

//-------------------------- Template implementations ------------------------------------
//!@cond
namespace MPIUtilPrivate
{
#ifdef MPI_ENABLED
	template<typename T> struct DataType;
	#define DECLARE_DataType(cName, mpiName) template<> struct DataType<cName> { static MPI_Datatype get() { return MPI_##mpiName; } };
	DECLARE_DataType(char, CHAR)
	DECLARE_DataType(unsigned char, UNSIGNED_CHAR)
	DECLARE_DataType(short, SHORT)
	DECLARE_DataType(int, INT)
	DECLARE_DataType(long, LONG)
	DECLARE_DataType(unsigned long, UNSIGNED_LONG)
	DECLARE_DataType(float, FLOAT)
	DECLARE_DataType(double, DOUBLE)
	#undef DECLARE_DataType
	
	static MPI_Op mpiOp(MPIUtil::ReduceOp op)
	{	switch(op)
		{	case MPIUtil::ReduceMax: return MPI_MAX;
			case MPIUtil::ReduceMin: return MPI_MIN;
			case MPIUtil::ReduceSum: return MPI_SUM;
			case MPIUtil::ReduceProd: return MPI_PROD;
			case MPIUtil::ReduceLAnd: return MPI_LAND;
			case MPIUtil::ReduceBAnd: return MPI_BAND;
			case MPIUtil::ReduceLOr: return MPI_LOR;
			case MPIUtil::ReduceBOr: return MPI_BOR;
			case MPIUtil::ReduceLXor: return MPI_LXOR;
			case MPIUtil::ReduceBXor: return MPI_BXOR;
		}
		return 0;
	}
#endif
}

template<typename T> void MPIUtil::send(const T* data, size_t nData, int dest, int tag) const
{	using namespace MPIUtilPrivate;
	#ifdef MPI_ENABLED
	if(nProcs>1) MPI_Send((T*)data, nData, DataType<T>::get(), dest, tag, MPI_COMM_WORLD);
	#endif
}

template<typename T> void MPIUtil::recv(T* data, size_t nData, int src, int tag) const
{	using namespace MPIUtilPrivate;
	#ifdef MPI_ENABLED
	if(nProcs>1) MPI_Recv(data, nData, DataType<T>::get(), src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	#endif
}

template<typename T> void MPIUtil::send(const T& data, int dest, int tag) const
{	send(&data, 1, dest, tag);
}

template<typename T> void MPIUtil::recv(T& data, int src, int tag) const
{	recv(&data, 1, src, tag);
}

template<typename T> void MPIUtil::bcast(T* data, size_t nData, int root) const
{	using namespace MPIUtilPrivate;
	#ifdef MPI_ENABLED
	if(nProcs>1) MPI_Bcast(data, nData, DataType<T>::get(), root, MPI_COMM_WORLD);
	#endif
}

template<typename T> void MPIUtil::bcast(T& data, int root) const
{	bcast(&data, 1, root);
}

template<typename T> void MPIUtil::allReduce(T* data, size_t nData, MPIUtil::ReduceOp op, bool safeMode) const
{	using namespace MPIUtilPrivate;
	#ifdef MPI_ENABLED
	if(nProcs>1)
	{	if(safeMode) //Reduce to root node and then broadcast result (to ensure identical values)
		{	MPI_Reduce(isHead()?MPI_IN_PLACE:data, data, nData, DataType<T>::get(), mpiOp(op), 0, MPI_COMM_WORLD);
			bcast(data, nData, 0);
		}
		else //standard Allreduce
			MPI_Allreduce(MPI_IN_PLACE, data, nData, DataType<T>::get(), mpiOp(op), MPI_COMM_WORLD);
	}
	#endif
}

template<typename T> void MPIUtil::allReduce(T& data, MPIUtil::ReduceOp op, bool safeMode) const
{	allReduce(&data, 1, op, safeMode);
}

//!@endcond
#endif // JDFTX_CORE_MPIUTIL_H
