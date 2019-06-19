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

#include <core/MPIUtil.h>
#include <core/Util.h>
#include <core/scalar.h>
#include <algorithm>
#include <climits>
#include <core/Random.h>

//---------- class MPIUtil::ProcDivision ----------

MPIUtil::ProcDivision::ProcDivision(const MPIUtil *mpiUtil, size_t nGroups, size_t iGroup)
: mpiUtil(mpiUtil), nGroups(nGroups), iGroup(
	(nGroups and mpiUtil)
	?  mpiUtil->iProcess() * nGroups / mpiUtil->nProcesses()
	: iGroup )
{
}


//---------- class MPIUtil ----------

MPIUtil::MPIUtil(int argc, char** argv, ProcDivision procDivision)
: procDivision(procDivision)
{
	#ifdef MPI_ENABLED

	if(procDivision)
	{	//Split parent communicator in procDivision using iGroup:
		MPI_Comm_split(procDivision.mpiUtil->comm, procDivision.iGroup,
			procDivision.mpiUtil->iProcess(), &comm);
	}
	else
	{	//Initialize MPI (use COMM_WORLD)
		int rc = MPI_Init(&argc, &argv);
		if(rc != MPI_SUCCESS) { printf("Error starting MPI program. Terminating.\n"); MPI_Abort(MPI_COMM_WORLD, rc); }
		comm = MPI_COMM_WORLD;
	}

	MPI_Comm_size(comm, &nProcs);
	MPI_Comm_rank(comm, &iProc);

	#else
	//No MPI:
	nProcs = 1;
	iProc = 0;
	#endif
	
	if(!procDivision) Random::seed(iProc); //Reproducible random seed per process in mpiWorld
}

MPIUtil::MPIUtil(const MPIUtil* mpiUtil, std::vector<int> ranks)
{
	#ifdef MPI_ENABLED
	#if MPI_VERSION < 3
		#define MPI_Comm_create_group MPIX_Comm_create_group //For older MPICH compatibility
	#endif
	//Create sub-communicator:
	MPI_Group parent, subset;
	MPI_Comm_group(mpiUtil->comm, &parent); //create group asociated with parent communicator
	MPI_Group_incl(parent, int(ranks.size()), ranks.data(), &subset); //create subgroup
	MPI_Comm_create_group(mpiUtil->comm, subset, 0, &comm);
	MPI_Group_free(&subset);
	MPI_Group_free(&parent);
	//Get rank and count within it:
	MPI_Comm_size(comm, &nProcs);
	MPI_Comm_rank(comm, &iProc);
	#else
	//No MPI:
	assert(ranks.size()==1 && ranks[0]==0);
	nProcs = 1;
	iProc = 0;
	#endif
}

MPIUtil::~MPIUtil()
{
	#ifdef MPI_ENABLED
	//Check if MPI has already finalized:
	int finalized;
	MPI_Finalized(&finalized);
	if(finalized) return; //to prevent double-free type errors
	//Finalize communicators or MPI as appropriate:
	if(comm == MPI_COMM_WORLD)
		MPI_Finalize();
	else
		MPI_Comm_free(&comm);
	#endif
}

void MPIUtil::exit(int errCode) const
{
	#ifdef MPI_ENABLED
	MPI_Abort(comm, errCode);
	#else
	::exit(errCode);
	#endif
}

void MPIUtil::checkErrors(const ostringstream& oss) const
{	const string buf = oss.str();
	int nChars = buf.length();
	//Find total length:
	int nCharsTot = nChars; allReduce(nCharsTot, ReduceSum);
	if(!nCharsTot) return; //no errors
	//Put together all error messages:
	string bufTot(nCharsTot, ' ');
	char* bufTotData = &bufTot[0];
	for(int jProcess=0; jProcess<nProcesses(); jProcess++)
	{	int nChars_j = nChars; bcast(nChars_j, jProcess);
		if(jProcess==iProcess())
			memcpy(bufTotData, &buf[0], nChars_j);
		bcast(bufTotData, nChars_j, jProcess);
		bufTotData += nChars_j;
	}
	bufTot += '\n';
	//Mimic the behaviour of die with collected error message:
	fputs(bufTot.c_str(), globalLog);
	if(isHead() && globalLog != stdout)
		fputs(bufTot.c_str(), stderr);
	finalizeSystem(false);
	::exit(1);
}

//-------------------- Asynchronous support functions ---------------------------

void MPIUtil::wait(MPIUtil::Request request)
{
#ifdef MPI_ENABLED
	MPI_Wait(&request, MPI_STATUS_IGNORE);
#endif
}

void MPIUtil::waitAll(const std::vector<Request>& requests)
{
#ifdef MPI_ENABLED
	MPI_Waitall(requests.size(), (Request*)requests.data(), MPI_STATUS_IGNORE);
#endif
}


//----------------------- Point-to-point routines -------------------------------

void MPIUtil::send(const bool* data, size_t nData, int dest, int tag, Request* request) const
{	if(request) throw string("Asynchronous send not supported for bool");
	std::vector<int> intCopy(nData);
	std::copy(data, data+nData, intCopy.begin());  //Copy data into an integer version
	send(&intCopy[0], nData, dest, tag); //Send integer version
}

void MPIUtil::recv(bool* data, size_t nData, int src, int tag, Request* request) const
{	if(request) throw string("Asynchronous recv not supported for bool");
	std::vector<int> intCopy(nData);
	recv(&intCopy[0], nData, src, tag); //Receive integer version of data
	std::copy(intCopy.begin(), intCopy.end(), data);  //Copy data to the target bool version
}

void MPIUtil::send(const string& s, int dest, int tag, Request* request) const
{	if(request) throw string("Asynchronous send not supported for string");
	unsigned long len = s.length();
	send(len, dest, tag);
	send(&s[0], len, dest, tag);
}

void MPIUtil::recv(string& s, int src, int tag, Request* request) const
{	if(request) throw string("Asynchronous recv not supported for string");
	unsigned long len = 0;
	recv(len, src, tag); s.resize(len);
	recv(&s[0], len, src, tag);
}


//----------------------- Broadcast routines -------------------------------

void MPIUtil::bcast(bool* data, size_t nData, int root, Request* request) const
{	if(nProcs>1)
	{	if(request) throw string("Asynchronous bcast not supported for bool");
		std::vector<int> intCopy(nData); //Copy data into an integer version (bool is not natively supported by MPI)
		std::copy(data, data+nData, intCopy.begin());
		bcast(&intCopy[0], nData, root);
		std::copy(intCopy.begin(), intCopy.end(), data);
	}
}

void MPIUtil::bcast(string& s, int root, Request* request) const
{	if(nProcs>1)
	{	if(request) throw string("Asynchronous bcast not supported for string");
		//Synchronize length of string:
		unsigned long len = s.length();
		bcast(len, root);
		if(iProc!=root) s.resize(len);
		//Bcast content:
		bcast(&s[0], len, root);
	}
}

//----------------------- Reduction routines -------------------------------

void MPIUtil::allReduce(bool* data, size_t nData, MPIUtil::ReduceOp op, bool safeMode, Request* request) const
{	if(nProcs>1)
	{	if(request) throw string("Asynchronous allReduce not supported for bool");
		std::vector<int> intCopy(nData); //Copy data into an integer version (bool is not natively supported by MPI)
		std::copy(data, data+nData, intCopy.begin());
		allReduce(&intCopy[0], nData, op);
		std::copy(intCopy.begin(), intCopy.end(), data);
	}
}

void MPIUtil::reduce(bool* data, size_t nData, MPIUtil::ReduceOp op, int root, Request* request) const
{	if(nProcs>1)
	{	if(request) throw string("Asynchronous reduce not supported for bool");
		std::vector<int> intCopy(nData); //Copy data into an integer version (bool is not natively supported by MPI)
		std::copy(data, data+nData, intCopy.begin());
		reduce(&intCopy[0], nData, op, root);
		std::copy(intCopy.begin(), intCopy.end(), data);
	}
}


//----------------------- File I/O routines -------------------------------

void MPIUtil::fopenRead(File& fp, const char* fname, size_t fsizeExpected, const char* fsizeErrMsg) const
{	if(fsizeExpected)
	{	intptr_t fsize = fileSize(fname);
		if(fsize < 0)
			die("Error opening file '%s' for reading.\n", fname);
		if(fsize != intptr_t(fsizeExpected))
			die("Length of '%s' was %" PRIdPTR " instead of the expected %zu bytes.\n%s\n", fname, fsize, fsizeExpected, fsizeErrMsg ? fsizeErrMsg : "");
	}
	#ifdef MPI_ENABLED
	if(MPI_File_open(comm, (char*)fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp) != MPI_SUCCESS)
	#else
	fp = ::fopen(fname, "rb");
	if(!fp)
	#endif
		die("Error opening file '%s' for reading.\n", fname);
}

void MPIUtil::fopenWrite(File& fp, const char* fname) const
{
	#ifdef MPI_ENABLED
	if(isHead()) MPI_File_delete((char*)fname, MPI_INFO_NULL); //delete existing file, if any
	MPI_Barrier(comm);
	if(MPI_File_open(comm, (char*)fname, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fp) != MPI_SUCCESS)
	#else
	fp = ::fopen(fname, "wb");
	if(!fp)
	#endif
		 die("Error opening file '%s' for writing.\n", fname);
}

void MPIUtil::fopenAppend(File& fp, const char* fname) const
{
	#ifdef MPI_ENABLED
	if(MPI_File_open(comm, (char*)fname, MPI_MODE_APPEND|MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fp) != MPI_SUCCESS)
	#else
	fp = ::fopen(fname, "a");
	if(!fp)
	#endif
		 die("Error opening file '%s' for writing.\n", fname);
	#ifdef MPI_ENABLED
	MPI_Barrier(comm);
	#endif
}

void MPIUtil::fclose(File& fp) const
{
	#ifdef MPI_ENABLED
	MPI_File_close(&fp);
	#else
	::fclose(fp);
	#endif
}

void MPIUtil::fseek(File fp, long offset, int whence) const
{
	#ifdef MPI_ENABLED
	int mpi_whence = 0;
	switch(whence)
	{	case SEEK_CUR: mpi_whence = MPI_SEEK_CUR; break;
		case SEEK_SET: mpi_whence = MPI_SEEK_SET; break;
		case SEEK_END: mpi_whence = MPI_SEEK_END; break;
		default: assert(!"Invalid seek offset mode");
	}
	if(MPI_File_seek(fp, offset, mpi_whence) != MPI_SUCCESS)
	#else
	if(::fseek(fp, offset, whence) != 0)
	#endif
		die("Error in file seek.\n");
}

void MPIUtil::fread(void *ptr, size_t size, size_t nmemb, File fp) const
{
	#ifdef MPI_ENABLED
	size_t blockSize = size_t(INT_MAX)/(2*size);
	size_t nBlocks = ceildiv(nmemb, blockSize);
	for(size_t iBlock=0; iBlock<nBlocks; iBlock++)
	{	size_t offset = iBlock * blockSize;
		int length = int(std::min(nmemb, (iBlock+1)*blockSize) - offset);
		MPI_Status status;
		char* ptrOffset = ((char*)ptr)+offset*size;
		MPI_File_read(fp, ptrOffset, length*size, MPI_BYTE, &status);
		convertFromLE(ptrOffset, size, length); //modify byte-order in place, if necessary
		int count; MPI_Get_count(&status, MPI_BYTE, &count);
		if(count != int(length*size)) die("Error in file read.\n");
	}
	#else
	if(freadLE(ptr, size, nmemb, fp) != nmemb)
		die("Error in file read.\n");
	#endif
}

void MPIUtil::fwrite(const void *ptr, size_t size, size_t nmemb, File fp) const
{
	#ifdef MPI_ENABLED
	size_t blockSize = size_t(INT_MAX)/(2*size);
	size_t nBlocks = ceildiv(nmemb, blockSize);
	for(size_t iBlock=0; iBlock<nBlocks; iBlock++)
	{	size_t offset = iBlock * blockSize;
		int length = int(std::min(nmemb, (iBlock+1)*blockSize) - offset);
		MPI_Status status;
		char* ptrOffset = ((char*)ptr)+offset*size;
		convertToLE(ptrOffset, size, length); //modify byte-order in place, if necessary
		MPI_File_write(fp, ptrOffset, length*size, MPI_BYTE, &status);
		convertFromLE(ptrOffset, size, length); //restore byte-order (since data should not be logically modified)
		int count; MPI_Get_count(&status, MPI_BYTE, &count);
		if(count != int(length*size)) die("Error in file write.\n");
	}
	#else
	if(fwriteLE(ptr, size, nmemb, fp) != nmemb)
		die("Error in file write.\n");
	#endif
}

//------- class TaskDivision ---------

TaskDivision::TaskDivision(size_t nTasks, const MPIUtil* mpiUtil)
: startMine(0), stopMine(nTasks), stopArr(1, nTasks) //initialize for no MPI division
{
	if(mpiUtil) init(nTasks, mpiUtil); //actually initialize if MPIUtil pointer provided
}

void TaskDivision::init(size_t nTasks, const MPIUtil* mpiUtil)
{	stopArr.resize(mpiUtil->nProcesses());
	for(int iProc=0; iProc<mpiUtil->nProcesses(); iProc++)
		stopArr[iProc] = (nTasks * (iProc+1)) / mpiUtil->nProcesses();
	startMine = start(mpiUtil->iProcess());
	stopMine = stop(mpiUtil->iProcess());
}

int TaskDivision::whose(size_t q) const
{	if(stopArr.size()>1)
		return std::upper_bound(stopArr.begin(),stopArr.end(), q) - stopArr.begin();
	else return 0;
}
