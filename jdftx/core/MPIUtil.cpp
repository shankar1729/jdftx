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
#include <vector>
#include <algorithm>

MPIUtil::MPIUtil(int argc, char** argv)
{
	#ifdef MPI_ENABLED
	int rc = MPI_Init(&argc, &argv);
	if(rc != MPI_SUCCESS) { printf("Error starting MPI program. Terminating.\n"); MPI_Abort(MPI_COMM_WORLD, rc); }
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &iProc);
	#else
	//No MPI:
	nProcs = 1;
	iProc = 0;
	#endif
}

MPIUtil::~MPIUtil()
{
	#ifdef MPI_ENABLED
	MPI_Finalize();
	#endif
}

void MPIUtil::exit(int errCode) const
{
	#ifdef MPI_ENABLED
	MPI_Abort(MPI_COMM_WORLD, errCode);
	#else
	exit(errCode);
	#endif
}

//----------------------- Point-to-point routines -------------------------------

void MPIUtil::send(const bool* data, size_t nData, int dest, int tag) const
{	std::vector<int> intCopy(nData);
	std::copy(data, data+nData, intCopy.begin());  //Copy data into an integer version
	send(&intCopy[0], nData, dest, tag); //Send integer version
}

void MPIUtil::recv(bool* data, size_t nData, int src, int tag) const
{	std::vector<int> intCopy(nData);
	recv(&intCopy[0], nData, src, tag); //Receive integer version of data
	std::copy(intCopy.begin(), intCopy.end(), data);  //Copy data to the target bool version
}

void MPIUtil::send(const string& s, int dest, int tag) const
{	unsigned long len = s.length();
	send(len, dest, tag);
	send(&s[0], len, dest, tag);
}

void MPIUtil::recv(string& s, int src, int tag) const
{	unsigned long len = 0;
	recv(len, src, tag); s.resize(len);
	recv(&s[0], len, src, tag);
}


//----------------------- Broadcast routines -------------------------------

void MPIUtil::bcast(bool* data, size_t nData, int root) const
{	if(nProcs>1)
	{	std::vector<int> intCopy(nData); //Copy data into an integer version (bool is not natively supported by MPI)
		std::copy(data, data+nData, intCopy.begin());
		bcast(&intCopy[0], nData, root);
		std::copy(intCopy.begin(), intCopy.end(), data);
	}
}

void MPIUtil::bcast(string& s, int root) const
{	if(nProcs>1)
	{	//Synchronize length of string:
		unsigned long len = s.length();
		bcast(len, root);
		if(iProc!=root) s.resize(len);
		//Bcast content:
		bcast(&s[0], len, root);
	}
}

//----------------------- Reduction routines -------------------------------

void MPIUtil::allReduce(bool* data, size_t nData, MPIUtil::ReduceOp op, bool safeMode) const
{	if(nProcs>1)
	{	std::vector<int> intCopy(nData); //Copy data into an integer version (bool is not natively supported by MPI)
		std::copy(data, data+nData, intCopy.begin());
		allReduce(&intCopy[0], nData, op);
		std::copy(intCopy.begin(), intCopy.end(), data);
	}
}


//----------------------- File I/O routines -------------------------------

void MPIUtil::fopenRead(File& fp, const char* fname, size_t fsizeExpected, const char* fsizeErrMsg) const
{	if(fsizeExpected)
	{	off_t fsize = fileSize(fname);
		if(fsize != off_t(fsizeExpected))
			die("Length of '%s' was %ld instead of the expected %ld bytes.\n%s\n", fname, fsize, fsizeExpected, fsizeErrMsg ? fsizeErrMsg : "");
	}
	#ifdef MPI_ENABLED
	if(MPI_File_open(MPI_COMM_WORLD, (char*)fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp) != MPI_SUCCESS)
	#else
	fp = ::fopen(fname, "rb");
	if(!fp)
	#endif
		die("Error opening file '%s' for reading.\n", fname);
}

void MPIUtil::fopenWrite(File& fp, const char* fname) const
{
	#ifdef MPI_ENABLED
	if(MPI_File_open(MPI_COMM_WORLD, (char*)fname, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fp) != MPI_SUCCESS)
	#else
	fp = ::fopen(fname, "wb");
	if(!fp)
	#endif
		 die("Error opening file '%s' for writing.\n", fname);
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
	MPI_Status status;
	MPI_File_read(fp, ptr, size*nmemb, MPI_BYTE, &status);
	int count; MPI_Get_count(&status, MPI_BYTE, &count);
	if(size_t(count) != size*nmemb)
	#else
	if(::fread(ptr, size, nmemb, fp) != nmemb)
	#endif
		die("Error in file read.\n");
}

void MPIUtil::fwrite(const void *ptr, size_t size, size_t nmemb, File fp) const
{
	#ifdef MPI_ENABLED
	MPI_Status status;
	MPI_File_write(fp, (void*)ptr, size*nmemb, MPI_BYTE, &status);
	int count; MPI_Get_count(&status, MPI_BYTE, &count);
	if(size_t(count) != size*nmemb)
	#else
	if(::fwrite(ptr, size, nmemb, fp) != nmemb)
	#endif
		die("Error in file write.\n");
}
