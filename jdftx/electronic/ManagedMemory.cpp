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

#include <electronic/ManagedMemory.h>
#include <core/BlasExtra.h>
#include <core/GpuUtil.h>
#include <fftw3.h>

// Construct, optionally with data allocation
ManagedMemory::ManagedMemory(size_t nElements, bool onGpu)
: nElements(0),c(0),onGpu(false)
{	memInit(nElements,onGpu);
}

ManagedMemory::~ManagedMemory()
{	memFree();
}

// Do the job of destructor
void ManagedMemory::memFree()
{	if(!nElements) return; //nothing to free
	if(onGpu)
	{
		#ifdef GPU_ENABLED
		assert(isGpuMine());
		cudaFree(c);
		gpuErrorCheck();
		#else
		assert(!"onGpu=true without GPU_ENABLED"); //Should never get here!
		#endif
	}
	else
	{	fftw_free(c);
	}
	c = 0;
	nElements = 0;
}

// Do the job of the constructor
void ManagedMemory::memInit(size_t nElem, bool ontoGpu)
{	if(nElements==nElem && onGpu==ontoGpu) return; //already in required state
	memFree();
	nElements = nElem;
	onGpu = ontoGpu;
	if(onGpu)
	{
		#ifdef GPU_ENABLED
		assert(isGpuMine());
		cudaMalloc(&c, sizeof(complex)*nElements);
		gpuErrorCheck();
		#else
		assert(!"onGpu=true without GPU_ENABLED");
		#endif
	}
	else
	{	c = (complex*)fftw_malloc(sizeof(complex)*nElements);
		if(!c) die("Memory allocation failed (out of memory)\n");
	}
}

void ManagedMemory::memMove(ManagedMemory&& mOther)
{	std::swap(nElements, mOther.nElements);
	std::swap(onGpu, mOther.onGpu);
	std::swap(c, mOther.c);
	//Now mOther will be empty, while *this will have all its contents
}

complex* ManagedMemory::data()
{
	#ifdef GPU_ENABLED
	toCpu();
	#endif
	return c;
}

const complex* ManagedMemory::data() const
{
	#ifdef GPU_ENABLED
	((ManagedMemory*)this)->toCpu(); //logically const, but may change data location
	#endif
	return c;
}



#ifdef GPU_ENABLED

complex* ManagedMemory::dataGpu()
{	toGpu();
	return (complex*)c;
}

const complex* ManagedMemory::dataGpu() const
{	((ManagedMemory*)this)->toGpu(); //logically const, but may change data location
	return (complex*)c;
}

//Move data to CPU
void ManagedMemory::toCpu()
{	if(!onGpu || !c) return; //already on cpu, or no data
	assert(isGpuMine());
	complex* cCpu = (complex*)fftw_malloc(sizeof(complex)*nElements);
	if(!cCpu) die("Memory allocation failed (out of memory)\n");
	cudaMemcpy(cCpu, c, sizeof(complex)*nElements, cudaMemcpyDeviceToHost); gpuErrorCheck();
	cudaFree(c); gpuErrorCheck(); //Free GPU mem
	c = cCpu; //Make c a cpu pointer
	onGpu = false;
	//printf("\nDevice -> Host (%d elements at ptr %lX):\n", nElements, (unsigned long)c); printStack();
}

// Move data to GPU
void ManagedMemory::toGpu()
{	if(onGpu || !c) return; //already on gpu, or no data
	assert(isGpuMine());
	complex* cGpu;
	cudaMalloc(&cGpu, sizeof(complex)*nElements); gpuErrorCheck();
	cudaMemcpy(cGpu, c, sizeof(complex)*nElements, cudaMemcpyHostToDevice);
	fftw_free(c); //Free CPU mem
	c = cGpu; //Make c a gpu pointer
	onGpu = true;
	//printf("\nHost -> Device (%d elements at ptr %lX)\n", nElements, (unsigned long)c); printStack();
}

#endif

void ManagedMemory::send(int dest, int tag) const
{	assert(mpiUtil->nProcesses()>1);
	mpiUtil->send((const double*)data(), 2*nData(), dest, tag);
}
void ManagedMemory::recv(int src, int tag)
{	assert(mpiUtil->nProcesses()>1);
	mpiUtil->recv((double*)data(), 2*nData(), src, tag);
}
void ManagedMemory::bcast(int root)
{	if(mpiUtil->nProcesses()>1)
		mpiUtil->bcast((double*)data(), 2*nData(), root);
}
void ManagedMemory::allReduce(MPIUtil::ReduceOp op, bool safeMode)
{	assert(op!=MPIUtil::ReduceProd && op!=MPIUtil::ReduceMax && op!=MPIUtil::ReduceMin); //not supported for complex
	if(mpiUtil->nProcesses()>1)
		mpiUtil->allReduce((double*)data(), 2*nData(), op, safeMode);
}


void ManagedMemory::write(const char *fname) const
{	FILE *fp = fopen(fname,"wb");
	if(!fp) die("Error opening %s for writing.\n", fname);
	write(fp);
	fclose(fp);
}
void ManagedMemory::writea(const char *fname) const
{	FILE *fp = fopen(fname,"ab");
	if(!fp) die("Error opening %s for appending.\n", fname);
	write(fp);
	fclose(fp);
}
void ManagedMemory::write(FILE *fp) const
{	size_t nDone = fwrite(data(), sizeof(complex), nData(), fp);
	if(nDone<nData()) die("Error after processing %lu of %lu records.\n", nDone, nData());
}


void ManagedMemory::read(const char *fname)
{	FILE *fp = fopen(fname, "rb");
	if(!fp) die("Error opening %s for reading.\n", fname);
	read(fp);
	fclose(fp);
}
void ManagedMemory::read(FILE *fp)
{	size_t nDone = fread(data(), sizeof(complex), nData(), fp);
	if(nDone<nData()) die("Error after processing %lu of %lu records.\n", nDone, nData());
}


void ManagedMemory::write_real(const char *fname) const
{	FILE *fp = fopen(fname,"wb");
	if(!fp) die("Error opening %s for writing.\n", fname);
	write_real(fp);
	fclose(fp);
}
void ManagedMemory::write_real(FILE *fp) const
{	const complex* thisData = this->data();
	double *dataReal = new double[nData()];
	for(size_t i=0; i<nData(); i++) dataReal[i] = thisData[i].real();
	fwrite(dataReal, sizeof(double), nData(), fp);
	delete[] dataReal;
}

void ManagedMemory::read_real(const char *fname)
{	FILE *fp = fopen(fname,"rb");
	read_real(fp);
	fclose(fp);
}
void ManagedMemory::read_real(FILE *fp)
{	double *dataReal = new double[nData()];
	fread(dataReal, sizeof(double), nData(), fp);
	complex* thisData = this->data();
	for (size_t i=0; i<nData(); i++) thisData[i] = dataReal[i];
	delete[] dataReal;
}

void ManagedMemory::zero()
{	callPref(eblas_zero)(nData(), dataPref());
}

void memcpy(ManagedMemory& a, const ManagedMemory& b)
{	assert(a.nData() == b.nData());
	if(!a.nData()) return; //no data to copy
	#ifdef GPU_ENABLED
	cudaMemcpy(a.dataGpu(), b.dataGpu(), a.nData()*sizeof(complex), cudaMemcpyDeviceToDevice);
	#else
	memcpy(a.data(), b.data(), a.nData()*sizeof(complex));
	#endif
}

void scale(double alpha, ManagedMemory& y)
{	callPref(eblas_zdscal)(y.nData(), alpha, y.dataPref(), 1);
}
void scale(complex alpha, ManagedMemory& y)
{	callPref(eblas_zscal)(y.nData(), alpha, y.dataPref(), 1);
}

void axpy(complex alpha, const ManagedMemory& x, ManagedMemory& y)
{	assert(x.nData() == y.nData());
	callPref(eblas_zaxpy)(x.nData(), alpha, x.dataPref(), 1, y.dataPref(), 1);
}

double nrm2(const ManagedMemory& a)
{	return callPref(eblas_dznrm2)(a.nData(), a.dataPref(), 1);
}

complex dotc(const ManagedMemory& a, const ManagedMemory& b)
{	assert(a.nData() == b.nData());
	return callPref(eblas_zdotc)(a.nData(), a.dataPref(), 1, b.dataPref(), 1);
}
