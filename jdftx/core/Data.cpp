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

#include <core/Data.h>
#include <core/GridInfo.h>
#include <core/GpuUtil.h>
#include <core/BlasExtra.h>
#include <string.h>

Data::Data(const GridInfo& gInfo, int nElem, int nDoublesPerElem, bool onGpu)
:nElem(nElem),scale(1.0),gInfo(gInfo),nDoubles(nElem*nDoublesPerElem),onGpu(onGpu)
{	if(onGpu)
	{
		#ifdef GPU_ENABLED
		assert(isGpuMine()); //Cannot allocate GPU memory from non-gpu-owner thread
		cudaMalloc(&pData, nDoubles*sizeof(double));
		gpuErrorCheck();
		#else
		assert(!"In Data(), wound up with onGpu=true with no GPU_ENABLED!!!\n");
		#endif
	}
	else
	{	pData = fftw_malloc(nDoubles*sizeof(double));
		if(!pData) die("Memory allocation failed (out of memory)\n");
	}
}
Data::~Data()
{
	if(onGpu)
	{
		#ifdef GPU_ENABLED
		assert(isGpuMine()); //Cannot free GPU memory from non-gpu-owner thread
		cudaFree(pData);
		gpuErrorCheck();
		#else
		assert(!"In ~Data(), wound up with onGpu=true with no GPU_ENABLED!!!\n");
		#endif
	}
	else fftw_free(pData);
}
void Data::copyData(const Data& other)
{	scale = other.scale;
	#ifdef GPU_ENABLED
	cudaMemcpy(dataGpu(false), other.dataGpu(false), nDoubles*sizeof(double), cudaMemcpyDeviceToDevice);
	#else
	memcpy(data(false), other.data(false), nDoubles*sizeof(double));
	#endif
}

void Data::send(int dest, int tag) const
{	assert(mpiUtil->nProcesses()>1);
	mpiUtil->send((const double*)data(), nDoubles, dest, tag);
}
void Data::recv(int src, int tag)
{	assert(mpiUtil->nProcesses()>1);
	mpiUtil->recv((double*)data(), nDoubles, src, tag);
}
void Data::bcast(int root)
{	if(mpiUtil->nProcesses()>1)
		mpiUtil->bcast((double*)data(), nDoubles, root);
}
void Data::allReduce(MPIUtil::ReduceOp op, bool safeMode)
{	if(mpiUtil->nProcesses()>1)
		mpiUtil->allReduce((double*)data(), nDoubles, op, safeMode);
}

void Data::absorbScale() const
{	if(scale != 1.0)
	{	Data* X = (Data*)this; //cast to non-const (this function modifies data, but is logically constant)
		callPref(eblas_dscal)(nDoubles, scale, (double*)X->dataPref(false), 1);
		X->scale=1.0;
	}
}
void Data::zero()
{	scale=1.0;
	callPref(eblas_zero)(nDoubles, (double*)dataPref(false));
}
void* Data::data(bool shouldAbsorbScale)
{	if(shouldAbsorbScale) absorbScale();
	#ifdef GPU_ENABLED
	toCpu();
	#endif
	return pData;
}
const void* Data::data(bool shouldAbsorbScale) const
{	if(shouldAbsorbScale) absorbScale();
	#ifdef GPU_ENABLED
	((Data*)this)->toCpu(); //logically const, but may change data location
	#endif
	return pData;
}
#ifdef GPU_ENABLED
void* Data::dataGpu(bool shouldAbsorbScale)
{	toGpu();
	if(shouldAbsorbScale) absorbScale();
	return pData;
}
const void* Data::dataGpu(bool shouldAbsorbScale) const
{	((Data*)this)->toGpu(); //logically const, but may change data location
	if(shouldAbsorbScale) absorbScale();
	return pData;
}
//Move data to CPU
void Data::toCpu()
{	if(!onGpu) return; //already on cpu
	assert(isGpuMine()); //Cannot initiate GPU->CPU transfer from non-gpu-owner thread
	complex* pDataCpu = (complex*)fftw_malloc(nDoubles*sizeof(double));
	if(!pDataCpu) die("Memory allocation failed (out of memory)\n");
	cudaMemcpy(pDataCpu, pData, nDoubles*sizeof(double), cudaMemcpyDeviceToHost); gpuErrorCheck();
	cudaFree(pData); gpuErrorCheck(); //Free GPU mem
	pData = pDataCpu; //Make pData a cpu pointer
	onGpu = false;
}
// Move data to GPU
void Data::toGpu()
{	if(onGpu) return; //already on gpu
	assert(isGpuMine()); //Cannot initiate CPU->GPU transfer from non-gpu-owner thread
	void* pDataGpu; cudaMalloc(&pDataGpu, nDoubles*sizeof(double)); gpuErrorCheck();
	cudaMemcpy(pDataGpu, pData, nDoubles*sizeof(double), cudaMemcpyHostToDevice);
	fftw_free(pData);
	pData = pDataGpu; //Make pData a gpu pointer
	onGpu = true;
}
#endif



DataR::DataR(const GridInfo& gInfo, bool onGpu) : Data(gInfo, gInfo.nr, 1, onGpu)
{
}
DataRptr DataR::clone() const
{	DataRptr copy(DataR::alloc(gInfo, isOnGpu()));
	copy->copyData(*this);
	return copy;
}
DataRptr DataR::alloc(const GridInfo& gInfo, bool onGpu) { return DataRptr(new DataR(gInfo, onGpu)); }



DataG::DataG(const GridInfo& gInfo, bool onGpu) : Data(gInfo, gInfo.nG, 2, onGpu)
{
}
DataGptr DataG::clone() const
{	DataGptr copy(DataG::alloc(gInfo, isOnGpu()));
	copy->copyData(*this);
	return copy;
}
DataGptr DataG::alloc(const GridInfo& gInfo, bool onGpu) { return DataGptr(new DataG(gInfo, onGpu)); }

double DataG::getGzero() const
{	
	#ifdef GPU_ENABLED
	if(isOnGpu())
	{	double ret;
		cudaMemcpy(&ret, dataGpu(false), sizeof(double), cudaMemcpyDeviceToHost);
		return ret * scale;
	}
	#endif
	return data(false)[0].real() * scale;
}
void DataG::setGzero(double Gzero)
{	if(!scale) absorbScale(); //otherwise division by zero below
	double scaledGzero = Gzero / scale;
	#ifdef GPU_ENABLED
	if(isOnGpu())
	{	cudaMemcpy(dataGpu(false), &scaledGzero, sizeof(double), cudaMemcpyHostToDevice);
		return;
	}
	#endif
	data(false)[0].real() = scaledGzero;
}



complexDataR::complexDataR(const GridInfo& gInfo, bool onGpu) : Data(gInfo, gInfo.nr, 2, onGpu)
{
}
complexDataRptr complexDataR::clone() const
{	complexDataRptr copy(complexDataR::alloc(gInfo, isOnGpu()));
	copy->copyData(*this);
	return copy;
}
complexDataRptr complexDataR::alloc(const GridInfo& gInfo, bool onGpu) { return complexDataRptr(new complexDataR(gInfo, onGpu)); }


complexDataG::complexDataG(const GridInfo& gInfo, bool onGpu) : Data(gInfo, gInfo.nr, 2, onGpu)
{
}
complexDataGptr complexDataG::clone() const
{	complexDataGptr copy(complexDataG::alloc(gInfo, isOnGpu()));
	copy->copyData(*this);
	return copy;
}
complexDataGptr complexDataG::alloc(const GridInfo& gInfo, bool onGpu) { return complexDataGptr(new complexDataG(gInfo, onGpu)); }




RealKernel::RealKernel(const GridInfo& gInfo) : gInfo(gInfo), nElem(gInfo.nG)
{	data = new double[nElem];
	#ifdef GPU_ENABLED
	cudaMalloc(&dataGpu, nElem*sizeof(double));
	dataPref = dataGpu;
	#else
	dataPref = data;
	#endif
}
RealKernel::~RealKernel()
{	delete[] data;
	#ifdef GPU_ENABLED
	cudaFree(dataGpu);
	#endif
}
void RealKernel::set()
{
	#ifdef GPU_ENABLED
	cudaMemcpy(dataGpu, data, nElem*sizeof(double), cudaMemcpyHostToDevice);
	#endif
}
