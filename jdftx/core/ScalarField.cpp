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

#include <core/ScalarField.h>
#include <core/GridInfo.h>
#include <core/GpuUtil.h>
#include <core/BlasExtra.h>
#include <string.h>

//--------------- Base class FieldData -------------------------

FieldData::FieldData(const GridInfo& gInfo, string category, int nElem, int nDoublesPerElem, bool onGpu)
: nElem(nElem), scale(1.), gInfo(gInfo), isReal(nDoublesPerElem==1)
{
	int nDoubles = nElem * nDoublesPerElem;
	int nComplex = ceildiv(nDoubles,2);
	memInit(category, nComplex, onGpu);
}

void FieldData::copyData(const FieldData& other)
{	scale = other.scale;
	memcpy((ManagedMemory&)(*this), other);
}

void FieldData::absorbScale() const
{	if(scale != 1.)
	{	FieldData& X = (FieldData&)(*this); //cast to non-const (this function modifies data, but is logically constant)
		::scale(scale, (ManagedMemory&)X);
		X.scale = 1.;
	}
}

void* FieldData::data(bool shouldAbsorbScale)
{	if(shouldAbsorbScale) absorbScale();
	return ManagedMemory::data();
}
const void* FieldData::data(bool shouldAbsorbScale) const
{	if(shouldAbsorbScale) absorbScale();
	return ManagedMemory::data();
}
#ifdef GPU_ENABLED
void* FieldData::dataGpu(bool shouldAbsorbScale)
{	if(shouldAbsorbScale) absorbScale();
	return ManagedMemory::dataGpu();
}
const void* FieldData::dataGpu(bool shouldAbsorbScale) const
{	if(shouldAbsorbScale) absorbScale();
	return ManagedMemory::dataGpu();
}
#endif

//------------ class ScalarFieldData ---------------

ScalarFieldData::ScalarFieldData(const GridInfo& gInfo, bool onGpu, PrivateTag) : FieldData(gInfo, "ScalarField", gInfo.nr, 1, onGpu)
{
}
ScalarField ScalarFieldData::clone() const
{	ScalarField copy(ScalarFieldData::alloc(gInfo, isOnGpu()));
	copy->copyData(*this);
	return copy;
}
ScalarField ScalarFieldData::alloc(const GridInfo& gInfo, bool onGpu) { return std::make_shared<ScalarFieldData>(gInfo, onGpu, PrivateTag()); }


//------------ class ScalarFieldTildeData ---------------

ScalarFieldTildeData::ScalarFieldTildeData(const GridInfo& gInfo, bool onGpu, PrivateTag) : FieldData(gInfo, "ScalarFieldTilde", gInfo.nG, 2, onGpu)
{
}
ScalarFieldTilde ScalarFieldTildeData::clone() const
{	ScalarFieldTilde copy(ScalarFieldTildeData::alloc(gInfo, isOnGpu()));
	copy->copyData(*this);
	return copy;
}
ScalarFieldTilde ScalarFieldTildeData::alloc(const GridInfo& gInfo, bool onGpu) { return std::make_shared<ScalarFieldTildeData>(gInfo, onGpu, PrivateTag()); }

double ScalarFieldTildeData::getGzero() const
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
void ScalarFieldTildeData::setGzero(double Gzero)
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


//------------ class complexScalarFieldData ---------------

complexScalarFieldData::complexScalarFieldData(const GridInfo& gInfo, bool onGpu, PrivateTag) : FieldData(gInfo, "complexScalarField", gInfo.nr, 2, onGpu)
{
}
complexScalarField complexScalarFieldData::clone() const
{	complexScalarField copy(complexScalarFieldData::alloc(gInfo, isOnGpu()));
	copy->copyData(*this);
	return copy;
}
complexScalarField complexScalarFieldData::alloc(const GridInfo& gInfo, bool onGpu) { return std::make_shared<complexScalarFieldData>(gInfo, onGpu, PrivateTag()); }

//------------ class complexScalarFieldTildeData ---------------

complexScalarFieldTildeData::complexScalarFieldTildeData(const GridInfo& gInfo, bool onGpu, PrivateTag) : FieldData(gInfo, "complexScalarFieldTilde", gInfo.nr, 2, onGpu)
{
}
complexScalarFieldTilde complexScalarFieldTildeData::clone() const
{	complexScalarFieldTilde copy(complexScalarFieldTildeData::alloc(gInfo, isOnGpu()));
	copy->copyData(*this);
	return copy;
}
complexScalarFieldTilde complexScalarFieldTildeData::alloc(const GridInfo& gInfo, bool onGpu) { return std::make_shared<complexScalarFieldTildeData>(gInfo, onGpu, PrivateTag()); }


//------------ class RealKernel ---------------

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
