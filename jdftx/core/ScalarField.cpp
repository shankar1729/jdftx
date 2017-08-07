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


//------------ class ScalarFieldData ---------------

ScalarFieldData::ScalarFieldData(const GridInfo& gInfo, bool onGpu, PrivateTag)
: FieldData<double>(gInfo, "ScalarField", gInfo.nr, onGpu)
{
}
ScalarField ScalarFieldData::clone() const
{	ScalarField copy(ScalarFieldData::alloc(gInfo, isOnGpu()));
	copy->copyData(*this);
	return copy;
}
ScalarField ScalarFieldData::alloc(const GridInfo& gInfo, bool onGpu) { return std::make_shared<ScalarFieldData>(gInfo, onGpu, PrivateTag()); }
 
matrix ScalarFieldData::toMatrix() const
{
  matrix mat = zeroes(gInfo.nr,1);
  callPref(eblas_daxpy)(gInfo.nr, 1.0, dataPref(), 1, (double*) mat.dataPref(), 2);
  return mat;
}



//------------ class ScalarFieldTildeData ---------------

ScalarFieldTildeData::ScalarFieldTildeData(const GridInfo& gInfo, bool onGpu, PrivateTag)
: FieldData<complex>(gInfo, "ScalarFieldTilde", gInfo.nG, onGpu)
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

complexScalarFieldData::complexScalarFieldData(const GridInfo& gInfo, bool onGpu, PrivateTag)
: FieldData<complex>(gInfo, "complexScalarField", gInfo.nr, onGpu)
{
}
complexScalarField complexScalarFieldData::clone() const
{	complexScalarField copy(complexScalarFieldData::alloc(gInfo, isOnGpu()));
	copy->copyData(*this);
	return copy;
}
complexScalarField complexScalarFieldData::alloc(const GridInfo& gInfo, bool onGpu) { return std::make_shared<complexScalarFieldData>(gInfo, onGpu, PrivateTag()); }

matrix complexScalarFieldData::toMatrix() const
{
  matrix mat(gInfo.nr,1,isOnGpu());
  callPref(eblas_copy)(mat.dataPref(), dataPref(), gInfo.nr);
  return mat;
}


//------------ class complexScalarFieldTildeData ---------------

complexScalarFieldTildeData::complexScalarFieldTildeData(const GridInfo& gInfo, bool onGpu, PrivateTag)
: FieldData<complex>(gInfo, "complexScalarFieldTilde", gInfo.nr, onGpu)
{
}
complexScalarFieldTilde complexScalarFieldTildeData::clone() const
{	complexScalarFieldTilde copy(complexScalarFieldTildeData::alloc(gInfo, isOnGpu()));
	copy->copyData(*this);
	return copy;
}
complexScalarFieldTilde complexScalarFieldTildeData::alloc(const GridInfo& gInfo, bool onGpu) { return std::make_shared<complexScalarFieldTildeData>(gInfo, onGpu, PrivateTag()); }


//------------ class RealKernel ---------------

RealKernel::RealKernel(const GridInfo& gInfo)
 : FieldData<double>(gInfo, "RealKernel", gInfo.nG, false) //basically real version of scalarFieldTilde
{
}

