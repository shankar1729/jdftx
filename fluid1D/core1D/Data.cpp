/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

This file is part of Fluid1D.

Fluid1D is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Fluid1D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Fluid1D.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#include <core1D/Data.h>
#include <core/Util.h>
#include <core/BlasExtra.h>

//-------- class ManagedMemory --------------

//No management at this point - simple container
double* ManagedMemory::data()
{	return pData;
}
const double* ManagedMemory::data() const
{	return pData;
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
{	size_t nDone = fwrite(data(), sizeof(double), nData(), fp);
	if(nDone<nData()) die("Error after processing %lu of %lu records.\n", nDone, nData());
}
void ManagedMemory::read(const char *fname)
{	FILE *fp = fopen(fname, "rb");
	if(!fp) die("Error opening %s for reading.\n", fname);
	read(fp);
	fclose(fp);
}
void ManagedMemory::read(FILE *fp)
{	size_t nDone = fread(data(), sizeof(double), nData(), fp);
	if(nDone<nData()) die("Error after processing %lu of %lu records.\n", nDone, nData());
}
void ManagedMemory::zero()
{	eblas_zero(nData(), data());
}

ManagedMemory::ManagedMemory(size_t nElements) : nElements(0)
{	memInit(nElements);
}
ManagedMemory::~ManagedMemory()
{	memFree();
}
void ManagedMemory::memFree()
{	if(!nElements) return; //nothing to free
	fftw_free(pData);
	pData = 0;
	nElements = 0;
}
void ManagedMemory::memInit(size_t nElem)
{	if(nElements==nElem) return; //already in required state
	memFree();
	nElements = nElem;
	pData = (double*)fftw_malloc(sizeof(double)*nElements);
	if(!pData) assert(!"Memory allocation failed (out of memory)\n");
}
void ManagedMemory::memMove(ManagedMemory&& other)
{	std::swap(nElements, other.nElements);
	std::swap(pData, other.pData);
	//Now other will be empty, while *this will have all its contents
}

void memcpy(ManagedMemory& a, const ManagedMemory& b)
{	assert(a.nData() == b.nData());
	if(!a.nData()) return; //no data to copy
	memcpy(a.data(), b.data(), a.nData()*sizeof(double));
}

void scale(double alpha, ManagedMemory& y)
{	eblas_dscal(y.nData(), alpha, y.data(), 1);
}

void axpy(double alpha, const ManagedMemory& x, ManagedMemory& y)
{	assert(x.nData() == y.nData());
	eblas_daxpy(x.nData(), alpha, x.data(), 1, y.data(), 1);
}

double nrm2(const ManagedMemory& a)
{	return eblas_dnrm2(a.nData(), a.data(), 1);
}

double dot(const ManagedMemory& a, const ManagedMemory& b)
{	assert(a.nData() == b.nData());
	return eblas_ddot(a.nData(), a.data(), 1, b.data(), 1);
}



//-------- class ScalarField ---------------

void ScalarField::init(const GridInfo *gInfo)
{	this->gInfo = gInfo;
	if(gInfo) memInit(gInfo->S);
}
ScalarField::ScalarField(const GridInfo *gInfo)
{	init(gInfo);
}
ScalarField::ScalarField(const ScalarField& other)
{	init(other.gInfo);
	memcpy((ManagedMemory&)*this, (const ManagedMemory&)other); //copy data
}
ScalarField::ScalarField(ScalarField&& other)
{	std::swap(gInfo, other.gInfo);
	memMove((ManagedMemory&&)other); //move data
}
ScalarField& ScalarField::operator=(const ScalarField& other)
{	init(other.gInfo);
	memcpy((ManagedMemory&)*this, (const ManagedMemory&)other); //copy data
	return *this;
}
ScalarField& ScalarField::operator=(ScalarField&& other)
{	std::swap(gInfo, other.gInfo);
	memMove((ManagedMemory&&)other); //move data
	return *this;
}

//-------- class ScalarFieldTilde ---------------

void ScalarFieldTilde::init(const GridInfo *gInfo)
{	this->gInfo = gInfo;
	if(gInfo) memInit(gInfo->S);
}
ScalarFieldTilde::ScalarFieldTilde(const GridInfo *gInfo)
{	init(gInfo);
}
ScalarFieldTilde::ScalarFieldTilde(const ScalarFieldTilde& other)
{	init(other.gInfo);
	memcpy((ManagedMemory&)*this, (const ManagedMemory&)other); //copy data
}
ScalarFieldTilde::ScalarFieldTilde(ScalarFieldTilde&& other)
{	std::swap(gInfo, other.gInfo);
	memMove((ManagedMemory&&)other); //move data
}
ScalarFieldTilde& ScalarFieldTilde::operator=(const ScalarFieldTilde& other)
{	init(other.gInfo);
	memcpy((ManagedMemory&)*this, (const ManagedMemory&)other); //copy data
	return *this;
}
ScalarFieldTilde& ScalarFieldTilde::operator=(ScalarFieldTilde&& other)
{	std::swap(gInfo, other.gInfo);
	memMove((ManagedMemory&&)other); //move data
	return *this;
}
