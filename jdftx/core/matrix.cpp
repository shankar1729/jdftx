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

#include <core/matrix.h>
#include <core/GpuUtil.h>
#include <core/GridInfo.h>

//---------------------- class diagMatrix --------------------------

bool diagMatrix::isScalar(double absTol, double relTol) const
{	double mean = 0.0; for(double d: *this) mean += d; mean /= nRows();
	double errThresh = fabs(absTol) + fabs(mean * relTol);
	for(double d: *this) if(fabs(d-mean) > errThresh) return false;
	return true;
}

diagMatrix diagMatrix::operator()(int iStart, int iStop) const
{	assert(iStart>=0 && iStart<nRows());
	assert(iStop>iStart && iStop<=nRows());
	int iDelta = iStop-iStart;
	diagMatrix ret(iDelta);
	for(int i=0; i<iDelta; i++) ret[i] = at(i+iStart);
	return ret;
}
void diagMatrix::set(int iStart, int iStop, const diagMatrix& m)
{	assert(iStart>=0 && iStart<nRows());
	assert(iStop>iStart && iStop<=nRows());
	int iDelta = iStop-iStart;
	assert(iDelta==m.nRows());
	for(int i=0; i<iDelta; i++) at(i+iStart) = m[i];
}
diagMatrix diagMatrix::operator()(int iStart, int iStep,  int iStop) const
{	assert(iStart>=0 && iStart<nRows());
	assert(iStop>iStart && iStop<=nRows());
	assert(iStep>0);
	int iDelta = ceildiv(iStop-iStart, iStep);
	diagMatrix ret(iDelta);
	for(int i=0; i<iDelta; i++) ret[i] = at(i*iStep+iStart);
	return ret;
}
void diagMatrix::set(int iStart, int iStep, int iStop, const diagMatrix& m)
{	assert(iStart>=0 && iStart<nRows());
	assert(iStop>iStart && iStop<=nRows());
	assert(iStep>0);
	int iDelta = ceildiv(iStop-iStart, iStep);
	assert(iDelta==m.nRows());
	for(int i=0; i<iDelta; i++) at(i*iStep+iStart) = m[i];
}

void diagMatrix::scan(FILE* fp)
{	for(double& d: *this) fscanf(fp, "%lg", &d);
}

void diagMatrix::print(FILE* fp, const char* fmt) const
{	for(double d: *this) fprintf(fp, fmt, d);
	fprintf(fp,"\n");
}

void diagMatrix::send(int dest, int tag) const
{	assert(mpiWorld->nProcesses()>1);
	mpiWorld->send(data(), size(), dest, tag);
}
void diagMatrix::recv(int src, int tag)
{	assert(mpiWorld->nProcesses()>1);
	mpiWorld->recv(data(), size(), src, tag);
}
void diagMatrix::bcast(int root)
{	if(mpiWorld->nProcesses()>1)
		mpiWorld->bcast(data(), size(), root);
}
void diagMatrix::allReduce(MPIUtil::ReduceOp op, bool safeMode)
{	if(mpiWorld->nProcesses()>1)
		mpiWorld->allReduce(data(), size(), op, safeMode);
}

//----------------------- class matrix ---------------------------

//Initialization
void matrix::init(int nrows, int ncols, bool onGpu)
{
	nr = nrows;
	nc = ncols;
	
	if(nr*nc>0) memInit("matrix", nr*nc, onGpu);
}
//Reshaping
void matrix::reshape(int nrows, int ncols)
{	assert(nrows>=0);
	assert(ncols>=0);
	size_t nProd = nr * nc; //current size
	//Fill in missing dimensions if any:
	if(!nrows) { assert(ncols); nrows = nProd / ncols; }
	if(!ncols) { assert(nrows); ncols = nProd / nrows; }
	//Update dimensions:
	assert(nrows * ncols == int(nProd));
	nr = nrows;
	nc = ncols;
}
// Default constructor
matrix::matrix(int nrows, int ncols, bool onGpu)
{	init(nrows,ncols, onGpu);
}
// Copy constructor
matrix::matrix(const matrix& m1)
{	init(m1.nRows(), m1.nCols(), m1.isOnGpu());
	memcpy((ManagedMemory&)*this, (const ManagedMemory&)m1);
}
// Move constructor
matrix::matrix(matrix&& m1)
{	std::swap(nr, m1.nr);
	std::swap(nc, m1.nc);
	memMove((ManagedMemory&&)m1);
}
// Construct from a real diagonal matrix:
matrix::matrix(const diagMatrix& d)
{	nr = d.size();
	nc = d.size();
	if(d.size())
	{	memInit("matrix", nr*nc); zero();
		complex* thisData = data();
		for(int i=0; i<nRows(); i++) thisData[index(i,i)] = d[i];
	}
}
// Construct from a complex diagonal matrix:
matrix::matrix(const std::vector<complex>& d)
{	nr = d.size();
	nc = d.size();
	if(d.size())
	{	memInit("matrix", nr*nc); zero();
		complex* thisData = data();
		for(int i=0; i<nRows(); i++) thisData[index(i,i)] = d[i];
	}
}
// Construct from a complex diagonal matrix:
matrix::matrix(const matrix3<>& m)
{	nr = 3;
	nc = 3;
	memInit("matrix", nr*nc);
	for(int j=0; j<3; j++)
		for(int i=0; i<3; i++)
			set(i,j, m(i,j));
}
//Copy assignment
matrix& matrix::operator=(const matrix &m1)
{	init(m1.nRows(), m1.nCols(), m1.isOnGpu());
	memcpy((ManagedMemory&)*this, (const ManagedMemory&)m1);
	return *this;
}
// Move assignment
matrix& matrix::operator=(matrix&& m1)
{	std::swap(nr, m1.nr);
	std::swap(nc, m1.nc);
	memMove((ManagedMemory&&)m1);
	return *this;
}

//----------- Formatted (ascii) read/write ----------

void matrix::scan(FILE* fp, const char* fmt)
{	complex* thisData = this->data();
	for(int i=0; i<nRows(); i++)
	{	for(int j=0; j<nCols(); j++)
		{	complex& c = thisData[index(i,j)];
			fscanf(fp, fmt, &c.real(), &c.imag());
		}
	}
}
void matrix::scan_real(FILE* fp)
{	complex* thisData = this->data();
	for(int i=0; i<nRows(); i++)
	{	for(int j=0; j<nCols(); j++)
		{	complex& c = thisData[index(i,j)];
			fscanf(fp, "%lg", &c.real());
			c.imag() = 0;
		}
	}
}
void matrix::print(FILE* fp, const char* fmt) const
{	const complex* thisData = this->data();
	for(int i=0; i<nRows(); i++)
	{	for(int j=0; j<nCols(); j++)
		{	const complex& c = thisData[index(i,j)];
			fprintf(fp, fmt, c.real(), c.imag());
		}
		fprintf(fp,"\n");
	}
}
void matrix::print_real(FILE* fp, const char* fmt) const
{	const complex* thisData = this->data();
	for(int i=0; i<nRows(); i++)
	{	for(int j=0; j<nCols(); j++)
			fprintf(fp, fmt, thisData[index(i,j)].real());
		fprintf(fp,"\n");
	}
}


//---------octave-like slicing operators on scalar field matrices--------------

//get a particular element (from the gpu if needed)
complex matrix::getElement(vector3<int> index, GridInfo& gInfo)
{	matrix result(1,1); //use a 1x1 matrix to store the result locally (CPU/GPU)
	callPref(eblas_copy)(result.dataPref(), dataPref()+gInfo.fullRindex(index), 1);
	return trace(result).real(); //this takes care of GPU->CPU copy
}
