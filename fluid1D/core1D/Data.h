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

#ifndef FLUID1D_CORE1D_DATA_H
#define FLUID1D_CORE1D_DATA_H

#include <vector>
#include <core1D/GridInfo.h>

typedef std::vector<double> SphericalKernel;

//! Base class for managed-memory objects such as scalar fields
class ManagedMemory
{
public:
	double* data(); //!< Return a pointer to the actual data
	const double* data() const; //!< Return a const pointer to the actual data
	size_t nData() const { return nElements; } //!< number of data points

	void write(const char *fname) const; //!< binary-write to a file
	void writea(const char *fname) const; //!< binary-append to a file
	void write(FILE *fp) const; //!< binary-write toa stream
	void read(const char *fname); //!< binary read from a file
	void read(FILE *fp); //!< binary read from a stream
	void zero(); //!< set all elements to zero
	
	explicit operator bool() const { return nElements; } //!< Test nullness

protected:
	void memFree(); //!< Do the job of the destructor
	void memInit(size_t nElem); //!< Do the job of the constructor

	void memMove(ManagedMemory&&); //!< Steal the other object's data (used for move constructors/assignment)
	ManagedMemory(size_t nElements=0); //!< Construct, optionally with data allocation
	~ManagedMemory();

private:
	size_t nElements;
	double* pData; //!< Actual data storage
};

//Some common elementwise / vector-like operations
void memcpy(ManagedMemory&, const ManagedMemory&); //!< copy entire object over
void scale(double alpha, ManagedMemory& y); //! scale y *= alpha
void axpy(double alpha, const ManagedMemory& x, ManagedMemory& y); //!< standard blas y += alpha*x
double nrm2(const ManagedMemory&); //!< 2-norm, pretending it is a vector
double dot(const ManagedMemory& a, const ManagedMemory& b); //!< return a^T b

//! Real-space real scalar field
class ScalarField : public ManagedMemory
{
public:
	const GridInfo *gInfo;

	void init(const GridInfo *gInfo);
	ScalarField(const GridInfo *gInfo=0);
	ScalarField(const ScalarField&); //!< copy constructor
	ScalarField(ScalarField&&); //!< move constructor
	ScalarField& operator=(const ScalarField&); //!< copy-assignment
	ScalarField& operator=(ScalarField&&); //!< move-assignment
};

//! Basis-space real scalar field
class ScalarFieldTilde : public ManagedMemory
{
public:
	const GridInfo *gInfo;

	void init(const GridInfo *gInfo);
	ScalarFieldTilde(const GridInfo *gInfo=0);
	ScalarFieldTilde(const ScalarFieldTilde&); //!< copy constructor
	ScalarFieldTilde(ScalarFieldTilde&&); //!< move constructor
	ScalarFieldTilde& operator=(const ScalarFieldTilde&); //!< copy-assignment
	ScalarFieldTilde& operator=(ScalarFieldTilde&&); //!< move-assignment
};

#endif // FLUID1D_CORE1D_DATA_H
