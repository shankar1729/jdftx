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

#ifndef JDFTX_ELECTRONIC_COLUMNBUNDLE_H
#define JDFTX_ELECTRONIC_COLUMNBUNDLE_H

#include <electronic/common.h>
#include <core/ManagedMemory.h>
#include <electronic/Basis.h>
#include <electronic/matrix.h>
#include <core/ScalarField.h>
#include <core/scaled.h>

class ColumnBundle : public ManagedMemory
{
	int ncols;
	size_t col_length; //typically equal to either basis->nbasis (or 2*basis->nbasis for spinors in noncollinear mode)

public:
	int nCols() const { return ncols; } //!< number of columns accessor
	size_t colLength() const { return col_length; } //!< column length accessor
	explicit operator bool() const { return ncols*col_length; } //!< Test nullness of ColumnBundle

	//! Index of the i'th column j'th point into the data array
	size_t index(int i, size_t j) const { return i*col_length+j; }

	bool isSpinor() const { return basis && (col_length==2*basis->nbasis); }
	int spinorLength() const { return isSpinor() ? 2 : 1; }
	
	const QuantumNumber *qnum;
	const Basis *basis;

	void init(int nc, size_t len, const Basis* b, const QuantumNumber* q, bool onGpu=false); //!< onGpu MUST remain false if no GPU_ENABLED!
	void free(); //!< Force cleanup
	ColumnBundle(int nc=0, size_t len=0, const Basis* b=NULL, const QuantumNumber* q=NULL, bool onGpu=false); //!< onGpu MUST remain false if no GPU_ENABLED!
	ColumnBundle(const ColumnBundle&); //!< copy constructor
	ColumnBundle(ColumnBundle&&); //!< move constructor

	//! Create a column bundle like this (without copying data)
	//! optionally with different number of columns (if ncOverride>=0)
	ColumnBundle similar(int ncOverride=-1) const;
	
	ColumnBundle& operator=(const ColumnBundle&); //!< copy-assignment
	ColumnBundle& operator=(ColumnBundle&&); //!< move-assignment
	
	// Get/set columns
	ColumnBundle getSub(int colStart, int colStop) const; //!< get a range of columns as a ColumnBundle 
	void setSub(int colStart, const ColumnBundle&); //!< set columns (starting at colStart) from a ColumnBundle, ignoring columns that would go beyond nCols()
	
	complexScalarFieldTilde getColumn(int i, int s) const; //!< Expand the i'th column and s'th spinor component from reduced to full G-space
	void setColumn(int i, int s, const complexScalarFieldTilde&); //!< Redeuce a full G-space vector and store it as the i'th column and s'th spinor component
	void accumColumn(int i, int s, const complexScalarFieldTilde&); //!< Redeuce a full G-space vector and accumulate onto the i'th column and s'th spinor component
	
	void randomize(int colStart, int colStop); //!< randomize a selected range of columns
};

//! Initialize an array of column bundles (with appropriate wavefunction sizes if ncols, basis, qnum and eInfo are all non-zero)
void init(std::vector<ColumnBundle>&, int nbundles, int ncols=0, const Basis* basis=0, const ElecInfo* eInfo=0);

void randomize(std::vector<ColumnBundle>&, const ElecInfo& eInfo);
void write(const std::vector<ColumnBundle>&, const char *fname, const ElecInfo& eInfo);

struct ColumnBundleReadConversion
{	bool realSpace; //!< whether to read realspace wavefunctions
	int nBandsOld; //!< nBands for the input wavefunction
	double Ecut, EcutOld; //!< Ecut for the current calculation and input wavefunction in fourier space
	vector3<int> S_old; //!< fftbox size for the input wavefunction in double space
	
	ColumnBundleReadConversion();
};
void read(std::vector<ColumnBundle>&, const char *fname, const ElecInfo& eInfo, const ColumnBundleReadConversion* conversion=0);

// Used in the CG template Minimize.h
ColumnBundle clone(const ColumnBundle&);  //! Copies the input
void randomize(ColumnBundle& x); //!< Initialize to random numbers
double dot(const ColumnBundle& x, const ColumnBundle& y); //!< inner product


//----------- Arithmetic ------------

ColumnBundle& operator+=(ColumnBundle& Y, const ColumnBundle &X);
ColumnBundle& operator-=(ColumnBundle& Y, const ColumnBundle &X);
ColumnBundle operator+(const ColumnBundle &Y1, const ColumnBundle &Y2);
ColumnBundle operator-(const ColumnBundle &Y1,const ColumnBundle &Y2);

ColumnBundle& operator*=(ColumnBundle& X, double s);
scaled<ColumnBundle> operator*(double s, const ColumnBundle &Y);
scaled<ColumnBundle> operator*(const ColumnBundle &Y, double s);
scaled<ColumnBundle> operator-(const ColumnBundle &Y);
ColumnBundle& operator*=(ColumnBundle& X, complex s);
ColumnBundle operator*(complex s, const ColumnBundle &Y);
ColumnBundle operator*(const ColumnBundle &Y, complex s);

ColumnBundle operator*(const scaled<ColumnBundle>&, const matrixScaledTransOp&);
ColumnBundle operator*(const scaled<ColumnBundle>&, const diagMatrix&);
matrix operator^(const scaled<ColumnBundle>&, const scaled<ColumnBundle>&); //!< inner product

#endif // JDFTX_ELECTRONIC_COLUMNBUNDLE_H
