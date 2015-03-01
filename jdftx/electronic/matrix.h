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

#ifndef JDFTX_ELECTRONIC_MATRIX_H
#define JDFTX_ELECTRONIC_MATRIX_H

#include <electronic/common.h>
#include <core/ManagedMemory.h>
#include <core/matrix3.h>
#include <core/scaled.h>
#include <gsl/gsl_cblas.h>

//! Real diagonal matrix
class diagMatrix : public std::vector<double>
{	
public:
	diagMatrix(int N=0, double d=0.) : std::vector<double>(N,d) {}
	int nRows() const { return size(); }
	int nCols() const { return size(); }
	bool isScalar(double absTol=1e-14, double relTol=1e-14) const; //!< Check if all entries are identical within threshold
	
	diagMatrix& operator*=(double s) { for(double& d: *this) d*=s; return *this; }
	
	//Splicing operations:
	diagMatrix operator()(int iStart, int iStop) const; //! get submatrix of elements (iStart \<= i \< iStop)
	diagMatrix operator()(int iStart, int iStep, int iStop) const; //! get submatrix of elements (iStart \<= i \< iStop) with arbitrary increments
	void set(int iStart, int iStop, const diagMatrix& m); //! set submatrix to m
	void set(int iStart, int iStep, int iStop, const diagMatrix& m); //! set submatrix to m at arbitrary increments

	void scan(FILE* fp); //!< read (ascii) from stream
	void print(FILE* fp, const char* fmt="%lg\t") const; //!< print (ascii) to stream

	//Inter-process communication:
	void send(int dest, int tag=0) const; //send to another process
	void recv(int src, int tag=0); //receive from another process
	void bcast(int root=0); //synchronize across processes (using value on specified root process)
	void allReduce(MPIUtil::ReduceOp op, bool safeMode=false); //apply all-to-all reduction (see MPIUtil::allReduce)
};

//! General complex matrix
class matrix : public ManagedMemory
{
	int nr; //!< number of rows
	int nc; //!< number of columns

public:
	int nRows() const { return nr; }
	int nCols() const { return nc; }
	int index(int i, int j) const { return nr*j+i; } //!< Index into data()
	explicit operator bool() const { return nr*nc; } //!< Test null-ness of matrix

	void init(int nrows,int ncols, bool onGpu=false); //!< called by the constructors
	void reshape(int nrows, int ncols); //!< change the dimensions without altering the data (zero dimensions are filled in to match size)
	matrix(int nrows=0, int ncols=0, bool onGpu=false);
	matrix(const matrix& m1); //!< copy constructor
	matrix(matrix&& m1); //!< move constructor
	matrix(const diagMatrix&); //!< convert from a real diagonal matrix
	matrix(const std::vector<complex>&); //!< convert from a complex diagonal matrix
	explicit matrix(const matrix3<>&); //!< convert from a 3x3 matrix
	
	matrix& operator=(const matrix& m1); //!< copy-assignment
	matrix& operator=(matrix&& m1); //!< move-assignment

	//! get a specific element of the matrix
	complex operator()(int i, int j) const;
	
	//! get submatrix of elements (iStart \<= i \< iStop, jStart \<= j \< jStop)
	matrix operator()(int iStart, int iStop, int jStart, int jStop) const { return (*this)(iStart,1,iStop, jStart,1,jStop); }
	
	//! get submatrix of elements (iStart \<= i \< iStop, jStart \<= j \< jStop) with arbitrary increments
	matrix operator()(int iStart, int iStep, int iStop, int jStart, int jStep, int jStop) const;

	//! set element to m
	void set(int i, int j, complex m);
	
	//! set submatrix to m
	void set(int iStart, int iStop, int jStart, int jStop, const matrix& m) { set(iStart,1,iStop, jStart,1,jStop, m); }
	
	//! set submatrix to m at arbitrary increments
	void set(int iStart, int iStep, int iStop, int jStart, int jStep, int jStop, const matrix& m);
	
	void scan(FILE* fp, const char* fmt="%lg%+lgi"); //!< read (ascii) from stream
	void scan_real(FILE* fp); //!< read (ascii) real parts from stream, setting imaginary parts to 0
	void print(FILE* fp, const char* fmt="%lg%+lgi\t") const; //!< print (ascii) to stream
	void print_real(FILE* fp, const char* fmt="%lg\t") const; //!< print (ascii) real parts to stream
	
	void diagonalize(matrix& evecs, diagMatrix& eigs) const; //!< diagonalize a hermitian matrix
	void diagonalize(matrix& levecs, std::vector<complex>& eigs, matrix& revecs) const; //!< diagonalize an arbitrary matrix
	void svd(matrix& U, diagMatrix& S, matrix& Vdag) const; //!< singular value decomposition (for dimensions of this: MxN, on output U: MxM, S: min(M,N), Vdag: NxN)
};

//! Matrix with a pending scale and transpose operation
struct matrixScaledTransOp
{	const matrix& mat; //!< matrix
	double scale; //!< pending scale factor
	CBLAS_TRANSPOSE op; //!< pending operation (none, transpose or dagger)
	
	int nRows() const { return op==CblasNoTrans ? mat.nRows() : mat.nCols(); }
	int nCols() const { return op==CblasNoTrans ? mat.nCols() : mat.nRows(); }
	complex conjOp(complex a) const { return (op==CblasConjTrans) ? a.conj() : a; } //!< return conjugate if op requires it
	int index(int i, int j) const { return op==CblasNoTrans ? mat.index(i,j) : mat.index(j,i); } //!< Index into data() with possible transpose
	
	//!Create from a matrix with an optional scale and op:
	matrixScaledTransOp(const matrix& mat, double scale=1.0, CBLAS_TRANSPOSE op=CblasNoTrans)
	: mat(mat), scale(scale), op(op) {}
	
	//! Create from a scaled matrix, with an optional op
	matrixScaledTransOp(const scaled<matrix>& smat, CBLAS_TRANSPOSE op=CblasNoTrans)
	: mat(smat.data), scale(smat.scale), op(op) {}
	
	operator matrix() const; //!< convert to matrix
	
	//Scaling:
	matrixScaledTransOp& operator*=(double s) { scale *= s; return *this; }
	matrixScaledTransOp operator*(double s) const { return matrixScaledTransOp(mat,scale*s,op); }
	friend matrixScaledTransOp operator*(double s, const matrixScaledTransOp& A) { return A * s; }
};
matrix conj(const scaled<matrix>& A); //!< return element-wise complex conjugate of A
matrixScaledTransOp dagger(const scaled<matrix>& A); //!< return hermitian adjoint of A
matrixScaledTransOp transpose(const scaled<matrix>& A); //!< return transpose of A
matrix dagger_symmetrize(const scaled<matrix>& A); //! return adjoint symmetric part: (A + Adag)/2
matrix transpose_symmetrize(const scaled<matrix>& A); //! return transpose symmetric part: (A + AT)/2


//------------ Arithmetic ------------------

inline matrix& operator*=(matrix& m, double s) { scale(s, m); return m; }
inline scaled<matrix> operator*(double s, const matrix &m) { return scaled<matrix>(m, s); }
inline scaled<matrix> operator*(const matrix &m, double s) { return scaled<matrix>(m, s); }
inline scaled<matrix> operator-(const matrix &m) { return scaled<matrix>(m, -1); }
inline matrix& operator*=(matrix& m, complex s) { scale(s, m); return m; }
inline matrix operator*(complex s, const matrix &m) { matrix sm(m); sm *= s; return sm; }
inline matrix operator*(const matrix &m, complex s) { matrix sm(m); sm *= s; return sm; }
inline diagMatrix operator*(double s, const diagMatrix &m) { diagMatrix ret(m); return ret*=s; }
inline diagMatrix operator*(const diagMatrix &m, double s) { diagMatrix ret(m); return ret*=s; }
inline diagMatrix operator-(const diagMatrix &m) { return m * (-1); }

matrix operator*(const matrixScaledTransOp&, const matrixScaledTransOp&);
matrix operator*(const matrix&, const diagMatrix&);
matrix operator*(const diagMatrix&, const matrix&);
diagMatrix operator*(const diagMatrix&, const diagMatrix&);

inline matrix& operator+=(matrix &m1, const matrix &m2) { if(m1) axpy(1.0, m2, m1); else m1 = m2; return m1; }
inline matrix& operator-=(matrix &m1, const matrix &m2) { if(m1) axpy(-1.0, m2, m1); else m1 = -m2; return m1; }
inline matrix operator+(const matrix &m1, const matrix &m2) { matrix tm(m1); tm += m2; return tm; }
inline matrix operator-(const matrix &m1, const matrix &m2) { matrix tm(m1); tm -= m2; return tm; }
void axpy(double alpha, const diagMatrix& x, matrix& y);
inline matrix& operator+=(matrix& m, const diagMatrix& d) { if(m) axpy(1., d, m); else m = d; return m; }
inline matrix& operator-=(matrix& m, const diagMatrix& d) { if(m) axpy(-1., d, m); else m=-d; return m; }
inline matrix operator+(const matrix& m, const diagMatrix& d) { matrix ret(m); ret+=d; return ret; }
inline matrix operator+(const diagMatrix& d, const matrix& m) { matrix ret(m); ret+=d; return ret; }
inline matrix operator-(const matrix& m, const diagMatrix& d) { matrix ret(m); ret-=d; return ret; }
inline matrix operator-(const diagMatrix& d, const matrix& m) { matrix ret(-m); ret+=d; return ret; }
void axpy(double alpha, const diagMatrix& x, diagMatrix& y);
inline diagMatrix& operator+=(diagMatrix& d1, const diagMatrix& d2) { axpy(1., d2, d1); return d1; }
inline diagMatrix& operator-=(diagMatrix& d1, const diagMatrix& d2) { axpy(-1., d2, d1); return d1; }
inline diagMatrix operator+(const diagMatrix& d1, const diagMatrix& d2) { diagMatrix ret(d1); return ret+=d2; }
inline diagMatrix operator-(const diagMatrix& d1, const diagMatrix& d2) { diagMatrix ret(d1); return ret-=d2; }

//Functions required for using the Minimizable template on matrices:
diagMatrix clone(const diagMatrix& x);
matrix clone(const matrix& x);
double dot(const diagMatrix& x, const diagMatrix& y);
double dot(const matrix& x, const matrix& y);
void randomize(diagMatrix& x);
void randomize(matrix& x);


//------- Nonlinear matrix functions and their gradients ---------

//! Compute inverse of an arbitrary matrix A (via LU decomposition)
matrix inv(const matrix& A);
diagMatrix inv(const diagMatrix& A);

//! Compute the LU decomposition of the matrix
matrix LU(const matrix& A);

//! Compute the determinant of an arbitrary matrix A (via LU decomposition)
//! If skipZeros is true, skips diagonal entries in the diagonal of the LU decomposition that are below a tolerance
complex det(const matrix& A);

//! Compute the determinant of an diagonal matrix A
double det(const diagMatrix& A);

//! Compute matrix A^exponent, and optionally the eigensystem of A (if non-null)
matrix pow(const matrix& A, double exponent, matrix* Aevecs=0, diagMatrix* Aeigs=0);

//! Compute matrix A^-0.5 and optionally the eigensystem of A (if non-null)
matrix invsqrt(const matrix& A, matrix* Aevecs=0, diagMatrix* Aeigs=0);

//! Compute cis(A) = exp(iota A) and optionally the eigensystem of A (if non-null)
matrix cis(const matrix& A, matrix* Aevecs=0, diagMatrix* Aeigs=0);

//! Inverse function of cis: find matrix B such that A = exp(iota B). A must be unitary.
//! Optionally retrieve the eigensystem of the result (note the real eigenvalues are of the input, not the result)
matrix cisInv(const matrix& A, matrix* Bevecs=0, diagMatrix* Beigs=0);

//! Return gradient w.r.t A given gradient w.r.t sqrt(A) and A's eigensystem
matrix sqrt_grad(const matrix& grad_sqrtA, const matrix& Aevecs, const diagMatrix& Aeigs);

//! Return gradient w.r.t A given gradient w.r.t cis(A) and A's eigensystem
matrix cis_grad(const matrix& grad_cisA, const matrix& Aevecs, const diagMatrix& Aeigs);


//------------ Misc matrix functions --------------
complex trace(const matrix &m);
double trace(const diagMatrix& m);
double nrm2(const diagMatrix& m);
diagMatrix diag(const matrix &m); //!< obtain the real diagonal part of a hermitian matrix
diagMatrix eye(int N); //!< identity
matrix zeroes(int nRows, int nCols); //!< a dense-matrix of zeroes

//! A block matrix formed by repeating (tiling) a dense matrix along the diagonal
class tiledBlockMatrix
{
	const matrix& mBlock; //!< dense matrix for each block
	int nBlocks; //!< number of blocks
	const std::vector<complex>* phaseArr; //!< optional phase for each block
public:
	tiledBlockMatrix(const matrix& mBlock, int nBlocks, const std::vector<complex>* phaseArr=0);
	
	matrix operator*(const matrix&) const; //!< multiply block matrix by dense matrix
};

#endif  // JDFTX_ELECTRONIC_MATRIX_H
