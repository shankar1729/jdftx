/*-------------------------------------------------------------------
Copyright 2018 Ravishankar Sundararaman

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
#include <core/Random.h>

//------------- Sub-matrices ------------------------

complex matrix::operator()(int i, int j) const
{	assert(i<nr and i>=0);
	assert(j<nc and j>=0);
	if(isOnGpu())
	{	complex ret;
		#ifdef GPU_ENABLED
		cudaMemcpy(&ret, dataGpu()+index(i,j), sizeof(complex), cudaMemcpyDeviceToHost);
		#else
		assert(!"onGpu=true without GPU_ENABLED");
		#endif
		return ret;
	}
	else return data()[index(i,j)];
}

#ifdef GPU_ENABLED
void matrixSubGet_gpu(int nr, int iStart, int iStep, int iDelta, int jStart, int jStep, int jDelta, const complex* in, complex* out); //implemented in operators.cu
#endif
matrix matrix::operator()(int iStart, int iStep,  int iStop, int jStart, int jStep, int jStop) const
{	if(iStart==0 && iStep==1 && iStop==nr && jStart==0 && jStep==1 && jStop==nc)
		return *this; //faster to copy matrix for this special case
	
	assert(iStart>=0 && iStart<nr);
	assert(iStop>iStart && iStop<=nr);
	assert(iStep>0);
	assert(jStart>=0 && jStart<nc);
	assert(jStop>jStart && jStop<=nc);
	assert(jStep>0);
	
	int iDelta = ceildiv(iStop-iStart, iStep), jDelta = ceildiv(jStop-jStart, jStep);
	matrix ret(iDelta,jDelta, isGpuEnabled()); complex* retData = ret.dataPref(); const complex* thisData = this->dataPref();
	#ifdef GPU_ENABLED
	matrixSubGet_gpu(nr, iStart,iStep,iDelta, jStart,jStep,jDelta, thisData, retData);
	#else
	for(int i=0; i<iDelta; i++)
		for(int j=0; j<jDelta; j++)
			retData[ret.index(i,j)] = thisData[this->index(i*iStep+iStart, j*jStep+jStart)];
	#endif
	return ret;
}


void matrix::set(int i, int j, complex m)
{	assert(i<nr and i>=0);
	assert(j<nc and j>=0);
	if(isOnGpu())
	{
		#ifdef GPU_ENABLED
		cudaMemcpy(dataGpu()+index(i,j), &m, sizeof(complex), cudaMemcpyHostToDevice);
		#else
		assert(!"onGpu=true without GPU_ENABLED");
		#endif
	}
	else data()[index(i,j)] = m;
}

#ifdef GPU_ENABLED
void matrixSubSet_gpu(int nr, int iStart, int iStep, int iDelta, int jStart, int jStep, int jDelta, const complex* in, complex* out);
#endif
void matrix::set(int iStart, int iStep, int iStop, int jStart, int jStep, int jStop, const matrix& m)
{	assert(iStart>=0 && iStart<nr);
	assert(iStop>iStart && iStop<=nr);
	assert(iStep>0);
	assert(jStart>=0 || jStart<nc);
	assert(jStop>jStart || jStop<=nc);
	assert(jStep>0);
	int iDelta = ceildiv(iStop-iStart, iStep), jDelta = ceildiv(jStop-jStart, jStep);
	assert(iDelta==m.nr);
	assert(jDelta==m.nc);

	const complex* mData = m.dataPref(); complex* thisData = this->dataPref();
	#ifdef GPU_ENABLED
	matrixSubSet_gpu(nr, iStart,iStep,iDelta, jStart,jStep,jDelta, mData, thisData);
	#else
	for(int i=0; i<iDelta; i++)
		for(int j=0; j<jDelta; j++)
			thisData[this->index(i*iStep+iStart, j*jStep+jStart)] = mData[m.index(i,j)];
	#endif
}

//----------------------- Arithmetic ---------------------

matrix operator*(const matrixScaledTransOp &m1st, const matrixScaledTransOp &m2st)
{	assert(m1st.nCols() == m2st.nRows());
	const matrix& m1 = m1st.mat;
	const matrix& m2 = m2st.mat;
	double scaleFac = m1st.scale * m2st.scale;
	matrix ret(m1st.nRows(), m2st.nCols(), isGpuEnabled());
	callPref(eblas_zgemm)(m1st.op, m2st.op, ret.nRows(), ret.nCols(), m1st.nCols(),
		scaleFac, m1.dataPref(), m1.nRows(), m2.dataPref(), m2.nRows(),
		0.0, ret.dataPref(), ret.nRows());
	return ret;
}

matrix operator*(const diagMatrix& d, const matrix& m)
{	assert(d.nCols()==m.nRows());
	matrix ret(m); //copy input to output
	ManagedArray<double> dManaged(d); //transfer diagonal matrix to managed memory (for potential GPU access)
	//Elementwise-multiply each column by the scale factors (better contiguous access pattern than row-wise)
	for(int j=0; j<ret.nCols(); j++)
		callPref(eblas_zmuld)(ret.nRows(), dManaged.dataPref(), 1, ret.dataPref()+ret.index(0,j), 1);
	return ret;
}

matrix operator*(const matrix& m, const diagMatrix& d)
{	assert(m.nCols()==d.nRows());
	matrix ret(m); //copy input to out
	//Scale each column:
	for(int j=0; j<ret.nCols(); j++)
		callPref(eblas_zdscal)(ret.nRows(), d[j], ret.dataPref()+ret.index(0,j), 1);
	return ret;
}

matrix operator*(const matrix& m, const std::vector<complex>& d)
{	assert(m.nCols()==int(d.size()));
	matrix ret(m); //copy input to out
	//Scale each column:
	for(int j=0; j<ret.nCols(); j++)
		callPref(eblas_zscal)(ret.nRows(), d[j], ret.dataPref()+ret.index(0,j), 1);
	return ret;
}

diagMatrix operator*(const diagMatrix& d1, const diagMatrix& d2)
{	assert(d1.nCols()==d2.nRows());
	diagMatrix ret(d1);
	for(int i=0; i<ret.nRows(); i++) ret[i] *= d2[i]; //elementwise multiply
	return ret;
}

void axpy(double alpha, const diagMatrix& x, matrix& y)
{	assert(x.nRows()==y.nRows());
	assert(x.nCols()==y.nCols());
	complex* yData = y.data();
	for(int i=0; i<y.nRows(); i++) yData[y.index(i,i)] += alpha * x[i];
}

void axpy(double alpha, const diagMatrix& x, diagMatrix& y)
{	assert(x.nRows()==y.nRows());
	for(int i=0; i<y.nRows(); i++) y[i] += alpha * x[i];
}

diagMatrix clone(const diagMatrix& x) { return x; }
matrix clone(const matrix& x) { return x; }

double dot(const diagMatrix& x, const diagMatrix& y)
{	assert(x.size()==y.size());
	double ret = 0.;
	for(size_t i=0; i<x.size(); i++)
		ret += x[i]*y[i];
	return ret;
}
double dot(const matrix& x, const matrix& y) { return dotc(x,y).real(); }

void randomize(diagMatrix& x) { for(size_t i=0; i<x.size(); i++) x[i] = Random::normal(); }
void randomize(matrix& x)
{	complex* xData = x.data();
	for(size_t i=0; i<x.nData(); i++)
		xData[i] = Random::normalComplex();
}

//---------------- Misc matrix ops --------------------

complex det(const matrix& A)
{
	matrix decomposition = LU(A);
	int N = A.nRows();
	
	// Multiplies the diagonal entries to get the determinant up to a sign
	complex determinant(1., 0.);
	for(int i=0; i<N; i++)
		determinant *= decomposition(i,i);

	return determinant;

}

double det(const diagMatrix& A)
{	double determinant = 1.;
	for(int i=0; i<A.nCols(); i++)
		determinant *= A[i];
	return determinant;
}

complex trace(const matrix &A)
{	assert(A.nRows() == A.nCols());
	matrix one(eye(1));
	return callPref(eblas_zdotc)(A.nRows(), one.dataPref(),0, A.dataPref(),A.nCols()+1);
}

double trace(const diagMatrix& A)
{	double ret=0.0;
	for(double d: A) ret += d;
	return ret;
}

double nrm2(const diagMatrix& A)
{	double ret=0.0;
	for(double d: A) ret += d*d;
	return sqrt(ret);
}

diagMatrix diag(const matrix &A)
{	assert(A.nRows()==A.nCols());
	diagMatrix ret(A.nRows());
	const complex* Adata = A.data();
	for(int i=0; i<A.nRows(); i++) ret[i] = Adata[A.index(i,i)].real();
	return ret;
}

diagMatrix eye(int N)
{	diagMatrix ret(N, 1.);
	return ret;
}

matrix zeroes(int nRows, int nCols)
{	matrix ret(nRows, nCols, isGpuEnabled());
	ret.zero();
	return ret;
}


//Apply pending transpose / dagger operations:
matrixScaledTransOp::operator matrix() const
{	if(op==CblasNoTrans) return scale * mat;
	else
	{	const complex* matData = mat.data();
		matrix ret(mat.nCols(), mat.nRows()); complex* retData = ret.data();
		for(int i=0; i < mat.nCols(); i++)
			for(int j=0; j < mat.nRows(); j++)
				retData[ret.index(i,j)] = scale * conjOp(matData[mat.index(j,i)]);
		return ret;
	}
}

//Complex conjugate:
matrix conj(const scaled<matrix>& A)
{	matrix B = A.data;
	double* Bdata = (double*)B.dataPref();
	callPref(eblas_dscal)(B.nData(), A.scale, Bdata, 2); //real parts
	callPref(eblas_dscal)(B.nData(), -A.scale, Bdata+1, 2); //imag parts
	return B;
}

// Hermitian adjoint
matrixScaledTransOp dagger(const scaled<matrix> &A)
{	return matrixScaledTransOp(A.data, A.scale, CblasConjTrans);
}
matrix dagger_symmetrize(const scaled<matrix> &A)
{	return 0.5*(dagger(A) + A);
}
// Transpose
matrixScaledTransOp transpose(const scaled<matrix> &A)
{	return matrixScaledTransOp(A.data, A.scale, CblasTrans);
}
matrix transpose_symmetrize(const scaled<matrix>& A)
{	return 0.5*(transpose(A) + A);
}

double relativeHermiticityError(int N, const complex* data)
{	double errNum=0.0, errDen=1e-20; //avoid irrelevant error for zero matrix
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
		{	int index = N*i + j;
			int indexT = N*j + i;
			errNum += norm(data[index]-data[indexT].conj());
			errDen += norm(data[index]);
		}
	return sqrt(errNum / (errDen*N));
}

//------------ Block matrices ------------

tiledBlockMatrix::tiledBlockMatrix(const matrix& mBlock, int nBlocks, const std::vector<complex>* phaseArr) : mBlock(mBlock), nBlocks(nBlocks), phaseArr(phaseArr)
{	if(phaseArr) assert(nBlocks==int(phaseArr->size()));
}

matrix tiledBlockMatrix::operator*(const matrix& other) const
{	assert(mBlock.nCols()*nBlocks == other.nRows());
	matrix result(mBlock.nRows()*nBlocks, other.nCols(), isGpuEnabled());
	//Dense matrix multiply for each block:
	for(int iBlock=0; iBlock<nBlocks; iBlock++)
	{	int offs = iBlock * mBlock.nCols();
		complex phase = phaseArr ? phaseArr->at(iBlock) : 1.;
		callPref(eblas_zgemm)(CblasNoTrans, CblasNoTrans, mBlock.nRows(), other.nCols(), mBlock.nCols(),
			phase, mBlock.dataPref(), mBlock.nRows(), other.dataPref()+offs, other.nRows(),
			0.0, result.dataPref()+offs, result.nRows());
	}
	return result;
}

matrix operator*(const matrix& m, const tiledBlockMatrix& tbm)
{	assert(m.nCols() == tbm.mBlock.nRows()*tbm.nBlocks);
	matrix result(m.nRows(), tbm.mBlock.nCols()*tbm.nBlocks, isGpuEnabled());
	//Dense matrix multiply for each block:
	for(int iBlock=0; iBlock<tbm.nBlocks; iBlock++)
	{	int offsIn = iBlock * tbm.mBlock.nRows() * m.nRows();
		int offsOut = iBlock * tbm.mBlock.nCols() * m.nRows();
		complex phase = tbm.phaseArr ? tbm.phaseArr->at(iBlock) : 1.;
		callPref(eblas_zgemm)(CblasNoTrans, CblasNoTrans, m.nRows(), tbm.mBlock.nCols(), tbm.mBlock.nRows(),
			phase, m.dataPref()+offsIn, m.nRows(), tbm.mBlock.dataPref(), tbm.mBlock.nRows(),
			0.0, result.dataPref()+offsOut, result.nRows());
	}
	return result;
}
