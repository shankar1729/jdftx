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

#if defined(GPU_ENABLED) and defined(CUSOLVER_ENABLED)
	#define USE_CUSOLVER
	#define NcutCuSolver 32  //minimum matrix dimension for which to use CuSolver (CPU LAPACK faster for small matrices)
	#include <cusolverDn.h>
#endif

//Lapack forward declarations
extern "C"
{	void zheevr_(char* JOBZ, char* RANGE, char* UPLO, int * N, complex* A, int * LDA,
		double* VL, double* VU, int* IL, int* IU, double* ABSTOL, int* M,
		double* W, complex* Z, int* LDZ, int* ISUPPZ, complex* WORK, int* LWORK,
		double* RWORK, int* LRWORK, int* IWORK, int* LIWORK, int* INFO);
	void zgeev_(char* JOBVL, char* JOBVR, int* N, complex* A, int* LDA,
		complex* W, complex* VL, int* LDVL, complex* VR, int* LDVR,
		complex* WORK, int* LWORK, double* RWORK, int* INFO);
	void zgesdd_(char* JOBZ, int* M, int* N, complex* A, int* LDA,
		double* S, complex* U, int* LDU, complex* VT, int* LDVT,
		complex* WORK, int* LWORK, double* RWORK, int* IWORK, int* INFO);
	void zgesvd_(char* JOBU, char* JOBVT, int* M, int* N, complex* A, int* LDA,
		double* S, complex* U, int* LDU, complex* VT, int* LDVT,
		complex* WORK, int* LWORK, double* RWORK, int* INFO);
	void zgetrf_(int* M, int* N, complex* A, int* LDA, int* IPIV, int* INFO);
	void zgetri_(int* N, complex* A, int* LDA, int* IPIV, complex* WORK, int* LWORK, int* INFO);
	void zposv_(char* UPLO, int* N, int* NRHS, complex* A, int* LDA, complex* B, int* LDB, int* INFO);
}

//------------------------- Eigensystem -----------------------------------

#ifdef GPU_ENABLED
double relativeHermiticityError_gpu(int N, const complex* data); //implemented in matrixOperators.cu
#endif
double relativeHermiticityError(int N, const complex* data); //implemented in matrixOperators.cpp

void matrix::diagonalize(matrix& evecs, diagMatrix& eigs) const
{	static StopWatch watch("matrix::diagonalize");
	watch.start();
	
	assert(nCols()==nRows());
	int N = nRows();
	assert(N > 0);
	
	//Check hermiticity
	const double hermErr = callPref(relativeHermiticityError)(N, dataPref());
	if(hermErr > 1e-10)
	{	logPrintf("Relative hermiticity error of %le (>1e-10) encountered in diagonalize\n", hermErr);
		stackTraceExit(1);
	}
	
#ifdef USE_CUSOLVER
	if(N >= NcutCuSolver)
	{	cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR;
		cublasFillMode_t uplo = CUBLAS_FILL_MODE_UPPER;
		ManagedArray<double> eigsManaged; eigsManaged.init(N, true);
		evecs = *this;
		//Determine buffer size:
		int lwork = 0;
		syevjInfo_t params; cusolverDnCreateSyevjInfo(&params);
		cusolverDnZheevj_bufferSize(cusolverHandle, jobz, uplo, N, (const double2*)evecs.dataPref(), N, eigsManaged.dataPref(), &lwork, params);
		//Main call:
		ManagedArray<double2> work; work.init(lwork, true);
		ManagedArray<int> infoArr; infoArr.init(1, true);
		cusolverDnZheevj(cusolverHandle, jobz, uplo, N, (double2*)evecs.dataPref(), N, eigsManaged.dataPref(), work.dataPref(), lwork, infoArr.dataPref(), params);
		cusolverDnDestroySyevjInfo(params);
		gpuErrorCheck();
		int info = infoArr.data()[0];
		if(!info) //Success
		{	eigs.resize(N);
			eblas_copy(eigs.data(), eigsManaged.data(), N); //eigenvectors generated in place in evecs above
			watch.stop();
			return;
		}
		if(info<0) { logPrintf("Argument# %d to cusolverDn eigenvalue routine Zheevj is invalid.\n", -info); stackTraceExit(1); }
		if(info>0) logPrintf("WARNING: %d elements failed to converge in cusolverDn eigenvalue routine Zheevj; falling back to CPU LAPACK.\n", info);
	}
#endif
	char jobz = 'V'; //compute eigenvectors and eigenvalues
	char range = 'A'; //compute all eigenvalues
	char uplo = 'U'; //use upper-triangular part
	matrix A = *this; //copy input matrix (zheevr destroys input matrix)
	double eigMin = 0., eigMax = 0.; //eigenvalue range (not used for range-type 'A')
	int indexMin = 0, indexMax = 0; //eignevalue index range (not used for range-type 'A')
	double absTol = 0.; int nEigsFound;
	eigs.resize(N);
	evecs.init(N, N);
	std::vector<int> iSuppz(2*N);
	int lwork = (64+1)*N; std::vector<complex> work(lwork); //Magic number 64 obtained by running ILAENV as suggested in doc of zheevr (and taking the max over all N)
	int lrwork = 24*N; std::vector<double> rwork(lrwork); //from doc of zheevr
	int liwork = 10*N; std::vector<int> iwork(liwork); //from doc of zheevr
	int info=0;
	zheevr_(&jobz, &range, &uplo, &N, A.data(), &N,
		&eigMin, &eigMax, &indexMin, &indexMax, &absTol, &nEigsFound,
		eigs.data(), evecs.data(), &N, iSuppz.data(), work.data(), &lwork,
		rwork.data(), &lrwork, iwork.data(), &liwork, &info);
	if(info<0) { logPrintf("Argument# %d to LAPACK eigenvalue routine ZHEEVR is invalid.\n", -info); stackTraceExit(1); }
	if(info>0) { logPrintf("Error code %d in LAPACK eigenvalue routine ZHEEVR.\n", info); stackTraceExit(1); }
	watch.stop();
}

void matrix::diagonalize(matrix& levecs, std::vector<complex>& eigs, matrix& revecs) const
{	static StopWatch watch("matrix::diagonalizeNH");
	watch.start();
	int N = nRows();
	assert(N > 0);
	assert(nCols()==N);
#ifdef USE_CUSOLVER
	//No general-matrix eigenvalue solver in CuSolver yet (as of 9.1)
#endif
	//Prepare inputs and outputs:
	matrix A = *this; //destructible copy
	eigs.resize(N);
	levecs.init(N, N);
	revecs.init(N, N);
	//Prepare temporaries:
	char jobz = 'V'; //compute eigenvectors and eigenvalues
	int lwork = (64+1)*N; std::vector<complex> work(lwork); //Magic number 64 obtained by running ILAENV as suggested in doc of zheevr (and taking the max over all N)
	std::vector<double> rwork(2*N);
	//Call LAPACK and check errors:
	int info=0;
	zgeev_(&jobz, &jobz, &N, A.data(), &N, eigs.data(), levecs.data(), &N, revecs.data(), &N, work.data(), &lwork, rwork.data(), &info);
	if(info<0) { logPrintf("Argument# %d to LAPACK eigenvalue routine ZGEEV is invalid.\n", -info); stackTraceExit(1); }
	if(info>0) { logPrintf("Error code %d in LAPACK eigenvalue routine ZGEEV.\n", info); stackTraceExit(1); }
	watch.stop();
}

void matrix::svd(matrix& U, diagMatrix& S, matrix& Vdag) const
{	static StopWatch watch("matrix::svd");
	watch.start();
	//Initialize input and outputs:
	int M = nRows();
	int N = nCols();
	U.init(M,M, isGpuEnabled());
	Vdag.init(N,N);
	S.resize(std::min(M,N));
#ifdef USE_CUSOLVER
	if(M>=N && M>NcutCuSolver)
	{	//Determine buffer size and allocate buffers:
		matrix A = *this; //destructible copy
		ManagedArray<double> Smanaged; Smanaged.init(S.size(), true);
		matrix V(N,N, true);
		cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR;
		int lwork = 0;
		int econ = 0; //return full matrices
		gesvdjInfo_t params; cusolverDnCreateGesvdjInfo(&params);
		cusolverDnZgesvdj_bufferSize(cusolverHandle, jobz, econ, M, N, (double2*)A.dataPref(), M,
			Smanaged.dataPref(), (double2*)U.dataPref(), M, (double2*)V.dataPref(), N,
			&lwork, params);
		ManagedArray<double2> work; work.init(lwork, true);
		ManagedArray<int> infoArr; infoArr.init(1, true);
		//Main call:
		cusolverDnZgesvdj(cusolverHandle, jobz, econ, M, N, (double2*)A.dataPref(), M,
			Smanaged.dataPref(), (double2*)U.dataPref(), M, (double2*)V.dataPref(), N,
			work.dataPref(), lwork, infoArr.dataPref(), params);
		cusolverDnDestroyGesvdjInfo(params);
		gpuErrorCheck();
		int info = infoArr.data()[0];
		if(!info) //Success
		{	eblas_copy(S.data(), Smanaged.data(), S.size());
			Vdag = dagger(V);
			watch.stop();
			return;
		}
		if(info<0) { logPrintf("Argument# %d to CuSolver SVD routine Zgesvd is invalid.\n", -info); stackTraceExit(1); }
		if(info>0) logPrintf("WARNING: %d elements did not converge in CuSolver SVD routine Zgesvd; falling back to CPU LAPACK.\n", info);
	}
#endif
	//Initialize temporaries:
	matrix A = *this; //destructible copy
	char jobz = 'A'; //full SVD (return complete unitary matrices)
	int lwork = 2*(M*N + M + N);
	std::vector<complex> work(lwork);
	std::vector<double> rwork(S.nRows() * std::max(5*S.nRows()+7, 2*(M+N)+1));
	std::vector<int> iwork(8*S.nRows());
	//Call LAPACK and check errors:
	int info=0;
	zgesdd_(&jobz, &M, &N, A.data(), &M, S.data(), U.data(), &M, Vdag.data(), &N,
		work.data(), &lwork, rwork.data(), iwork.data(), &info);
	if(info>0) //convergence failure; try the slower stabler version
	{	int info=0;
		matrix A = *this; //destructible copy
		zgesvd_(&jobz, &jobz, &M, &N, A.data(), &M, S.data(), U.data(), &M, Vdag.data(), &N,
			work.data(), &lwork, rwork.data(), &info);
		if(info<0) { logPrintf("Argument# %d to LAPACK SVD routine ZGESVD is invalid.\n", -info); stackTraceExit(1); }
		if(info>0) { logPrintf("Error code %d in LAPACK SVD routine ZGESVD.\n", info); stackTraceExit(1); }
	}
	if(info<0) { logPrintf("Argument# %d to LAPACK SVD routine ZGESDD is invalid.\n", -info); stackTraceExit(1); }
	watch.stop();
}


//------------- Matrix nonlinear functions ---------------

//Common implementation for the matrix nonlinear functions:
#define MATRIX_FUNC(code) \
	assert(A.nRows()==A.nCols()); \
	matrix evecs; diagMatrix eigs(A.nRows()); \
	A.diagonalize(evecs, eigs); \
	std::vector<complex> eigOut(A.nRows()); \
	\
	for(int i=0; i<A.nRows(); i++) \
	{ \
		code \
	} \
	\
	if(Aevecs) *Aevecs = evecs; \
	if(Aeigs) *Aeigs = eigs; \
	return evecs * matrix(eigOut) * dagger(evecs);

// Compute matrix A^exponent, and optionally the eigensystem of A (if non-null)
matrix pow(const matrix& A, double exponent, matrix* Aevecs, diagMatrix* Aeigs, bool* isSingular)
{	if(isSingular) *isSingular = false;
	MATRIX_FUNC
	(	if(exponent<0. && eigs[i]<=0.0)
		{	if(isSingular)
				*isSingular = true; //project out current eigenvalue; calling function will handle if needed
			else //Unhandled: print stack-trace
			{	logPrintf("Eigenvalue# %d is non-positive (%le) in pow (exponent %lg)\n", i, eigs[i], exponent);
				stackTraceExit(1);
			}
		}
		else if(exponent>=0. && eigs[i]<0.0)
		{	logPrintf("WARNING: Eigenvalue# %d is negative (%le) in pow (exponent %lg); zeroing it out.\n", i, eigs[i], exponent);
			eigs[i] = 0.;
		}
		else eigOut[i] = pow(eigs[i], exponent);
	)
}

// Compute matrix A^-0.5 and optionally the eigensystem of A (if non-null)
matrix invsqrt(const matrix& A, matrix* Aevecs, diagMatrix* Aeigs, bool* isSingular)
{	return pow(A, -0.5, Aevecs, Aeigs, isSingular);
}

// Compute cis(A) = exp(iota A) and optionally the eigensystem of A (if non-null)
matrix cis(const matrix& A, matrix* Aevecs, diagMatrix* Aeigs)
{	MATRIX_FUNC
	(	eigOut[i] = cis(eigs[i]);
	)
}

#undef MATRIX_FUNC

//Inverse of cis: get the Hermitian arg() of a unitary matrix
matrix cisInv(const matrix& A, matrix* Bevecs, diagMatrix* Beigs)
{	//Make sure A is unitary:
	assert(A.nRows()==A.nCols());
	assert(nrm2(A*dagger(A) - eye(A.nRows())) < 1e-10*sqrt(A.nData()));
	//Diagonalize:
	matrix Alevecs, Arevecs; std::vector<complex> Aeigs;
	A.diagonalize(Alevecs, Aeigs, Arevecs);
	assert(nrm2(Alevecs-Arevecs) < 1e-10*sqrt(A.nData()));
	//Compute diagonal of result:
	diagMatrix B(A.nRows());
	for(int i=0; i<A.nRows(); i++)
		B[i] = Aeigs[i].arg();
	//Return results:
	if(Bevecs) *Bevecs = Arevecs;
	if(Beigs) *Beigs = B;
	matrix ret = Arevecs * B * dagger(Arevecs);
	return ret;
}



//Common implementation of the matrix nonlinear function gradients
#define MATRIX_FUNC_GRAD(code) \
	assert(gradIn.nRows()==gradIn.nCols()); \
	assert(Aevecs.nRows()==Aevecs.nCols()); \
	assert(Aevecs.nRows()==gradIn.nCols()); \
	matrix AevecsDag = dagger(Aevecs); \
	\
	matrix gradOut = AevecsDag * gradIn * Aevecs; \
	complex* gradOutData = gradOut.data(); \
	for(int i=0; i<gradOut.nRows(); i++) \
		for(int j=0; j<gradOut.nCols(); j++) \
		{	complex& elem = gradOutData[gradOut.index(i,j)]; \
			code \
		} \
	return Aevecs * gradOut * AevecsDag;

// Return gradient w.r.t A given gradient w.r.t sqrt(A) and A's eigensystem
matrix sqrt_grad(const matrix& gradIn, const matrix& Aevecs, const diagMatrix& Aeigs)
{	MATRIX_FUNC_GRAD
	(	elem /= (sqrt(Aeigs[i])+sqrt(Aeigs[j]));
	)
}
// Return gradient w.r.t A given gradient w.r.t cis(A) and A's eigensystem
matrix cis_grad(const matrix& gradIn, const matrix& Aevecs, const diagMatrix& Aeigs)
{	MATRIX_FUNC_GRAD
	(	double x = Aeigs[j] - Aeigs[i];
		elem *= fabs(x)<1.0e-13 ? complex(-0.5*x,1) : (cis(x)-1.0)/x;
	)
}
#undef MATRIX_FUNC_GRAD

//--------- LU, linear solve and inverse ----------

//Return LU decomposition if calcInv = false and inverse if calcInv = true
matrix invOrLU(const matrix& A, bool calcInv)
{	int N = A.nRows();
	assert(N > 0);
	assert(N == A.nCols());
	matrix LU(A); //destructible copy
#ifdef USE_CUSOLVER
	if(N > NcutCuSolver)
	{	//Get buffer size and allocate for LU decomposition:
		int lwork = 0;
		cusolverDnZgetrf_bufferSize(cusolverHandle, N, N, (double2*)LU.dataPref(), N, &lwork);
		ManagedArray<double2> work; work.init(lwork, true);
		ManagedArray<int> iPivot; iPivot.init(N, true); //pivot info
		ManagedArray<int> infoArr; infoArr.init(1, true);
		//Main call:
		cusolverDnZgetrf(cusolverHandle, N, N, (double2*)LU.dataPref(), N, work.dataPref(), iPivot.dataPref(), infoArr.dataPref());
		gpuErrorCheck();
		int info = infoArr.data()[0];
		if(info<0) { logPrintf("Argument# %d to CuSolver LU decomposition routine Zgetrf is invalid.\n", -info); stackTraceExit(1); }
		if(!calcInv) return LU; //rest only needed to calc inv() from LU
		if(info>0) { logPrintf("CuSolver LU decomposition routine Zgetrf found input matrix to be singular at the %d'th step.\n", info); stackTraceExit(1); }
		//Calculate inverse:
		matrix result(eye(N)); //will contain inv(A) on output
		cusolverDnZgetrs(cusolverHandle, CUBLAS_OP_N, N, N, (double2*)LU.dataPref(), N,
           iPivot.dataPref(), (double2*)result.dataPref(), N, infoArr.dataPref());
		info = infoArr.data()[0];
		if(info<0) { logPrintf("Argument# %d to CuSolver linear solve routine Zgetrs is invalid.\n", -info); stackTraceExit(1); }
		return result;
	}
#endif
	std::vector<int> iPivot(N); //pivot info
	int info; //error code in return
	//LU decomposition (in place):
	zgetrf_(&N, &N, LU.data(), &N, iPivot.data(), &info);
	if(info<0) { logPrintf("Argument# %d to LAPACK LU decomposition routine ZGETRF is invalid.\n", -info); stackTraceExit(1); }
	if(!calcInv) return LU; //rest only needed to calc inv() from LU
	if(info>0) { logPrintf("LAPACK LU decomposition routine ZGETRF found input matrix to be singular at the %d'th step.\n", info); stackTraceExit(1); }
	//Compute inverse in place:
	int lWork = (64+1)*N;
	std::vector<complex> work(lWork);
	zgetri_(&N, LU.data(), &N, iPivot.data(), work.data(), &lWork, &info);
	if(info<0) { logPrintf("Argument# %d to LAPACK matrix inversion routine ZGETRI is invalid.\n", -info); stackTraceExit(1); }
	if(info>0) { logPrintf("LAPACK matrix inversion routine ZGETRI found input matrix to be singular at the %d'th step.\n", info); stackTraceExit(1); }
	return LU;
}

matrix LU(const matrix& A)
{	static StopWatch watch("LU(matrix)");
	watch.start();
	matrix result = invOrLU(A, false);
	watch.stop();
	return result;
}

matrix inv(const matrix& A)
{	static StopWatch watch("inv(matrix)");
	watch.start();
	matrix result = invOrLU(A, true);
	watch.stop();
	return result;
}

diagMatrix inv(const diagMatrix& A)
{	diagMatrix invA = A;
	for(double& x: invA) x = 1./x;
	return invA;
}

matrix invApply(const matrix& A, const matrix& b)
{	static StopWatch watch("invApply(matrix)");
	watch.start();
	
	//Check dimensions:
	assert(A.nCols()==A.nRows());
	int N = A.nRows();
	assert(N > 0);
	assert(N == b.nRows());
	int Nrhs = b.nCols();
	assert(Nrhs > 0);
	
	//Check hermiticity
	const double hermErr = callPref(relativeHermiticityError)(N, A.dataPref());
	if(hermErr > 1e-10)
	{	logPrintf("Relative hermiticity error of %le (>1e-10) encountered in invApply\n", hermErr);
		stackTraceExit(1);
	}

#ifdef USE_CUSOLVER
	if(N > NcutCuSolver)
	{	//Get buffer size and allocate for Cholesky factorization:
		matrix Acopy = A; //destructible copy; routine will factorize matrix in place
		cublasFillMode_t uplo = CUBLAS_FILL_MODE_UPPER;
		int lwork = 0;
		cusolverDnZpotrf_bufferSize(cusolverHandle, uplo, N, (double2*)Acopy.dataPref(), N, &lwork);
		ManagedArray<double2> work; work.init(lwork, true);
		ManagedArray<int> infoArr; infoArr.init(1, true);
		//Main call:
		cusolverDnZpotrf(cusolverHandle, uplo, N, (double2*)Acopy.dataPref(), N, work.dataPref(), lwork, infoArr.dataPref());
		int info = infoArr.data()[0];
		if(info<0) { logPrintf("Argument# %d to CuSolver Cholesky routine Zpotrf is invalid.\n", -info); stackTraceExit(1); }
		if(info>0) { logPrintf("Matrix not positive-definite at leading minor# %d in CuSolver Cholesky routine Zpotrf.\n", info); stackTraceExit(1); }
		//Solve:
		matrix x = b; //solution will happen in place
		cusolverDnZpotrs(cusolverHandle, uplo, N, Nrhs, (double2*)Acopy.dataPref(), N, (double2*)x.dataPref(), N, infoArr.dataPref());
		info = infoArr.data()[0];
		if(info<0) { logPrintf("Argument# %d to CuSolver solver routine Zpotrs is invalid.\n", -info); stackTraceExit(1); }
		return x;
	}
#endif
	//Apply inverse using LAPACK routine:
	char uplo = 'U';
	matrix Acopy = A; //destructible copy; routine will factorize matrix in place
	matrix x = b; //solution will happen in place
	int info = 0;
	zposv_(&uplo, &N, &Nrhs, Acopy.data(), &N, x.data(), &N, &info);
	if(info<0) { logPrintf("Argument# %d to LAPACK linear solve routine ZPOSV is invalid.\n", -info); stackTraceExit(1); }
	if(info>0) { logPrintf("Matrix not positive-definite at leading minor# %d in LAPACK linear solve routine ZPOSV.\n", info); stackTraceExit(1); }
	watch.stop();
	return x;
}
