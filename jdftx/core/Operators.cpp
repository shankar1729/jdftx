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

#include <core/Operators.h>
#include <core/Operators_internal.h>
#include <core/GridInfo.h>
#include <core/DataMultiplet.h>
#include <core/Random.h>
#include <string.h>

//------------------------------ Conversion operators ------------------------------


DataRptr Real(const complexDataRptr& C)
{	const GridInfo& gInfo = C->gInfo;
	DataRptr R; nullToZero(R, gInfo);
	callPref(eblas_daxpy)(gInfo.nr, C->scale, (const double*)C->dataPref(false), 2, R->dataPref(false), 1);
	R->scale = 1;
	return R;
}

void RealG_sub(size_t iStart, size_t iStop, const vector3<int> S, const complex* vFull, complex* vHalf, double scaleFac)
{	THREAD_halfGspaceLoop( RealG_calc(i, iG, S, vFull, vHalf, scaleFac); )
}
#ifdef GPU_ENABLED
void RealG_gpu(const vector3<int> S, const complex* vFull, complex* vHalf, double scaleFac);
#endif
DataGptr Real(const complexDataGptr& full)
{	const GridInfo& gInfo = full->gInfo;
	DataGptr half = DataG::alloc(gInfo, isGpuEnabled());
	#ifdef GPU_ENABLED
	RealG_gpu(gInfo.S, full->dataGpu(false), half->dataGpu(false), full->scale);
	#else
	threadLaunch(shouldThreadOperators() ? 0 : 1, //0 => max threads
		RealG_sub, gInfo.nG, gInfo.S, full->data(false), half->data(false), full->scale);
	#endif
	half->scale = 1;
	return half;
}

DataRptr Imag(const complexDataRptr& C)
{	const GridInfo& gInfo = C->gInfo;
	DataRptr I; nullToZero(I, gInfo);
	callPref(eblas_daxpy)(gInfo.nr, C->scale, ((const double*)C->dataPref(false))+1, 2, I->dataPref(false), 1);
	I->scale = 1;
	return I;
}

void ImagG_sub(size_t iStart, size_t iStop, const vector3<int> S, const complex* vFull, complex* vHalf, double scaleFac)
{	THREAD_halfGspaceLoop( ImagG_calc(i, iG, S, vFull, vHalf, scaleFac); )
}
#ifdef GPU_ENABLED
void ImagG_gpu(const vector3<int> S, const complex* vFull, complex* vHalf, double scaleFac);
#endif
DataGptr Imag(const complexDataGptr& full)
{	const GridInfo& gInfo = full->gInfo;
	DataGptr half = DataG::alloc(gInfo, isGpuEnabled());
	#ifdef GPU_ENABLED
	ImagG_gpu(gInfo.S, full->dataGpu(false), half->dataGpu(false), full->scale);
	#else
	threadLaunch(shouldThreadOperators() ? 0 : 1, //0 => max threads
		ImagG_sub, gInfo.nG, gInfo.S, full->data(false), half->data(false), full->scale);
	#endif
	half->scale = 1;
	return half;
}

complexDataRptr Complex(const DataRptr& R)
{	const GridInfo& gInfo = R->gInfo;
	complexDataRptr C; nullToZero(C, gInfo);
	callPref(eblas_daxpy)(gInfo.nr, R->scale, R->dataPref(false), 1, (double*)C->dataPref(false), 2);
	C->scale = 1;
	return C;
}

void ComplexG_sub(size_t iStart, size_t iStop, const vector3<int> S, const complex* vHalf, complex* vFull, double scaleFac)
{	THREAD_halfGspaceLoop( ComplexG_calc(i, iG, S, vHalf, vFull, scaleFac); )
}
#ifdef GPU_ENABLED
void ComplexG_gpu(const vector3<int> S, const complex* vHalf, complex* vFull, double scaleFac);
#endif
complexDataGptr Complex(const DataGptr& half)
{	const GridInfo& gInfo = half->gInfo;
	complexDataGptr full = complexDataG::alloc(gInfo, isGpuEnabled());
	#ifdef GPU_ENABLED
	ComplexG_gpu(gInfo.S, half->dataGpu(false), full->dataGpu(false), half->scale);
	#else
	threadLaunch(shouldThreadOperators() ? 0 : 1, //0 => max threads
		ComplexG_sub, gInfo.nG, gInfo.S, half->data(false), full->data(false), half->scale);
	#endif
	full->scale = 1;
	return full;
}


//------------------------------ Linear Unary operators ------------------------------

DataGptr O(const DataGptr& in) { return in * in->gInfo.detR; }
DataGptr O(DataGptr&& in) { return in *= in->gInfo.detR; }
complexDataGptr O(const complexDataGptr& in) { return in * in->gInfo.detR; }
complexDataGptr O(complexDataGptr&& in) { return in *= in->gInfo.detR; }


//Forward transform
DataRptr I(DataGptr&& in, bool compat)
{	//CPU c2r transforms destroy input, but this input can be destroyed
	DataRptr out(DataR::alloc(in->gInfo, isGpuEnabled()));
	#ifdef GPU_ENABLED
	cufftExecZ2D(compat ? in->gInfo.planZ2Dcompat : in->gInfo.planZ2D, (double2*)in->dataGpu(false), out->dataGpu(false));
	#else
	fftw_execute_dft_c2r(shouldThreadOperators() ? in->gInfo.planCtoRmulti : in->gInfo.planCtoRsingle,
		(fftw_complex*)in->data(false), out->data(false));
	#endif
	out->scale = in->scale;
	return out;
}
DataRptr I(const DataGptr& in, bool compat)
{	//CPU c2r transforms destroy input, hence copy input in that case (and let that copy be destroyed above)
	#ifdef GPU_ENABLED
	return I((DataGptr&&)in, compat);
	#else
	return  I(in->clone(), compat);
	#endif
}
complexDataRptr I(const complexDataGptr& in)
{	complexDataRptr out(complexDataR::alloc(in->gInfo, isGpuEnabled()));
	#ifdef GPU_ENABLED
	cufftExecZ2Z(in->gInfo.planZ2Z, (double2*)in->dataGpu(false), (double2*)out->dataGpu(false), CUFFT_INVERSE);
	#else
	fftw_execute_dft(shouldThreadOperators() ? in->gInfo.planInverseMulti : in->gInfo.planInverseSingle,
		(fftw_complex*)in->data(false), (fftw_complex*)out->data(false));
	#endif
	out->scale = in->scale;
	return out;
}
complexDataRptr I(complexDataGptr&& in)
{	//Destructible input (transform in place):
	#ifdef GPU_ENABLED
	cufftExecZ2Z(in->gInfo.planZ2Z, (double2*)in->dataGpu(false), (double2*)in->dataGpu(false), CUFFT_INVERSE);
	#else
	fftw_execute_dft(shouldThreadOperators() ? in->gInfo.planInverseInPlaceMulti : in->gInfo.planInverseInPlaceSingle,
		(fftw_complex*)in->data(false), (fftw_complex*)in->data(false));
	#endif
	return (complexDataRptr&&)in;
}

//Forward transform h.c.
DataGptr Idag(const DataRptr& in)
{	//r2c transform does not destroy input (no backing up needed)
	DataGptr out(DataG::alloc(in->gInfo, isGpuEnabled()));
	#ifdef GPU_ENABLED
	cufftExecD2Z(in->gInfo.planD2Z, in->dataGpu(false), (double2*)out->dataGpu(false));
	#else
	fftw_execute_dft_r2c(shouldThreadOperators() ? in->gInfo.planRtoCmulti : in->gInfo.planRtoCsingle,
		in->data(false), (fftw_complex*)out->data(false));
	#endif
	out->scale = in->scale;
	return out;
}
complexDataGptr Idag(const complexDataRptr& in)
{	complexDataGptr out(complexDataG::alloc(in->gInfo, isGpuEnabled()));
	#ifdef GPU_ENABLED
	cufftExecZ2Z(in->gInfo.planZ2Z, (double2*)in->dataGpu(false), (double2*)out->dataGpu(false), CUFFT_FORWARD);
	#else
	fftw_execute_dft(shouldThreadOperators() ? in->gInfo.planForwardMulti : in->gInfo.planForwardSingle,
		(fftw_complex*)in->data(false), (fftw_complex*)out->data(false));
	#endif
	out->scale = in->scale;
	return out;
}
complexDataGptr Idag(complexDataRptr&& in)
{	//Destructible input (transform in place):
	#ifdef GPU_ENABLED
	cufftExecZ2Z(in->gInfo.planZ2Z, (double2*)in->dataGpu(false), (double2*)in->dataGpu(false), CUFFT_FORWARD);
	#else
	fftw_execute_dft(shouldThreadOperators() ? in->gInfo.planForwardInPlaceMulti : in->gInfo.planForwardInPlaceSingle,
		(fftw_complex*)in->data(false), (fftw_complex*)in->data(false));
	#endif
	return (complexDataGptr&&)in;
}

//Reverse transform (same as Idag upto the normalization factor)
DataGptr J(const DataRptr& in) { return (1.0/in->gInfo.nr)*Idag(in); }
complexDataGptr J(const complexDataRptr& in) { return (1.0/in->gInfo.nr)*Idag(in); }
complexDataGptr J(complexDataRptr&& in) { return Idag((complexDataRptr&&)(in *= 1.0/in->gInfo.nr)); }

//Reverse transform h.c. (same as I upto the normalization factor)
DataRptr Jdag(const DataGptr& in, bool compat) { return (1.0/in->gInfo.nr)*I(in, compat); }
DataRptr Jdag(DataGptr&& in, bool compat) { return (1.0/in->gInfo.nr)*I(in, compat); }
complexDataRptr Jdag(const complexDataGptr& in) { return (1.0/in->gInfo.nr)*I(in); }
complexDataRptr Jdag(complexDataGptr&& in) { return I((complexDataGptr&&)(in *= 1.0/in->gInfo.nr)); }

DataRptr JdagOJ(const DataRptr& in) { return in * in->gInfo.dV; }
DataRptr JdagOJ(DataRptr&& in) { return in *= in->gInfo.dV; }
complexDataRptr JdagOJ(const complexDataRptr& in) { return in * in->gInfo.dV; }
complexDataRptr JdagOJ(complexDataRptr&& in) { return in *= in->gInfo.dV; }



inline void L_sub(int i, double Gsq, complex* v)
{	v[i] *= Gsq;
}
#ifdef GPU_ENABLED //implemented in Operators.cu
void L_gpu(const vector3<int> S, const matrix3<> GGT, complex* v);
#endif
DataGptr L(DataGptr&& in)
{	const GridInfo& gInfo = in->gInfo;
	in *= -gInfo.detR;
	#ifdef GPU_ENABLED
	L_gpu(gInfo.S, gInfo.GGT, in->dataGpu(false));
	#else
	applyFuncGsq(gInfo, L_sub, in->data(false));
	#endif
	return in;
}
DataGptr L(const DataGptr& in) { return L(in->clone()); }

inline void Linv_sub(int i, double Gsq, complex* v)
{	if(i==0) v[i]=0.0;
	else v[i] /= Gsq;
}
#ifdef GPU_ENABLED //implemented in Operators.cu
void Linv_gpu(const vector3<int> S, const matrix3<> GGT, complex* v);
#endif
DataGptr Linv(DataGptr&& in)
{	const GridInfo& gInfo = in->gInfo;
	in *= (-1.0/gInfo.detR);
	#ifdef GPU_ENABLED
	Linv_gpu(gInfo.S, gInfo.GGT, in->dataGpu(false));
	#else
	applyFuncGsq(gInfo, Linv_sub, in->data(false));
	#endif
	return in;
}
DataGptr Linv(const DataGptr& in) { return Linv(in->clone()); }


void fullL_sub(size_t iStart, size_t iStop, const vector3<int> S, const matrix3<> GGT, complex* v)
{	THREAD_fullGspaceLoop( v[i] *= GGT.metric_length_squared(iG); )
}
#ifdef GPU_ENABLED //implemented in Operators.cu
void fullL_gpu(const vector3<int> S, const matrix3<> GGT, complex* v);
#endif
complexDataGptr L(complexDataGptr&& in)
{	const GridInfo& gInfo = in->gInfo;
	in *= -gInfo.detR;
	#ifdef GPU_ENABLED
	fullL_gpu(gInfo.S, gInfo.GGT, in->dataGpu(false));
	#else
	threadLaunch(shouldThreadOperators() ? 0 : 1, //0 => max threads
		fullL_sub, gInfo.nr, gInfo.S, gInfo.GGT, in->data(false));
	#endif
	return in;
}
complexDataGptr L(const complexDataGptr& in) { return L(in->clone()); }

void fullLinv_sub(size_t iStart, size_t iStop, const vector3<int> S, const matrix3<> GGT, complex* v)
{	THREAD_fullGspaceLoop( v[i] *= i ? (1.0/GGT.metric_length_squared(iG)) : 0.0; )
}
#ifdef GPU_ENABLED //implemented in Operators.cu
void fullLinv_gpu(const vector3<int> S, const matrix3<> GGT, complex* v);
#endif
complexDataGptr Linv(complexDataGptr&& in)
{	const GridInfo& gInfo = in->gInfo;
	in *= (-1.0/gInfo.detR);
	#ifdef GPU_ENABLED
	fullLinv_gpu(gInfo.S, gInfo.GGT, in->dataGpu(false));
	#else
	threadLaunch(shouldThreadOperators() ? 0 : 1, //0 => max threads
		fullLinv_sub, gInfo.nr, gInfo.S, gInfo.GGT, in->data(false));
	#endif
	return in;
}
complexDataGptr Linv(const complexDataGptr& in) { return Linv(in->clone()); }



template<typename Scalar> void zeroNyquist_sub(size_t iStart, size_t iStop, const vector3<int> S, Scalar* data)
{	THREAD_halfGspaceLoop( if(IS_NYQUIST) data[i] = Scalar(0.0); )
}
void zeroNyquist(RealKernel& K)
{	threadLaunch(shouldThreadOperators() ? 0 : 1, //0 => max threads
		zeroNyquist_sub<double>, K.gInfo.nG, K.gInfo.S, K.data);
}
void zeroNyquist(DataGptr& Gptr)
{	const GridInfo& gInfo = Gptr->gInfo;
	threadLaunch(shouldThreadOperators() ? 0 : 1, //0 => max threads
		zeroNyquist_sub<complex>, gInfo.nG, gInfo.S, Gptr->data());
}
void zeroNyquist(DataRptr& Rptr) { DataGptr Rtilde=J(Rptr); zeroNyquist(Rtilde); Rptr = I((DataGptr&&)Rtilde); }

//------------------------------ Nonlinear Unary operators ------------------------------

void exp_sub(size_t i, double* X, double prefac) { X[i] = exp(prefac*X[i]); }
#ifdef GPU_ENABLED
void exp_gpu(int N, double* X, double prefac); //in operators.cu
#endif
DataRptr exp(DataRptr&& X)
{
	#ifdef GPU_ENABLED
	exp_gpu(X->nElem, X->dataGpu(false), X->scale);
	#else
	threadedLoop(exp_sub, X->nElem, X->data(false), X->scale);
	#endif
	X->scale = 1.0;
	return X;
}
DataRptr exp(const DataRptr& X) { return exp(X->clone()); }

void log_sub(size_t i, double* X, double prefac) { X[i] = log(prefac*X[i]); }
#ifdef GPU_ENABLED
void log_gpu(int N, double* X, double prefac); //in operators.cu
#endif
DataRptr log(DataRptr&& X)
{
	#ifdef GPU_ENABLED
	log_gpu(X->nElem, X->dataGpu(false), X->scale);
	#else
	threadedLoop(log_sub, X->nElem, X->data(false), X->scale);
	#endif
	X->scale = 1.0;
	return X;
}
DataRptr log(const DataRptr& X) { return log(X->clone()); }


void sqrt_sub(size_t i, double* X, double prefac) { X[i] = sqrt(prefac*X[i]); }
#ifdef GPU_ENABLED
void sqrt_gpu(int N, double* X, double prefac); //in operators.cu
#endif
DataRptr sqrt(DataRptr&& X)
{
	#ifdef GPU_ENABLED
	sqrt_gpu(X->nElem, X->dataGpu(false), X->scale);
	#else
	threadedLoop(sqrt_sub, X->nElem, X->data(false), X->scale);
	#endif
	X->scale = 1.0;
	return X;
}
DataRptr sqrt(const DataRptr& X) { return sqrt(X->clone()); }


void inv_sub(size_t i, double* X, double prefac) { X[i] = prefac/X[i]; }
#ifdef GPU_ENABLED
void inv_gpu(int N, double* X, double prefac); //in operators.cu
#endif
DataRptr inv(DataRptr&& X)
{
	#ifdef GPU_ENABLED
	inv_gpu(X->nElem, X->dataGpu(false), 1.0/X->scale);
	#else
	threadedLoop(inv_sub, X->nElem, X->data(false), 1.0/X->scale);
	#endif
	X->scale = 1.0;
	return X;
}
DataRptr inv(const DataRptr& X) { return inv(X->clone()); }

void pow_sub(size_t i, double* X, double scale, double alpha) { X[i] = pow(scale*X[i],alpha); }
#ifdef GPU_ENABLED
void pow_gpu(int N, double* X, double scale, double alpha); //in operators.cu
#endif
DataRptr pow(DataRptr&& X, double alpha)
{
	#ifdef GPU_ENABLED
	pow_gpu(X->nElem, X->dataGpu(false), X->scale, alpha);
	#else
	threadedLoop(pow_sub, X->nElem, X->data(false), X->scale, alpha);
	#endif
	X->scale = 1.0;
	return X;
}
DataRptr pow(const DataRptr& X, double alpha) { return pow(X->clone(), alpha); }


//------------------------------ Multiplication operators------------------------------

DataRptr& operator*=(DataRptr& in, const DataRptr& other)
{	in->scale *= other->scale;
	callPref(eblas_dmul)(in->nElem, other->dataPref(false), 1, in->dataPref(false), 1);
	return in;
}

complexDataRptr& operator*=(complexDataRptr& in, const DataRptr& other)
{	in->scale *= other->scale;
	callPref(eblas_zmuld)(in->nElem, other->dataPref(false), 1, in->dataPref(false), 1);
	return in;
}
complexDataRptr operator*(const complexDataRptr& inC, const DataRptr& inR) { complexDataRptr out(inC->clone()); return out *= inR; }
complexDataRptr operator*(const DataRptr& inR, const complexDataRptr& inC) { complexDataRptr out(inC->clone()); return out *= inR; }
complexDataRptr operator*(complexDataRptr&& inC, const DataRptr& inR) { return inC *= inR; }
complexDataRptr operator*(const DataRptr& inR, complexDataRptr&& inC) { return inC *= inR; }

DataGptr& operator*=(DataGptr& inG, const RealKernel& inR)
{	callPref(eblas_zmuld)(inG->nElem, inR.dataPref, 1, inG->dataPref(false), 1);
	return inG;
}
DataGptr operator*(const RealKernel& inR, const DataGptr& inG) { DataGptr out(inG->clone()); return out *= inR; }
DataGptr operator*(const DataGptr& inG, const RealKernel& inR) { DataGptr out(inG->clone()); return out *= inR; }
DataGptr operator*(const RealKernel& inR, DataGptr&& inG) { return inG *= inR; }
DataGptr operator*(DataGptr&& inG, const RealKernel& inR) { return inG *= inR; }


//------------------------------ Linear combine operators ------------------------------

void axpy(double alpha, const DataRptr& X, DataRptr& Y)
{	if(X)
	{	if(Y)
		{	if(Y->scale == 0.0) { Y = X * alpha; }
			else callPref(eblas_daxpy)(X->nElem, alpha*X->scale/Y->scale, X->dataPref(false), 1, Y->dataPref(false), 1);
		}
		else Y = X * alpha;
	}
	//if X is null, nothing needs to be done, Y remains unchanged
}
DataRptr& operator+=(DataRptr& in, double scalar)
{	Data dataScalar(in->gInfo, 1, 1, false); *((double*)dataScalar.data()) = scalar;
	callPref(eblas_daxpy)(in->nElem, 1.0, (double*)dataScalar.dataPref(), 0, in->dataPref(), 1);
	return in;
}
DataRptr operator+(double scalar, const DataRptr& in) { DataRptr out(in->clone()); return out += scalar; }
DataRptr operator+(const DataRptr& in, double scalar) { DataRptr out(in->clone()); return out += scalar; }
DataRptr operator+(double scalar, DataRptr&& in) { return in += scalar; }
DataRptr operator+(DataRptr&& in, double scalar) { return in += scalar; }
DataRptr& operator-=(DataRptr& in, double scalar)
{	return (in += -scalar);
}
DataRptr operator-(double scalar, const DataRptr& in) { DataRptr out(in->clone()); return (out *= -1.0) += scalar; }
DataRptr operator-(const DataRptr& in, double scalar) { DataRptr out(in->clone()); return out -= scalar; }
DataRptr operator-(double scalar, DataRptr&& in) { return (in *= -1.0) += scalar; }
DataRptr operator-(DataRptr&& in, double scalar) { return in -= scalar; }


//------------------------------ Dot products and 2-norms ------------------------------

double dot(const DataRptr& X, const DataRptr& Y)
{	return X->scale * Y->scale * callPref(eblas_ddot)(X->nElem, X->dataPref(false), 1, Y->dataPref(false), 1);
}

double dot(const DataGptr& X, const DataGptr& Y)
{	int N = X->nElem;
	int S2 = X->gInfo.S[2]/2 + 1; //inner dimension
	int S01 = X->gInfo.S[0] * X->gInfo.S[1]; //number of inner dimension slices
	complex complexDot = callPref(eblas_zdotc)(N, X->dataPref(false), 1, Y->dataPref(false), 1);
	complex correction1 = callPref(eblas_zdotc)(S01, X->dataPref(false), S2, Y->dataPref(false), S2);
	complex correction2 = callPref(eblas_zdotc)(S01, X->dataPref(false)+S2-1, S2, Y->dataPref(false)+S2-1, S2);
	if(S2==1) correction2=complex(0,0); //because slices 1 and 2 are the same
	return X->scale * Y->scale * (2.0*complexDot - correction1 - correction2).real();
}

double nrm2(const DataRptr& X)
{	return fabs(X->scale) * callPref(eblas_dnrm2)(X->nElem, X->dataPref(false), 1);
}

double nrm2(const DataGptr& X)
{	int N = X->nElem;
	int S2 = X->gInfo.S[2]/2 + 1; //inner dimension
	int S01 = X->gInfo.S[0] * X->gInfo.S[1]; //number of inner dimension slices
	double complexNorm = callPref(eblas_dznrm2)(N, X->dataPref(false), 1);
	double correction1 = callPref(eblas_dznrm2)(S01, X->dataPref(false), S2);
	double correction2 = callPref(eblas_dznrm2)(S01, X->dataPref(false)+S2-1, S2);
	if(S2==1) correction2=0; //because slices 1 and 2 are the same
	return fabs(X->scale) * sqrt(2*pow(complexNorm,2) - pow(correction1,2) - pow(correction2,2));
}

double sum(const DataRptr& X)
{	Data dataScale(X->gInfo, 1, 1, false); *((double*)dataScale.data()) = X->scale;
	return callPref(eblas_ddot)(X->nElem, X->dataPref(false), 1, (double*)dataScale.dataPref(false), 0);
}

double sum(const DataGptr& X)
{	int N = X->nElem;
	int S2 = X->gInfo.S[2]/2 + 1; //inner dimension
	int S01 = X->gInfo.S[0] * X->gInfo.S[1]; //number of inner dimension slices
	Data dataOne(X->gInfo, 1, 2, false); *((complex*)dataOne.data()) = 1.0;
	complex complexSum = callPref(eblas_zdotc)(N, X->dataPref(false), 1, (complex*)dataOne.dataPref(false), 0);
	complex correction1 = callPref(eblas_zdotc)(S01, X->dataPref(false), S2, (complex*)dataOne.dataPref(false), 0);
	complex correction2 = callPref(eblas_zdotc)(S01, X->dataPref(false)+S2-1, S2, (complex*)dataOne.dataPref(false), 0);
	if(S2==1) correction2 = 0; //because slices 1 and 2 are the same
	return X->scale * (2.0*complexSum - correction1 - correction2).real();
}



double integral(const DataRptr& X)
{	return X->gInfo.dV * sum(X);
}
double integral(const DataGptr& X)
{	double XdataZero;
	#ifdef GPU_ENABLED
	cudaMemcpy(&XdataZero, X->dataGpu(false), sizeof(double), cudaMemcpyDeviceToHost);
	#else
	XdataZero = X->data(false)[0].real();
	#endif
	return XdataZero * X->gInfo.detR * X->scale;
}
complex integral(const complexDataRptr& X)
{	return X->gInfo.dV * sum(X);
}
complex integral(const complexDataGptr& X)
{	complex XdataZero;
	#ifdef GPU_ENABLED
	cudaMemcpy(&XdataZero, X->dataGpu(false), sizeof(complex), cudaMemcpyDeviceToHost);
	#else
	XdataZero = X->data(false)[0];
	#endif
	return XdataZero * X->gInfo.detR * X->scale;
}

//------------------------------ Grid conversion utilities ------------------------------

void changeGrid_sub(size_t iStart, size_t iStop, const vector3<int>& S, const vector3<int>& Sin, const vector3<int>& Sout, const complex* in, complex* out)
{	THREAD_halfGspaceLoop( changeGrid_calc(iG, Sin, Sout, in, out); )
}
void changeGrid(const vector3<int>& S, const vector3<int>& Sin, const vector3<int>& Sout, const complex* in, complex* out)
{	threadLaunch(changeGrid_sub, S[0]*S[1]*(1+S[2]/2), S, Sin, Sout, in, out);
}
#ifdef GPU_ENABLED
void changeGrid_gpu(const vector3<int>& S, const vector3<int>& Sin, const vector3<int>& Sout, const complex* in, complex* out);
#endif
DataGptr changeGrid(const DataGptr& in, const GridInfo& gInfoNew)
{	static StopWatch watch("changeGrid"); watch.start();
	DataGptr out; nullToZero(out, gInfoNew);
	assert(gInfoNew.R == in->gInfo.R);
	const vector3<int>& Sin = in->gInfo.S;
	const vector3<int>& Sout = gInfoNew.S;
	vector3<int> Smax; for(int k=0; k<3; k++) Smax[k] = std::max(Sin[k],Sout[k]);
	callPref(changeGrid)(Smax, Sin, Sout, in->dataPref(), out->dataPref());
	watch.stop();
	return out;
}

DataRptr changeGrid(const DataRptr& in, const GridInfo& gInfoNew)
{	return I(changeGrid(J(in), gInfoNew), true);
}


void changeGridFull_sub(size_t iStart, size_t iStop, const vector3<int>& S, const vector3<int>& Sin, const vector3<int>& Sout, const complex* in, complex* out)
{	THREAD_fullGspaceLoop( changeGridFull_calc(iG, Sin, Sout, in, out); )
}
void changeGridFull(const vector3<int>& S, const vector3<int>& Sin, const vector3<int>& Sout, const complex* in, complex* out)
{	threadLaunch(changeGridFull_sub, S[0]*S[1]*S[2], S, Sin, Sout, in, out);
}
#ifdef GPU_ENABLED
void changeGridFull_gpu(const vector3<int>& S, const vector3<int>& Sin, const vector3<int>& Sout, const complex* in, complex* out);
#endif
complexDataGptr changeGrid(const complexDataGptr& in, const GridInfo& gInfoNew)
{	static StopWatch watch("changeGridFull"); watch.start();
	complexDataGptr out; nullToZero(out, gInfoNew);
	assert(gInfoNew.R == in->gInfo.R);
	const vector3<int>& Sin = in->gInfo.S;
	const vector3<int>& Sout = gInfoNew.S;
	vector3<int> Smax; for(int k=0; k<3; k++) Smax[k] = std::max(Sin[k],Sout[k]);
	callPref(changeGridFull)(Smax, Sin, Sout, in->dataPref(), out->dataPref());
	watch.stop();
	return out;
}

complexDataRptr changeGrid(const complexDataRptr& in, const GridInfo& gInfoNew)
{	return I(changeGrid(J(in), gInfoNew));
}


//------------------------------ Initialization utilities ------------------------------

void initRandom(DataRptr& X, double cap)
{	double* Xdata = X->data();
	for(int i=0; i<X->nElem; i++)
		Xdata[i] = Random::normal(0, 1, cap);
}

void initRandomFlat(DataRptr& X)
{	double* Xdata = X->data();
	for(int i=0; i<X->nElem; i++)
		Xdata[i] = Random::uniform();
}

void initGaussianKernel_sub(int i, double Gsq, double expfac, double* X)
{	X[i] = exp(expfac*Gsq);
}
void initGaussianKernel(RealKernel& X, double x0)
{	applyFuncGsq(X.gInfo, initGaussianKernel_sub, -pow(0.5*x0,2), X.data);
	X.set();
}

void initTranslation_sub(size_t iStart, size_t iStop, const vector3<int> S, const vector3<> Gr, complex* X)
{	THREAD_halfGspaceLoop( X[i] = cis(-dot(iG,Gr)); )
}
void initTranslation(DataGptr& X, const vector3<>& r)
{	const GridInfo& gInfo = X->gInfo;
	threadLaunch(shouldThreadOperators() ? 0 : 1, //0 => max threads
		initTranslation_sub, gInfo.nG, gInfo.S, gInfo.G*r, X->data());
}


void gaussConvolve_sub(size_t iStart, size_t iStop, const vector3<int>& S, const matrix3<>& GGT, complex* data, double sigma)
{	THREAD_halfGspaceLoop( data[i] *= exp(-0.5*sigma*sigma*GGT.metric_length_squared(iG)); )
}
void gaussConvolve(const vector3<int>& S, const matrix3<>& GGT, complex* data, double sigma)
{	threadLaunch(shouldThreadOperators() ? 0 : 1, //0 => max threads
		gaussConvolve_sub, S[0]*S[1]*(1+S[2]/2), S, GGT, data, sigma);
}
#ifdef GPU_ENABLED
void gaussConvolve_gpu(const vector3<int>& S, const matrix3<>& GGT, complex* data, double sigma);
#endif
DataGptr gaussConvolve(DataGptr&& in, double sigma)
{	assert(in);
	callPref(gaussConvolve)(in->gInfo.S, in->gInfo.GGT, in->dataPref(false), sigma);
	return in;
}
DataGptr gaussConvolve(const DataGptr& in, double sigma)
{	DataGptr out(in->clone());
	return gaussConvolve((DataGptr&&)out, sigma);
}

//------------------------------ Debug utilities ------------------------------

void printStats(const DataRptr& X, const char* name, FILE* fp)
{	int N = X->nElem;
	double mean = sum(X)/N;
	double stdDev = sqrt(fabs(dot(X,X)/N - mean*mean));
	fprintf(fp, "vector3 %s\t= %.15le +/- %.15le\n", name, mean, stdDev);
}

//------------------------------ From DataMultiplet.h ------------------------------

inline void gradient_sub(size_t iStart, size_t iStop, const vector3<int> S,
	const matrix3<> G, const complex* Xtilde, vector3<complex*> gradTilde)
{	THREAD_halfGspaceLoop( gradient_calc(i, iG, IS_NYQUIST, G, Xtilde, gradTilde); )
}
#ifdef GPU_ENABLED //implemented in Operators.cu
void gradient_gpu(const vector3<int> S, const matrix3<> G, const complex* Xtilde, vector3<complex*> gradTilde);
#endif
DataGptrVec gradient(const DataGptr& Xtilde)
{	const GridInfo& gInfo = Xtilde->gInfo;
	DataGptrVec gradTilde(gInfo, isGpuEnabled());
	#ifdef GPU_ENABLED
	gradient_gpu(gInfo.S, gInfo.G, Xtilde->dataGpu(), gradTilde.dataGpu());
	#else
	threadLaunch(shouldThreadOperators() ? 0 : 1, //0 => max threads
		gradient_sub, gInfo.nG, gInfo.S, gInfo.G, Xtilde->data(), gradTilde.data());
	#endif
	return gradTilde;
}
DataRptrVec gradient(const DataRptr& X) { return I(gradient(J(X))); }


inline void divergence_sub(size_t iStart, size_t iStop, const vector3<int> S,
	const matrix3<> G, vector3<const complex*> Vtilde, complex* divTilde)
{	THREAD_halfGspaceLoop( divergence_calc(i, iG, IS_NYQUIST, G, Vtilde, divTilde); )
}
#ifdef GPU_ENABLED
void divergence_gpu(const vector3<int> S, const matrix3<> G, vector3<const complex*> Vtilde, complex* divTilde);
#endif
DataGptr divergence(const DataGptrVec& Vtilde)
{	const GridInfo& gInfo = Vtilde[0]->gInfo;
	DataGptr divTilde(DataG::alloc(gInfo, isGpuEnabled()));
	#ifdef GPU_ENABLED
	divergence_gpu(gInfo.S, gInfo.G, Vtilde.dataGpu(), divTilde->dataGpu());
	#else
	threadLaunch(shouldThreadOperators() ? 0 : 1, //0 => max threads
		divergence_sub, gInfo.nG, gInfo.S, gInfo.G, Vtilde.data(), divTilde->data());
	#endif
	return divTilde;
}
DataRptr divergence(const DataRptrVec& V) { return I(divergence(J(V))); }


inline void tensorGradient_sub(size_t iStart, size_t iStop, const vector3<int> S,
	const matrix3<> G, const complex* Xtilde, tensor3<complex*> gradTilde)
{	THREAD_halfGspaceLoop( tensorGradient_calc(i, iG, IS_NYQUIST, G, Xtilde, gradTilde); )
}
#ifdef GPU_ENABLED //implemented in Operators.cu
void tensorGradient_gpu(const vector3<int> S, const matrix3<> G, const complex* Xtilde, tensor3<complex*> gradTilde);
#endif
DataGptrTensor tensorGradient(const DataGptr& Xtilde)
{	const GridInfo& gInfo = Xtilde->gInfo;
	DataGptrTensor gradTilde(gInfo, isGpuEnabled());
	#ifdef GPU_ENABLED
	tensorGradient_gpu(gInfo.S, gInfo.G, Xtilde->dataGpu(), gradTilde.dataGpu());
	#else
	threadLaunch(shouldThreadOperators() ? 0 : 1, //0 => max threads
		tensorGradient_sub, gInfo.nG, gInfo.S, gInfo.G, Xtilde->data(), gradTilde.data());
	#endif
	return gradTilde;
}


inline void tensorDivergence_sub(size_t iStart, size_t iStop, const vector3<int> S,
	const matrix3<> G, tensor3<const complex*> Vtilde, complex* divTilde)
{	THREAD_halfGspaceLoop( tensorDivergence_calc(i, iG, IS_NYQUIST, G, Vtilde, divTilde); )
}
#ifdef GPU_ENABLED
void tensorDivergence_gpu(const vector3<int> S, const matrix3<> G, tensor3<const complex*> Vtilde, complex* divTilde);
#endif
DataGptr tensorDivergence(const DataGptrTensor& Vtilde)
{	const GridInfo& gInfo = Vtilde[0]->gInfo;
	DataGptr divTilde(DataG::alloc(gInfo, isGpuEnabled()));
	#ifdef GPU_ENABLED
	tensorDivergence_gpu(gInfo.S, gInfo.G, Vtilde.dataGpu(), divTilde->dataGpu());
	#else
	threadLaunch(shouldThreadOperators() ? 0 : 1, //0 => max threads
		tensorDivergence_sub, gInfo.nG, gInfo.S, gInfo.G, Vtilde.data(), divTilde->data());
	#endif
	return divTilde;
}

