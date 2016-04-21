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
#include <core/VectorField.h>
#include <core/Random.h>
#include <string.h>

//------------------------------ Conversion operators ------------------------------


ScalarField Real(const complexScalarField& C)
{	const GridInfo& gInfo = C->gInfo;
	ScalarField R; nullToZero(R, gInfo);
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
ScalarFieldTilde Real(const complexScalarFieldTilde& full)
{	const GridInfo& gInfo = full->gInfo;
	ScalarFieldTilde half = ScalarFieldTildeData::alloc(gInfo, isGpuEnabled());
	#ifdef GPU_ENABLED
	RealG_gpu(gInfo.S, full->dataGpu(false), half->dataGpu(false), full->scale);
	#else
	threadLaunch(RealG_sub, gInfo.nG, gInfo.S, full->data(false), half->data(false), full->scale);
	#endif
	half->scale = 1;
	return half;
}

ScalarField Imag(const complexScalarField& C)
{	const GridInfo& gInfo = C->gInfo;
	ScalarField I; nullToZero(I, gInfo);
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
ScalarFieldTilde Imag(const complexScalarFieldTilde& full)
{	const GridInfo& gInfo = full->gInfo;
	ScalarFieldTilde half = ScalarFieldTildeData::alloc(gInfo, isGpuEnabled());
	#ifdef GPU_ENABLED
	ImagG_gpu(gInfo.S, full->dataGpu(false), half->dataGpu(false), full->scale);
	#else
	threadLaunch(ImagG_sub, gInfo.nG, gInfo.S, full->data(false), half->data(false), full->scale);
	#endif
	half->scale = 1;
	return half;
}

complexScalarField Complex(const ScalarField& R)
{	const GridInfo& gInfo = R->gInfo;
	complexScalarField C; nullToZero(C, gInfo);
	callPref(eblas_daxpy)(gInfo.nr, R->scale, R->dataPref(false), 1, (double*)C->dataPref(false), 2);
	C->scale = 1;
	return C;
}

complexScalarField Complex(const ScalarField& re, const ScalarField& im)
{	const GridInfo& gInfo = re->gInfo;
	complexScalarField C; nullToZero(C, gInfo);
	callPref(eblas_daxpy)(gInfo.nr, re->scale, re->dataPref(false), 1, (double*)(C->dataPref(false))+0, 2);
	callPref(eblas_daxpy)(gInfo.nr, im->scale, im->dataPref(false), 1, (double*)(C->dataPref(false))+1, 2);
	C->scale = 1;
	return C;
}

void ComplexG_sub(size_t iStart, size_t iStop, const vector3<int> S, const complex* vHalf, complex* vFull, double scaleFac)
{	THREAD_halfGspaceLoop( ComplexG_calc(i, iG, S, vHalf, vFull, scaleFac); )
}
#ifdef GPU_ENABLED
void ComplexG_gpu(const vector3<int> S, const complex* vHalf, complex* vFull, double scaleFac);
#endif
complexScalarFieldTilde Complex(const ScalarFieldTilde& half)
{	const GridInfo& gInfo = half->gInfo;
	complexScalarFieldTilde full = complexScalarFieldTildeData::alloc(gInfo, isGpuEnabled());
	#ifdef GPU_ENABLED
	ComplexG_gpu(gInfo.S, half->dataGpu(false), full->dataGpu(false), half->scale);
	#else
	threadLaunch(ComplexG_sub, gInfo.nG, gInfo.S, half->data(false), full->data(false), half->scale);
	#endif
	full->scale = 1;
	return full;
}


//------------------------------ Linear Unary operators ------------------------------

ScalarFieldTilde O(const ScalarFieldTilde& in) { return in * in->gInfo.detR; }
ScalarFieldTilde O(ScalarFieldTilde&& in) { return in *= in->gInfo.detR; }
complexScalarFieldTilde O(const complexScalarFieldTilde& in) { return in * in->gInfo.detR; }
complexScalarFieldTilde O(complexScalarFieldTilde&& in) { return in *= in->gInfo.detR; }


//Forward transform
ScalarField I(ScalarFieldTilde&& in, bool compat, int nThreads)
{	//CPU c2r transforms destroy input, but this input can be destroyed
	ScalarField out(ScalarFieldData::alloc(in->gInfo, isGpuEnabled()));
	#ifdef GPU_ENABLED
	cufftExecZ2D(compat ? in->gInfo.planZ2Dcompat : in->gInfo.planZ2D, (double2*)in->dataGpu(false), out->dataGpu(false));
	#else
	if(!nThreads) nThreads = shouldThreadOperators() ? nProcsAvailable : 1;
	fftw_execute_dft_c2r(in->gInfo.getPlan(GridInfo::PlanCtoR, nThreads),
		(fftw_complex*)in->data(false), out->data(false));
	#endif
	out->scale = in->scale;
	return out;
}
ScalarField I(const ScalarFieldTilde& in, bool compat, int nThreads)
{	//CPU c2r transforms destroy input, hence copy input in that case (and let that copy be destroyed above)
	#ifdef GPU_ENABLED
	return I((ScalarFieldTilde&&)in, compat, nThreads);
	#else
	return  I(in->clone(), compat, nThreads);
	#endif
}
complexScalarField I(const complexScalarFieldTilde& in, int nThreads)
{	complexScalarField out(complexScalarFieldData::alloc(in->gInfo, isGpuEnabled()));
	#ifdef GPU_ENABLED
	cufftExecZ2Z(in->gInfo.planZ2Z, (double2*)in->dataGpu(false), (double2*)out->dataGpu(false), CUFFT_INVERSE);
	#else
	if(!nThreads) nThreads = shouldThreadOperators() ? nProcsAvailable : 1;
	fftw_execute_dft(in->gInfo.getPlan(GridInfo::PlanInverse, nThreads),
		(fftw_complex*)in->data(false), (fftw_complex*)out->data(false));
	#endif
	out->scale = in->scale;
	return out;
}
complexScalarField I(complexScalarFieldTilde&& in, int nThreads)
{	//Destructible input (transform in place):
	#ifdef GPU_ENABLED
	cufftExecZ2Z(in->gInfo.planZ2Z, (double2*)in->dataGpu(false), (double2*)in->dataGpu(false), CUFFT_INVERSE);
	#else
	if(!nThreads) nThreads = shouldThreadOperators() ? nProcsAvailable : 1;
	fftw_execute_dft(in->gInfo.getPlan(GridInfo::PlanInverseInPlace, nThreads),
		(fftw_complex*)in->data(false), (fftw_complex*)in->data(false));
	#endif
	return std::static_pointer_cast<complexScalarFieldData>(std::static_pointer_cast<FieldData>(in));
}

//Forward transform h.c.
ScalarFieldTilde Idag(const ScalarField& in, int nThreads)
{	//r2c transform does not destroy input (no backing up needed)
	ScalarFieldTilde out(ScalarFieldTildeData::alloc(in->gInfo, isGpuEnabled()));
	#ifdef GPU_ENABLED
	cufftExecD2Z(in->gInfo.planD2Z, in->dataGpu(false), (double2*)out->dataGpu(false));
	#else
	if(!nThreads) nThreads = shouldThreadOperators() ? nProcsAvailable : 1;
	fftw_execute_dft_r2c(in->gInfo.getPlan(GridInfo::PlanRtoC, nThreads),
		in->data(false), (fftw_complex*)out->data(false));
	#endif
	out->scale = in->scale;
	return out;
}
complexScalarFieldTilde Idag(const complexScalarField& in, int nThreads)
{	complexScalarFieldTilde out(complexScalarFieldTildeData::alloc(in->gInfo, isGpuEnabled()));
	#ifdef GPU_ENABLED
	cufftExecZ2Z(in->gInfo.planZ2Z, (double2*)in->dataGpu(false), (double2*)out->dataGpu(false), CUFFT_FORWARD);
	#else
	if(!nThreads) nThreads = shouldThreadOperators() ? nProcsAvailable : 1;
	fftw_execute_dft(in->gInfo.getPlan(GridInfo::PlanForward, nThreads),
		(fftw_complex*)in->data(false), (fftw_complex*)out->data(false));
	#endif
	out->scale = in->scale;
	return out;
}
complexScalarFieldTilde Idag(complexScalarField&& in, int nThreads)
{	//Destructible input (transform in place):
	#ifdef GPU_ENABLED
	cufftExecZ2Z(in->gInfo.planZ2Z, (double2*)in->dataGpu(false), (double2*)in->dataGpu(false), CUFFT_FORWARD);
	#else
	if(!nThreads) nThreads = shouldThreadOperators() ? nProcsAvailable : 1;
	fftw_execute_dft(in->gInfo.getPlan(GridInfo::PlanForwardInPlace, nThreads),
		(fftw_complex*)in->data(false), (fftw_complex*)in->data(false));
	#endif
	return std::static_pointer_cast<complexScalarFieldTildeData>(std::static_pointer_cast<FieldData>(in));
}

//Reverse transform (same as Idag upto the normalization factor)
ScalarFieldTilde J(const ScalarField& in, int nThreads) { return (1.0/in->gInfo.nr)*Idag(in, nThreads); }
complexScalarFieldTilde J(const complexScalarField& in, int nThreads) { return (1.0/in->gInfo.nr)*Idag(in, nThreads); }
complexScalarFieldTilde J(complexScalarField&& in, int nThreads) { return Idag((complexScalarField&&)(in *= 1.0/in->gInfo.nr), nThreads); }

//Reverse transform h.c. (same as I upto the normalization factor)
ScalarField Jdag(const ScalarFieldTilde& in, bool compat, int nThreads) { return (1.0/in->gInfo.nr)*I(in, compat, nThreads); }
ScalarField Jdag(ScalarFieldTilde&& in, bool compat, int nThreads) { return (1.0/in->gInfo.nr)*I(in, compat, nThreads); }
complexScalarField Jdag(const complexScalarFieldTilde& in, int nThreads) { return (1.0/in->gInfo.nr)*I(in, nThreads); }
complexScalarField Jdag(complexScalarFieldTilde&& in, int nThreads) { return I((complexScalarFieldTilde&&)(in *= 1.0/in->gInfo.nr), nThreads); }

ScalarField JdagOJ(const ScalarField& in) { return in * in->gInfo.dV; }
ScalarField JdagOJ(ScalarField&& in) { return in *= in->gInfo.dV; }
complexScalarField JdagOJ(const complexScalarField& in) { return in * in->gInfo.dV; }
complexScalarField JdagOJ(complexScalarField&& in) { return in *= in->gInfo.dV; }



inline void L_sub(int i, double Gsq, complex* v)
{	v[i] *= Gsq;
}
#ifdef GPU_ENABLED //implemented in Operators.cu
void L_gpu(const vector3<int> S, const matrix3<> GGT, complex* v);
#endif
ScalarFieldTilde L(ScalarFieldTilde&& in)
{	const GridInfo& gInfo = in->gInfo;
	in *= -gInfo.detR;
	#ifdef GPU_ENABLED
	L_gpu(gInfo.S, gInfo.GGT, in->dataGpu(false));
	#else
	applyFuncGsq(gInfo, L_sub, in->data(false));
	#endif
	return in;
}
ScalarFieldTilde L(const ScalarFieldTilde& in) { return L(in->clone()); }

inline void Linv_sub(int i, double Gsq, complex* v)
{	if(i==0) v[i]=0.0;
	else v[i] /= Gsq;
}
#ifdef GPU_ENABLED //implemented in Operators.cu
void Linv_gpu(const vector3<int> S, const matrix3<> GGT, complex* v);
#endif
ScalarFieldTilde Linv(ScalarFieldTilde&& in)
{	const GridInfo& gInfo = in->gInfo;
	in *= (-1.0/gInfo.detR);
	#ifdef GPU_ENABLED
	Linv_gpu(gInfo.S, gInfo.GGT, in->dataGpu(false));
	#else
	applyFuncGsq(gInfo, Linv_sub, in->data(false));
	#endif
	return in;
}
ScalarFieldTilde Linv(const ScalarFieldTilde& in) { return Linv(in->clone()); }


void fullL_sub(size_t iStart, size_t iStop, const vector3<int> S, const matrix3<> GGT, complex* v)
{	THREAD_fullGspaceLoop( v[i] *= GGT.metric_length_squared(iG); )
}
#ifdef GPU_ENABLED //implemented in Operators.cu
void fullL_gpu(const vector3<int> S, const matrix3<> GGT, complex* v);
#endif
complexScalarFieldTilde L(complexScalarFieldTilde&& in)
{	const GridInfo& gInfo = in->gInfo;
	in *= -gInfo.detR;
	#ifdef GPU_ENABLED
	fullL_gpu(gInfo.S, gInfo.GGT, in->dataGpu(false));
	#else
	threadLaunch(fullL_sub, gInfo.nr, gInfo.S, gInfo.GGT, in->data(false));
	#endif
	return in;
}
complexScalarFieldTilde L(const complexScalarFieldTilde& in) { return L(in->clone()); }

void fullLinv_sub(size_t iStart, size_t iStop, const vector3<int> S, const matrix3<> GGT, complex* v)
{	THREAD_fullGspaceLoop( v[i] *= i ? (1.0/GGT.metric_length_squared(iG)) : 0.0; )
}
#ifdef GPU_ENABLED //implemented in Operators.cu
void fullLinv_gpu(const vector3<int> S, const matrix3<> GGT, complex* v);
#endif
complexScalarFieldTilde Linv(complexScalarFieldTilde&& in)
{	const GridInfo& gInfo = in->gInfo;
	in *= (-1.0/gInfo.detR);
	#ifdef GPU_ENABLED
	fullLinv_gpu(gInfo.S, gInfo.GGT, in->dataGpu(false));
	#else
	threadLaunch(fullLinv_sub, gInfo.nr, gInfo.S, gInfo.GGT, in->data(false));
	#endif
	return in;
}
complexScalarFieldTilde Linv(const complexScalarFieldTilde& in) { return Linv(in->clone()); }



template<typename Scalar> void zeroNyquist_sub(size_t iStart, size_t iStop, const vector3<int> S, Scalar* data)
{	THREAD_halfGspaceLoop( if(IS_NYQUIST) data[i] = Scalar(0.0); )
}
void zeroNyquist(RealKernel& K)
{	threadLaunch(zeroNyquist_sub<double>, K.gInfo.nG, K.gInfo.S, K.data);
}
void zeroNyquist(ScalarFieldTilde& Gptr)
{	const GridInfo& gInfo = Gptr->gInfo;
	threadLaunch(zeroNyquist_sub<complex>, gInfo.nG, gInfo.S, Gptr->data());
}
void zeroNyquist(ScalarField& Rptr) { ScalarFieldTilde Rtilde=J(Rptr); zeroNyquist(Rtilde); Rptr = I((ScalarFieldTilde&&)Rtilde); }

//------------------------------ Nonlinear Unary operators ------------------------------

void exp_sub(size_t i, double* X, double prefac) { X[i] = exp(prefac*X[i]); }
#ifdef GPU_ENABLED
void exp_gpu(int N, double* X, double prefac); //in operators.cu
#endif
ScalarField exp(ScalarField&& X)
{
	#ifdef GPU_ENABLED
	exp_gpu(X->nElem, X->dataGpu(false), X->scale);
	#else
	threadedLoop(exp_sub, X->nElem, X->data(false), X->scale);
	#endif
	X->scale = 1.0;
	return X;
}
ScalarField exp(const ScalarField& X) { return exp(X->clone()); }

void log_sub(size_t i, double* X, double prefac) { X[i] = log(prefac*X[i]); }
#ifdef GPU_ENABLED
void log_gpu(int N, double* X, double prefac); //in operators.cu
#endif
ScalarField log(ScalarField&& X)
{
	#ifdef GPU_ENABLED
	log_gpu(X->nElem, X->dataGpu(false), X->scale);
	#else
	threadedLoop(log_sub, X->nElem, X->data(false), X->scale);
	#endif
	X->scale = 1.0;
	return X;
}
ScalarField log(const ScalarField& X) { return log(X->clone()); }


void sqrt_sub(size_t i, double* X, double prefac) { X[i] = sqrt(prefac*X[i]); }
#ifdef GPU_ENABLED
void sqrt_gpu(int N, double* X, double prefac); //in operators.cu
#endif
ScalarField sqrt(ScalarField&& X)
{
	#ifdef GPU_ENABLED
	sqrt_gpu(X->nElem, X->dataGpu(false), X->scale);
	#else
	threadedLoop(sqrt_sub, X->nElem, X->data(false), X->scale);
	#endif
	X->scale = 1.0;
	return X;
}
ScalarField sqrt(const ScalarField& X) { return sqrt(X->clone()); }


void inv_sub(size_t i, double* X, double prefac) { X[i] = prefac/X[i]; }
#ifdef GPU_ENABLED
void inv_gpu(int N, double* X, double prefac); //in operators.cu
#endif
ScalarField inv(ScalarField&& X)
{
	#ifdef GPU_ENABLED
	inv_gpu(X->nElem, X->dataGpu(false), 1.0/X->scale);
	#else
	threadedLoop(inv_sub, X->nElem, X->data(false), 1.0/X->scale);
	#endif
	X->scale = 1.0;
	return X;
}
ScalarField inv(const ScalarField& X) { return inv(X->clone()); }

void pow_sub(size_t i, double* X, double scale, double alpha) { X[i] = pow(scale*X[i],alpha); }
#ifdef GPU_ENABLED
void pow_gpu(int N, double* X, double scale, double alpha); //in operators.cu
#endif
ScalarField pow(ScalarField&& X, double alpha)
{
	#ifdef GPU_ENABLED
	pow_gpu(X->nElem, X->dataGpu(false), X->scale, alpha);
	#else
	threadedLoop(pow_sub, X->nElem, X->data(false), X->scale, alpha);
	#endif
	X->scale = 1.0;
	return X;
}
ScalarField pow(const ScalarField& X, double alpha) { return pow(X->clone(), alpha); }


//------------------------------ Multiplication operators------------------------------

ScalarField& operator*=(ScalarField& in, const ScalarField& other)
{	in->scale *= other->scale;
	callPref(eblas_dmul)(in->nElem, other->dataPref(false), 1, in->dataPref(false), 1);
	return in;
}

complexScalarField& operator*=(complexScalarField& in, const ScalarField& other)
{	in->scale *= other->scale;
	callPref(eblas_zmuld)(in->nElem, other->dataPref(false), 1, in->dataPref(false), 1);
	return in;
}
complexScalarField operator*(const complexScalarField& inC, const ScalarField& inR) { complexScalarField out(inC->clone()); return out *= inR; }
complexScalarField operator*(const ScalarField& inR, const complexScalarField& inC) { complexScalarField out(inC->clone()); return out *= inR; }
complexScalarField operator*(complexScalarField&& inC, const ScalarField& inR) { return inC *= inR; }
complexScalarField operator*(const ScalarField& inR, complexScalarField&& inC) { return inC *= inR; }

ScalarFieldTilde& operator*=(ScalarFieldTilde& inG, const RealKernel& inR)
{	callPref(eblas_zmuld)(inG->nElem, inR.dataPref, 1, inG->dataPref(false), 1);
	return inG;
}
ScalarFieldTilde operator*(const RealKernel& inR, const ScalarFieldTilde& inG) { ScalarFieldTilde out(inG->clone()); return out *= inR; }
ScalarFieldTilde operator*(const ScalarFieldTilde& inG, const RealKernel& inR) { ScalarFieldTilde out(inG->clone()); return out *= inR; }
ScalarFieldTilde operator*(const RealKernel& inR, ScalarFieldTilde&& inG) { return inG *= inR; }
ScalarFieldTilde operator*(ScalarFieldTilde&& inG, const RealKernel& inR) { return inG *= inR; }


//------------------------------ Linear combine operators ------------------------------

void axpy(double alpha, const ScalarField& X, ScalarField& Y)
{	if(X)
	{	if(Y)
		{	if(Y->scale == 0.0) { Y = X * alpha; }
			else callPref(eblas_daxpy)(X->nElem, alpha*X->scale/Y->scale, X->dataPref(false), 1, Y->dataPref(false), 1);
		}
		else Y = X * alpha;
	}
	//if X is null, nothing needs to be done, Y remains unchanged
}
ScalarField& operator+=(ScalarField& in, double scalar)
{	FieldData dataScalar(in->gInfo, "ScalarField", 1, 1, false); *((double*)dataScalar.data()) = scalar;
	callPref(eblas_daxpy)(in->nElem, 1.0, (double*)dataScalar.dataPref(), 0, in->dataPref(), 1);
	return in;
}
ScalarField operator+(double scalar, const ScalarField& in) { ScalarField out(in->clone()); return out += scalar; }
ScalarField operator+(const ScalarField& in, double scalar) { ScalarField out(in->clone()); return out += scalar; }
ScalarField operator+(double scalar, ScalarField&& in) { return in += scalar; }
ScalarField operator+(ScalarField&& in, double scalar) { return in += scalar; }
ScalarField& operator-=(ScalarField& in, double scalar)
{	return (in += -scalar);
}
ScalarField operator-(double scalar, const ScalarField& in) { ScalarField out(in->clone()); return (out *= -1.0) += scalar; }
ScalarField operator-(const ScalarField& in, double scalar) { ScalarField out(in->clone()); return out -= scalar; }
ScalarField operator-(double scalar, ScalarField&& in) { return (in *= -1.0) += scalar; }
ScalarField operator-(ScalarField&& in, double scalar) { return in -= scalar; }


//------------------------------ Dot products and 2-norms ------------------------------

double dot(const ScalarField& X, const ScalarField& Y)
{	return X->scale * Y->scale * callPref(eblas_ddot)(X->nElem, X->dataPref(false), 1, Y->dataPref(false), 1);
}

double dot(const ScalarFieldTilde& X, const ScalarFieldTilde& Y)
{	int N = X->nElem;
	int S2 = X->gInfo.S[2]/2 + 1; //inner dimension
	int S01 = X->gInfo.S[0] * X->gInfo.S[1]; //number of inner dimension slices
	complex complexDot = callPref(eblas_zdotc)(N, X->dataPref(false), 1, Y->dataPref(false), 1);
	complex correction1 = callPref(eblas_zdotc)(S01, X->dataPref(false), S2, Y->dataPref(false), S2);
	complex correction2 = callPref(eblas_zdotc)(S01, X->dataPref(false)+S2-1, S2, Y->dataPref(false)+S2-1, S2);
	if(S2==1) correction2=complex(0,0); //because slices 1 and 2 are the same
	return X->scale * Y->scale * (2.0*complexDot - correction1 - correction2).real();
}

double nrm2(const ScalarField& X)
{	return fabs(X->scale) * callPref(eblas_dnrm2)(X->nElem, X->dataPref(false), 1);
}

double nrm2(const ScalarFieldTilde& X)
{	int N = X->nElem;
	int S2 = X->gInfo.S[2]/2 + 1; //inner dimension
	int S01 = X->gInfo.S[0] * X->gInfo.S[1]; //number of inner dimension slices
	double complexNorm = callPref(eblas_dznrm2)(N, X->dataPref(false), 1);
	double correction1 = callPref(eblas_dznrm2)(S01, X->dataPref(false), S2);
	double correction2 = callPref(eblas_dznrm2)(S01, X->dataPref(false)+S2-1, S2);
	if(S2==1) correction2=0; //because slices 1 and 2 are the same
	return fabs(X->scale) * sqrt(2*pow(complexNorm,2) - pow(correction1,2) - pow(correction2,2));
}

double sum(const ScalarField& X)
{	FieldData dataScale(X->gInfo, "ScalarField", 1, 1, false); *((double*)dataScale.data()) = X->scale;
	return callPref(eblas_ddot)(X->nElem, X->dataPref(false), 1, (double*)dataScale.dataPref(false), 0);
}

double sum(const ScalarFieldTilde& X)
{	int N = X->nElem;
	int S2 = X->gInfo.S[2]/2 + 1; //inner dimension
	int S01 = X->gInfo.S[0] * X->gInfo.S[1]; //number of inner dimension slices
	FieldData dataOne(X->gInfo, "ScalarFieldTilde", 1, 2, false); *((complex*)dataOne.data()) = 1.0;
	complex complexSum = callPref(eblas_zdotc)(N, X->dataPref(false), 1, (complex*)dataOne.dataPref(false), 0);
	complex correction1 = callPref(eblas_zdotc)(S01, X->dataPref(false), S2, (complex*)dataOne.dataPref(false), 0);
	complex correction2 = callPref(eblas_zdotc)(S01, X->dataPref(false)+S2-1, S2, (complex*)dataOne.dataPref(false), 0);
	if(S2==1) correction2 = 0; //because slices 1 and 2 are the same
	return X->scale * (2.0*complexSum - correction1 - correction2).real();
}



double integral(const ScalarField& X)
{	return X->gInfo.dV * sum(X);
}
double integral(const ScalarFieldTilde& X)
{	double XdataZero;
	#ifdef GPU_ENABLED
	cudaMemcpy(&XdataZero, X->dataGpu(false), sizeof(double), cudaMemcpyDeviceToHost);
	#else
	XdataZero = X->data(false)[0].real();
	#endif
	return XdataZero * X->gInfo.detR * X->scale;
}
complex integral(const complexScalarField& X)
{	return X->gInfo.dV * sum(X);
}
complex integral(const complexScalarFieldTilde& X)
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
ScalarFieldTilde changeGrid(const ScalarFieldTilde& in, const GridInfo& gInfoNew)
{	static StopWatch watch("changeGrid"); watch.start();
	ScalarFieldTilde out; nullToZero(out, gInfoNew);
	assert(gInfoNew.R == in->gInfo.R);
	const vector3<int>& Sin = in->gInfo.S;
	const vector3<int>& Sout = gInfoNew.S;
	vector3<int> Smax; for(int k=0; k<3; k++) Smax[k] = std::max(Sin[k],Sout[k]);
	callPref(changeGrid)(Smax, Sin, Sout, in->dataPref(), out->dataPref());
	watch.stop();
	return out;
}

ScalarField changeGrid(const ScalarField& in, const GridInfo& gInfoNew)
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
complexScalarFieldTilde changeGrid(const complexScalarFieldTilde& in, const GridInfo& gInfoNew)
{	static StopWatch watch("changeGridFull"); watch.start();
	complexScalarFieldTilde out; nullToZero(out, gInfoNew);
	assert(gInfoNew.R == in->gInfo.R);
	const vector3<int>& Sin = in->gInfo.S;
	const vector3<int>& Sout = gInfoNew.S;
	vector3<int> Smax; for(int k=0; k<3; k++) Smax[k] = std::max(Sin[k],Sout[k]);
	callPref(changeGridFull)(Smax, Sin, Sout, in->dataPref(), out->dataPref());
	watch.stop();
	return out;
}

complexScalarField changeGrid(const complexScalarField& in, const GridInfo& gInfoNew)
{	return I(changeGrid(J(in), gInfoNew));
}


//------------------------------ Initialization utilities ------------------------------

void initRandom(ScalarField& X, double cap)
{	double* Xdata = X->data();
	for(int i=0; i<X->nElem; i++)
		Xdata[i] = Random::normal(0, 1, cap);
}

void initRandomFlat(ScalarField& X)
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
void initTranslation(ScalarFieldTilde& X, const vector3<>& r)
{	const GridInfo& gInfo = X->gInfo;
	threadLaunch(initTranslation_sub, gInfo.nG, gInfo.S, gInfo.G*r, X->data());
}


void gaussConvolve_sub(size_t iStart, size_t iStop, const vector3<int>& S, const matrix3<>& GGT, complex* data, double sigma)
{	THREAD_halfGspaceLoop( data[i] *= exp(-0.5*sigma*sigma*GGT.metric_length_squared(iG)); )
}
void gaussConvolve(const vector3<int>& S, const matrix3<>& GGT, complex* data, double sigma)
{	threadLaunch(gaussConvolve_sub, S[0]*S[1]*(1+S[2]/2), S, GGT, data, sigma);
}
#ifdef GPU_ENABLED
void gaussConvolve_gpu(const vector3<int>& S, const matrix3<>& GGT, complex* data, double sigma);
#endif
ScalarFieldTilde gaussConvolve(ScalarFieldTilde&& in, double sigma)
{	assert(in);
	callPref(gaussConvolve)(in->gInfo.S, in->gInfo.GGT, in->dataPref(false), sigma);
	return in;
}
ScalarFieldTilde gaussConvolve(const ScalarFieldTilde& in, double sigma)
{	ScalarFieldTilde out(in->clone());
	return gaussConvolve((ScalarFieldTilde&&)out, sigma);
}

//------------------------------ Debug utilities ------------------------------

void printStats(const ScalarField& X, const char* name, FILE* fp)
{	int N = X->nElem;
	double mean = sum(X)/N;
	double stdDev = sqrt(fabs(dot(X,X)/N - mean*mean));
	double minVal, maxVal; callPref(eblas_capMinMax)(N, X->dataPref(), minVal, maxVal);
	fprintf(fp, "vector %s\t= %.15le +/- %.15le  min: %le  max: %le\n", name, mean, stdDev, minVal, maxVal);
}

//------------------------------ From VectorField.h ------------------------------

inline void gradient_sub(size_t iStart, size_t iStop, const vector3<int> S,
	const matrix3<> G, const complex* Xtilde, vector3<complex*> gradTilde)
{	THREAD_halfGspaceLoop( gradient_calc(i, iG, IS_NYQUIST, G, Xtilde, gradTilde); )
}
#ifdef GPU_ENABLED //implemented in Operators.cu
void gradient_gpu(const vector3<int> S, const matrix3<> G, const complex* Xtilde, vector3<complex*> gradTilde);
#endif
VectorFieldTilde gradient(const ScalarFieldTilde& Xtilde)
{	const GridInfo& gInfo = Xtilde->gInfo;
	VectorFieldTilde gradTilde(gInfo, isGpuEnabled());
	#ifdef GPU_ENABLED
	gradient_gpu(gInfo.S, gInfo.G, Xtilde->dataGpu(), gradTilde.dataGpu());
	#else
	threadLaunch(gradient_sub, gInfo.nG, gInfo.S, gInfo.G, Xtilde->data(), gradTilde.data());
	#endif
	return gradTilde;
}
VectorField gradient(const ScalarField& X) { return I(gradient(J(X))); }


inline void divergence_sub(size_t iStart, size_t iStop, const vector3<int> S,
	const matrix3<> G, vector3<const complex*> Vtilde, complex* divTilde)
{	THREAD_halfGspaceLoop( divergence_calc(i, iG, IS_NYQUIST, G, Vtilde, divTilde); )
}
#ifdef GPU_ENABLED
void divergence_gpu(const vector3<int> S, const matrix3<> G, vector3<const complex*> Vtilde, complex* divTilde);
#endif
ScalarFieldTilde divergence(const VectorFieldTilde& Vtilde)
{	const GridInfo& gInfo = Vtilde[0]->gInfo;
	ScalarFieldTilde divTilde(ScalarFieldTildeData::alloc(gInfo, isGpuEnabled()));
	#ifdef GPU_ENABLED
	divergence_gpu(gInfo.S, gInfo.G, Vtilde.dataGpu(), divTilde->dataGpu());
	#else
	threadLaunch(divergence_sub, gInfo.nG, gInfo.S, gInfo.G, Vtilde.data(), divTilde->data());
	#endif
	return divTilde;
}
ScalarField divergence(const VectorField& V) { return I(divergence(J(V))); }


inline void tensorGradient_sub(size_t iStart, size_t iStop, const vector3<int> S,
	const matrix3<> G, const complex* Xtilde, tensor3<complex*> gradTilde)
{	THREAD_halfGspaceLoop( tensorGradient_calc(i, iG, IS_NYQUIST, G, Xtilde, gradTilde); )
}
#ifdef GPU_ENABLED //implemented in Operators.cu
void tensorGradient_gpu(const vector3<int> S, const matrix3<> G, const complex* Xtilde, tensor3<complex*> gradTilde);
#endif
TensorFieldTilde tensorGradient(const ScalarFieldTilde& Xtilde)
{	const GridInfo& gInfo = Xtilde->gInfo;
	TensorFieldTilde gradTilde(gInfo, isGpuEnabled());
	#ifdef GPU_ENABLED
	tensorGradient_gpu(gInfo.S, gInfo.G, Xtilde->dataGpu(), gradTilde.dataGpu());
	#else
	threadLaunch(tensorGradient_sub, gInfo.nG, gInfo.S, gInfo.G, Xtilde->data(), gradTilde.data());
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
ScalarFieldTilde tensorDivergence(const TensorFieldTilde& Vtilde)
{	const GridInfo& gInfo = Vtilde[0]->gInfo;
	ScalarFieldTilde divTilde(ScalarFieldTildeData::alloc(gInfo, isGpuEnabled()));
	#ifdef GPU_ENABLED
	tensorDivergence_gpu(gInfo.S, gInfo.G, Vtilde.dataGpu(), divTilde->dataGpu());
	#else
	threadLaunch(tensorDivergence_sub, gInfo.nG, gInfo.S, gInfo.G, Vtilde.data(), divTilde->data());
	#endif
	return divTilde;
}

