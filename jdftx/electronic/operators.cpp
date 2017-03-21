/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman
Copyright 1996-2003 Sohrab Ismail-Beigi

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

#include <electronic/operators.h>
#include <electronic/operators_internal.h>
#include <electronic/ColumnBundle.h>
#include <electronic/matrix.h>
#include <electronic/Basis.h>
#include <electronic/ElecInfo.h>
#include <electronic/IonInfo.h>
#include <core/Thread.h>
#include <core/BlasExtra.h>
#include <core/GpuUtil.h>
#include <core/GridInfo.h>
#include <core/LoopMacros.h>
#include <core/Operators.h>

void removePhase(size_t N, complex* data, double& meanPhase, double& sigmaPhase, double& rmsImagErr)
{	//Find mean phase
	double w0=0.0, r1=0.0, r2=0.0, i1=0.0, i2=0.0; //moments of normalized real and imaginary parts
	for(size_t i=0; i<N; i++)
	{	complex c = data[i]*data[i];
		double w = abs(c); //weight
		if(w > 1e-300)
		{	w0 += w;
			r1 += c.real();
			i1 += c.imag();
			r2 += c.real() * c.real() / w;
			i2 += c.imag() * c.imag() / w;
		}
	}
	double rMean=r1/w0, rSigma=sqrt(std::max(0.,r2/w0-pow(rMean,2)));
	double iMean=i1/w0, iSigma=sqrt(std::max(0.,i2/w0-pow(iMean,2)));
	meanPhase = 0.5*atan2(iMean, rMean);
	sigmaPhase = 0.5*hypot(iMean*rSigma, rMean*iSigma)/(pow(rMean,2) + pow(iMean,2));
	
	//Remove phase:
	complex phaseCompensate = cis(-meanPhase);
	double imSqSum=0.0, reSqSum=0.0;
	for(size_t i=0; i<N; i++)
	{	complex& c = data[i];
		c *= phaseCompensate;
		reSqSum += pow(c.real(),2);
		imSqSum += pow(c.imag(),2);
		c = c.real();
	}
	rmsImagErr = sqrt(imSqSum/(imSqSum+reSqSum));
}


//------------------------------ Spatial gradients of vectors ---------------------------------

//first cartesian derivative
void D_sub(size_t iStart, size_t iStop, const vector3<int> S, const complex* in, complex* out, vector3<> Ge)
{	THREAD_halfGspaceLoop( D_calc(i, iG, in, out, Ge); )
}
#ifdef GPU_ENABLED
void D_gpu(const vector3<int> S, const complex* in, complex* out, vector3<> Ge);
#endif
ScalarFieldTilde D(const ScalarFieldTilde& in, int iDir)
{	const GridInfo& gInfo = in->gInfo;
	ScalarFieldTilde out(ScalarFieldTildeData::alloc(gInfo, isGpuEnabled()));
	#ifdef GPU_ENABLED
	D_gpu(gInfo.S, in->dataGpu(), out->dataGpu(), gInfo.G.column(iDir));
	#else
	threadLaunch(D_sub, gInfo.nG, gInfo.S, in->data(), out->data(), gInfo.G.column(iDir));
	#endif
	return out;
}

//second cartesian derivative
void DD_sub(size_t iStart, size_t iStop, const vector3<int> S, const complex* in, complex* out, vector3<> Ge1, vector3<> Ge2)
{	THREAD_halfGspaceLoop( DD_calc(i, iG, in, out, Ge1, Ge2); )
}
#ifdef GPU_ENABLED
void DD_gpu(const vector3<int> S, const complex* in, complex* out, vector3<> g, vector3<> q);
#endif
ScalarFieldTilde DD(const ScalarFieldTilde& in, int iDir, int jDir)
{	const GridInfo& gInfo = in->gInfo;
	ScalarFieldTilde out(ScalarFieldTildeData::alloc(gInfo, isGpuEnabled()));
	#ifdef GPU_ENABLED
	DD_gpu(gInfo.S, in->dataGpu(), out->dataGpu(), gInfo.G.column(iDir), gInfo.G.column(jDir));
	#else
	threadLaunch(DD_sub, gInfo.nG, gInfo.S, in->data(), out->data(), gInfo.G.column(iDir), gInfo.G.column(jDir));
	#endif
	return out;
}


template<int l> void lGradient_sub(size_t iStart, size_t iStop, const vector3<int>& S, const complex* in, const array<complex*, 2*l+1>& out, const matrix3<>& G)
{	THREAD_halfGspaceLoop( lGradient_calc<l>(i, iG, IS_NYQUIST, in, out, G); )
}
template<int l> void lGradient(const vector3<int>& S, const complex* in, array<complex*, 2*l+1> out, const matrix3<>& G)
{	threadLaunch(lGradient_sub<l>, S[0]*S[1]*(S[2]/2+1), S, in, out, G);
}
void lGradient(const vector3<int>& S, const complex* in, std::vector<complex*> out, int l, const matrix3<>& G)
{	SwitchTemplate_l(l, lGradient, (S, in, out, G))
}
#ifdef GPU_ENABLED
void lGradient_gpu(const vector3<int>& S, const complex* in, std::vector<complex*> out, int l, const matrix3<>& G);
#endif

ScalarFieldTildeArray lGradient(const ScalarFieldTilde& in, int l)
{	ScalarFieldTildeArray out; nullToZero(out, in->gInfo, 2*l+1);
	callPref(lGradient)(in->gInfo.S, in->dataPref(), dataPref(out), l, in->gInfo.G);
	return out;
}

template<int l> void lDivergence_sub(size_t iStart, size_t iStop, const vector3<int>& S, const array<const complex*,2*l+1>& in, complex* out, const matrix3<>& G)
{	THREAD_halfGspaceLoop( lDivergence_calc<l>(i, iG, IS_NYQUIST, in, out, G); )
}
template<int l> void lDivergence(const vector3<int>& S, array<const complex*,2*l+1> in, complex* out, const matrix3<>& G)
{	threadLaunch(lDivergence_sub<l>, S[0]*S[1]*(S[2]/2+1), S, in, out, G);
}
void lDivergence(const vector3<int>& S, const std::vector<const complex*>& in, complex* out, int l, const matrix3<>& G)
{	SwitchTemplate_l(l, lDivergence, (S, in, out, G))
}
#ifdef GPU_ENABLED
void lDivergence_gpu(const vector3<int>& S, const std::vector<const complex*>& in, complex* out, int l, const matrix3<>& G);
#endif

ScalarFieldTilde lDivergence(const ScalarFieldTildeArray& in, int l)
{	assert(int(in.size()) == 2*l+1);
	ScalarFieldTilde out; nullToZero(out, in[0]->gInfo);
	callPref(lDivergence)(in[0]->gInfo.S, constDataPref(in), out->dataPref(), l, in[0]->gInfo.G);
	return out;
}



void multiplyBlochPhase_sub(size_t iStart, size_t iStop,
	const vector3<int>& S, const vector3<>& invS, complex* v, const vector3<>& k)
{	THREAD_rLoop( v[i] *= blochPhase_calc(iv, invS, k); )
}
#ifdef GPU_ENABLED
void multiplyBlochPhase_gpu(const vector3<int>& S, const vector3<>& invS, complex* v, const vector3<>& k);
#endif
void multiplyBlochPhase(complexScalarField& v, const vector3<>& k)
{	const GridInfo& gInfo = v->gInfo;
	vector3<> invS(1./gInfo.S[0], 1./gInfo.S[1], 1./gInfo.S[2]);
	#ifdef GPU_ENABLED
	multiplyBlochPhase_gpu(gInfo.S, invS, v->dataGpu(), k);
	#else
	threadLaunch(multiplyBlochPhase_sub, gInfo.nr, gInfo.S, invS, v->data(), k);
	#endif
}


void radialFunction_sub(size_t iStart, size_t iStop, const vector3<int> S, const matrix3<>& GGT,
	complex* F, const RadialFunctionG& f, vector3<> r0 )
{	THREAD_halfGspaceLoop( F[i] = radialFunction_calc(iG, GGT, f, r0); )
}
#ifdef GPU_ENABLED
void radialFunction_gpu(const vector3<int> S, const matrix3<>& GGT,
	complex* F, const RadialFunctionG& f, vector3<> r0);
#endif
ScalarFieldTilde radialFunctionG(const GridInfo& gInfo, const RadialFunctionG& f, vector3<> r0)
{	
	ScalarFieldTilde F(ScalarFieldTildeData::alloc(gInfo,isGpuEnabled()));
	#ifdef GPU_ENABLED
	radialFunction_gpu(gInfo.S, gInfo.GGT, F->dataGpu(), f, r0);
	#else
	threadLaunch(radialFunction_sub, gInfo.nG, gInfo.S, gInfo.GGT, F->data(), f, r0);
	#endif
	return F;
}

ScalarField radialFunction(const GridInfo& gInfo, const RadialFunctionG& f, vector3<> r0)
{	
	ScalarFieldTilde F = radialFunctionG(gInfo, f, r0);
	return (1.0/gInfo.detR) * I(F);
}

void radialFunctionG(const RadialFunctionG& f, RealKernel& Kernel)
{	
	ScalarFieldTilde F = radialFunctionG(Kernel.gInfo, f, vector3<>(0,0,0));
	const complex* FData = F->data(); //put F into Kernel
	for(int i=0; i<Kernel.gInfo.nG; i++)
		Kernel.data[i] = FData[i].real();
	Kernel.set();
}


void radialFunctionMultiply_sub(size_t iStart, size_t iStop, const vector3<int> S, const matrix3<>& GGT,
	complex* in, const RadialFunctionG& f)
{	THREAD_halfGspaceLoop( in[i] *= f(sqrt(GGT.metric_length_squared(iG))); )
}
#ifdef GPU_ENABLED
void radialFunctionMultiply_gpu(const vector3<int> S, const matrix3<>& GGT, complex* in, const RadialFunctionG& f);
#endif

ScalarFieldTilde operator*(const RadialFunctionG& f, ScalarFieldTilde&& in)
{	const GridInfo& gInfo = in->gInfo;
	#ifdef GPU_ENABLED
	radialFunctionMultiply_gpu(gInfo.S, gInfo.GGT, in->dataGpu(), f);
	#else
	threadLaunch(radialFunctionMultiply_sub, gInfo.nG, gInfo.S, gInfo.GGT, in->data(), f);
	#endif
	return in;
}

ScalarFieldTilde operator*(const RadialFunctionG& f, const ScalarFieldTilde& in)
{	ScalarFieldTilde out(in->clone()); //destructible copy
	return f * ((ScalarFieldTilde&&)out);
}

VectorFieldTilde operator*(const RadialFunctionG& f, VectorFieldTilde&& in)
{	for(int k=0; k<3; k++) in[k] = f * (ScalarFieldTilde&&)in[k];
	return in;
}

VectorFieldTilde operator*(const RadialFunctionG& f, const VectorFieldTilde& in)
{	VectorFieldTilde out;
	for(int k=0; k<3; k++) out[k] = f * in[k];
	return out;
}


//------------------------------ ColumnBundle operators ---------------------------------

void Idag_DiagV_I_sub(int colStart, int colEnd, const ColumnBundle* C, const ScalarFieldArray* V, ColumnBundle* VC)
{	const ScalarField& Vs = V->at(V->size()==1 ? 0 : C->qnum->index());
	int nSpinor = VC->spinorLength();
	for(int col=colStart; col<colEnd; col++)
		for(int s=0; s<nSpinor; s++)
			VC->setColumn(col,s, Idag(Vs * I(C->getColumn(col,s))));
}

//Noncollinear version of above (with the preprocessing of complex off-diagonal potentials done in calling function)
void Idag_DiagVmat_I_sub(int colStart, int colEnd, const ColumnBundle* C, const ScalarField* Vup, const ScalarField* Vdn,
	const complexScalarField* VupDn, const complexScalarField* VdnUp, ColumnBundle* VC)
{	for(int col=colStart; col<colEnd; col++)
	{	complexScalarField ICup = I(C->getColumn(col,0));
		complexScalarField ICdn = I(C->getColumn(col,1));
		VC->setColumn(col,0, Idag((*Vup)*ICup + (*VupDn)*ICdn));
		VC->setColumn(col,1, Idag((*Vdn)*ICdn + (*VdnUp)*ICup));
	}
	
}

ColumnBundle Idag_DiagV_I(const ColumnBundle& C, const ScalarFieldArray& V)
{	static StopWatch watch("Idag_DiagV_I"); watch.start();
	ColumnBundle VC = C.similar();
	//Convert V to wfns grid if necessary:
	const GridInfo& gInfoWfns = *(C.basis->gInfo);
	ScalarFieldArray Vtmp;
	if(&(V[0]->gInfo) != &gInfoWfns)
		for(const ScalarField& Vs: V)
			Vtmp.push_back(Jdag(changeGrid(Idag(Vs), gInfoWfns), true));
	const ScalarFieldArray& Vwfns = Vtmp.size() ? Vtmp : V;
	assert(Vwfns.size()==1 || Vwfns.size()==2 || Vwfns.size()==4);
	if(Vwfns.size()==2) assert(!C.isSpinor());
	if(Vwfns.size()==1 || Vwfns.size()==2)
	{	threadLaunch(isGpuEnabled()?1:0, Idag_DiagV_I_sub, C.nCols(), &C, &Vwfns, &VC);
	}
	else //Vwfns.size()==4
	{	assert(C.isSpinor());
		complexScalarField VupDn = 0.5*Complex(Vwfns[2], Vwfns[3]);
		complexScalarField VdnUp = conj(VupDn);
		threadLaunch(isGpuEnabled()?1:0, Idag_DiagVmat_I_sub, C.nCols(), &C, &Vwfns[0], &Vwfns[1], &VupDn, &VdnUp, &VC);
	}
	watch.stop();
	return VC;
}


//Laplacian of a column bundle
#ifdef GPU_ENABLED
void reducedL_gpu(int nbasis, int ncols, const complex* Y, complex* LY,
	const matrix3<> GGT, const vector3<int>* iGarr, const vector3<> k, double detR);
#endif
ColumnBundle L(const ColumnBundle &Y)
{	ColumnBundle LY = Y.similar();
	assert(Y.basis);
	const Basis& basis = *(Y.basis);
	const matrix3<>& GGT = basis.gInfo->GGT;
	int nSpinors = Y.spinorLength();
	#ifdef GPU_ENABLED
	reducedL_gpu(basis.nbasis, Y.nCols()*nSpinors, Y.dataGpu(), LY.dataGpu(), GGT, basis.iGarrGpu, Y.qnum->k, basis.gInfo->detR);
	#else
	threadedLoop(reducedL_calc, basis.nbasis,
		basis.nbasis, Y.nCols()*nSpinors, Y.data(), LY.data(), GGT, basis.iGarr, Y.qnum->k, basis.gInfo->detR);
	#endif
	return LY;
}

//Inverse-Laplacian of a column bundle
#ifdef GPU_ENABLED
void reducedLinv_gpu(int nbasis, int ncols, const complex* Y, complex* LinvY,
	const matrix3<> GGT, const vector3<int>* iGarr, const vector3<> k, double detR);
#endif
ColumnBundle Linv(const ColumnBundle &Y)
{	ColumnBundle LinvY = Y.similar();
	assert(Y.basis);
	const Basis& basis = *(Y.basis);
	const matrix3<>& GGT = basis.gInfo->GGT;
	int nSpinors = Y.spinorLength();
	#ifdef GPU_ENABLED
	reducedLinv_gpu(basis.nbasis, Y.nCols()*nSpinors, Y.dataGpu(), LinvY.dataGpu(), GGT, basis.iGarrGpu, Y.qnum->k, basis.gInfo->detR);
	#else
	threadedLoop(reducedLinv_calc, basis.nbasis,
		basis.nbasis, Y.nCols()*nSpinors, Y.data(), LinvY.data(), GGT, basis.iGarr, Y.qnum->k, basis.gInfo->detR);
	#endif
	return LinvY;
}


// Overlap operator (scale by unit cell volume in PW basis)
ColumnBundle O(const ColumnBundle &Y, std::vector<matrix>* VdagY)
{	ColumnBundle OY = Y * Y.basis->gInfo->detR; //basic planewave overlap
	Y.basis->iInfo->augmentOverlap(Y, OY, VdagY); //pseudopotential augmentation
	return OY;
}

//Compute cartesian gradient of column bundle in direction #iDir
#ifdef GPU_ENABLED
void reducedD_gpu(int nbasis, int ncols, const complex* Ydata, complex* DYdata,
	const vector3<int>* iGarr, double kdotGe, const vector3<> Ge);
#endif
ColumnBundle D(const ColumnBundle &Y, int iDir)
{	assert(Y.basis);
	const Basis& basis = *(Y.basis);
	ColumnBundle DY = Y.similar();
	int nSpinors = Y.spinorLength();
	const vector3<> Ge = basis.gInfo->G.column(iDir);
	double kdotGe = dot(Y.qnum->k, Ge);
	#ifdef GPU_ENABLED
	reducedD_gpu(basis.nbasis, Y.nCols()*nSpinors, Y.dataGpu(), DY.dataGpu(), basis.iGarrGpu, kdotGe, Ge);
	#else
	threadedLoop(reducedD_calc, basis.nbasis,
		basis.nbasis, Y.nCols()*nSpinors, Y.data(), DY.data(), basis.iGarr, kdotGe, Ge);
	#endif
	return DY;
}


//Compute cartesian gradient of column bundle in direction #iDir
#ifdef GPU_ENABLED
void reducedDD_gpu(int nbasis, int ncols, const complex* Ydata, complex* DDYdata,
	const vector3<int>* iGarr, double kdotGe1, double kdotGe2, const vector3<> Ge1, const vector3<> Ge2);
#endif
ColumnBundle DD(const ColumnBundle &Y, int iDir, int jDir)
{	assert(Y.basis);
	const Basis& basis = *(Y.basis);
	ColumnBundle DDY = Y.similar();
	int nSpinors = Y.spinorLength();
	const vector3<> Ge1 = basis.gInfo->G.column(iDir);
	const vector3<> Ge2 = basis.gInfo->G.column(jDir);
	double kdotGe1 = dot(Y.qnum->k, Ge1);
	double kdotGe2 = dot(Y.qnum->k, Ge2);
	#ifdef GPU_ENABLED
	reducedDD_gpu(basis.nbasis, Y.nCols()*nSpinors, Y.dataGpu(), DDY.dataGpu(), basis.iGarrGpu, kdotGe1, kdotGe2, Ge1, Ge2);
	#else
	threadedLoop(reducedDD_calc, basis.nbasis,
		basis.nbasis, Y.nCols()*nSpinors, Y.data(), DDY.data(), basis.iGarr, kdotGe1, kdotGe2, Ge1, Ge2);
	#endif
	return DDY;
}


// Multiply each column by f(0.5*|k+G|^2/KErollover)
// with f(x) = (1+x+x^2+x^3+...+x^8)/(1+x+x^2+...+x^9) = (1-x^N)/(1-x^(N+1))
void precond_inv_kinetic(int nbasis, int ncols, complex* Ydata,
	double KErollover, const matrix3<>& GGT, const vector3<int>* iGarr, const vector3<> k, double invdetR)
{
	threadedLoop(precond_inv_kinetic_calc, nbasis, nbasis, ncols, Ydata, KErollover, GGT, iGarr, k, invdetR);
}
#ifdef GPU_ENABLED
void precond_inv_kinetic_gpu(int nbasis, int ncols, complex* Ydata,
	double KErollover, const matrix3<> GGT, const vector3<int>* iGarr, const vector3<> k, double invdetR);
#endif
void precond_inv_kinetic(ColumnBundle &Y, double KErollover)
{	assert(Y.basis);
	const Basis& basis = *Y.basis;
	const matrix3<>& GGT = basis.gInfo->GGT;
	int  nSpinors = Y.spinorLength();
	callPref(precond_inv_kinetic)(basis.nbasis, Y.nCols()*nSpinors, Y.dataPref(),
		KErollover, GGT, basis.iGarrPref, Y.qnum->k, 1./basis.gInfo->detR);
}

diagMatrix diagDot(const ColumnBundle& X, const ColumnBundle& Y)
{	assert(X.nCols()==Y.nCols());
	assert(X.basis==Y.basis);
	diagMatrix ret(X.nCols());
	const complex* Xdata = X.dataPref();
	const complex* Ydata = Y.dataPref();
	for(size_t b=0; b<ret.size(); b++)
		ret[b] = callPref(eblas_zdotc)(X.colLength(), Xdata+X.index(b,0),1, Ydata+Y.index(b,0),1).real();
	return ret;
}

void precond_inv_kinetic_band(int nbasis, int ncols, complex* Ydata, const double* KEref,
	const matrix3<>& GGT, const vector3<int>* iGarr, const vector3<>& k)
{	threadedLoop(precond_inv_kinetic_band_calc, nbasis, nbasis, ncols, Ydata, KEref, GGT, iGarr, k);
}
#ifdef GPU_ENABLED
void precond_inv_kinetic_band_gpu(int nbasis, int ncols, complex* Ydata, const double* KEref,
	const matrix3<>& GGT, const vector3<int>* iGarr, const vector3<>& k);
#endif
void precond_inv_kinetic_band(ColumnBundle& Y, const diagMatrix& KErefIn)
{	assert(Y.basis);
	const Basis& basis = *Y.basis;
	assert(Y.nCols()==KErefIn.nCols());
	int nSpinors = Y.spinorLength();
	//Adapt KEref array for spinors:
	diagMatrix KEtmp;
	if(nSpinors > 1)
	{	KEtmp.reserve(Y.nCols()*nSpinors);
		for(const double& KE: KErefIn)
			KEtmp.insert(KEtmp.end(), nSpinors, KE);
	}
	const diagMatrix& KEref = KEtmp.size() ? KEtmp : KErefIn;
	#ifdef GPU_ENABLED
	matrix KErefCopy(KEref.nCols(), 1); //used just a a dummy ManagedMemory object
	eblas_copy((double*)KErefCopy.data(), KEref.data(), KEref.nCols());
	const double* KErefData = (const double*)KErefCopy.dataGpu();
	#else
	const double* KErefData = KEref.data();
	#endif
	callPref(precond_inv_kinetic_band)(basis.nbasis, Y.nCols()*nSpinors, Y.dataPref(), KErefData,
		basis.gInfo->GGT, basis.iGarrPref, Y.qnum->k);
}


#ifdef GPU_ENABLED
void translate_gpu(int nbasis, int ncols, complex* Y, const vector3<int>* iGarr, const vector3<>& k, const vector3<>& dr);
#endif
ColumnBundle translate(ColumnBundle&& Y, vector3<> dr)
{	assert(Y.basis);
	const Basis& basis = *Y.basis;
	int nSpinors = Y.spinorLength();
	#ifdef GPU_ENABLED
	translate_gpu(basis.nbasis, Y.nCols()*nSpinors, Y.dataGpu(), basis.iGarrGpu, Y.qnum->k, dr);
	#else
	threadedLoop(translate_calc, basis.nbasis, basis.nbasis, Y.nCols()*nSpinors, Y.data(), basis.iGarr, Y.qnum->k, dr);
	#endif
	return Y;
}
ColumnBundle translate(const ColumnBundle& Y, vector3<> dr)
{	return translate((ColumnBundle&&)ColumnBundle(Y), dr); //call above function on a destructible copy
}

void translateColumns(int nbasis, int ncols, complex* Y, const vector3<int>* iGarr, const vector3<>& k, const vector3<>* dr)
{	threadedLoop(translateColumns_calc, nbasis, nbasis, ncols, Y, iGarr, k, dr);
}
#ifdef GPU_ENABLED
void translateColumns_gpu(int nbasis, int ncols, complex* Y, const vector3<int>* iGarr, const vector3<>& k, const vector3<>* dr);
#endif
void translateColumns(ColumnBundle& Y, const vector3<>* dr)
{	assert(Y.basis);
	const Basis& basis = *Y.basis;
	int nSpinors = Y.spinorLength();
	callPref(translateColumns)(basis.nbasis, Y.nCols()*nSpinors, Y.dataPref(), basis.iGarrPref, Y.qnum->k, dr);
}


ColumnBundle switchBasis(const ColumnBundle& in, const Basis& basisOut)
{	if(in.basis == &basisOut) return in; //no basis change required
	int nSpinors = in.spinorLength();
	ColumnBundle out(in.nCols(), basisOut.nbasis*nSpinors, &basisOut, 0, isGpuEnabled());
	for(int b=0; b<in.nCols(); b++)
		for(int s=0; s<nSpinors; s++)
			out.setColumn(b,s, in.getColumn(b,s)); //convert using the full G-space as an intermediate
	return out;
}

//------------------------------ ColumnBundle reductions ---------------------------------

// Returns trace(F*X^Y)
complex traceinner(const diagMatrix &F, const ColumnBundle &X, const ColumnBundle &Y)
{	assert(X.colLength()==Y.colLength());
	assert(X.nCols()==Y.nCols());
	assert(X.nCols()==F.nRows());
	complex result = 0.0;
	for (int i=0; i < X.nCols(); i++)
		result += F[i] * callPref(eblas_zdotc)(X.colLength(), X.dataPref()+X.index(i,0), 1, Y.dataPref()+Y.index(i,0), 1);
	return result;
}

// Compute the density from a subset of columns of a ColumnBundle
void diagouterI_sub(int iThread, int nThreads, const diagMatrix *F, const ColumnBundle *X, std::vector<ScalarFieldArray>* nSub)
{
	//Determine column range:
	int colStart = (( iThread ) * X->nCols())/nThreads;
	int colStop  = ((iThread+1) * X->nCols())/nThreads;
	
	ScalarFieldArray& nLocal = (*nSub)[iThread];
	nullToZero(nLocal, *(X->basis->gInfo)); //sets to zero
	int nDensities = nLocal.size();
	if(nDensities==1) //Note that nDensities==2 below will also enter this branch sinc eonly one component is non-zero
	{	int nSpinor = X->spinorLength();
		for(int i=colStart; i<colStop; i++)
			for(int s=0; s<nSpinor; s++)
				callPref(eblas_accumNorm)(X->basis->gInfo->nr, (*F)[i], I(X->getColumn(i,s))->dataPref(), nLocal[0]->dataPref());
	}
	else //nDensities==4 (ensured by assertions in launching function below)
	{	for(int i=colStart; i<colStop; i++)
		{	complexScalarField psiUp = I(X->getColumn(i,0));
			complexScalarField psiDn = I(X->getColumn(i,1));
			callPref(eblas_accumNorm)(X->basis->gInfo->nr, (*F)[i], psiUp->dataPref(), nLocal[0]->dataPref()); //UpUp
			callPref(eblas_accumNorm)(X->basis->gInfo->nr, (*F)[i], psiDn->dataPref(), nLocal[1]->dataPref()); //DnDn
			callPref(eblas_accumProd)(X->basis->gInfo->nr, (*F)[i], psiUp->dataPref(), psiDn->dataPref(), nLocal[2]->dataPref(), nLocal[3]->dataPref()); //Re and Im parts of UpDn
		}
	}
}

// Collect all contributions from nSub into the first entry
void diagouterI_collect(size_t iStart, size_t iStop, std::vector<ScalarFieldArray>* nSub)
{	assert(!isGpuEnabled()); // this is needed and should be called only in CPU mode
	for(size_t s=0; s<(*nSub)[0].size(); s++)
	{	//Get the data pointers for each piece in nSub:
		int nThreads = nSub->size();
		std::vector<double*> nSubData(nThreads);
		for(int j=0; j<nThreads; j++) nSubData[j] = (*nSub)[j][s]->data();

		//Accumulate pointwise into the first piece:
		for(size_t i=iStart; i<iStop; i++)
			for(int j=1; j<nThreads; j++)
				nSubData[0][i] += nSubData[j][i];
	}
}

// Returns diag((I*X)*F*(I*X)^) where X^ is the hermetian adjoint of X.
ScalarFieldArray diagouterI(const diagMatrix &F,const ColumnBundle &X,  int nDensities, const GridInfo* gInfoOut)
{	static StopWatch watch("diagouterI"); watch.start();
	//Check sizes:
	assert(F.nRows()==X.nCols());
	assert(nDensities==1 || nDensities==2 || nDensities==4);
	if(nDensities==2) assert(!X.isSpinor());
	if(nDensities==4) assert(X.isSpinor());
	
	//Collect the contributions for different sets of columns in separate scalar fields (one per thread):
	int nThreads = isGpuEnabled() ? 1: nProcsAvailable;
	std::vector<ScalarFieldArray> nSub(nThreads, ScalarFieldArray(nDensities==2 ? 1 : nDensities)); //collinear spin-polarized will have only one non-zero output channel
	threadLaunch(nThreads, diagouterI_sub, 0, &F, &X, &nSub);

	//If more than one thread, accumulate all vectors in nSub into the first:
	if(nThreads>1) threadLaunch(diagouterI_collect, X.basis->gInfo->nr, &nSub);
	watch.stop();
	
	//Change grid if necessary:
	if(gInfoOut && (X.basis->gInfo!=gInfoOut))
		for(ScalarField& nSub0s: nSub[0])
			nSub0s = changeGrid(nSub0s, *gInfoOut);
	
	//Correct the location of the single non-zero channel of collinear spin-polarized densities:
	if(nDensities==2)
	{	nSub[0].resize(2);
		if(X.qnum->index()==1) std::swap(nSub[0][0], nSub[0][1]);
	}
	return nSub[0]; //rest cleaned up destructor
}
