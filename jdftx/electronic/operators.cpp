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
	register double w0=0.0, r1=0.0, r2=0.0, i1=0.0, i2=0.0; //moments of normalized real and imaginary parts
	for(size_t i=0; i<N; i++)
	{	register complex c = data[i]*data[i];
		register double w = abs(c); //weight
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
DataGptr D(const DataGptr& in, int iDir)
{	const GridInfo& gInfo = in->gInfo;
	DataGptr out(DataG::alloc(gInfo, isGpuEnabled()));
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
DataGptr DD(const DataGptr& in, int iDir, int jDir)
{	const GridInfo& gInfo = in->gInfo;
	DataGptr out(DataG::alloc(gInfo, isGpuEnabled()));
	#ifdef GPU_ENABLED
	DD_gpu(gInfo.S, in->dataGpu(), out->dataGpu(), gInfo.G.column(iDir), gInfo.G.column(jDir));
	#else
	threadLaunch(DD_sub, gInfo.nG, gInfo.S, in->data(), out->data(), gInfo.G.column(iDir), gInfo.G.column(jDir));
	#endif
	return out;
}

void multiplyBlochPhase_sub(size_t iStart, size_t iStop,
	const vector3<int>& S, const vector3<>& invS, complex* v, const vector3<>& k)
{	THREAD_rLoop( v[i] *= blochPhase_calc(iv, invS, k); )
}
#ifdef GPU_ENABLED
void multiplyBlochPhase_gpu(const vector3<int>& S, const vector3<>& invS, complex* v, const vector3<>& k);
#endif
void multiplyBlochPhase(complexDataRptr& v, const vector3<>& k)
{	const GridInfo& gInfo = v->gInfo;
	vector3<> invS(1./gInfo.S[0], 1./gInfo.S[1], 1./gInfo.S[2]);
	#ifdef GPU_ENABLED
	multiplyBlochPhase_gpu(gInfo.S, invS, v->dataGpu(), k);
	#else
	threadLaunch(multiplyBlochPhase_sub, gInfo.nr, gInfo.S, invS, v->data(), k);
	#endif
}


//point group scatter
template<typename scalar> void pointGroupScatter_sub(size_t iStart, size_t iStop,
	const vector3<int>& S, const scalar* in, scalar* out, const matrix3<int>& mMesh)
{	THREAD_rLoop( pointGroupScatter_calc(i, iv, S, in, out, mMesh); )
}
#ifdef GPU_ENABLED
void pointGroupScatter_gpu(const vector3<int>& S, const double* in, double* out, const matrix3<int>& mMesh);
void pointGroupScatter_gpu(const vector3<int>& S, const complex* in, complex* out, const matrix3<int>& mMesh);
#endif
template<typename T> std::shared_ptr<T> pointGroupScatter(const std::shared_ptr<T>& in, const matrix3<int>& mMesh)
{	if(mMesh == matrix3<int>(1,1,1)) return in; //shortcut for identity
	const GridInfo& gInfo = in->gInfo;
	std::shared_ptr<T> out(T::alloc(gInfo, isGpuEnabled()));
	#ifdef GPU_ENABLED
	pointGroupScatter_gpu(gInfo.S, in->dataGpu(), out->dataGpu(), mMesh);
	#else
	threadLaunch(shouldThreadOperators() ? 0 : 1,
		pointGroupScatter_sub<typename T::DataType>, gInfo.nr, gInfo.S, in->data(), out->data(), mMesh);
	#endif
	return out;
}
DataRptr pointGroupScatter(const DataRptr& in, const matrix3<int>& mMesh)
{	return pointGroupScatter<DataR>(in, mMesh);
}
complexDataRptr pointGroupScatter(const complexDataRptr& in, const matrix3<int>& mMesh)
{	return pointGroupScatter<complexDataR>(in, mMesh);
}


//point group gather
template<typename Tptr> Tptr pointGroupGather(const Tptr& in, const matrix3<int>& mMesh)
{	if(mMesh == matrix3<int>(1,1,1)) return in; //shortcut for identity
	//Gathering is equivalent to gathering with inverse rotation
	//Scattered stores are faster than scattered loads (so implement scatter in terms of gather)
	int mMeshDet = det(mMesh);
	assert(abs(mMeshDet)==1);
	matrix3<int> mMeshInv = adjugate(mMesh)*mMeshDet; //inverse = adjugate*det since |det|=1
	return pointGroupScatter(in, mMeshInv);
}
DataRptr pointGroupGather(const DataRptr& in, const matrix3<int>& mMesh)
{	return pointGroupGather<DataRptr>(in, mMesh);
}
complexDataRptr pointGroupGather(const complexDataRptr& in, const matrix3<int>& mMesh)
{	return pointGroupGather<complexDataRptr>(in, mMesh);
}


void radialFunction_sub(size_t iStart, size_t iStop, const vector3<int> S, const matrix3<>& GGT,
	complex* F, const RadialFunctionG& f, vector3<> r0 )
{	THREAD_halfGspaceLoop( F[i] = radialFunction_calc(iG, GGT, f, r0); )
}
#ifdef GPU_ENABLED
void radialFunction_gpu(const vector3<int> S, const matrix3<>& GGT,
	complex* F, const RadialFunctionG& f, vector3<> r0);
#endif
DataGptr radialFunctionG(const GridInfo& gInfo, const RadialFunctionG& f, vector3<> r0)
{	
	DataGptr F(DataG::alloc(gInfo,isGpuEnabled()));
	#ifdef GPU_ENABLED
	radialFunction_gpu(gInfo.S, gInfo.GGT, F->dataGpu(), f, r0);
	#else
	threadLaunch(radialFunction_sub, gInfo.nG, gInfo.S, gInfo.GGT, F->data(), f, r0);
	#endif
	return F;
}

DataRptr radialFunction(const GridInfo& gInfo, const RadialFunctionG& f, vector3<> r0)
{	
	DataGptr F = radialFunctionG(gInfo, f, r0);
	return (1.0/gInfo.detR) * I(F, true);
}

void radialFunctionG(const RadialFunctionG& f, RealKernel& Kernel)
{	
	DataGptr F = radialFunctionG(Kernel.gInfo, f, vector3<>(0,0,0));
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

DataGptr operator*(const RadialFunctionG& f, DataGptr&& in)
{	const GridInfo& gInfo = in->gInfo;
	#ifdef GPU_ENABLED
	radialFunctionMultiply_gpu(gInfo.S, gInfo.GGT, in->dataGpu(), f);
	#else
	threadLaunch(radialFunctionMultiply_sub, gInfo.nG, gInfo.S, gInfo.GGT, in->data(), f);
	#endif
	return in;
}

DataGptr operator*(const RadialFunctionG& f, const DataGptr& in)
{	DataGptr out(in->clone()); //destructible copy
	return f * ((DataGptr&&)out);
}

DataGptrVec operator*(const RadialFunctionG& f, DataGptrVec&& in)
{	for(int k=0; k<3; k++) in[k] = f * (DataGptr&&)in[k];
	return in;
}

DataGptrVec operator*(const RadialFunctionG& f, const DataGptrVec& in)
{	DataGptrVec out;
	for(int k=0; k<3; k++) out[k] = f * in[k];
	return out;
}


//------------------------------ ColumnBundle operators ---------------------------------

void Idag_DiagV_I_sub(int colStart, int colEnd, const ColumnBundle* C, const DataRptr& V, ColumnBundle* VC)
{	for(int col=colStart; col<colEnd; col++)
		VC->setColumn(col, Idag(V * I(C->getColumn(col))));
}
ColumnBundle Idag_DiagV_I(const ColumnBundle& C, const DataRptr& V)
{	static StopWatch watch("Idag_DiagV_I"); watch.start();
	ColumnBundle VC = C.similar();
	suspendOperatorThreading();
	threadLaunch(isGpuEnabled()?1:0, Idag_DiagV_I_sub, C.nCols(), &C, V, &VC);
	resumeOperatorThreading();
	watch.stop();
	return VC;
}


ColumnBundle Pbar(const ColumnBundle &C, const ColumnBundle &Y)
{	return Y - O(C*(C^Y));
}


//Laplacian of a column bundle
#ifdef GPU_ENABLED
void reducedL_gpu(int nbasis, int ncols, const complex* Y, complex* LY,
	const matrix3<> GGT, const vector3<int>* iGarr, const vector3<> k, double detR);
#endif
ColumnBundle L(const ColumnBundle &Y)
{	ColumnBundle LY = Y.similar();
	const Basis& basis = *(Y.basis);
	const matrix3<>& GGT = basis.gInfo->GGT;
	#ifdef GPU_ENABLED
	reducedL_gpu(Y.colLength(), Y.nCols(), Y.dataGpu(), LY.dataGpu(), GGT, basis.iGarrGpu, Y.qnum->k, basis.gInfo->detR);
	#else
	threadedLoop(reducedL_calc, Y.colLength(),
		Y.colLength(), Y.nCols(), Y.data(), LY.data(), GGT, basis.iGarr, Y.qnum->k, basis.gInfo->detR);
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
	const Basis& basis = *(Y.basis);
	const matrix3<>& GGT = basis.gInfo->GGT;
	#ifdef GPU_ENABLED
	reducedLinv_gpu(Y.colLength(), Y.nCols(), Y.dataGpu(), LinvY.dataGpu(), GGT, basis.iGarrGpu, Y.qnum->k, basis.gInfo->detR);
	#else
	threadedLoop(reducedLinv_calc, Y.colLength(),
		Y.colLength(), Y.nCols(), Y.data(), LinvY.data(), GGT, basis.iGarr, Y.qnum->k, basis.gInfo->detR);
	#endif
	return LinvY;
}


// Overlap operator (scale by unit cell volume in PW basis)
ColumnBundle O(const ColumnBundle &Y)
{	ColumnBundle OY = Y * Y.basis->gInfo->detR; //basic planewave overlap
	Y.basis->iInfo->augmentOverlap(Y, OY); //pseudopotential augmentation
	return OY;
}

//Compute cartesian gradient of column bundle in direction #iDir
#ifdef GPU_ENABLED
void reducedD_gpu(int nbasis, int ncols, const complex* Ydata, complex* DYdata,
	const vector3<int>* iGarr, double kdotGe, const vector3<> Ge);
#endif
ColumnBundle D(const ColumnBundle &Y, int iDir)
{	const Basis& basis = *(Y.basis);
	ColumnBundle DY = Y.similar();
	const vector3<> Ge = basis.gInfo->G.column(iDir);
	double kdotGe = dot(Y.qnum->k, Ge);
	#ifdef GPU_ENABLED
	reducedD_gpu(Y.colLength(), Y.nCols(), Y.dataGpu(), DY.dataGpu(), basis.iGarrGpu, kdotGe, Ge);
	#else
	threadedLoop(reducedD_calc, Y.colLength(),
		Y.colLength(), Y.nCols(), Y.data(), DY.data(), basis.iGarr, kdotGe, Ge);
	#endif
	return DY;
}


// Multiply each column by f(0.5*|k+G|^2/KErollover)
// with f(x) = (1+x+x^2+x^3+...+x^8)/(1+x+x^2+...+x^9) = (1-x^N)/(1-x^(N+1))
#ifdef GPU_ENABLED
void precond_inv_kinetic_gpu(int nbasis, int ncols, const complex* Ydata, complex* KYdata,
	double KErollover, const matrix3<> GGT, const vector3<int>* iGarr, const vector3<> k, double invdetR);
#endif
ColumnBundle precond_inv_kinetic(const ColumnBundle &Y, double KErollover)
{	const Basis& basis = *Y.basis;
	const matrix3<>& GGT = basis.gInfo->GGT;
	ColumnBundle KY = Y.similar();
	#ifdef GPU_ENABLED
	precond_inv_kinetic_gpu(Y.colLength(), Y.nCols(), Y.dataGpu(), KY.dataGpu(),
		KErollover, GGT, basis.iGarrGpu, Y.qnum->k, 1/basis.gInfo->detR);
	#else
	threadedLoop(precond_inv_kinetic_calc, Y.colLength(),
		Y.colLength(), Y.nCols(), Y.data(), KY.data(),
		KErollover, GGT, basis.iGarr, Y.qnum->k, 1/basis.gInfo->detR);
	#endif
	return KY;
}

#ifdef GPU_ENABLED
void translate_gpu(int nbasis, int ncols, complex* Y, const vector3<int>* iGarr, const vector3<>& k, const vector3<>& dr);
#endif
ColumnBundle translate(ColumnBundle&& Y, vector3<> dr)
{	const Basis& basis = *Y.basis;
	#ifdef GPU_ENABLED
	translate_gpu(Y.colLength(), Y.nCols(), Y.dataGpu(), basis.iGarrGpu, Y.qnum->k, dr);
	#else
	threadedLoop(translate_calc, Y.colLength(), Y.colLength(), Y.nCols(), Y.data(), basis.iGarr, Y.qnum->k, dr);
	#endif
	return Y;
}
ColumnBundle translate(const ColumnBundle& Y, vector3<> dr)
{	return translate((ColumnBundle&&)ColumnBundle(Y), dr); //call above function on a destructible copy
}


ColumnBundle switchBasis(const ColumnBundle& in, const Basis& basisOut)
{	if(in.basis == &basisOut) return in; //no basis change required
	ColumnBundle out(in.nCols(), basisOut.nbasis, &basisOut, 0, isGpuEnabled());
	for(int i=0; i<in.nCols(); i++)
		out.setColumn(i, in.getColumn(i)); //convert using the full G-space as an intermediate
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
void diagouterI_sub(int iThread, int nThreads, const diagMatrix *F, const ColumnBundle *X, std::vector<DataRptr>* nSub)
{
	//Determine column range:
	int colStart = (( iThread ) * X->nCols())/nThreads;
	int colStop  = ((iThread+1) * X->nCols())/nThreads;
	
	DataRptr& nLocal = (*nSub)[iThread];
	nullToZero(nLocal, *(X->basis->gInfo)); //sets to zero
	// loop over subset of columns of X and accumulate density into nLocal
	for(int i=colStart; i<colStop; i++)
		callPref(eblas_accumNorm)(X->basis->gInfo->nr, (*F)[i], I(X->getColumn(i))->dataPref(), nLocal->dataPref());
}

// Collect all contributions from nSub into the first entry
void diagouterI_collect(size_t iStart, size_t iStop, std::vector<DataRptr>* nSub)
{	assert(!isGpuEnabled()); // this is needed and should be called only in CPU mode
	//Get the data pointers for each piece in nSub:
	int nThreads = nSub->size();
	std::vector<double*> nSubData(nThreads);
	for(int j=0; j<nThreads; j++) nSubData[j] = (*nSub)[j]->data();

	//Accumulate pointwise into the first piece:
	for(size_t i=iStart; i<iStop; i++)
		for(int j=1; j<nThreads; j++)
			nSubData[0][i] += nSubData[j][i];
}

// Returns diag((I*X)*F*(I*X)^) where X^ is the hermetian adjoint of X.
DataRptr diagouterI(const diagMatrix &F,const ColumnBundle &X)
{	static StopWatch watch("diagouterI"); watch.start();
	assert(F.nRows()==X.nCols());

	//Collect the contributions for different sets of columns in separate vectors (one per thread):
	int nThreads = isGpuEnabled() ? 1: nProcsAvailable;
	std::vector<DataRptr> nSub(nThreads);
	threadLaunch(nThreads, diagouterI_sub, 0, &F, &X, &nSub);

	//If more than one thread, accumulate all vectors in nSub into the first:
	if(nThreads>1) threadLaunch(diagouterI_collect, X.basis->gInfo->nr, &nSub);
	watch.stop();
	return nSub[0]; //rest will be cleaned up by the destructors
}

