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

#include <electronic/ColumnBundle.h>
#include <electronic/ColumnBundleOperators_internal.h>
#include <electronic/Basis.h>
#include <electronic/ElecInfo.h>
#include <electronic/IonInfo.h>
#include <core/matrix.h>
#include <core/Thread.h>
#include <core/BlasExtra.h>
#include <core/GpuUtil.h>
#include <core/GridInfo.h>
#include <core/LoopMacros.h>
#include <core/Operators.h>

//------------------------ Arithmetic operators --------------------

ColumnBundle& operator+=(ColumnBundle& Y, const scaled<ColumnBundle> &X) { if(Y) axpy(+X.scale, X.data, Y); else Y=X; return Y; }
ColumnBundle& operator-=(ColumnBundle& Y, const scaled<ColumnBundle> &X) { if(Y) axpy(-X.scale, X.data, Y); else Y=-X; return Y; }
ColumnBundle operator+(const scaled<ColumnBundle> &Y1, const scaled<ColumnBundle> &Y2) { ColumnBundle Ysum(Y1); Ysum += Y2; return Ysum; }
ColumnBundle operator-(const scaled<ColumnBundle> &Y1, const scaled<ColumnBundle> &Y2) { ColumnBundle Ydiff(Y1); Ydiff -= Y2; return Ydiff; }

ColumnBundle& operator*=(ColumnBundle& X, double s) { scale(s, X); return X; }
ColumnBundle operator*(double s, ColumnBundle&& Y) { scale(s, Y); return Y; }
ColumnBundle operator*(ColumnBundle&& Y, double s) { scale(s, Y); return Y; }
scaled<ColumnBundle> operator*(double s, const ColumnBundle &Y) { return scaled<ColumnBundle>(Y, s); }
scaled<ColumnBundle> operator*(const ColumnBundle &Y, double s) { return scaled<ColumnBundle>(Y, s); }
scaled<ColumnBundle> operator-(const ColumnBundle &Y) { return scaled<ColumnBundle>(Y, -1); }
ColumnBundle& operator*=(ColumnBundle& X, complex s) { scale(s, X); return X; }
ColumnBundle operator*(complex s, const ColumnBundle &Y) { ColumnBundle sY(Y); sY *= s; return sY; }
ColumnBundle operator*(const ColumnBundle &Y, complex s) { ColumnBundle sY(Y); sY *= s; return sY; }

ColumnBundleMatrixProduct::operator ColumnBundle() const
{	ColumnBundle YM;
	scaleAccumulate(1., 0., YM);
	return YM;
}

void ColumnBundleMatrixProduct::scaleAccumulate(double alpha, double beta, ColumnBundle& YM) const
{	static StopWatch watch("Y*M");
	watch.start();
	double scaleFac = alpha * scale * Mst.scale;
	bool spinorMode = (2*Y.nCols() == Mst.nRows()); //treat each column of non-spinor Y as two identical consecutive spinor ones with opposite spins
	assert(spinorMode || Y.nCols()==Mst.nRows());
	matrix Mtmp;
	CBLAS_TRANSPOSE Mop;
	const complex* Mdata; //pointer to start of M's data
	int ldM; //leading dimension of M
	int nColsOut; //number of columns in the result
	if(spinorMode)
	{	matrix mIn(Mst); Mop=CblasNoTrans; //pre-apply the op in this case
		Mtmp.init(Y.nCols(), 2*mIn.nCols(), isGpuEnabled());
		Mtmp.set(0,1,Y.nCols(), 0,2,Mtmp.nCols(), mIn(0,2,mIn.nRows(), 0,1,mIn.nCols()));
		Mtmp.set(0,1,Y.nCols(), 1,2,Mtmp.nCols(), mIn(1,2,mIn.nRows(), 0,1,mIn.nCols()));
		Mdata = Mtmp.dataPref();
		ldM = Mtmp.nRows();
		nColsOut = Mtmp.nCols();
		assert(!Y.isSpinor());
		if(beta) { assert(YM); assert(YM.nCols()==mIn.nCols()); assert(YM.colLength()==Y.colLength()*2); }
		else YM.init(mIn.nCols(), Y.colLength()*2, Y.basis, Y.qnum, isGpuEnabled());
	}
	else
	{	Mop = Mst.op;
		Mdata = Mst.mat.dataPref() + Mst.index(0,0);
		ldM = Mst.mat.nRows();
		nColsOut = Mst.nCols();
		if(beta) { assert(YM); assert(YM.nCols()==nColsOut); assert(YM.colLength()==Y.colLength()); }
		else YM = Y.similar(nColsOut);
	}
	callPref(eblas_zgemm)(CblasNoTrans, Mop, Y.colLength(), nColsOut, Y.nCols(),
		scaleFac, Y.dataPref(), Y.colLength(), Mdata, ldM,
		beta, YM.dataPref(), Y.colLength());
	watch.stop();
}

ColumnBundleMatrixProduct operator*(const scaled<ColumnBundle>& sY, const matrixScaledTransOp& Mst)
{	return ColumnBundleMatrixProduct(sY.data, Mst, sY.scale);
}
ColumnBundle& operator+=(ColumnBundle& Y, const ColumnBundleMatrixProduct &XM)
{	XM.scaleAccumulate(+1.,1.,Y);
	return Y;
}
ColumnBundle& operator-=(ColumnBundle& Y, const ColumnBundleMatrixProduct &XM)
{	XM.scaleAccumulate(-1.,1.,Y);
	return Y;
}
ColumnBundle operator+(const ColumnBundleMatrixProduct &XM1, const ColumnBundleMatrixProduct &XM2)
{	ColumnBundle result(XM1);
	result += XM2;
	return result;
}
ColumnBundle operator-(const ColumnBundleMatrixProduct &XM1, const ColumnBundleMatrixProduct &XM2)
{	ColumnBundle result(XM1);
	result -= XM2;
	return result;
}

ColumnBundle operator*(const scaled<ColumnBundle> &sY, const diagMatrix& d)
{	const ColumnBundle& Y = sY.data;
	assert(Y.nCols()==d.nRows());
	ColumnBundle Yd = Y; complex* YdData = Yd.dataPref();
	for(int i=0; i<d.nCols(); i++)
		callPref(eblas_zscal)(Yd.colLength(), sY.scale*d[i], YdData+Yd.index(i,0), 1);
	return Yd;
}

matrix operator^(const scaled<ColumnBundle> &sY1, const scaled<ColumnBundle> &sY2)
{	static StopWatch watch("Y1^Y2");
	watch.start();
	const ColumnBundle& Y1 = sY1.data;
	const ColumnBundle& Y2 = sY2.data;
	double scaleFac = sY1.scale * sY2.scale;
	int nCols1, nCols2, colLength;
	if(Y1.colLength() == Y2.colLength()) //standard mode
	{	nCols1 = Y1.nCols();
		nCols2 = Y2.nCols();
		colLength = Y1.colLength();
	}
	else //exactly one of the two columnbundles is a spinor (but they have a common basis)
	{	assert(Y1.basis);
		assert(Y2.basis);
		assert(Y1.basis->nbasis == Y2.basis->nbasis);
		assert(Y1.isSpinor() xor Y2.isSpinor());
		nCols1 = Y1.nCols() * Y1.spinorLength();
		nCols2 = Y2.nCols() * Y2.spinorLength();
		colLength = Y1.basis->nbasis;
	}
	matrix Y1dY2(nCols1, nCols2, isGpuEnabled());
	callPref(eblas_zgemm)(CblasConjTrans, CblasNoTrans, nCols1, nCols2, colLength,
		scaleFac, Y1.dataPref(), colLength, Y2.dataPref(), colLength,
		0.0, Y1dY2.dataPref(), Y1dY2.nRows());
	watch.stop();
	//If one of the columnbundles was spinor, shape the matrix as if the non-spinor columnbundle had consecutive spinor columns with identical pure up and down spinors
	if(Y1.nCols() != nCols1) //Y1 is spinor, so double the dimension of output along Y2
	{	matrix out(Y1.nCols(), 2*nCols2);
		out.set(0,1,Y1.nCols(), 0,2,2*nCols2, Y1dY2(0,2,nCols1, 0,1,nCols2));
		out.set(0,1,Y1.nCols(), 1,2,2*nCols2, Y1dY2(1,2,nCols1, 0,1,nCols2));
		return out;
	}
	else if(Y2.nCols() != nCols2) //Y2 is spinor, so double the dimension of output along Y1
	{	matrix out(2*nCols1, Y2.nCols());
		out.set(0,2,2*nCols1, 0,1,Y2.nCols(), Y1dY2(0,1,nCols1, 0,2,nCols2));
		out.set(1,2,2*nCols1, 0,1,Y2.nCols(), Y1dY2(0,1,nCols1, 1,2,nCols2));
		return out;
	}
	else return Y1dY2; //normal mode (neither is a spinor)
}

vector3<matrix> spinOverlap(const scaled<ColumnBundle> &sY)
{	const ColumnBundle& Y = sY.data;
	double scaleFac = std::pow(sY.scale, 2) * Y.basis->gInfo->detR; //norm conserving part of O
	assert(Y.isSpinor());
	//Matrix multiply with an effective reshape to 2*nCols and 1/2 colLength:
	int nCols = 2*Y.nCols();
	int colLength = Y.basis->nbasis; //= Y.colLength()/2
	matrix YdotY(nCols, nCols, isGpuEnabled());
	callPref(eblas_zgemm)(CblasConjTrans, CblasNoTrans, nCols, nCols, colLength,
		scaleFac, Y.dataPref(), colLength, Y.dataPref(), colLength,
		0.0, YdotY.dataPref(), YdotY.nRows());
	//Extract 2x2 components:
	matrix S[2][2];
	for(int s1=0; s1<2; s1++)
		for(int s2=0; s2<2; s2++)
			S[s1][s2] = YdotY(s1,2,nCols, s2,2,nCols);
	//Recombine into spin vector components:
	vector3<matrix> Svec(
		S[0][1] + S[1][0], //Sx
		(S[0][1] - S[1][0]) * complex(0,-1), //Sy
		S[0][0] - S[1][1] //Sz
	);
	//Ultrasoft augmentation (if necessary):
	for(const auto& sp: Y.basis->iInfo->species)
		sp->augmentSpinOverlap(Y, Svec);
	return Svec;
}


//------------------------------ Other operators ---------------------------------

void Idag_DiagV_I_sub(int colStart, int colEnd, const ColumnBundle* C, const ScalarFieldArray* V, ColumnBundle* VC)
{	const ScalarField& Vs = V->at(V->size()==1 ? 0 : C->qnum->index());
	int nSpinor = VC->spinorLength();
	for(int col=colStart; col<colEnd; col++)
		for(int s=0; s<nSpinor; s++)
			VC->accumColumn(col,s, Idag(Vs * I(C->getColumn(col,s)))); //note VC is zero'd just before
}

//Noncollinear version of above (with the preprocessing of complex off-diagonal potentials done in calling function)
void Idag_DiagVmat_I_sub(int colStart, int colEnd, const ColumnBundle* C, const ScalarField* Vup, const ScalarField* Vdn,
	const complexScalarField* VupDn, const complexScalarField* VdnUp, ColumnBundle* VC)
{	for(int col=colStart; col<colEnd; col++)
	{	complexScalarField ICup = I(C->getColumn(col,0));
		complexScalarField ICdn = I(C->getColumn(col,1));
		VC->accumColumn(col,0, Idag((*Vup)*ICup + (*VupDn)*ICdn));
		VC->accumColumn(col,1, Idag((*Vdn)*ICdn + (*VdnUp)*ICup));
	}
	
}

ColumnBundle Idag_DiagV_I(const ColumnBundle& C, const ScalarFieldArray& V)
{	static StopWatch watch("Idag_DiagV_I"); watch.start();
	ColumnBundle VC = C.similar(); VC.zero();
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
	reducedL_gpu(basis.nbasis, Y.nCols()*nSpinors, Y.dataGpu(), LY.dataGpu(), GGT, basis.iGarr.dataGpu(), Y.qnum->k, basis.gInfo->detR);
	#else
	threadedLoop(reducedL_calc, basis.nbasis,
		basis.nbasis, Y.nCols()*nSpinors, Y.data(), LY.data(), GGT, basis.iGarr.data(), Y.qnum->k, basis.gInfo->detR);
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
	reducedLinv_gpu(basis.nbasis, Y.nCols()*nSpinors, Y.dataGpu(), LinvY.dataGpu(), GGT, basis.iGarr.dataGpu(), Y.qnum->k, basis.gInfo->detR);
	#else
	threadedLoop(reducedLinv_calc, basis.nbasis,
		basis.nbasis, Y.nCols()*nSpinors, Y.data(), LinvY.data(), GGT, basis.iGarr.data(), Y.qnum->k, basis.gInfo->detR);
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
	reducedD_gpu(basis.nbasis, Y.nCols()*nSpinors, Y.dataGpu(), DY.dataGpu(), basis.iGarr.dataGpu(), kdotGe, Ge);
	#else
	threadedLoop(reducedD_calc, basis.nbasis,
		basis.nbasis, Y.nCols()*nSpinors, Y.data(), DY.data(), basis.iGarr.data(), kdotGe, Ge);
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
	reducedDD_gpu(basis.nbasis, Y.nCols()*nSpinors, Y.dataGpu(), DDY.dataGpu(), basis.iGarr.dataGpu(), kdotGe1, kdotGe2, Ge1, Ge2);
	#else
	threadedLoop(reducedDD_calc, basis.nbasis,
		basis.nbasis, Y.nCols()*nSpinors, Y.data(), DDY.data(), basis.iGarr.data(), kdotGe1, kdotGe2, Ge1, Ge2);
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
		KErollover, GGT, basis.iGarr.dataPref(), Y.qnum->k, 1./basis.gInfo->detR);
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
		basis.gInfo->GGT, basis.iGarr.dataPref(), Y.qnum->k);
}


#ifdef GPU_ENABLED
void translate_gpu(int nbasis, int ncols, complex* Y, const vector3<int>* iGarr, const vector3<>& k, const vector3<>& dr);
#endif
ColumnBundle translate(ColumnBundle&& Y, vector3<> dr)
{	assert(Y.basis);
	const Basis& basis = *Y.basis;
	int nSpinors = Y.spinorLength();
	#ifdef GPU_ENABLED
	translate_gpu(basis.nbasis, Y.nCols()*nSpinors, Y.dataGpu(), basis.iGarr.dataGpu(), Y.qnum->k, dr);
	#else
	threadedLoop(translate_calc, basis.nbasis, basis.nbasis, Y.nCols()*nSpinors, Y.data(), basis.iGarr.data(), Y.qnum->k, dr);
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
	int nColsTot = Y.nCols()*nSpinors;
	ManagedArray<vector3<>> drManaged(dr, nColsTot);
	callPref(translateColumns)(basis.nbasis, nColsTot, Y.dataPref(), basis.iGarr.dataPref(), Y.qnum->k, drManaged.dataPref());
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
