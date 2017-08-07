/*-------------------------------------------------------------------
Copyright 2014 Ravishankar Sundararaman

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

#include <wannier/WannierMinimizer.h>
#include <electronic/SpeciesInfo_internal.h>
#include <core/BlasExtra.h>
#include <core/Random.h>

void WannierGradient::init(const WannierMinimizer* wmin)
{	this->wmin = wmin;
	resize(wmin->kMesh.size());
}
size_t WannierGradient::ikStart() const { return wmin->ikStart; }
size_t WannierGradient::ikStop() const { return wmin->ikStop; }

//---- linear algebra functions required by Minimizable<WannierGradient> -----

WannierGradient clone(const WannierGradient& grad) { return grad; }
double dot(const WannierGradient& x, const WannierGradient& y)
{	assert(x.size()==y.size());
	double result = 0.;
	for(unsigned ik=x.ikStart(); ik<x.ikStop(); ik++)
	{	result += dotc(x[ik], y[ik]).real();
		//For rectangular matrices, account for the fact that we are actually working with the hermitian completion
		if(x[ik].nCols() != x[ik].nRows())
		{	int rStart=0, rStop=x[ik].nRows();
			int cStart=0, cStop=x[ik].nCols();
			if(rStop>cStop) rStart=cStop; else cStart=rStop;
			result += dotc(x[ik](rStart,rStop,cStart,cStop), y[ik](rStart,rStop,cStart,cStop)).real();
		}
	}
	mpiUtil->allReduce(result, MPIUtil::ReduceSum);
	return result;
}
WannierGradient& operator*=(WannierGradient& x, double alpha)
{	for(unsigned ik=x.ikStart(); ik<x.ikStop(); ik++)
		x[ik] *= alpha;
	return x;
}
void axpy(double alpha, const WannierGradient& x, WannierGradient& y)
{	assert(x.size()==y.size());
	for(unsigned ik=x.ikStart(); ik<x.ikStop(); ik++)
		axpy(alpha, x[ik], y[ik]);
}

matrix randomMatrix(int nRows, int nCols)
{	matrix ret(nRows, nCols, false);
	complex* retData = ret.data();
	for(unsigned j=0; j<ret.nData(); j++)
		retData[j] = Random::normalComplex();
	return ret;
}
void randomize(WannierGradient& x)
{	for(unsigned ik=x.ikStart(); ik<x.ikStop(); ik++) if(x[ik].nData())
	{	int minDim = std::min(x[ik].nRows(), x[ik].nCols());
		x[ik].set(0,minDim, 0,minDim, dagger_symmetrize(randomMatrix(minDim,minDim)));
		if(x[ik].nRows()>minDim) x[ik].set(minDim,x[ik].nRows(), 0,minDim, randomMatrix(x[ik].nRows()-minDim,minDim));
		if(x[ik].nCols()>minDim) x[ik].set(0,minDim, minDim,x[ik].nCols(), randomMatrix(minDim,x[ik].nCols()-minDim));
	}
}

//---- energy/gradient functions required by Minimizable<WannierGradient> -----

void WannierMinimizer::step(const WannierGradient& grad, double alpha)
{	assert(grad.wmin == this);
	for(unsigned ik=ikStart; ik<ikStop; ik++)
	{	KmeshEntry& ki = kMesh[ik];
		matrix B = alpha * grad[ik];
		//Stage 1:
		if(ki.nIn > nCenters)
		{	matrix B1block = B(ki.nFixed,nCenters, nCenters,ki.nIn);
			matrix B1 = zeroes(ki.nIn, ki.nIn);
			B1.set(ki.nFixed,nCenters, nCenters,ki.nIn, B1block);
			B1.set(nCenters,ki.nIn, ki.nFixed,nCenters, dagger(B1block));
			ki.U1 = ki.U1 * cis(B1);
		}
		//Stage 2:
		matrix B2 = zeroes(nCenters, nCenters);
		B2.set(nFrozen,nCenters, nFrozen,nCenters, B(nFrozen,nCenters, nFrozen,nCenters));
		ki.U2 = ki.U2 * cis(B2);
		//Net rotation:
		ki.U = ki.U1(0,nBands, 0,nCenters) * ki.U2;
	}
	for(size_t ik=0; ik<kMesh.size(); ik++) kMesh[ik].U.bcast(whose(ik)); //Make U available on all processes
}


double WannierMinimizer::compute(WannierGradient* grad, WannierGradient* Kgrad)
{
	if(grad) for(KmeshEntry& ki: kMesh) ki.Omega_U = zeroes(nCenters, nBands); //Clear Omega_U
	
	double Omega = getOmega(grad);
	
	//Collect Omega_U and propagate to Omega_B if necessary:
	if(grad)
	{	for(KmeshEntry& ki: kMesh) ki.Omega_U.allReduce(MPIUtil::ReduceSum); //Collect Omega_U
		grad->init(this);
		for(size_t ik=ikStart; ik<ikStop; ik++)
		{	KmeshEntry& ki = kMesh[ik];
			(*grad)[ik] = zeroes(nCenters, ki.nIn);
			if(ki.nIn > nCenters) //Stage 1:
			{	matrix Omega_B1 = complex(0,0.5)*(ki.U2 * ki.Omega_U * ki.U1);
				(*grad)[ik].set(ki.nFixed,nCenters, nCenters,ki.nIn, Omega_B1(ki.nFixed,nCenters, nCenters,ki.nIn));
			}
			matrix Omega_B2 = dagger_symmetrize(complex(0,1)*(ki.Omega_U * ki.U));
			(*grad)[ik].set(nFrozen,nCenters, nFrozen,nCenters, Omega_B2(nFrozen,nCenters, nFrozen,nCenters));
		}
		
		if(Kgrad) *Kgrad = precondition(*grad);
	}
	return Omega;
}

void WannierMinimizer::constrain(WannierGradient& grad)
{	for(size_t ik=ikStart; ik<ikStop; ik++)
	{	const KmeshEntry& ki = kMesh[ik];
		//Pick out the parts that are allowed to be non-zero:
		matrix gradB1;
		if(ki.nIn > nCenters)
			gradB1 = grad[ik](ki.nFixed,nCenters, nCenters,ki.nIn);
		matrix gradB2 = grad[ik](nFrozen,nCenters, nFrozen,nCenters);
		//Zero everything and only set back the non-zero parts:
		grad[ik].zero();
		if(gradB1) grad[ik].set(ki.nFixed,nCenters, nCenters,ki.nIn, gradB1);
		grad[ik].set(nFrozen,nCenters, nFrozen,nCenters, dagger_symmetrize(gradB2));
	}
}

WannierGradient WannierMinimizer::precondition(const WannierGradient& grad)
{	WannierGradient Kgrad = grad;
	constrain(Kgrad);
	return Kgrad;
}


matrix WannierMinimizer::fixUnitary(const matrix& U)
{	return U * invsqrt(dagger(U) * U);
}

bool WannierMinimizer::report(int iter)
{	//Check unitarity:
	bool needRestart = false;
	for(size_t ik=ikStart; ik<ikStop; ik++)
		if(nrm2(dagger(kMesh[ik].U) * kMesh[ik].U - eye(nCenters)) > 1e-6)
		{	needRestart = true;
			break;
		}
	mpiUtil->allReduce(needRestart, MPIUtil::ReduceLOr);
	if(needRestart)
	{	logPrintf("%s\tUpdating rotations to enforce unitarity\n", wannier.minParams.linePrefix);
		for(size_t ik=ikStart; ik<ikStop; ik++)
		{	KmeshEntry& ki = kMesh[ik];
			ki.U1 = fixUnitary(ki.U1);
			ki.U2 = fixUnitary(ki.U2);
			ki.U = ki.U1(0,nBands, 0,nCenters) * ki.U2;
		}
		return true;
	}
    return false;
}

double WannierMinimizer::sync(double x) const
{	mpiUtil->bcast(x);
	return x;
}

//---------------- kpoint and wavefunction handling -------------------

bool WannierMinimizer::Kpoint::operator<(const WannierMinimizer::Kpoint& other) const
{	if(iReduced!=other.iReduced) return iReduced<other.iReduced;
	if(iSym!=other.iSym) return iSym<other.iSym;
	if(invert!=other.invert) return invert<other.invert;
	if(!(offset==other.offset)) return offset<other.offset;
	return false; //all equal
}

bool WannierMinimizer::Kpoint::operator==(const WannierMinimizer::Kpoint& other) const
{	if(iReduced!=other.iReduced) return false;
	if(iSym!=other.iSym) return false;
	if(invert!=other.invert) return false;
	if(!(offset==other.offset)) return false;
	return true;
}

ColumnBundle WannierMinimizer::getWfns(const WannierMinimizer::Kpoint& kpoint, int iSpin) const
{	ColumnBundle ret(nBands, basis.nbasis*nSpinor, &basis, &kpoint, isGpuEnabled());
	ret.zero();
	axpyWfns(1., matrix(), kpoint, iSpin, ret);
	return ret;
}

#define axpyWfns_COMMON(result) \
	/* Pick transform: */ \
	const ColumnBundleTransform& transform = *(((result.basis==&basisSuper) ? transformMapSuper : transformMap).find(kpoint)->second); \
	/* Pick source ColumnBundle: */ \
	int q = kpoint.iReduced + iSpin*qCount; \
	const ColumnBundle& Cin = e.eInfo.isMine(q) ? e.eVars.C[q] : Cother[q]; \
	assert(Cin); \
	const ColumnBundle* C = &Cin; \

void WannierMinimizer::axpyWfns(double alpha, const matrix& A, const WannierMinimizer::Kpoint& kpoint, int iSpin, ColumnBundle& result) const
{	static StopWatch watch("WannierMinimizer::axpyWfns"); watch.start();
	axpyWfns_COMMON(result)
	//Apply transformation if provided:
	ColumnBundle Cout;
	if(A)
	{	matrix Astar = (kpoint.invert<0 ? conj(A) : A);
		Cout = (*C) * Astar;
		C = &Cout;
	}
	//Scatter from reduced basis to common basis with transformations:
	assert(C->nCols() == result.nCols());
	transform.scatterAxpy(alpha, *C, result,0,1);
	watch.stop();
}

void WannierMinimizer::axpyWfns_grad(double alpha, matrix& Omega_A, const WannierMinimizer::Kpoint& kpoint, int iSpin, const ColumnBundle& Omega_result) const
{	static StopWatch watch("WannierMinimizer::axpyWfns_grad"); watch.start();
	axpyWfns_COMMON(Omega_result)
	//Gather from common basis to reduced basis (=> conjugate transformations):
	ColumnBundle Omega_C = C->similar(Omega_result.nCols());
	Omega_C.zero();
	transform.gatherAxpy(alpha, Omega_result,0,1, Omega_C);
	//Propagate gradient to rotation matrix:
	matrix Omega_Astar = Omega_C ^ *C;
	Omega_A += (kpoint.invert<0 ? conj(Omega_Astar) : Omega_Astar);
	watch.stop();
}

#undef axpyWfns_COMMON

//Fourier transform of hydrogenic orbitals
inline double hydrogenicTilde(double G, double a, int nIn, int l, double normPrefac)
{	int n = nIn+1 + l; //conventional principal quantum number
	double nG = n*G*a/(l+1), nGsq = nG*nG;
	double prefac = normPrefac / pow(1.+nGsq, n+1);
	switch(l)
	{	case 0:
			switch(n)
			{	case 1: return prefac;
				case 2: return prefac*8.*(-1.+nGsq);
				case 3: return prefac*9.*(3.+nGsq*(-10.+nGsq*3.));
				case 4: return prefac*64.*(-1.+nGsq*(7.+nGsq*(-7.+nGsq)));
			}
		case 1:
			switch(n)
			{	case 2: return prefac*16.*nG;
				case 3: return prefac*144.*nG*(-1.+nGsq);
				case 4: return prefac*128.*nG*(5.+nGsq*(-14.+nGsq*5.));
			}
		case 2:
			switch(n)
			{	case 3: return prefac*288.*nGsq;
				case 4: return prefac*3072.*nGsq*(-1.+nGsq);
			}
		case 3:
			switch(n)
			{	case 4: return prefac*6144.*nG*nGsq;
			}
	}
	return 0.;
}

ColumnBundle WannierMinimizer::trialWfns(const WannierMinimizer::Kpoint& kpoint) const
{	ColumnBundle ret(nCenters-nFrozen, basis.nbasis*nSpinor, &basis, &kpoint, isGpuEnabled());
	ColumnBundle temp = ret.similar(1); //single column for intermediate computations
	//Generate atomic orbitals if necessary:
	std::vector<ColumnBundle> psiAtomic;
	if(wannier.needAtomicOrbitals)
	{	psiAtomic.resize(e.iInfo.species.size());
		for(unsigned sp=0; sp<e.iInfo.species.size(); sp++)
		{	psiAtomic[sp].init(e.iInfo.species[sp]->nAtomicOrbitals(), basis.nbasis*nSpinor, &basis, &kpoint, isGpuEnabled());
			e.iInfo.species[sp]->setAtomicOrbitals(psiAtomic[sp], false);
		}
	}
	ret.zero();
	complex* retData = ret.dataPref();
	for(const Wannier::TrialOrbital& t: wannier.trialOrbitals)
	{	for(const Wannier::AtomicOrbital& ao: t)
		{	//Handle numerical orbitals:
			if(ao.numericalOrbIndex >= 0)
			{	const ColumnBundle& Cnum = *(numericalOrbitals.find(kpoint)->second);
				//Apply offset to selected column:
				assert(ao.numericalOrbIndex < Cnum.nCols());
				temp = translate(Cnum.getSub(ao.numericalOrbIndex,ao.numericalOrbIndex+1), ao.r);
				//Accumulate to result
				callPref(eblas_zaxpy)(ret.colLength(), ao.coeff, temp.dataPref(),1, retData,1);
				continue;
			}
			//Handle atomic orbitals that are actually atom-centered:
			const DOS::Weight::OrbitalDesc& od = ao.orbitalDesc;
			complex lPhase =  cis(0.5*M_PI*od.l); //including this phase ensures odd l projectors are real (i^l term in spherical wave expansion)
			if(ao.atom >= 0)
			{	int iCol = e.iInfo.species[ao.sp]->atomicOrbitalOffset(ao.atom, od.n, od.l, od.m, od.s);
				callPref(eblas_zaxpy)(ret.colLength(), ao.coeff*lPhase, psiAtomic[ao.sp].dataPref()+iCol*ret.colLength(),1, retData,1);
				continue;
			}
			//--- Copy the center to managed memory:
			ManagedArray<vector3<>> pos(&ao.r, 1);
			//--- Get / create the radial part:
			RadialFunctionG hRadial;
			if(ao.sp < 0)
			{	double normPrefac = pow((od.l+1)/ao.a,3);
				for(unsigned p=od.n+1; p<=od.n+1+2*od.l; p++)
					normPrefac *= p;
				normPrefac = 16*M_PI/sqrt(normPrefac);
				hRadial.init(od.l, 0.02, e.gInfo.GmaxSphere, hydrogenicTilde, ao.a, od.n, od.l, normPrefac);
			}
			const RadialFunctionG& atRadial = (ao.sp<0) ? hRadial : e.iInfo.species[ao.sp]->OpsiRadial->at(od.l)[od.n];
			//--- Initialize the projector:
			assert(od.s < nSpinor);
			if(nSpinor > 1) { temp.zero(); assert(od.spinType==SpinZ); } //The relativistic orbitals must be handled above via atom-centered orbitals
			callPref(Vnl)(basis.nbasis, basis.nbasis, 1, od.l, od.m, kpoint.k, basis.iGarr.dataPref(), e.gInfo.G, pos.dataPref(), atRadial, temp.dataPref()+od.s*basis.nbasis);
			if(ao.sp < 0) hRadial.free();
			//--- Accumulate to trial orbital:
			callPref(eblas_zaxpy)(ret.colLength(), ao.coeff*lPhase/e.gInfo.detR, temp.dataPref(),1, retData,1);
		}
		retData += ret.colLength();
	}
	return ret;
}

matrix WannierMinimizer::overlap(const ColumnBundle& C1, const ColumnBundle& C2) const
{	const GridInfo& gInfo = *(C1.basis->gInfo);
	const IonInfo& iInfo = *(C1.basis->iInfo);
	matrix ret = gInfo.detR * (C1 ^ C2);
	//k-point difference:
	vector3<> dkVec = C2.qnum->k - C1.qnum->k;
	double dk = sqrt(gInfo.GGT.metric_length_squared(dkVec));
	vector3<> dkHat = gInfo.GT * dkVec * (dk ? 1.0/dk : 0.0); //the unit Vector along dkVec (set dkHat to 0 for dk=0 (doesn't matter))
	//Augment at each species:
	for(const auto& sp: iInfo.species) if(sp->Qint.size())
	{	//Create the Q matrix appropriate for current k-point difference:
		matrix Qk = zeroes(sp->QintAll.nRows(), sp->QintAll.nCols());
		complex* QkData = Qk.data();
		int i1 = 0;
		for(int l1=0; l1<int(sp->VnlRadial.size()); l1++)
		for(int p1=0; p1<int(sp->VnlRadial[l1].size()); p1++)
		for(int m1=-l1; m1<=l1; m1++)
		{	//Triple loop over second projector:
			int i2 = 0;
			for(int l2=0; l2<int(sp->VnlRadial.size()); l2++)
			for(int p2=0; p2<int(sp->VnlRadial[l2].size()); p2++)
			for(int m2=-l2; m2<=l2; m2++)
			{	if(i2<=i1) //rest handled by i1<->i2 symmetry
				{	std::vector<YlmProdTerm> terms = expandYlmProd(l1,m1, l2,m2);
					complex q12 = 0.;
					for(const YlmProdTerm& term: terms)
					{	SpeciesInfo::QijIndex qIndex = { l1, p1, l2, p2, term.l };
						auto Qijl = sp->Qradial.find(qIndex);
						if(Qijl==sp->Qradial.end()) continue; //no entry at this l
						q12 += term.coeff * cis(-0.5*M_PI*term.l) * Ylm(term.l,term.m, dkHat) * Qijl->second(dk);
					}
					QkData[Qk.index(i1,i2)] = q12;
					QkData[Qk.index(i2,i1)] = q12.conj();
				}
				i2++;
			}
			i1++;
		}
		//Phases for each atom:
		std::vector<complex> phaseArr;
		for(vector3<> x: sp->atpos)
			phaseArr.push_back(cis(-2*M_PI*dot(dkVec,x)));
		//Augment the overlap
		matrix VdagC1 = (*sp->getV(C1)) ^ C1;
		matrix VdagC2 = (*sp->getV(C2)) ^ C2;
		ret += dagger(VdagC1) * (tiledBlockMatrix(Qk, sp->atpos.size(),&phaseArr) * VdagC2);
	}
	return ret;
}
