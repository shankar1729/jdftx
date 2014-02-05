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
#include <electronic/operators.h>
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
		axpy(alpha, grad[ik], kMesh[ik].B);
}

double WannierMinimizer::compute(WannierGradient* grad)
{	//Compute the unitary matrices:
	for(size_t ik=ikStart; ik<ikStop; ik++)
	{	KmeshEntry& ki = kMesh[ik];
		//Stage 1:
		if(ki.nIn > nCenters)
		{	matrix B1block = ki.B(ki.nFixed,nCenters, nCenters,ki.nIn);
			matrix B1 = zeroes(ki.nIn, ki.nIn);
			B1.set(ki.nFixed,nCenters, nCenters,ki.nIn, B1block);
			B1.set(nCenters,ki.nIn, ki.nFixed,nCenters, dagger(B1block));
			ki.V1 = cis(B1, &ki.B1evecs, &ki.B1eigs)(0,ki.nIn, 0,nCenters);
		}
		else ki.V1 = eye(nCenters);
		//Stage 2:
		ki.V2 = cis(ki.B(0,nCenters, 0,nCenters), &ki.B2evecs, &ki.B2eigs);
		//Net rotation:
		ki.U = ki.U1 * ki.V1 * ki.U2 * ki.V2;
	}
	for(size_t ik=0; ik<kMesh.size(); ik++) kMesh[ik].U.bcast(whose(ik));
	
	//Compute the expectation values of r and rSq for each center (split over processes)
	rSqExpect.assign(nCenters, 0.);
	rExpect.assign(nCenters, vector3<>(0,0,0));
	OmegaI = 0.;
	for(size_t ik=ikStart; ik<ikStop; ik++)
	{	const KmeshEntry& ki = kMesh[ik];
		for(const EdgeFD& edge: ki.edge)
		{	const KmeshEntry& kj = kMesh[edge.ik];
			const matrix M = dagger(ki.U) * edge.M0 * kj.U;
			OmegaI += ki.point.weight * edge.wb * (nCenters - trace(M * dagger(M)).real());
			const complex* Mdata = M.data();
			for(int n=0; n<nCenters; n++)
			{	complex Mnn = Mdata[M.index(n,n)];
				double argMnn = Mnn.arg();
				rExpect[n] -= (ki.point.weight * edge.wb * argMnn) * edge.b;
				rSqExpect[n] += ki.point.weight * edge.wb * (argMnn*argMnn + 1. - Mnn.norm());
			}
		}
	}
	mpiUtil->allReduce(OmegaI, MPIUtil::ReduceSum);
	mpiUtil->allReduce(rSqExpect.data(), nCenters, MPIUtil::ReduceSum);
	mpiUtil->allReduce((double*)rExpect.data(), 3*nCenters, MPIUtil::ReduceSum);
	
	//Compute the total variance of the Wannier centers
	double Omega = 0.;
	for(int n=0; n<nCenters; n++)
		Omega += (rSqExpect[n] - rExpect[n].length_squared());
	
	//Compute the gradients of the mean variance (if required)
	if(grad)
	{	//Accumulate gradients from each edge (split over processes):
		for(KmeshEntry& ki: kMesh)
			ki.Omega_U = zeroes(nCenters, nBands);
		for(size_t ik=ikStart; ik<ikStop; ik++)
		{	KmeshEntry& ki = kMesh[ik];
			for(EdgeFD& edge: ki.edge)
			{	KmeshEntry& kj = kMesh[edge.ik];
				const matrix M = dagger(ki.U) * edge.M0 * kj.U;
				//Compute dOmega/dM:
				matrix Omega_M = zeroes(nCenters, nCenters);
				const complex* Mdata = M.data();
				complex* Omega_Mdata = Omega_M.data();
				for(int n=0; n<nCenters; n++)
				{	complex Mnn = Mdata[M.index(n,n)];
					double argMnn = atan2(Mnn.imag(), Mnn.real());
					Omega_Mdata[Omega_M.index(n,n)] =
						2. * ki.point.weight * edge.wb
						* ((argMnn + dot(rExpect[n],edge.b))*complex(0,-1)/Mnn - Mnn.conj());
				}
				//Propagate Omega_M to Omega_U:
				ki.Omega_U += dagger(edge.M0 * kj.U * Omega_M);
				kj.Omega_U += Omega_M * dagger(ki.U) * edge.M0;
			}
		}
		for(KmeshEntry& ki: kMesh)
			ki.Omega_U.allReduce(MPIUtil::ReduceSum);
		//Propagate to gradients w.r.t B:
		grad->init(this);
		for(size_t ik=ikStart; ik<ikStop; ik++)
		{	KmeshEntry& ki = kMesh[ik];
			(*grad)[ik] = zeroes(nCenters, ki.nIn);
			if(ki.nIn > nCenters)
			{	matrix Omega_B1 = dagger_symmetrize(cis_grad(ki.V1 * ki.U2 * ki.V2 * ki.Omega_U * ki.U1, ki.B1evecs, ki.B1eigs));
				(*grad)[ik].set(ki.nFixed,nCenters, nCenters,ki.nIn, Omega_B1(ki.nFixed,nCenters, nCenters,ki.nIn));
			}
			(*grad)[ik].set(0,nCenters, 0,nCenters, dagger_symmetrize(cis_grad(ki.V2 * ki.Omega_U * ki.U1 * ki.V1 * ki.U2, ki.B2evecs, ki.B2eigs)));
		}
	}
	return Omega;
}

//---------------- kpoint and wavefunction handling -------------------

int WannierMinimizer::whose(size_t ik) const
{	if(mpiUtil->nProcesses()>1) return std::upper_bound(ikStopArr.begin(),ikStopArr.end(), ik) - ikStopArr.begin();
	else return 0;
}

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


WannierMinimizer::Index::Index(int nIndices, bool needSuper) : nIndices(nIndices), dataPref(0), dataSuperPref(0)
{	data = new int[nIndices];
	dataSuper = needSuper ? new int[nIndices] : 0;
	#ifdef GPU_ENABLED
	dataGpu = 0;
	dataSuperGpu = 0;
	#endif
}
WannierMinimizer::Index::~Index()
{	delete[] data;
	if(dataSuper) delete[] dataSuper;
	#ifdef GPU_ENABLED
	if(dataGpu) cudaFree(dataGpu);
	if(dataSuperGpu) cudaFree(dataSuperGpu);
	#endif
}
void WannierMinimizer::Index::set()
{
	#ifdef GPU_ENABLED
	cudaMalloc(&dataGpu, sizeof(int)*nIndices); gpuErrorCheck();
	cudaMemcpy(dataGpu, data, sizeof(int)*nIndices, cudaMemcpyHostToDevice); gpuErrorCheck();
	dataPref = dataGpu;
	if(dataSuper)
	{	cudaMalloc(&dataSuperGpu, sizeof(int)*nIndices); gpuErrorCheck();
		cudaMemcpy(dataSuperGpu, dataSuper, sizeof(int)*nIndices, cudaMemcpyHostToDevice); gpuErrorCheck();
	}
	else dataSuperGpu = 0;
	dataSuperPref = dataSuperGpu;
	#else
	dataPref = data;
	dataSuperPref = dataSuper;
	#endif
}

void WannierMinimizer::addIndex(const WannierMinimizer::Kpoint& kpoint)
{	if(indexMap.find(kpoint)!=indexMap.end()) return; //previously computed
	//Determine integer offset due to k-point in supercell basis:
	const matrix3<int>& super = e.coulombParams.supercell->super;
	vector3<int> ksuper; //integer version of above
	if(needSuper)
	{	vector3<> ksuperTemp = kpoint.k * super - qnumSuper.k; //note reciprocal lattice vectors transform on the right (or on the left by the transpose)
		for(int l=0; l<3; l++)
		{	ksuper[l] = int(round(ksuperTemp[l]));
			assert(fabs(ksuper[l]-ksuperTemp[l]) < symmThreshold);
		}
	}
	//Compute transformed index array (mapping to full G-space)
	const Basis& basis = e.basis[kpoint.iReduced];
	std::shared_ptr<Index> index(new Index(basis.nbasis, needSuper));
	const matrix3<int> mRot = (~sym[kpoint.iSym]) * kpoint.invert;
	for(int j=0; j<index->nIndices; j++)
	{	vector3<int> iGrot = mRot * basis.iGarr[j] - kpoint.offset;
		index->data[j] = e.gInfo.fullGindex(iGrot);
		if(needSuper)
			index->dataSuper[j] = gInfoSuper.fullGindex(ksuper + iGrot*super);
	}
	//Save to map:
	indexMap[kpoint] = index;
}

ColumnBundle WannierMinimizer::getWfns(const WannierMinimizer::Kpoint& kpoint, int iSpin) const
{	ColumnBundle ret(nBands, basis.nbasis, &basis, &kpoint, isGpuEnabled());
	ret.zero();
	axpyWfns(1., matrix(), kpoint, iSpin, ret);
	return ret;
}

void WannierMinimizer::axpyWfns(double alpha, const matrix& A, const WannierMinimizer::Kpoint& kpoint, int iSpin, ColumnBundle& result) const
{	static StopWatch watch("WannierMinimizer::accumWfns"); watch.start();
	//Figure out basis:
	const Index& index = *(indexMap.find(kpoint)->second);
	const int* indexData = (result.basis==&basisSuper) ? index.dataSuperPref : index.dataPref;
	//Pick source ColumnBundle:
	int q = kpoint.iReduced + iSpin*qCount;
	const ColumnBundle& Cin = e.eInfo.isMine(q) ? e.eVars.C[q] : Cother[q];
	assert(Cin);
	const ColumnBundle* C = &Cin;
	//Complex conjugate if inversion symmetry employed:
	ColumnBundle Cout;
	if(kpoint.invert < 0)
	{	Cout = *C;
		callPref(eblas_dscal)(Cout.nData(), -1., ((double*)Cout.dataPref())+1, 2); //negate the imaginary parts
		C = &Cout;
	}
	//Apply transformation if provided:
	if(A)
	{	Cout = (*C) * A;
		C = &Cout;
	}
	//Scatter from reduced basis to common basis with transformations:
	assert(C->nCols() == result.nCols());
	for(int b=0; b<C->nCols(); b++)
		callPref(eblas_scatter_zdaxpy)(index.nIndices, alpha, indexData, C->dataPref()+C->index(b,0), result.dataPref()+result.index(b,0));
	watch.stop();
}

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
{	ColumnBundle ret(nCenters, basis.nbasis, &basis, &kpoint, isGpuEnabled());
	ColumnBundle temp = ret.similar(1); //single column for intermediate computations
	#ifdef GPU_ENABLED
	vector3<>* pos; cudaMalloc(&pos, sizeof(vector3<>));
	#endif
	ret.zero();
	complex* retData = ret.dataPref();
	complex* tempData = temp.dataPref();
	for(const Wannier::TrialOrbital& t: wannier.trialOrbitals)
	{	for(const Wannier::AtomicOrbital& ao: t)
		{	const DOS::Weight::OrbitalDesc& od = ao.orbitalDesc;
			//--- Copy the center to GPU if necessary:
			#ifdef GPU_ENABLED
			cudaMemcpy(pos, &ao.r, sizeof(vector3<>), cudaMemcpyHostToDevice);
			#else
			const vector3<>* pos = &ao.r;
			#endif
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
			callPref(Vnl)(basis.nbasis, basis.nbasis, 1, od.l, od.m, kpoint.k, basis.iGarrPref, e.gInfo.G, pos, atRadial, tempData);
			if(ao.sp < 0) hRadial.free();
			//--- Accumulate to trial orbital:
			callPref(eblas_zaxpy)(basis.nbasis, ao.coeff * cis(0.5*M_PI*od.l)/e.gInfo.detR, tempData,1, retData,1);  //phase ensures odd l projectors are real
		}
		retData += basis.nbasis;
	}
	#ifdef GPU_ENABLED
	cudaFree(pos);
	#endif
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
