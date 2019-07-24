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
	B1.resize(wmin->kMesh.size());
	B2.resize(wmin->kMesh.size());
}
size_t WannierGradient::ikStart() const { return wmin->ikStart; }
size_t WannierGradient::ikStop() const { return wmin->ikStop; }

//---- linear algebra functions required by Minimizable<WannierGradient> -----

WannierGradient clone(const WannierGradient& grad) { return grad; }
double dot(const WannierGradient& x, const WannierGradient& y)
{	assert(x.wmin == y.wmin);
	double result = 0.;
	for(unsigned ik=x.ikStart(); ik<x.ikStop(); ik++)
	{	if(x.B1[ik]) result += dotc(x.B1[ik], y.B1[ik]).real() * 2.; //off-diagonal blocks of a Hermitian matrix
		result += dotc(x.B2[ik], y.B2[ik]).real(); //diagonal block; always present
	}
	mpiWorld->allReduce(result, MPIUtil::ReduceSum);
	return result;
}
WannierGradient& operator*=(WannierGradient& x, double alpha)
{	for(unsigned ik=x.ikStart(); ik<x.ikStop(); ik++)
	{	if(x.B1[ik]) x.B1[ik] *= alpha;
		x.B2[ik] *= alpha;
	}
	return x;
}
void axpy(double alpha, const WannierGradient& x, WannierGradient& y)
{	assert(x.wmin == y.wmin);
	for(unsigned ik=x.ikStart(); ik<x.ikStop(); ik++)
	{	if(x.B1[ik]) axpy(alpha, x.B1[ik], y.B1[ik]);
		axpy(alpha, x.B2[ik], y.B2[ik]);
	}
}
void randomize(WannierGradient& x)
{	for(unsigned ik=x.ikStart(); ik<x.ikStop(); ik++)
	{	if(x.B1[ik]) randomize(x.B1[ik]); //off-diagonal block of Hermitian
		if(x.B2[ik]) { randomize(x.B2[ik]); x.B2[ik] = dagger_symmetrize(x.B2[ik]); } //diagonal block
	}
}

//---- energy/gradient functions required by Minimizable<WannierGradient> -----

void WannierMinimizer::step(const WannierGradient& grad, double alpha)
{	static StopWatch watch("WannierMinimizer::step"); watch.start();
	assert(grad.wmin == this);
	for(unsigned ik=ikStart; ik<ikStop; ik++)
	{	KmeshEntry& ki = kMesh[ik];
		//Stage 1:
		if(ki.nIn > nCenters)
		{	int nFreeIn = ki.nIn - ki.nFixed;
			int nFree = nCenters - ki.nFixed;
			matrix B1 = zeroes(nFreeIn, nFreeIn);
			B1.set(0,nFree, nFree,nFreeIn, alpha*grad.B1[ik]);
			B1.set(nFree,nFreeIn, 0,nFree, alpha*dagger(grad.B1[ik]));
			ki.U1.set(0,nBands, ki.nFixed,ki.nIn, ki.U1(0,nBands, ki.nFixed,ki.nIn) * cis(B1));
		}
		//Stage 2:
		ki.U2.set(nFrozen,nCenters, nFrozen,nCenters, cis(alpha*grad.B2[ik]) * ki.U2(nFrozen,nCenters, nFrozen,nCenters));
		//Net rotation:
		ki.U = ki.U1 * ki.U2;
	}
	watch.stop();
	bcastU(); //Make U available on all processes that need it
}

void WannierMinimizer::bcastU()
{	static StopWatch watch("WannierMinimizer::bcastU"); watch.start();
	std::vector<MPIUtil::Request> requests;
	for(KmeshEntry& ki: kMesh)
		if(ki.mpi)
		{	MPIUtil::Request request;
			ki.mpi->bcastData(ki.U, 0, &request); 
			requests.push_back(request);
		}
	MPIUtil::waitAll(requests);
	watch.stop();
}


double WannierMinimizer::compute(WannierGradient* grad, WannierGradient* Kgrad)
{	static StopWatch watch("WannierMinimizer::compute");
	if(grad) for(KmeshEntry& ki: kMesh) if(ki.U) ki.Omega_UdotU = zeroes(nCenters, ki.nIn); //Clear gradient
	
	double Omega = getOmega(grad);
	
	//Collect Omega_U and propagate to Omega_B if necessary:
	if(grad)
	{	watch.start();
		std::vector<MPIUtil::Request> requests;
		for(KmeshEntry& ki: kMesh)
			if(ki.mpi)
			{	MPIUtil::Request request;
				ki.mpi->reduceData(ki.Omega_UdotU, MPIUtil::ReduceSum, 0, &request); //Collect Omega_U
				requests.push_back(request);
			}
		MPIUtil::waitAll(requests);
		grad->init(this);
		for(size_t ik=ikStart; ik<ikStop; ik++)
		{	KmeshEntry& ki = kMesh[ik];
			matrix Omega_B = complex(0,0.5)*(ki.U2(0,ki.nIn, 0,nCenters) * ki.Omega_UdotU * dagger(ki.U2));
			if(ki.nIn > nCenters) //Stage 1:
				grad->B1[ik] = Omega_B(ki.nFixed,nCenters, nCenters,ki.nIn);
			matrix Omega_B2_hlf = Omega_B(nFrozen,nCenters, nFrozen,nCenters);
			grad->B2[ik] = Omega_B2_hlf + dagger(Omega_B2_hlf);
		}
		if(Kgrad) *Kgrad = precondition(*grad);
		watch.stop();
	}
	return Omega;
}

void WannierMinimizer::constrain(WannierGradient& grad)
{	for(size_t ik=ikStart; ik<ikStop; ik++)
		grad.B2[ik] = dagger_symmetrize(grad.B2[ik]);
}

WannierGradient WannierMinimizer::precondition(const WannierGradient& grad)
{	WannierGradient Kgrad = grad;
	constrain(Kgrad);
	return Kgrad;
}


matrix WannierMinimizer::fixUnitary(const matrix& U, bool* isSingular)
{	return U * invsqrt(dagger(U) * U, 0, 0, isSingular);
}

bool WannierMinimizer::report(int iter)
{	//Check unitarity:
	bool needRestart = false;
	for(size_t ik=ikStart; ik<ikStop; ik++)
	{	KmeshEntry& ki = kMesh[ik];
		if(nrm2(dagger(ki.U) * ki.U - eye(ki.nIn)) > 1e-6)
		{	needRestart = true;
			break;
		}
	}
	mpiWorld->allReduce(needRestart, MPIUtil::ReduceLOr);
	if(needRestart)
	{	logPrintf("%s\tUpdating rotations to enforce unitarity\n", wannier.minParams.linePrefix);
		ostringstream ossErr;
		for(size_t ik=ikStart; ik<ikStop; ik++)
		{	KmeshEntry& ki = kMesh[ik];
			bool isSingular = false;
			ki.U1 = fixUnitary(ki.U1, &isSingular);
			ki.U2 = fixUnitary(ki.U2, &isSingular);
			if(isSingular)
			{	ossErr << "Unitary rotations are singular at k = [ "
					<< ki.point.k[0] << ' ' << ki.point.k[1] << ' ' << ki.point.k[2] << " ]\n";
				break;
			}
			ki.U = ki.U1 * ki.U2;
		}
		mpiWorld->checkErrors(ossErr);
		bcastU();
		return true;
	}
    return false;
}

double WannierMinimizer::sync(double x) const
{	mpiWorld->bcast(x);
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

ColumnBundle WannierMinimizer::getWfns(const WannierMinimizer::Kpoint& kpoint, int iSpin, std::vector<matrix>* VdagResult) const
{	ColumnBundle ret(nBands, basis.nbasis*nSpinor, &basis, &kpoint, isGpuEnabled());
	ret.zero();
	axpyWfns(1., matrix(), kpoint, iSpin, ret, VdagResult);
	return ret;
}

#define axpyWfns_COMMON(result) \
	/* Pick transform: */ \
	const ColumnBundleTransform& transform = *(((result.basis==&basisSuper) ? transformMapSuper : transformMap).find(kpoint)->second); \
	/* Pick source ColumnBundle: */ \
	int q = kpoint.iReduced + iSpin*qCount; \
	const ColumnBundle* C = e.eInfo.isMine(q) ? &e.eVars.C[q] : &Cother[q]; \
	assert(*C);

void WannierMinimizer::axpyWfns(double alpha, const matrix& A, const WannierMinimizer::Kpoint& kpoint, int iSpin, ColumnBundle& result, std::vector<matrix>* VdagResult) const
{	static StopWatch watch("WannierMinimizer::axpyWfns"); watch.start();
	axpyWfns_COMMON(result)
	const std::vector<matrix>* VdagC = VdagResult ? (e.eInfo.isMine(q) ? &e.eVars.VdagC[q] : &VdagCother[q]) : 0;
	//Apply transformation if provided:
	ColumnBundle Cout; std::vector<matrix> VdagCout;
	if(A)
	{	matrix Astar = (kpoint.invert<0 ? conj(A) : A);
		Cout = (*C) * Astar;
		C = &Cout;
		//Similarly for projections, if needed:
		if(VdagResult)
		{	VdagCout.resize(VdagC->size());
			for(size_t iSp=0; iSp<VdagC->size(); iSp++)
				if(VdagC->at(iSp))
					VdagCout[iSp] = VdagC->at(iSp) * Astar;
			VdagC = &VdagCout;
		}
	}
	//Scatter from reduced basis to common basis with transformations:
	assert(C->nCols() == result.nCols());
	transform.scatterAxpy(alpha, *C, result,0,1);
	//Corresponding transformation in projections:
	if(VdagResult)
		*VdagResult = transform.transformVdagC(*VdagC, kpoint.iSym);
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

//Gaussian orbital of specified width and angular momentum
inline double gaussTilde(double G, double sigma, int l, double normPrefac)
{	double Gsigma = G*sigma;
	return normPrefac * std::pow(Gsigma,l) * exp(-0.5*Gsigma*Gsigma);
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
			//Handle atomic orbitals:
			const DOS::Weight::OrbitalDesc& od = ao.orbitalDesc;
			complex lPhase =  cis(0.5*M_PI*od.l); //including this phase ensures odd l projectors are real (i^l term in spherical wave expansion)
			if(ao.sp >= 0)
			{	int iCol = e.iInfo.species[ao.sp]->atomicOrbitalOffset(ao.atom, od.n, od.l, od.m, od.s);
				callPref(eblas_zaxpy)(ret.colLength(), ao.coeff*lPhase, psiAtomic[ao.sp].dataPref()+iCol*ret.colLength(),1, retData,1);
				continue;
			}
			//Gaussian orbital:
			if(ao.sigma > 0.)
			{	//--- Copy the center to managed memory:
				ManagedArray<vector3<>> pos(&ao.r, 1);
				//--- Get / create the radial part:
				RadialFunctionG hRadial;
				double Al = 0.25*sqrt(M_PI);
				for(int p=1; p<=od.l; p++)
					Al *= (p+0.5);
				double sigma = (od.n+1) * ao.sigma;
				double normPrefac = sqrt(std::pow(2*M_PI*sigma,3) / Al);
				hRadial.init(od.l, 0.02, e.gInfo.GmaxSphere, gaussTilde, sigma, od.l, normPrefac);
				//--- Initialize the projector:
				assert(od.s < nSpinor);
				if(nSpinor > 1) { temp.zero(); assert(od.spinType==SpinZ); }
				callPref(Vnl)(basis.nbasis, basis.nbasis, 1, od.l, od.m, kpoint.k, basis.iGarr.dataPref(), e.gInfo.G, pos.dataPref(), hRadial, temp.dataPref()+od.s*basis.nbasis);
				hRadial.free();
				//--- Accumulate to trial orbital:
				callPref(eblas_zaxpy)(ret.colLength(), ao.coeff*lPhase/e.gInfo.detR, temp.dataPref(),1, retData,1);
				continue;
			}
			assert(!"Orbital was neither Numerical, Gaussian nor Atomic.");
		}
		retData += ret.colLength();
	}
	return ret;
}

matrix WannierMinimizer::overlap(const ColumnBundle& C1, const ColumnBundle& C2, const std::vector<matrix>* VdagC1ptr, const std::vector<matrix>* VdagC2ptr) const
{	static StopWatch watch("WannierMinimizer::overlap"); watch.start();
	const GridInfo& gInfo = *(C1.basis->gInfo);
	const IonInfo& iInfo = *(C1.basis->iInfo);
	matrix ret = gInfo.detR * (C1 ^ C2);
	//k-point difference:
	vector3<> dkVec = C2.qnum->k - C1.qnum->k;
	double dk = sqrt(gInfo.GGT.metric_length_squared(dkVec));
	vector3<> dkHat = gInfo.GT * dkVec * (dk ? 1.0/dk : 0.0); //the unit Vector along dkVec (set dkHat to 0 for dk=0 (doesn't matter))
	//Augment at each species:
	for(size_t iSp=0; iSp<iInfo.species.size(); iSp++)
	{	const SpeciesInfo& sp = *(iInfo.species[iSp]);
		if(!sp.isUltrasoft()) continue; //no augmentation
		//Create the Q matrix appropriate for current k-point difference:
		matrix Qk = zeroes(sp.QintAll.nRows(), sp.QintAll.nCols());
		complex* QkData = Qk.data();
		int i1 = 0;
		for(int l1=0; l1<int(sp.VnlRadial.size()); l1++)
		for(int p1=0; p1<int(sp.VnlRadial[l1].size()); p1++)
		for(int m1=-l1; m1<=l1; m1++)
		{	//Triple loop over second projector:
			int i2 = 0;
			for(int l2=0; l2<int(sp.VnlRadial.size()); l2++)
			for(int p2=0; p2<int(sp.VnlRadial[l2].size()); p2++)
			for(int m2=-l2; m2<=l2; m2++)
			{	std::vector<YlmProdTerm> terms = expandYlmProd(l1,m1, l2,m2);
				complex q12 = 0.;
				for(const YlmProdTerm& term: terms)
				{	SpeciesInfo::QijIndex qIndex = { l1, p1, l2, p2, term.l };
					auto Qijl = sp.Qradial.find(qIndex);
					if(Qijl==sp.Qradial.end()) continue; //no entry at this l
					q12 += term.coeff * cis(0.5*M_PI*(l2-l1-term.l)) * Ylm(term.l,term.m, dkHat) * Qijl->second(dk);
				}
				for(int s=0; s<nSpinor; s++)
					QkData[Qk.index(i1+s,i2+s)] = q12;
				i2 += nSpinor;
			}
			i1 += nSpinor;
		}
		if(sp.isRelativistic()) Qk = sp.fljAll * Qk * sp.fljAll;
		//Phases for each atom:
		std::vector<complex> phaseArr;
		for(vector3<> x: sp.atpos)
			phaseArr.push_back(cis(-2*M_PI*dot(dkVec,x)));
		//Augment the overlap
		matrix VdagC1 = VdagC1ptr ? VdagC1ptr->at(iSp) : (*sp.getV(C1)) ^ C1;
		matrix VdagC2 = VdagC2ptr ? VdagC2ptr->at(iSp) : (*sp.getV(C2)) ^ C2;
		ret += dagger(VdagC1) * (tiledBlockMatrix(Qk, sp.atpos.size(), &phaseArr) * VdagC2);
	}
	watch.stop();
	return ret;
}

void WannierMinimizer::dumpWannierized(const matrix& Htilde, const matrix& phase, string varName, bool realPartOnly, int iSpin) const
{
	string fname = wannier.getFilename(Wannier::FilenameDump, varName, &iSpin);
	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
	FILE* fp = 0;
	if(mpiWorld->isHead())
	{	fp = fopen(fname.c_str(), "w");
		if(!fp) die_alone("could not open file for writing.\n");
	}
	//Determine block size:
	int nCells = phase.nCols();
	int blockSize = ceildiv(nCells, mpiWorld->nProcesses()); //so that memory before and after FT roughly similar
	int nBlocks = ceildiv(nCells, blockSize);
	//Loop over blocks:
	int iCellStart = 0;
	double nrm2totSq = 0., nrm2imSq = 0.;
	for(int iBlock=0; iBlock<nBlocks; iBlock++)
	{	int iCellStop = std::min(iCellStart+blockSize, nCells);
		matrix Hblock = Htilde * phase(0,phase.nRows(), iCellStart,iCellStop);
		mpiWorld->reduceData(Hblock, MPIUtil::ReduceSum);
		if(mpiWorld->isHead())
		{	//Write to file:
			if(realPartOnly)
			{	nrm2totSq += std::pow(nrm2(Hblock), 2); 
				nrm2imSq += std::pow(callPref(eblas_dnrm2)(Hblock.nData(), ((double*)Hblock.dataPref())+1, 2), 2); //imaginary parts with a stride of 2
				Hblock.write_real(fp);
			}
			else Hblock.write(fp);
		}
		iCellStart = iCellStop;
	}
	if(mpiWorld->isHead()) fclose(fp);
	if(realPartOnly)
	{	mpiWorld->bcast(nrm2totSq);
		mpiWorld->bcast(nrm2imSq);
		logPrintf("done. Relative discarded imaginary part: %le\n", sqrt(nrm2imSq / nrm2totSq));
	}
	else
		logPrintf("done.\n");
}
