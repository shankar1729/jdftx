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
	mpiWorld->allReduce(result, MPIUtil::ReduceSum);
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


matrix WannierMinimizer::fixUnitary(const matrix& U, bool* isSingular)
{	return U * invsqrt(dagger(U) * U, 0, 0, isSingular);
}

bool WannierMinimizer::report(int iter)
{	//Check unitarity:
	bool needRestart = false;
	for(size_t ik=ikStart; ik<ikStop; ik++)
		if(nrm2(dagger(kMesh[ik].U) * kMesh[ik].U - eye(nCenters)) > 1e-6)
		{	needRestart = true;
			break;
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
			ki.U = ki.U1(0,nBands, 0,nCenters) * ki.U2;
		}
		mpiWorld->checkErrors(ossErr);
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

#include <iomanip>
std::vector<complex> diagDotRatio(const ColumnBundle& X, const ColumnBundle& Y)
{	assert(X.nCols()==Y.nCols());
	assert(X.basis==Y.basis);
	diagMatrix Xnorm = diagDot(X,X);
	diagMatrix Ynorm = diagDot(Y,Y);
	std::vector<complex> ret(X.nCols());
	const complex* Xdata = X.dataPref();
	const complex* Ydata = Y.dataPref();
	for(size_t b=0; b<ret.size(); b++)
		ret[b] = callPref(eblas_zdotc)(X.colLength(), Xdata+X.index(b,0),1, Ydata+Y.index(b,0),1) / sqrt(Xnorm[b]*Ynorm[b]);
	return ret;
}

void testVtransform(const Everything& e, std::vector<std::vector<int>> nProjArr)
{	int nSpinor = e.eInfo.spinorLength();
	QuantumNumber qnumIn;
	for(int iDir=0; iDir<3; iDir++)
		qnumIn.k[iDir] = Random::uniform();
	Basis basisIn; logSuspend(); basisIn.setup(e.gInfo, e.iInfo, e.cntrl.Ecut, qnumIn.k); logResume();
	ColumnBundle Cin(1, basisIn.nbasis*nSpinor, &basisIn, &qnumIn);
	const std::vector<SpaceGroupOp>& sym = e.symm.getMatrices();
	const std::vector<std::vector<std::vector<int> > >& atomMap = e.symm.getAtomMap();
	//Species loop:
	for(size_t iSp=0; iSp<e.iInfo.species.size(); iSp++)
	{	const SpeciesInfo& sp = *(e.iInfo.species[iSp]);
		logPrintf("\n\nSpecies %s:\n", sp.name.c_str());
		//Projectors at kIn:
		std::shared_ptr<ColumnBundle> Vin = sp.getV(Cin);
		//Loop over symmetry operations and inversions:
		for(unsigned iSym=0; iSym<sym.size(); iSym++)
		for(int invert=-1; invert<=+1; invert+=2)
		{	QuantumNumber qnumOut;
			qnumOut.k = invert * (qnumIn.k * sym[iSym].rot);
			vector3<int> kOffset;
			for(int iDir=0; iDir<3; iDir++)
			{	kOffset[iDir] = int(round(Random::normal()));
				qnumOut.k[iDir] -= kOffset[iDir];
			}
			logPrintf("\niSym: %u  invert: %+d  |a|: %6.3lf  |off|^2: %2d\n",
				iSym, invert, sym[iSym].a.length(), kOffset.length_squared());
			Basis basisOut; logSuspend(); basisOut.setup(e.gInfo, e.iInfo, e.cntrl.Ecut, qnumOut.k); logResume();
			ColumnBundle Cout(1, basisOut.nbasis*nSpinor, &basisOut, &qnumOut);
			//Initialize transformation:
			ColumnBundleTransform transform(qnumIn.k, basisIn, qnumOut.k, basisOut, 1, sym[iSym], invert);
			int nAtoms = sp.atpos.size();
			std::vector<vector3<int>> offsets(nAtoms);
			for(int atom=0; atom<nAtoms; atom++)
			{	int atomOut = atomMap[iSp][atom][iSym];
				const SpaceGroupOp& op = sym[iSym];
				offsets[atom] = round((op.rot * sp.atpos[atom] + op.a) - sp.atpos[atomOut]);
			}
			//Set up projector transformation matrix:
			int nProj = 0;
			for(int l=0; l<int(nProjArr[iSp].size()); l++)
				nProj += nProjArr[iSp][l]*(2*l+1);
			int nProjTot = nProj * nAtoms;
			matrix rot = zeroes(nProjTot, nProjTot);
			int nProjPrev = 0;
			double lSign = 1.;
			for(int l=0; l<int(nProjArr[iSp].size()); l++)
			{	matrix sym_l = e.symm.getSphericalMatrices(l, false)[iSym](0,2,4*l+2, 0,2,4*l+2); //projectors done in (l,m) not (j,mj)
				int nm = sym_l.nRows(); //= (2l + 1)
				//Set for each atom, accounting for atom mapping under symmetry:
				for(int p=0; p<nProjArr[iSp][l]; p++)
				{	for(int atom=0; atom<nAtoms; atom++)
					{	int atomOut = atomMap[iSp][atom][iSym];
						int pStart = atom*nProj + nProjPrev;
						int pStartOut = atomOut*nProj + nProjPrev;
						rot.set(pStartOut,pStartOut+nm, pStart,pStart+nm, lSign*sym_l);
					}
					nProjPrev += nm;
				}
				if(invert<0.) lSign = -lSign; //(-1)^l due to effect of inversion on Ylm
			}
			assert(nProjPrev == nProj);
			//Projectors directly at kOut:
			std::shared_ptr<ColumnBundle> VoutD = sp.getV(Cout);
			//Projectors transformed from kIn:
			ColumnBundle VoutT = VoutD->similar();
			VoutT.zero();
			transform.scatterAxpy(1., (*Vin), VoutT,0,1);
			VoutT = VoutT * rot;
			//Compare:
			std::vector<complex> ratio = diagDotRatio(*VoutD, VoutT);
			//--- measured phases:
			logPrintf("\tPhases: [");
			for(complex r: ratio)
				logPrintf(" %+f", r.arg()/(2*M_PI));
			logPrintf(" ]\n");
			//--- offset.kIn:
			logPrintf("\toffset.kIn:  [");
			for(vector3<int> offset: offsets)
			{	double d = invert*dot(qnumIn.k, offset);
				d -= floor(d+0.5);
				logPrintf(" %+f", d);
			}
			logPrintf(" ]\n");
			//--- offset.kOut:
			logPrintf("\toffset.kOut: [");
			for(vector3<int> offset: offsets)
			{	double d = invert*dot(qnumOut.k, offset);
				d -= floor(d+0.5);
				logPrintf(" %+f", d);
			}
			logPrintf(" ]\n");
		}
	}
	die("\nTesting.\n");
}

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
	{	VdagResult->clear();
		VdagResult->resize(VdagC->size());
		const std::vector<std::vector<std::vector<int> > >& atomMap = e.symm.getAtomMap();
		for(size_t iSp=0; iSp<VdagC->size(); iSp++)
			if(VdagC->at(iSp))
			{	const SpeciesInfo& sp = *(e.iInfo.species[iSp]);
				
				//VdagResult->at(iSp) = *sp.getV(result) ^ result; continue; //Bypass HACK
				
				//Determine phases due to atom offsets:
				int nAtoms = sp.atpos.size();
				std::vector<complex> phase(nAtoms);
				for(int atom=0; atom<nAtoms; atom++)
				{	int atomOut = atomMap[iSp][atom][kpoint.iSym];
					const SpaceGroupOp& op = sym[kpoint.iSym];
					vector3<int> offset = round((op.rot * sp.atpos[atom] + op.a) - sp.atpos[atomOut]);
					phase[atom] = cis(-2*M_PI*dot(C->qnum->k, offset));
				}
				//Set up projector transformation matrix:
				int nProjTot = VdagC->at(iSp).nRows();
				int nProj = nProjTot / nAtoms; //projectors per atom
				matrix rot = zeroes(nProjTot, nProjTot);
				int nProjPrev = 0;
				double lSign = 1.;
				for(int l=0; l<int(sp.VnlRadial.size()); l++)
				{	const matrix& sym_l = e.symm.getSphericalMatrices(l, false)[kpoint.iSym]; //projectors done in (l,m) not (j,mj)
					int nms = sym_l.nRows(); //= (2l + 1) * nSpinor
					//Set for each atom, accounting for atom mapping under symmetry:
					for(size_t p=0; p<sp.VnlRadial[l].size(); p++)
					{	for(int atom=0; atom<nAtoms; atom++)
						{	int atomOut = atomMap[iSp][atom][kpoint.iSym];
							int pStart = atom*nProj + nProjPrev;
							int pStartOut = atomOut*nProj + nProjPrev;
							rot.set(pStartOut,pStartOut+nms, pStart,pStart+nms, (lSign*phase[atom])*sym_l);
						}
						nProjPrev += nms;
					}
					if(kpoint.invert<0.) lSign = -lSign; //(-1)^l due to effect of inversion on Ylm
				}
				assert(nProjPrev == nProj);
				//Account for spinor rotations:
				if(nSpinor > 1)
				{	matrix spinorRotDag = (kpoint.invert<0) ? transpose(transform.spinorRot) : dagger(transform.spinorRot);
					rot = tiledBlockMatrix(spinorRotDag, nProjTot/nSpinor) * rot;
				}
				//Apply transform
				VdagResult->at(iSp) = dagger(rot) * VdagC->at(iSp); //apply rotation
				if(kpoint.invert < 0) VdagResult->at(iSp) = conj(VdagResult->at(iSp));
				
				//HACK
				matrix directAll = (*sp.getV(result) ^ result);
				for(int atom=0; atom<nAtoms; atom++)
				{	matrix transformed = VdagResult->at(iSp)(atom*nProj,(atom+1)*nProj, 0, result.nCols());
					matrix direct = directAll(atom*nProj,(atom+1)*nProj, 0, result.nCols());
					double err = nrm2(transformed-direct)/nrm2(direct);
					double phaseSum=0., magSum=0., wSum = 0.;
					for(unsigned i=0; i<direct.nData(); i++)
					{	double w = direct.data()[i].abs();
						if(w > 1e-8)
						{	complex ratio = transformed.data()[i]/direct.data()[i];
							double mag=ratio.abs(); magSum+=w*mag;
							double phase=ratio.arg()/(2*M_PI); phaseSum+=w*phase;
							wSum += w;
						}
					}
					double magMean = magSum/wSum;
					double phaseMean = phaseSum/wSum;
					double phaseErr = nrm2(direct*cis(2*M_PI*phaseMean)-transformed)/nrm2(direct);
					if(err > 1e-14)
						printf("sp: %2s  atom: %d  iSym: %2d  invert: %+d  |a|: %6.3lf  |off|^2: %2d  err: %.2le  ratio: %lg  phase: %+.4lf  phaseErr: %.2le\n",
							sp.name.c_str(), atom, kpoint.iSym, kpoint.invert, sym[kpoint.iSym].a.length(), kpoint.offset.length_squared(), err, magMean, phaseMean, phaseErr);
				}
			}
		
		//HACK
		std::vector<std::vector<int>> nProjArr;
		for(const auto& sp: e.iInfo.species)
		{	std::vector<int> nProj_sp;
			for(const auto& Vl: sp->VnlRadial)
				nProj_sp.push_back(Vl.size());
			nProjArr.push_back(nProj_sp);
		}
		//testVtransform(e, nProjArr); 
	}
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
	//Augment at each species:
	for(size_t iSp=0; iSp<iInfo.species.size(); iSp++)
	{	const SpeciesInfo& sp = *(iInfo.species[iSp]);
		if(!sp.Qint.size()) continue; //no augmentation
		//Phases for each atom:
		std::vector<complex> phaseArr;
		for(vector3<> x: sp.atpos)
			phaseArr.push_back(cis(-2*M_PI*dot(dkVec,x)));
		//Augment the overlap
		matrix VdagC1 = VdagC1ptr ? VdagC1ptr->at(iSp) : (*sp.getV(C1)) ^ C1;
		matrix VdagC2 = VdagC2ptr ? VdagC2ptr->at(iSp) : (*sp.getV(C2)) ^ C2;
		ret += dagger(VdagC1) * (tiledBlockMatrix(sp.QintAll, sp.atpos.size(),&phaseArr) * VdagC2);
	}
	watch.stop();
	return ret;
}
