/*
 * PerturbationSolver.cpp
 *
 *  Created on: Jul 22, 2022
 *	  Author: brandon
 */

#include <electronic/PerturbationSolver.h>
#include <electronic/TestPerturbation.h>
#include <electronic/Everything.h>
#include <electronic/ElecMinimizer.h>
#include <electronic/ColumnBundle.h>
#include <electronic/ExCorr.h>
#include <core/matrix.h>
#include <core/Units.h>
#include <core/ScalarField.h>
#include <core/ScalarFieldIO.h>
#include <core/Operators.h>
#include <core/MPIUtil.h>
#include <cstdio>
#include <cmath>
#include <core/Operators.h>


void PerturbationGradient::init(const Everything& e)
{	eInfo = &e.eInfo;
	pInfo = &e.vptInfo;
	if (pInfo->incommensurate)
		dY.resize(eInfo->nStates*2);
	else
		dY.resize(eInfo->nStates);
}

PerturbationGradient& PerturbationGradient::operator*=(double alpha)
{
	int nStates = eInfo->nStates;

	for(int q=eInfo->qStart; q<eInfo->qStop; q++)
	{
		if(dY[q]) dY[q] *= alpha;
		if (pInfo->incommensurate) if(dY[q+nStates]) dY[q+nStates] *= alpha;
	}

	return *this;
}

void axpy(double alpha, const PerturbationGradient& x, PerturbationGradient& y)
{	assert(x.eInfo == y.eInfo);

	int nStates = x.eInfo->nStates;

	for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++)
	{
		if(x.dY[q]) { if(y.dY[q]) axpy(alpha, x.dY[q], y.dY[q]); else y.dY[q] = alpha*x.dY[q]; }
		if (x.pInfo->incommensurate) if(x.dY[q+nStates]) { if(y.dY[q+nStates]) axpy(alpha, x.dY[q+nStates], y.dY[q+nStates]); else y.dY[q+nStates] = alpha*x.dY[q+nStates]; }
	}
}

double dot(const PerturbationGradient& x, const PerturbationGradient& y)
{	assert(x.eInfo == y.eInfo);

	int nStates = x.eInfo->nStates;
	double result = 0.0;

	for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++)
	{
		if(x.dY[q] && y.dY[q]) result += dotc(x.dY[q], y.dY[q]).real()*2.0;
		if (x.pInfo->incommensurate) if(x.dY[q+nStates] && y.dY[q+nStates]) result += dotc(x.dY[q+nStates], y.dY[q+nStates]).real()*2.0;
	}

	mpiWorld->allReduce(result, MPIUtil::ReduceSum);
	return result;
}

double dot(const std::vector<ColumnBundle>& x, const std::vector<ColumnBundle>& y, PerturbationInfo* pInfo, ElecInfo *eInfo)
{
	int nStates = eInfo->nStates;
	double result = 0.0;

	for(int q=eInfo->qStart; q<eInfo->qStop; q++)
	{
		if(x[q] && y[q]) result += dotc(x[q], y[q]).real()*2.0;
		if (pInfo->incommensurate) if(x[q+nStates] && y[q+nStates]) result += dotc(x[q+nStates], y[q+nStates]).real()*2.0;
	}

	mpiWorld->allReduce(result, MPIUtil::ReduceSum);
	return result;
}

PerturbationGradient clone(const PerturbationGradient& x)
{	return x;
}

void randomize(PerturbationGradient& x)
{
	randomize(x.dY, *x.eInfo);
}

ScalarField randomRealSpaceVector(ColumnBundle& ref, const GridInfo* gInfo) {
	ColumnBundle res = ref;
	res.randomize(0, 1);
	ScalarField vec = Real(I(res.getColumn(0, 0)));
	ScalarFieldTilde vtilde = J(vec);
	vtilde->setGzero(0);
	return changeGrid(I(vtilde), *gInfo);
}

void PerturbationSolver::printCB(ColumnBundle C) {
	double *dat = (double*)(C.getColumn(0, 0)->dataPref());
	logPrintf("ColumnBundle values %g %g %g %g %g %g %g %g %g %g\n", dat[0], dat[1], dat[2], dat[3], dat[4], dat[5], dat[6], dat[7], dat[8], dat[9]);
}

void PerturbationSolver::printM(matrix C) {
	//logPrintf("Matrix values %g %g %g %g %g %g %g %g %g\n", C(0,0).x, C(0,1).x, C(0,2).x, C(1,0).x,C(1,1).x,C(1,2).x,C(2,0).x,C(2,1).x,C(2,2).x);
	logPrintf("Matrix values %g %g %g %g\n", C(0,0).x, C(0,1).x, C(1,0).x,C(1,1).x);
}

void PerturbationSolver::printV(ScalarField V) {
	double *dat = V->dataPref();
	logPrintf("Vector values %g %g %g %g %g %g %g %g %g %g\n", dat[0], dat[1], dat[2], dat[3], dat[4], dat[5], dat[6], dat[7], dat[8], dat[9]);
}

void PerturbationSolver::printV(ScalarFieldTilde V) {
	double *dat = (double*)(V->dataPref());
	logPrintf("Vector values %g %g %g %g %g %g %g %g %g %g\n", dat[0], dat[1], dat[2], dat[3], dat[4], dat[5], dat[6], dat[7], dat[8], dat[9]);
}

void printVvpt(ScalarField V) {
	double *dat = V->dataPref();
	logPrintf("Vector values %g %g %g %g %g %g %g %g %g %g\n", dat[0], dat[1], dat[2], dat[3], dat[4], dat[5], dat[6], dat[7], dat[8], dat[9]);
}

ScalarFieldArray PerturbationSolver::getdn(const std::vector<ColumnBundle>* dC, const std::vector<ColumnBundle>* C)
{	static StopWatch watch("getdn"); watch.start();
	ScalarFieldArray density(eVars.n.size());
	nullToZero(density, e.gInfo);

	//Runs over all states and accumulates density to the corresponding spin channel of the total density
	e.iInfo.augmentDensityInit();
	
	std::vector<matrix> VdagdCq;

	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{
		if (ultrasoftDensityAugRequired) {
			VdagdCq.resize(e.iInfo.species.size());
			e.iInfo.project((*dC)[q], VdagdCq);
		}
		
		if (C) {
			std::vector<matrix> VdagCq(e.iInfo.species.size());
			e.iInfo.project((*C)[q], VdagCq); //update the atomic projections
			
			density += 2.0*e.eInfo.qnums[q].weight * Real(diagouterI(eVars.F[q], (*C)[q], (*dC)[q], density.size(), &e.gInfo));
			e.iInfo.augmentDensitySpherical(e.eInfo.qnums[q], eVars.F[q], VdagCq, &VdagdCq);
		} else {
			density += 2.0*e.eInfo.qnums[q].weight * Real(diagouterI(eVars.F[q], eVars.C[q], (*dC)[q], density.size(), &e.gInfo));
			e.iInfo.augmentDensitySpherical(e.eInfo.qnums[q], eVars.F[q], eVars.VdagC[q], &VdagdCq);
		}
	}
	e.iInfo.augmentDensityGrid(density);
	for(ScalarField& ns: density)
	{	nullToZero(ns, e.gInfo);
		ns->allReduceData(mpiWorld, MPIUtil::ReduceSum);
	}
	e.symm.symmetrize(density);
	watch.stop();
	return density;
}


ScalarFieldArray PerturbationSolver::getdnatom()
{	
	/*auto m = pInfo.atomdisplacement;
	auto sp = iInfo.species[m.sp];
	ScalarFieldArray n1 = getn(eVars.C);
	vector3<> pos1 = sp->atpos[m.at];
	vector3<> pos2 = pos1 + e.gInfo.invR*m.dir;
	
	sp->atpos[m.at] = pos2;
	sp->sync_atpos();
	e.iInfo.update(e.ener);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		e.iInfo.project(eVars.C[q], eVars.VdagC[q]);
	}
	
	ScalarFieldArray n2 = getn(eVars.C);
	
	sp->atpos[m.at] = pos1;
	sp->sync_atpos();
	e.iInfo.update(e.ener);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		e.iInfo.project(eVars.C[q], eVars.VdagC[q]);
	}
	
	return (n2 - n1);*/
	
	static StopWatch watch("getdn"); watch.start();
	ScalarFieldArray density(eVars.n.size());
	
	nullToZero(density, e.gInfo);
	
	if (ultrasoftDensityAugRequired && pInfo.atposPerturbationExists) {
		PerturbationInfo::Mode m = pInfo.atomdisplacement;
		auto sp = iInfo.species[m.sp];
		//Runs over all states and accumulates density to the corresponding spin channel of the total density
		sp->augmentDensityInit(m.at);
		
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		{
			sp->augmentDensitySpherical(e.eInfo.qnums[q], eVars.F[q], pInfo.VdagCatom_cached[q], &pInfo.dVdagCatom_cached[q], 0, m.at);
		}
		sp->augmentDensityGrid(density, m.at);
		
		
		sp->augmentDensityInit(m.at);
		
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		{
			sp->augmentDensitySpherical(e.eInfo.qnums[q], eVars.F[q], pInfo.VdagCatom_cached[q], 0, 0, m.at);
		}
		sp->augmentDensityGrid(density, m.at, &m.dirLattice);
		
		for(ScalarField& ns: density)
		{	nullToZero(ns, e.gInfo);
			ns->allReduceData(mpiWorld, MPIUtil::ReduceSum);
		}
		e.symm.symmetrize(density);
	}
	watch.stop();
	return density;
}

//TODO Add nonlocal density aug terms in
void PerturbationSolver::getdnInc(const std::vector<ColumnBundle>* dC, const std::vector<ColumnBundle>* C, complexScalarFieldArray& dnpq, complexScalarFieldArray& dnmq)
{   static StopWatch watch("getdnInc"); watch.start();
	ScalarFieldArray dnpq_aug(dnpq.size()), dnmq_aug(dnmq.size());
	for (unsigned int s = 0; s < dnpq_aug.size(); s++) {
		nullToZero(dnpq_aug[s], e.gInfo);
		nullToZero(dnpq[s], e.gInfo);
		dnpq[s] = 0*dnpq[s];
	}
	for (unsigned int s = 0; s < dnmq_aug.size(); s++) {
		nullToZero(dnmq_aug[s], e.gInfo);
		nullToZero(dnmq[s], e.gInfo);
		dnmq[s] = 0*dnmq[s];
	}

	if (C)
		die("Density augmentation terms must be added to getdnInc.");

	
	std::vector<matrix> VdagdCTinvk_k, VdagdCTk_k;
	
	if (ultrasoftDensityAugRequired) {
		VdagdCTinvk_k.resize(e.iInfo.species.size());
		VdagdCTk_k.resize(e.iInfo.species.size());
		
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		{
			int Tk_k = q;
			int Tinvk_k = q+e.eInfo.nStates;
			
			e.iInfo.project((*dC)[Tinvk_k], VdagdCTinvk_k); //update the atomic projections
			e.iInfo.project((*dC)[Tk_k], VdagdCTk_k); //update the atomic projections
		}
	}
	
	//Runs over all states and accumulates density to the corresponding spin channel of the total density
	e.iInfo.augmentDensityInit();
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{
		int Tk_k = q;
		int Tinvk_k = q+e.eInfo.nStates;

		dnpq += e.eInfo.qnums[q].weight * diagouterI(eVars.F[q], eVars.C[q], (*dC)[Tinvk_k], dnpq.size(), &e.gInfo);
		dnpq += e.eInfo.qnums[q].weight * diagouterI(eVars.F[q], (*dC)[Tk_k], eVars.C[q], dnpq.size(), &e.gInfo);
		e.iInfo.augmentDensitySpherical(e.eInfo.qnums[q], eVars.F[q], eVars.VdagC[q], &VdagdCTk_k, &VdagdCTinvk_k); //pseudopotential contribution
	}

	e.iInfo.augmentDensityGrid(dnpq_aug);

	//logPrintf("testt");
	//printV(dnpq_aug[0]);

	//TODO ADD BACK
	//dnpq += Complex(dnpq_aug);

	for(complexScalarField& ns: dnpq)
	{	nullToZero(ns, e.gInfo);
		ns->allReduceData(mpiWorld, MPIUtil::ReduceSum);
	}

	//Redundant calculation since dnpq = dnmq* in theory
	/*
	e.iInfo.augmentDensityInit();
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{
		int Tk_k = q;
		int Tinvk_k = q+e.eInfo.nStates;

		dnmq += e.eInfo.qnums[q].weight * diagouterI(eVars.F[q], eVars.C[q], (*dC)[Tk_k], dnmq.size(), &e.gInfo);
		dnmq += e.eInfo.qnums[q].weight * diagouterI(eVars.F[q], (*dC)[Tinvk_k], eVars.C[q], dnmq.size(), &e.gInfo);
		e.iInfo.augmentDensitySpherical(e.eInfo.qnums[q], eVars.F[q], eVars.VdagC[q], &VdagdCTinvk_k, &VdagdCTk_k); //pseudopotential contribution
	}

	e.iInfo.augmentDensityGrid(dnmq_aug);

	//TODO ADD BACK
	//dnmq += Complex(dnmq_aug);
	for(complexScalarField& ns: dnmq)
	{	nullToZero(ns, e.gInfo);
		ns->allReduceData(mpiWorld, MPIUtil::ReduceSum);
	}*/

	e.symm.symmetrize(dnpq);
	//e.symm.symmetrize(dnmq);
	
	dnmq = conj(dnpq);
	
	watch.stop();
}


ScalarFieldArray PerturbationSolver::getn(const std::vector<ColumnBundle>& C)
{	static StopWatch watch("getn"); watch.start();
	ScalarFieldArray density(eVars.n.size());
	nullToZero(density, e.gInfo);

	//Runs over all states and accumulates density to the corresponding spin channel of the total density
	e.iInfo.augmentDensityInit();
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{
		std::vector<matrix> VdagCq(e.iInfo.species.size());
		e.iInfo.project(C[q], VdagCq); //update the atomic projections

		density += e.eInfo.qnums[q].weight * diagouterI(eVars.F[q], C[q], density.size(), &e.gInfo);
		e.iInfo.augmentDensitySpherical(e.eInfo.qnums[q], eVars.F[q], VdagCq); //pseudopotential contribution
	} //TODO get rid of eVars.VdagC and compute nl pp from C
	e.iInfo.augmentDensityGrid(density);
	for(ScalarField& ns: density)
	{	nullToZero(ns, e.gInfo);
		ns->allReduceData(mpiWorld, MPIUtil::ReduceSum);
	}
	e.symm.symmetrize(density);
	watch.stop();
	return density;
}


PerturbationSolver::PerturbationSolver(Everything& e) : e(e), eVars(e.eVars), eInfo(e.eInfo), iInfo(e.iInfo), pInfo(e.vptInfo) {}

double PerturbationSolver::compute(PerturbationGradient* grad, PerturbationGradient* Kgrad)
{
	if(grad) grad->init(e);
	if(Kgrad) Kgrad->init(e);

	if (!pInfo.incommensurate) {
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			//pInfo.dU[q] = (pInfo.dY[q]^O(eVars.C[q])) + (eVars.C[q]^O(pInfo.dY[q]));
			pInfo.dU[q] = 2*dagger_symmetrize(pInfo.dY[q]^pInfo.OC[q]);
			pInfo.dC[q] = pInfo.dY[q];
			pInfo.dC[q] -= 0.5*(eVars.C[q]*pInfo.dU[q]);
		}
	} else {
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			int Tk_k = q;
			int Tinvk_k = q+e.eInfo.nStates;
			int Tk = q ;
			int Tinvk = q + e.eInfo.nStates;

			pInfo.dC[Tk_k] = pInfo.dY[Tk_k] - O(pInfo.Cinc[Tk]*(pInfo.Cinc[Tk]^pInfo.dY[Tk_k]));
			pInfo.dC[Tinvk_k] = pInfo.dY[Tinvk_k] - O(pInfo.Cinc[Tinvk]*(pInfo.Cinc[Tinvk]^pInfo.dY[Tinvk_k]));
		}
	}

	if (!pInfo.incommensurate) {
		pInfo.dn = getdn(&pInfo.dC);
		getdVsclocPsi(pInfo.dn, pInfo.dVscloc);
		e.iInfo.augmentDensityGridGrad(pInfo.dVscloc);
		pInfo.E_nAug_dVsclocpsi = iInfo.getE_nAug();
	} else {
		getdnInc(&pInfo.dC, 0, pInfo.dnpq, pInfo.dnmq);
		getdVsclocPsi(pInfo.dnpq, pInfo.dVsclocpq, &pInfo.qvec);
		pInfo.dVsclocmq = conj(pInfo.dVsclocpq);
	}

	if (!pInfo.incommensurate) {
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			const QuantumNumber& qnum = e.eInfo.qnums[q];
			ColumnBundle dwGradEq;
			matrix dHtilde;
			dwGradEq.zero(); 

			diagMatrix F = eVars.F[q];
			applyH(e.eInfo.qnums[q], eVars.F[q], pInfo.HdC[q], pInfo.dC[q]);
			dHpsi(e.eInfo.qnums[q], pInfo.dHC[q], eVars.C[q], pInfo.dVscloc, &eVars.VdagC[q]);

			dHtilde = (pInfo.dC[q]^pInfo.HC[q]) + (eVars.C[q]^pInfo.dHC[q]) + (eVars.C[q]^pInfo.HdC[q]);

			dwGradEq = pInfo.dHC[q];
			dwGradEq += pInfo.HdC[q];
			//dwGradEq -= O(eVars.C[q]*dHtilde);
			//dwGradEq -= O(pInfo.dC[q]*eVars.Hsub[q]);
			dwGradEq -= O(eVars.C[q]*dHtilde+pInfo.dC[q]*eVars.Hsub[q]);
			//dU correction term?

			pInfo.dGradPsi[q] = dwGradEq*F*qnum.weight;
		}
	} else {
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {

			const QuantumNumber& qnum = e.eInfo.qnums[q];

			ColumnBundle dHC_Tk_k, dHC_Tinvk_k,
			HdC_Tk_k, HdC_Tinvk_k,
			dwGradEq_Tk_k, dwGradEq_Tinvk_k;

			diagMatrix F = eVars.F[q];
			int Tk_k = q;
			int Tinvk_k = q + e.eInfo.nStates;
			int Tk = q ;
			int Tinvk = q + e.eInfo.nStates;

			dHC_Tk_k = pInfo.Cinc[Tk].similar();
			dHC_Tinvk_k = pInfo.Cinc[Tinvk].similar();
			dHC_Tk_k.zero();
			dHC_Tinvk_k.zero();
			
			dHpsi(pInfo.Tk_vectors[q], dHC_Tk_k, eVars.C[q], pInfo.dVsclocpq);
			applyH(pInfo.Tk_vectors[q], F, HdC_Tk_k, pInfo.dC[Tk_k]);
			dwGradEq_Tk_k = dHC_Tk_k;
			dwGradEq_Tk_k += HdC_Tk_k;
			dwGradEq_Tk_k -= O(pInfo.Cinc[Tk]*(pInfo.Cinc[Tk]^dHC_Tk_k));
			dwGradEq_Tk_k -= O(pInfo.dC[Tk_k]*eVars.Hsub[q]);

			dHpsi(pInfo.Tinvk_vectors[q], dHC_Tinvk_k, eVars.C[q], pInfo.dVsclocmq);
			applyH(pInfo.Tinvk_vectors[q], F, HdC_Tinvk_k, pInfo.dC[Tinvk_k]);
			dwGradEq_Tinvk_k = dHC_Tinvk_k;
			dwGradEq_Tinvk_k += HdC_Tinvk_k;
			dwGradEq_Tinvk_k -= O(pInfo.Cinc[Tinvk]*(pInfo.Cinc[Tinvk]^dHC_Tinvk_k));
			dwGradEq_Tinvk_k -= O(pInfo.dC[Tinvk_k]*eVars.Hsub[q]);

			pInfo.dGradPsi[Tk_k] = (dwGradEq_Tk_k - O(pInfo.Cinc[Tk]*(pInfo.Cinc[Tk]^dwGradEq_Tk_k)))*F*qnum.weight;
			pInfo.dGradPsi[Tinvk_k] = (dwGradEq_Tinvk_k - O(pInfo.Cinc[Tinvk]*(pInfo.Cinc[Tinvk]^dwGradEq_Tinvk_k)))*F*qnum.weight;
		}
	}

	//Objective function is 1/2 x^T*A*x + b^T*x whose derivative is A*x+b
	double ener = 0.5*dot(pInfo.dY, pInfo.dGradPsi, &pInfo, &e.eInfo) + dot(pInfo.dY, pInfo.dGradTau, &pInfo, &e.eInfo);

	//logPrintf("x^Ax = %f\n", dot(pInfo.dY, pInfo.dGradPsi, &pInfo, &e.eInfo));
	//logPrintf("b^x = %f\n", dot(pInfo.dY, pInfo.dGradTau, &pInfo, &e.eInfo));

	//Apply preconditioner
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		if (grad) {
			grad->dY[q] = pInfo.dGradPsi[q] + pInfo.dGradTau[q];
			if (pInfo.incommensurate) grad->dY[q + eInfo.nStates] = pInfo.dGradPsi[q + eInfo.nStates] + pInfo.dGradTau[q + eInfo.nStates];

			if (Kgrad) {
				Kgrad->dY[q] = grad->dY[q];
				precond_inv_kinetic(Kgrad->dY[q], 1);

				if (pInfo.incommensurate) {
					Kgrad->dY[q + eInfo.nStates] = grad->dY[q + eInfo.nStates];
					precond_inv_kinetic(Kgrad->dY[q + eInfo.nStates], 1);
				}
			}
		}
	}

	return ener;
}

void PerturbationSolver::calcdGradTau() {

	if (!pInfo.incommensurate) {
		
		nullToZero(pInfo.dVsclocatom, e.gInfo);
		getdVsclocTau(pInfo.dVsclocatom);
		
		if (pInfo.atposPerturbationExists && ultrasoftDensityAugRequired) {
			for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
				ColumnBundle dOC = eVars.C[q].similar();
				dOC.zero();
				
				iInfo.species[pInfo.atomdisplacement.sp]->augmentOverlapDeriv(eVars.C[q], dOC, pInfo.Vatom_cached[q], pInfo.dVatom_cached[q]);
				pInfo.dUsqrtinvatom[q] = (-0.5)*eVars.C[q]^dOC;
				
				pInfo.dCatom[q] = eVars.C[q]*pInfo.dUsqrtinvatom[q];
			}
			
			pInfo.dnatom = getdn(&pInfo.dCatom);
			ScalarFieldArray dVsclocTmp;
			dVsclocTmp.resize(pInfo.dVsclocatom.size());
			getdVsclocPsi(pInfo.dnatom, dVsclocTmp);
			pInfo.dVsclocatom += dVsclocTmp;
		}
		
		iInfo.augmentDensityGridGrad(pInfo.dVsclocatom);
		pInfo.E_nAug_dVscloctau = iInfo.getE_nAug();
			
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			const QuantumNumber& qnum = e.eInfo.qnums[q];
			ColumnBundle dwGradEq;
			matrix dHtilde;
			diagMatrix F = eVars.F[q];

			dHtau(e.eInfo.qnums[q], pInfo.dHCatom[q], eVars.C[q], pInfo.dVsclocatom);
			dHtilde = eVars.C[q]^pInfo.dHCatom[q];

			dwGradEq = pInfo.dHCatom[q]*F;
			dwGradEq -= pInfo.OC[q]*(dHtilde*F);
			
			if(pInfo.atposPerturbationExists && ultrasoftDensityAugRequired) {
				//Atom perturbation causes wfns to become unnormalized, so dC needs to be factored in the gradient

				applyH(e.eInfo.qnums[q], eVars.F[q], pInfo.HdCatom[q], pInfo.dCatom[q]);
				dHtilde = (pInfo.dCatom[q]^pInfo.HC[q]) + (eVars.C[q]^pInfo.HdCatom[q]);

				dwGradEq += pInfo.HdCatom[q]*F;
				dwGradEq -= pInfo.OC[q]*(dHtilde*F);
				dwGradEq -= O(pInfo.dCatom[q]*(eVars.Hsub[q]*F));
				
				
				//Derivative of overlap operator w.r.t atpos
				
				ColumnBundle CCdagHCF, grad;
				
				CCdagHCF = (-1)*eVars.C[q]*(eVars.C[q]^(pInfo.HC[q]*F));
				grad = pInfo.HC[q]*F;
				grad -= pInfo.OC[q]*(eVars.Hsub[q]*F);
				
				iInfo.species[pInfo.atomdisplacement.sp]->augmentOverlapDeriv(CCdagHCF, dwGradEq, pInfo.Vatom_cached[q], pInfo.dVatom_cached[q]);
				dwGradEq += grad*pInfo.dUsqrtinvatom[q]; //Correction term 
			}
			
			pInfo.dGradTau[q] = dwGradEq*qnum.weight;
		}
	} else {
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			const QuantumNumber& qnum = e.eInfo.qnums[q];
			ColumnBundle dHC_Tk_k, dHC_Tinvk_k,
			dwGradEq_Tk_k, dwGradEq_Tinvk_k;
			diagMatrix F = eVars.F[q];

			int Tk_k = q;
			int Tinvk_k = q + e.eInfo.nStates;
			int Tk = q ;
			int Tinvk = q + e.eInfo.nStates;

			dHC_Tk_k = pInfo.Cinc[Tk].similar();
			dHC_Tinvk_k = pInfo.Cinc[Tinvk].similar();
			dHC_Tk_k.zero();
			dHC_Tinvk_k.zero();

			dHtau(pInfo.Tk_vectors[q], dHC_Tk_k, eVars.C[q], 1);
			dwGradEq_Tk_k = dHC_Tk_k - O(pInfo.Cinc[Tk]*(pInfo.Cinc[Tk]^dHC_Tk_k));

			dHtau(pInfo.Tinvk_vectors[q], dHC_Tinvk_k, eVars.C[q], -1);
			dwGradEq_Tinvk_k = dHC_Tinvk_k - O(pInfo.Cinc[Tinvk]*(pInfo.Cinc[Tinvk]^dHC_Tinvk_k));

			pInfo.dGradTau[Tk_k] = dwGradEq_Tk_k*F*qnum.weight;
			pInfo.dGradTau[Tinvk_k] = dwGradEq_Tinvk_k*F*qnum.weight;
		}
	}
}

void PerturbationSolver::getGrad(std::vector<ColumnBundle> *grad, std::vector<ColumnBundle> Y) {
	if (!pInfo.testing)
		die("This was not supposed to happen.");

	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		const QuantumNumber& qnum = e.eInfo.qnums[q];
		ColumnBundle HC, gradq, HCFUmhalf;
		matrix dHtilde, HtildeCommF;

		HCFUmhalf.zero();
		gradq.zero();

		diagMatrix F = eVars.F[q];

		applyH(qnum, F, HC, eVars.C[q]);

		HtildeCommF = eVars.Hsub[q]*F - F*eVars.Hsub[q];
		matrix Usqrtinv = invsqrt(Y[q]^O(Y[q]));
		ColumnBundle HCFUsqrtinv = HC*(F*Usqrtinv);
		gradq = HC;
		gradq -= O(eVars.C[q]*eVars.Hsub[q]);
		//gradq += 0.5*O(eVars.C[q]*HtildeCommF);
		(*grad)[q]  = gradq*F*Usqrtinv*qnum.weight;
	}
}

void PerturbationSolver::step(const PerturbationGradient& dir, double alpha)
{	assert(dir.eInfo == &eInfo);

	int nStates = eInfo.nStates;
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		axpy(alpha, dir.dY[q], e.vptInfo.dY[q]);
		if (pInfo.incommensurate) axpy(alpha, dir.dY[q+nStates], e.vptInfo.dY[q+nStates]);
	}
}

void PerturbationSolver::constrain(PerturbationGradient&)
{
	return;
}

//TODO Move checks to command section and use forbid
//TODO Implement GPU code for ultrasoft derivs
double PerturbationSolver::minimize(const MinimizeParams& params)
{
	logPrintf("VPT solver is starting.\n");
	if(!(eInfo.fillingsUpdate==ElecInfo::FillingsConst && eInfo.scalarFillings))
		die("Constant fillings only.");
	if(e.exCorr.exxFactor())
		die("Variational perturbation currently does not support exact exchange.");
	if(e.eInfo.hasU)
		die("Variational perturbation currently does not support DFT+U.");
	if (e.exCorr.needsKEdensity())
		die("Variational perturbation currently does not support KE density dependent functionals.");
	if(!e.exCorr.hasEnergy())
		die("Variational perturbation does not support potential functionals.");
	if(e.exCorr.orbitalDep)
		die("Variational perturbation currently does not support orbital dependent potential functionals.");
	if(eVars.fluidParams.fluidType != FluidNone)
		die("Variational perturbation does not support fluids.");
	if (e.symm.mode != SymmetriesNone)
		die("Variational perturbation has not been tested with symmetries yet.");
	if (eVars.n.size() > 1)
		die("Multiple spin channels not supported yet.");
	if (e.coulombParams.geometry != CoulombParams::Geometry::Periodic && pInfo.incommensurate)
		die("Periodic coloumb interaction required for incommensurate perturbations.")

	for (auto sp: e.iInfo.species) {
		if (sp->isRelativistic())
			die("Relativistic potentials are not supported yet.");
		if (sp->isUltrasoft() & pInfo.incommensurate)
			die("Ultrasoft potentials are compatible with commensurate perturbations only.");
		
		ultrasoftDensityAugRequired |= sp->isUltrasoft();
	}

	if (e.exCorr.needFiniteDifferencing())
	logPrintf("Excorr analytical derivative not available. Using finite differencing instead.\n");


	if (!pInfo.incommensurate) {
		if(eVars.wfnsFilename.empty()) {
			logPrintf("No input wavefunctions given. Performing electronic structure calculation.\n");
			elecFluidMinimize(e);
		}
	} else {
		if(eVars.wfnsFilename.empty() && pInfo.wfnsFilename.empty()) {
			die("Currently, incommensurate perturbations require prior specification of ground state wavefunctions.\n");
			//computeIncommensurateWfns();
		} else if (!eVars.wfnsFilename.empty() && pInfo.wfnsFilename.empty()) {
			die("Please specify offset wavefunctions")
		} else if (eVars.wfnsFilename.empty() && !pInfo.wfnsFilename.empty()) {
			die("Please specify ground state wavefunctions")
		}
	}
	
	//Initialize vars to random values
	if (!pInfo.incommensurate) {
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
			pInfo.dY[q].randomize(0, pInfo.dY[q].nCols());
	} else {
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			pInfo.dY[q].randomize(0, pInfo.dY[q].nCols());
			pInfo.dY[q+eInfo.nStates].randomize(0, pInfo.dY[q+eInfo.nStates].nCols());
		}
	}

	if (!eVars.wfnsFilename.empty()) {
		Energies ener;
		eVars.elecEnergyAndGrad(ener, 0, 0, true);
	}

	//e.iInfo.augmentDensityGridGrad(eVars.Vscloc);
	pInfo.updateExcorrCache(e.exCorr, e.gInfo, eVars.get_nXC()[0]);
	updateHC();
	updateNonlocalDerivs();
	
	calcdGradTau();

	if (pInfo.testing) {
		TestPerturbation t(e);
		t.ps = this;
		t.testVPT();
		return 0;
	}

	double E = Minimizable<PerturbationGradient>::minimize(params);
	logPrintf("VPT minimization completed.\n");

	return E;
}

void PerturbationSolver::computeIncommensurateWfns()
{
	die("Feature has not been implemented yet");
	//TODO: Implement
}

void PerturbationSolver::updateHC() {
	if (pInfo.incommensurate) return;
	
	iInfo.augmentDensityGridGrad(eVars.Vscloc);
	pInfo.E_nAug_cached = iInfo.getE_nAug();
	
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		QuantumNumber qnum = eInfo.qnums[q];
		applyH(qnum, eVars.F[q], pInfo.HC[q], eVars.C[q]);
		pInfo.OC[q] = O(eVars.C[q]);
	}
}

void PerturbationSolver::updateNonlocalDerivs()
{
	if (!pInfo.atposPerturbationExists || pInfo.incommensurate) return;
	
	PerturbationInfo::Mode m = pInfo.atomdisplacement;
    auto sp = iInfo.species[m.sp];

    for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
        pInfo.Vatom_cached[q] = *sp->getV(eVars.C[q], m.at);
        pInfo.dVatom_cached[q] = -DirectionalGradient(pInfo.Vatom_cached[q], m.dirCartesian);

        pInfo.VdagCatom_cached[q] = pInfo.Vatom_cached[q]^eVars.C[q];
        pInfo.dVdagCatom_cached[q] = pInfo.dVatom_cached[q]^eVars.C[q];
    }

    sp->augmentDensityGridGradDeriv(eVars.Vscloc, m.at, &m.dirLattice);
    pInfo.E_nAug_datom = sp->getE_nAug();
}


void PerturbationSolver::applyH(const QuantumNumber& qnum, const diagMatrix& Fq, ColumnBundle& HCq, const ColumnBundle& Cq)
{	static StopWatch watch("applyH"); watch.start();
	assert(Cq); //make sure wavefunction is available for this state
	std::vector<matrix> HVdagCq(e.iInfo.species.size());
	std::vector<matrix> VdagCq(e.iInfo.species.size());
	HCq.zero();

	//Propagate grad_n (Vscloc) to HCq (which is grad_Cq upto weights and fillings) if required
	//ScalarFieldArray Vscloc = getVscloc(eVars.n);
	//e.iInfo.augmentDensityGridGrad(Vscloc);
	
	iInfo.setE_nAug(pInfo.E_nAug_cached);

	//logPrintf("Vscloc\n");
	//printV(Vscloc[0]);

	HCq += Idag_DiagV_I(Cq, eVars.Vscloc);

	//TODO track how augment/project depends on q num
	e.iInfo.project(Cq, VdagCq); //update the atomic projections
	e.iInfo.augmentDensitySphericalGrad(qnum, VdagCq, HVdagCq);

	HCq += (-0.5) * L(Cq);

	//Nonlocal pseudopotentials:

	e.iInfo.EnlAndGrad(qnum, Fq, VdagCq, HVdagCq);
	e.iInfo.projectGrad(HVdagCq, Cq, HCq);
	
	watch.stop();
}

/*void PerturbationSolver::dAugatom(QuantumNumber qnum, ColumnBundle& HC, ColumnBundle& C) {
	
	matrix HVdagCq1;
	matrix HVdagCq2;
	matrix VdagCq;
	ColumnBundle tmp1, tmp2;

	
	auto m = pInfo.atomdisplacement;
	auto sp = iInfo.species[m.sp];
	vector3<> pos1 = sp->atpos[m.at];
	vector3<> pos2 = pos1 + m.dirLattice;
	
	sp->atpos[m.at] = pos2;
	sp->sync_atpos();
	e.iInfo.update(e.ener);
	//for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
	//	e.iInfo.project(eVars.C[q], eVars.VdagC[q]);
	//}
	
	e.iInfo.augmentDensityGridGrad(eVars.Vscloc);
	auto V =sp->getV(C);
	VdagCq = (*V) ^ C;
	sp->augmentDensitySphericalGrad(qnum, VdagCq, HVdagCq2);
	sp->EnlAndGrad(qnum, eVars.F[0], VdagCq, HVdagCq2);
	tmp2 = (*V)*(HVdagCq2);
	
	sp->atpos[m.at] = pos1;
	sp->sync_atpos();
	e.iInfo.update(e.ener);
	//for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
	//	e.iInfo.project(eVars.C[q], eVars.VdagC[q]);
	//}
	
	e.iInfo.augmentDensityGridGrad(eVars.Vscloc);
	V =sp->getV(C);
	VdagCq = (*V) ^ C;
	sp->augmentDensitySphericalGrad(qnum, VdagCq, HVdagCq1);
	sp->EnlAndGrad(qnum, eVars.F[0], VdagCq, HVdagCq1);
	tmp1 = (*V)*(HVdagCq1);
	
	
	
	HC += (tmp2-tmp1);
}*/

void PerturbationSolver::dH_Vscloc(const QuantumNumber& qnum, ColumnBundle& HCq, const ColumnBundle& Cq, const ScalarFieldArray& dVscloc, const std::vector<matrix>* VdagCq)
{	static StopWatch watch("dH_Vscloc"); watch.start();
	assert(Cq);

	std::vector<matrix> HVdagCq(e.iInfo.species.size());

	//e.iInfo.augmentDensityGridGrad(dVscloc);

	//TODO fix derivative w.r.t. Vscloc

	HCq.zero();

	HCq += Idag_DiagV_I(Cq, dVscloc, &HCq);

	if (VdagCq) {
		e.iInfo.augmentDensitySphericalGrad(qnum, *VdagCq, HVdagCq);
	} else {
		std::vector<matrix> VdagCq_fresh(e.iInfo.species.size());
		e.iInfo.project(Cq, VdagCq_fresh); //update the atomic projections
		e.iInfo.augmentDensitySphericalGrad(qnum, VdagCq_fresh, HVdagCq);
	}
	e.iInfo.projectGrad(HVdagCq, HCq, HCq);
	
	watch.stop();
}

void PerturbationSolver::dHpsi(const QuantumNumber& qnum, ColumnBundle& HCq, const ColumnBundle& Cq, const ScalarFieldArray& dVscloc, const std::vector<matrix>* VdagCq)
{   static StopWatch watch("dHpsi"); watch.start();
	
	iInfo.setE_nAug(pInfo.E_nAug_dVsclocpsi);
	
	dH_Vscloc(qnum, HCq, Cq, dVscloc, VdagCq);
	watch.stop();
}

void PerturbationSolver::dHpsi(const QuantumNumber& qnum, ColumnBundle& HCq, const ColumnBundle& Cq, const complexScalarFieldArray& dVscloc)
{   static StopWatch watch("dHpsi"); watch.start();
	//complexScalarFieldArray dVscloc(eVars.Vscloc.size());

	//getdVsclocPsi(dn, dVscloc, q);

	ColumnBundle HCr, HCi;
	HCr = HCq.similar();
	HCi = HCq.similar();
	HCr.zero(); HCi.zero();

	dH_Vscloc(qnum, HCr, Cq, Real(dVscloc));
	dH_Vscloc(qnum, HCi, Cq, Imag(dVscloc));

	HCq = HCr;

	HCq += complex(0,1)*HCi;
	
	watch.stop();
}

void PerturbationSolver::dHtau(const QuantumNumber& qnum, ColumnBundle& HCq, const ColumnBundle& Cq, const ScalarFieldArray& dVscloc)
{   static StopWatch watch("dHtau"); watch.start();
	assert(!pInfo.incommensurate);

	iInfo.setE_nAug(pInfo.E_nAug_dVscloctau);
	
	dH_Vscloc(qnum, HCq, Cq, dVscloc);
	
	if (pInfo.atposPerturbationExists) {
		//dAugatom(qnum, HC, C);
		
		matrix VatomdagCq, dVatomdagCq, HVatomdagCq, HdVatomdagCq;
		ColumnBundle dVHVdagCq, VHdVdagCq;
		
		PerturbationInfo::Mode m = pInfo.atomdisplacement;
		auto sp = iInfo.species[m.sp];
		auto Vatom = sp->getV(Cq, m.at);
		ColumnBundle dVatom = -DirectionalGradient(*Vatom, m.dirCartesian);
		
		VatomdagCq = (*Vatom)^Cq;
		dVatomdagCq = dVatom^Cq;
		
        e.iInfo.setE_nAug(pInfo.E_nAug_cached);
		sp->augmentDensitySphericalGrad(qnum, dVatomdagCq, HdVatomdagCq, m.at);
		sp->augmentDensitySphericalGrad(qnum, VatomdagCq, HVatomdagCq, m.at);
		
		sp->EnlAndGrad(qnum, eVars.F[qnum.index()], dVatomdagCq, HdVatomdagCq, m.at);
		sp->EnlAndGrad(qnum, eVars.F[qnum.index()], VatomdagCq, HVatomdagCq, m.at);
		
		dVHVdagCq = dVatom*HVatomdagCq;
		VHdVdagCq = (*Vatom)*HdVatomdagCq;
		
		HCq += VHdVdagCq;
		HCq += dVHVdagCq;
		
		if (ultrasoftDensityAugRequired) {
			matrix dHVatomdagCq;
			sp->setE_nAug(pInfo.E_nAug_datom);
			sp->augmentDensitySphericalGrad(qnum, VatomdagCq, dHVatomdagCq, m.at);
			HCq += (*Vatom)*dHVatomdagCq;
		}
    }
    
	watch.stop();
}

void PerturbationSolver::dHtau(const QuantumNumber& qnum, ColumnBundle& HCq, const ColumnBundle& Cq, int qsign)
{	static StopWatch watch("dHtau"); watch.start();
	assert(Cq);

	complexScalarFieldArray dVscloc(eVars.Vscloc.size());

	getdVsclocTau(dVscloc, qsign == -1);

	ColumnBundle HCr, HCi;
	HCr = HCq.similar();
	HCi = HCq.similar();
	HCr.zero(); HCi.zero();

	dH_Vscloc(qnum, HCr, Cq, Real(dVscloc));
	dH_Vscloc(qnum, HCi, Cq, Imag(dVscloc));
	
	assert(!pInfo.atposPerturbationExists);

	HCq = HCr;
	HCq += complex(0,1)*HCi;
	watch.stop();
}

//simplified version
ScalarFieldArray PerturbationSolver::getVscloc(ScalarFieldArray ntot)
{	static StopWatch watch("getVscloc"); watch.start();
	const IonInfo& iInfo = e.iInfo;

	ScalarFieldArray Vscloc(eVars.Vscloc.size());

	ScalarFieldTilde nTilde = J(ntot[0]);

	// Local part of pseudopotential:
	ScalarFieldTilde VsclocTilde = clone(iInfo.Vlocps);

	// Hartree term:
	ScalarFieldTilde dH = (*e.coulomb)(nTilde); //Note: external charge and nuclear charge contribute to d_vac as well (see below)
	VsclocTilde += dH;

	// External charge:
	if(eVars.rhoExternal)
	{	ScalarFieldTilde phiExternal = (*e.coulomb)(eVars.rhoExternal);
		VsclocTilde += phiExternal;
		logPrintf("EXTERNAL CHARGE");
	}

	// Exchange and correlation, and store the real space Vscloc with the odd historic normalization factor of JdagOJ:
	e.exCorr(eVars.get_nXC(&ntot), &eVars.Vxc, false, &eVars.tau, &eVars.Vtau);

	for(unsigned s=0; s<Vscloc.size(); s++)
	{	Vscloc[s] = JdagOJ(eVars.Vxc[s]);
		if(s<2) //Include all the spin-independent contributions along the diagonal alone
			Vscloc[s] += Jdag(O(VsclocTilde), true);

		//External potential contributions:
		if(eVars.Vexternal.size())
			Vscloc[s] += JdagOJ(eVars.Vexternal[s]);
	}
	e.symm.symmetrize(Vscloc);
	if(eVars.Vtau[0]) e.symm.symmetrize(eVars.Vtau);
	watch.stop();

	return Vscloc;
}


void PerturbationSolver::getdVsclocPsi(const ScalarFieldArray dn, ScalarFieldArray& dVscloc)
{	static StopWatch watch("getdVsclocPsi"); watch.start();
	assert(!pInfo.incommensurate);

	ScalarFieldTilde dnTilde = J(dn[0]);
	//TODO implement get_nTot

	// Hartree term:
	ScalarFieldTilde dVsclocTilde = (*e.coulomb)(dnTilde); //Note: external charge and nuclear charge contribute to d_vac as well (see below)

	//Second derivative of excorr energy
	ScalarFieldArray dVxc(eVars.Vxc.size());
	e.exCorr.getdVxc(eVars.get_nXC(), &dVxc, false, &eVars.tau, &eVars.Vtau, dn);

	for(unsigned s=0; s<eVars.Vscloc.size(); s++)
	{	dVscloc[s] = JdagOJ(dVxc[s]);

		if(s<2) //Include all the spin-independent contributions along the diagonal alone
			dVscloc[s] += Jdag(O(dVsclocTilde), true);
	}

	e.symm.symmetrize(dVscloc);
	watch.stop();
}


void PerturbationSolver::getdVsclocPsi(const complexScalarFieldArray dn, complexScalarFieldArray& dVscloc, vector3<>* q)
{	static StopWatch watch("getdVsclocPsi"); watch.start();

	assert(pInfo.incommensurate);

	complexScalarFieldTilde dnTilde = J(dn[0]);
	//TODO implement get_nTot

	// Hartree term:
	complexScalarFieldTilde dVsclocTilde;
	
	assert(e.coulombParams.geometry == CoulombParams::Geometry::Periodic);
	dVsclocTilde = -4*M_PI*Linv(O(dnTilde), q);

	//Second derivative of excorr energy, real and imaginary parts
	ScalarFieldArray dVxcRe(eVars.Vxc.size());
	ScalarFieldArray dVxcIm(eVars.Vxc.size());
	e.exCorr.getdVxc(eVars.get_nXC(), &dVxcRe, false, &eVars.tau, &eVars.Vtau, Real(dn));
	e.exCorr.getdVxc(eVars.get_nXC(), &dVxcIm, false, &eVars.tau, &eVars.Vtau, Imag(dn));

	for(unsigned s=0; s<eVars.Vscloc.size(); s++)
	{	dVscloc[s] = Complex(JdagOJ(dVxcRe[s]), JdagOJ(dVxcIm[s]));

		if(s<2) //Include all the spin-independent contributions along the diagonal alone
			dVscloc[s] += Jdag(O(dVsclocTilde), true);
	}

	e.symm.symmetrize(dVscloc);
	watch.stop();
}


void PerturbationSolver::getdVsclocTau(ScalarFieldArray& dVscloc)
{   static StopWatch watch("getdVsclocTau"); watch.start();
	
	nullToZero(dVscloc, e.gInfo);
	initZero(dVscloc);
	
	assert(!e.exCorr.orbitalDep);
		
	if (pInfo.atposPerturbationExists && ultrasoftDensityAugRequired) {
		ScalarFieldArray dn_tau;
		dn_tau = getdnatom();
		printV(dn_tau[0]);
		
		getdVsclocPsi(dn_tau, dVscloc);
	}
	
	for(unsigned s=0; s<eVars.Vscloc.size(); s++)
	{
		
		if (pInfo.VextPerturbationExists) 
			dVscloc[s] += JdagOJ(pInfo.dVext[s]);
		
		if (pInfo.atposPerturbationExists) {
			ScalarField dVlocps, drhoIon, dnCore;
			nullToZero(dVlocps, e.gInfo);
            nullToZero(drhoIon, e.gInfo);
			
            assert(!iInfo.nChargeball);
			
            ScalarFieldTilde Vlocps, rhoIon, nCoreTilde, tauCoreTilde, nChargeball;
			nullToZero(Vlocps, e.gInfo);
			nullToZero(rhoIon, e.gInfo);
			
			PerturbationInfo::Mode m = pInfo.atomdisplacement;
            auto sp = iInfo.species[m.sp];
			sp->updateLocal(Vlocps, rhoIon, nChargeball, nCoreTilde, tauCoreTilde, m.at);
			
            //Add long-range part to Vlocps and smoothen rhoIon:
            Vlocps += (*(e.coulomb))(rhoIon, Coulomb::PointChargeRight);
            assert(!(iInfo.computeStress and iInfo.ionWidth));
            //rhoIon = gaussConvolve(rhoIon, ionWidth);
			
			assert(!tauCoreTilde);
			
            //Process partial core density:
			if(nCoreTilde) {
				dnCore = -I(DirectionalGradient(nCoreTilde, m.dirCartesian));
				ScalarFieldArray dVxc(eVars.Vxc.size());
				e.exCorr.getdVxc(eVars.get_nXC(), &dVxc, false, &eVars.tau, &eVars.Vtau, getdnXC(dnCore));
				dVscloc[s] += JdagOJ(dVxc[0]);
				dnCoreA = dnCore;
				dVxcA = dVxc;
			}
			
			dVscloc[s] += -Jdag(O(DirectionalGradient(Vlocps, m.dirCartesian)), true);
			nullToZero(dVlocpsA, e.gInfo);
			initZero(dVlocpsA, e.gInfo);
			dVlocpsA = -DirectionalGradient(Vlocps, m.dirCartesian);
        }
    }

	e.symm.symmetrize(dVscloc);
	watch.stop();
}

void PerturbationSolver::getdVsclocTau(complexScalarFieldArray& dVscloc, int qsign)
{	static StopWatch watch("getdVsclocTau"); watch.start();

	assert(!pInfo.atposPerturbationExists);
	
	//debugBP("Vscloc size %d ", eVars.Vscloc.size());
	//debugBP("Dvext size %d ", pInfo.dVext.size());

	
	for(unsigned s=0; s<eVars.Vscloc.size(); s++)
	{
		initZero(dVscloc[s]);
		
		if (pInfo.VextPerturbationExists) {
			if (qsign == 1)
				dVscloc[s] = JdagOJ(pInfo.dVextpq[s]);
			else if (qsign == -1)
				dVscloc[s] = JdagOJ(pInfo.dVextmq[s]);
		}
	}

	e.symm.symmetrize(dVscloc);
	watch.stop();
}

void PerturbationSolver::orthonormalize(int q, ColumnBundle& Cq)
{	assert(e.eInfo.isMine(q));
	//VdagC[q].clear();
	//matrix rot = orthoMatrix(Cq^O(Cq)); //Compute matrix that orthonormalizes wavefunctions
	matrix rot = invsqrt(Cq^O(Cq));
	Cq = Cq * rot;
//	e->iInfo.project(C[q], VdagC[q], &rot); //update the atomic projections
}


ScalarFieldArray PerturbationSolver::addnXC(const ScalarFieldArray n)
{	if(e.iInfo.nCore)
	{	ScalarFieldArray nXC = clone(n);
		int nSpins = std::min(int(nXC.size()), 2); //1 for unpolarized and 2 for polarized
		for(int s=0; s<nSpins; s++) //note that off-diagonal components of spin-density matrix are excluded
			nXC[s] += (1./nSpins) * e.iInfo.nCore; //add core density
		return nXC;
	}
	else return n; //no cores
}

ScalarFieldArray PerturbationSolver::getdnXC(const ScalarField dnCore)
{
	ScalarFieldArray dnXC(eVars.n.size());
	nullToZero(dnXC, e.gInfo);
	if(e.iInfo.nCore)
	{
		int nSpins = std::min(int(dnXC.size()), 2); //1 for unpolarized and 2 for polarized
		for(int s=0; s<nSpins; s++) //note that off-diagonal components of spin-density matrix are excluded
		{	initZero(dnXC[s]);
			dnXC[s] += (1./nSpins) * dnCore; //add core density
		}
	}
	
	return dnXC;
}
