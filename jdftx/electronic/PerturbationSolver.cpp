/*
 * PerturbationSolver.cpp
 *
 *  Created on: Jul 22, 2022
 *      Author: brandon
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
		for(int q=eInfo->qStart; q<eInfo->qStop; q++)
		{
			if(dY[q]) dY[q] *= alpha;
		}

		if (pInfo->incommensurate) {
			for(int q=eInfo->qStart+eInfo->nStates; q<eInfo->qStop+eInfo->nStates; q++)
			{
				if(dY[q]) dY[q] *= alpha;
			}
		}

	return *this;
}

void axpy(double alpha, const PerturbationGradient& x, PerturbationGradient& y)
{	assert(x.eInfo == y.eInfo);
	for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++)
	{	if(x.dY[q]) { if(y.dY[q]) axpy(alpha, x.dY[q], y.dY[q]); else y.dY[q] = alpha*x.dY[q]; }
	}

	if (x.pInfo->incommensurate) {
		for(int q=x.eInfo->qStart+x.eInfo->nStates; q<x.eInfo->qStop+x.eInfo->nStates; q++)
		{
			if(x.dY[q]) { if(y.dY[q]) axpy(alpha, x.dY[q], y.dY[q]); else y.dY[q] = alpha*x.dY[q]; }
		}
	}
}

double dot(const PerturbationGradient& x, const PerturbationGradient& y)
{	assert(x.eInfo == y.eInfo);
	double result = 0.0;
	for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++)
	{	if(x.dY[q] && y.dY[q]) result += dotc(x.dY[q], y.dY[q]).real()*2.0;
	}

	if (x.pInfo->incommensurate) {
		for(int q=x.eInfo->qStart+x.eInfo->nStates; q<x.eInfo->qStop+x.eInfo->nStates; q++)
		{
			if(x.dY[q] && y.dY[q]) result += dotc(x.dY[q], y.dY[q]).real()*2.0;
		}
	}

	mpiWorld->allReduce(result, MPIUtil::ReduceSum);
	return result;
}

double dot(const std::vector<ColumnBundle>& x, const std::vector<ColumnBundle>& y, PerturbationInfo* pInfo, ElecInfo *eInfo)
{
	double result = 0.0;
	for(int q=eInfo->qStart; q<eInfo->qStop; q++)
	{	if(x[q] && y[q]) result += dotc(x[q], y[q]).real()*2.0;
	}

	if (pInfo->incommensurate) {
		for(int q=eInfo->qStart+eInfo->nStates; q<eInfo->qStop+eInfo->nStates; q++)
		{
			if(x[q] && y[q]) result += dotc(x[q], y[q]).real()*2.0;
		}
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
	logPrintf("Matrix values %g %g %g %g %g %g %g %g %g\n", C(0,0).x, C(0,1).x, C(0,2).x, C(1,0).x,C(1,1).x,C(1,2).x,C(2,0).x,C(2,1).x,C(2,2).x);
}

void PerturbationSolver::printV(ScalarField V) {
	double *dat = V->dataPref();
	logPrintf("Vector values %g %g %g %g %g %g %g %g %g %g\n", dat[0], dat[1], dat[2], dat[3], dat[4], dat[5], dat[6], dat[7], dat[8], dat[9]);
}

void printVvpt(ScalarField V) {
	double *dat = V->dataPref();
	logPrintf("Vector values %g %g %g %g %g %g %g %g %g %g\n", dat[0], dat[1], dat[2], dat[3], dat[4], dat[5], dat[6], dat[7], dat[8], dat[9]);
}

ScalarFieldArray PerturbationSolver::getdn(std::vector<ColumnBundle>* dC, std::vector<ColumnBundle>* C)
{	ScalarFieldArray density(eVars.n.size());

	//Runs over all states and accumulates density to the corresponding spin channel of the total density
	e.iInfo.augmentDensityInit();

	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{
		std::vector<matrix> VdagdCq(e.iInfo.species.size());
		e.iInfo.project((*dC)[q], VdagdCq); //update the atomic projections

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
	return density;
}

//TODO Add nonlocal density aug terms in
void PerturbationSolver::getdnInc(std::vector<ColumnBundle>* dC, std::vector<ColumnBundle>* C, complexScalarFieldArray& dnpq, complexScalarFieldArray& dnmq)
{

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

	//Runs over all states and accumulates density to the corresponding spin channel of the total density
	e.iInfo.augmentDensityInit();
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{
		int Tk_k = q;
		int Tinvk_k = q+e.eInfo.nStates;

		std::vector<matrix> VdagdCTinvk_k(e.iInfo.species.size());
		e.iInfo.project((*dC)[Tinvk_k], VdagdCTinvk_k); //update the atomic projections
		std::vector<matrix> VdagdCTk_k(e.iInfo.species.size());
		e.iInfo.project((*dC)[Tk_k], VdagdCTk_k); //update the atomic projections

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

	e.iInfo.augmentDensityInit();
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{
		int Tk_k = q;
		int Tinvk_k = q+e.eInfo.nStates;

		std::vector<matrix> VdagdCTinvk_k(e.iInfo.species.size());
		e.iInfo.project((*dC)[Tinvk_k], VdagdCTinvk_k); //update the atomic projections
		std::vector<matrix> VdagdCTk_k(e.iInfo.species.size());
		e.iInfo.project((*dC)[Tk_k], VdagdCTk_k); //update the atomic projections

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
	}

	e.symm.symmetrize(dnpq);
	e.symm.symmetrize(dnmq);
}


ScalarFieldArray PerturbationSolver::getn(std::vector<ColumnBundle>& C)
{	ScalarFieldArray density(eVars.n.size());
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
	return density;
}


PerturbationSolver::PerturbationSolver(Everything& e) : e(e), eVars(e.eVars), eInfo(e.eInfo), pInfo(e.vptInfo) {}

//TODO support USPPs
//TODO fix GGA derivative
//TODO M operator

double PerturbationSolver::compute(PerturbationGradient* grad, PerturbationGradient* Kgrad)
{
	if(grad) grad->init(e);
	if(Kgrad) Kgrad->init(e);

	//TODO Solve excorr derivative
	//TODO Augmented USPP derivative

	if (!pInfo.incommensurate) {
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			pInfo.dU[q] = (pInfo.dY[q]^O(eVars.C[q])) + (eVars.C[q]^O(pInfo.dY[q]));
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

	ScalarFieldArray n;

	if (!pInfo.incommensurate) {
		pInfo.dn = getdn(&pInfo.dC);
		n = eVars.n; //TODO change
	} else {
		//pInfo.dnpq
		getdnInc(&pInfo.dC, 0, pInfo.dnpq, pInfo.dnmq);
		n = eVars.n; //TODO change
	}

	if (!pInfo.incommensurate) {
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {

			const QuantumNumber& qnum = e.eInfo.qnums[q];
			ColumnBundle dHC, HC, HdC, dwGradEq, d_HCFUmhalf;
			matrix dHtilde, HtildeCommF, dHtildeCommF;

			dHC = eVars.C[q].similar();

			dwGradEq.zero(); d_HCFUmhalf.zero();

			diagMatrix F = eVars.F[q];

			dHpsi(qnum, dHC, eVars.C[q], pInfo.dn);
			//dHtau(dHC, eVars.C[q]);

			applyH(qnum, F, HC, eVars.C[q], n);
			applyH(qnum, F, HdC, pInfo.dC[q], n);

			HtildeCommF = eVars.Hsub[q]*F - F*eVars.Hsub[q];
			dHtilde = (pInfo.dC[q]^HC) + (eVars.C[q]^dHC) + (eVars.C[q]^HdC);
			dHtildeCommF = dHtilde*F - F*dHtilde;

			//TODO
			//d_HCFUmhalf = (dHC*F) + (HdC*F) - (0.5*(HC*(F*pInfo.dU[q])));
			d_HCFUmhalf = dHC*F;
			d_HCFUmhalf += HdC*F;
			d_HCFUmhalf -= 0.5*(HC*(F*pInfo.dU[q]));

			//dwGradEq = -O(pInfo.dC[q]*(eVars.C[q]^HC*F) + eVars.C[q]*(pInfo.dC[q]^HC*F)) + d_HCFUmhalf - O(eVars.C[q]*(eVars.C[q]^d_HCFUmhalf))
			//    + O(0.5*pInfo.dC[q]*HtildeCommF + 0.5*eVars.C[q]*dHtildeCommF - 0.125*eVars.C[q]*(HtildeCommF*pInfo.dU[q]+pInfo.dU[q]*HtildeCommF));
			dwGradEq = -O(pInfo.dC[q]*(eVars.C[q]^HC*F));
			dwGradEq -= O(eVars.C[q]*(pInfo.dC[q]^HC*F));
			dwGradEq += d_HCFUmhalf;

			dwGradEq -= O(eVars.C[q]*(eVars.C[q]^d_HCFUmhalf));
			dwGradEq += O(0.5*pInfo.dC[q]*HtildeCommF);
			dwGradEq += O(0.5*eVars.C[q]*dHtildeCommF);

			dwGradEq -= O(0.125*eVars.C[q]*(HtildeCommF*pInfo.dU[q]+pInfo.dU[q]*HtildeCommF));

			//printCB(dwGradEq);
			pInfo.dGradPsi[q] = dwGradEq*qnum.weight;
		}
	} else {
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {

			// Switched formulation
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

			dHpsi(pInfo.Tk_vectors[q], dHC_Tk_k, eVars.C[q], pInfo.dnpq, pInfo.qvec);
			applyH(pInfo.Tk_vectors[q], F, HdC_Tk_k, pInfo.dC[Tk_k], n);
			dwGradEq_Tk_k = dHC_Tk_k - O(pInfo.Cinc[Tk]*(pInfo.Cinc[Tk]^dHC_Tk_k)) + HdC_Tk_k - O(pInfo.dC[Tk_k]*eVars.Hsub[q]);

			dHpsi(pInfo.Tinvk_vectors[q], dHC_Tinvk_k, eVars.C[q], pInfo.dnmq, -pInfo.qvec);
			applyH(pInfo.Tinvk_vectors[q], F, HdC_Tinvk_k, pInfo.dC[Tinvk_k], n);
			dwGradEq_Tinvk_k = dHC_Tinvk_k - O(pInfo.Cinc[Tinvk]*(pInfo.Cinc[Tinvk]^dHC_Tinvk_k)) + HdC_Tinvk_k - O(pInfo.dC[Tinvk_k]*eVars.Hsub[q]);

			pInfo.dGradPsi[Tk_k] = (dwGradEq_Tk_k - O(pInfo.Cinc[Tk]*(pInfo.Cinc[Tk]^dwGradEq_Tk_k)))*F*qnum.weight;
			pInfo.dGradPsi[Tinvk_k] = (dwGradEq_Tinvk_k - O(pInfo.Cinc[Tinvk]*(pInfo.Cinc[Tinvk]^dwGradEq_Tinvk_k)))*F*qnum.weight;
		}
	}

	//Objective function is 1/2 x^T*A*x - b^T*x
	double ener = 0.5*dot(pInfo.dY, pInfo.dGradPsi, &pInfo, &e.eInfo) - dot(pInfo.dY, pInfo.dGradTau, &pInfo, &e.eInfo);
	//double ener = 0;

	logPrintf("x^Ax = %f", dot(pInfo.dY, pInfo.dGradPsi, &pInfo, &e.eInfo));
	logPrintf("b^x = %f", dot(pInfo.dY, pInfo.dGradTau, &pInfo, &e.eInfo));


	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		if (grad) {
			grad->dY[q] = pInfo.dGradPsi[q] - pInfo.dGradTau[q];
			if (pInfo.incommensurate) grad->dY[q + eInfo.nStates] = pInfo.dGradPsi[q + eInfo.nStates] - pInfo.dGradTau[q + eInfo.nStates];
			//grad->dY[q] = pInfo.dGradPsi[q];
			if (Kgrad) {
				//const QuantumNumber& qnum = eInfo.qnums[q];
				//double Nq = qnum.weight*trace(eVars.F[q]);
				//double KErollover = 2. * (Nq>1e-3 ? KEq/Nq : 1.);
				Kgrad->dY[q] = grad->dY[q];
				precond_inv_kinetic(Kgrad->dY[q], 1); //apply preconditioner

				if (pInfo.incommensurate) {
					Kgrad->dY[q + eInfo.nStates] = grad->dY[q + eInfo.nStates];
					precond_inv_kinetic(Kgrad->dY[q + eInfo.nStates], 1); //apply preconditioner
				}
			}
		}
	}

	//if (grad)
	//	logPrintf("Norm Ax-b: %g", nrm2(grad->dY[0]));getdVsclocPsi

	return ener;
}

void PerturbationSolver::calcdGradTau() {

	ScalarFieldArray n = eVars.n; //TODO change

	if (!pInfo.incommensurate) {
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			const QuantumNumber& qnum = e.eInfo.qnums[q];
			ColumnBundle dHC, HC, dwGradEq, d_HCFUmhalf;
			matrix dHtilde, HtildeCommF, dHtildeCommF;
			diagMatrix F = eVars.F[q];

			dHC = eVars.C[q].similar();

			//dHpsi(dHC, eVars.C[q], dn);
			dHtau(qnum, dHC, eVars.C[q]);
			applyH(qnum, F, HC, eVars.C[q], n);

			HtildeCommF = eVars.Hsub[q]*F - F*eVars.Hsub[q];
			dHtilde = eVars.C[q]^dHC;
			dHtildeCommF = dHtilde*F - F*dHtilde;

			//TODO
			//d_HCFUmhalf = (dHC*F) + (HdC*F) - (0.5*(HC*(F*pInfo.dU[q])));
			d_HCFUmhalf = dHC*F;

			//dwGradEq = -O(pInfo.dC[q]*(eVars.C[q]^HC*F) + eVars.C[q]*(pInfo.dC[q]^HC*F)) + d_HCFUmhalf - O(eVars.C[q]*(eVars.C[q]^d_HCFUmhalf))
			//    + O(0.5*pInfo.dC[q]*HtildeCommF + 0.5*eVars.C[q]*dHtildeCommF - 0.125*eVars.C[q]*(HtildeCommF*pInfo.dU[q]+pInfo.dU[q]*HtildeCommF));
			dwGradEq = d_HCFUmhalf;
			dwGradEq -= O(eVars.C[q]*(eVars.C[q]^d_HCFUmhalf));
			dwGradEq += O(0.5*eVars.C[q]*dHtildeCommF);
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

			//TODO is qnum correct?
			dHtau(pInfo.Tk_vectors[q], dHC_Tk_k, eVars.C[q], 1);
			dwGradEq_Tk_k = -(dHC_Tk_k - O(pInfo.Cinc[Tk]*(pInfo.Cinc[Tk]^dHC_Tk_k)));

			dHtau(pInfo.Tinvk_vectors[q], dHC_Tinvk_k, eVars.C[q], -1);
			dwGradEq_Tinvk_k = -(dHC_Tinvk_k - O(pInfo.Cinc[Tinvk]*(pInfo.Cinc[Tinvk]^dHC_Tinvk_k)));

			pInfo.dGradTau[Tk_k] = dwGradEq_Tk_k*F*qnum.weight;
			pInfo.dGradTau[Tinvk_k] = dwGradEq_Tinvk_k*F*qnum.weight;
		}
	}
}

void PerturbationSolver::getGrad(std::vector<ColumnBundle> *grad, std::vector<ColumnBundle> Y) {
	if (!pInfo.testing)
		die("This was not supposed to happen.");

	ScalarFieldArray n = eVars.n; //TODO change
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		const QuantumNumber& qnum = e.eInfo.qnums[q];
		ColumnBundle HC, gradq, HCFUmhalf;
		matrix dHtilde, HtildeCommF;

		HCFUmhalf.zero();
		gradq.zero();

		diagMatrix F = eVars.F[q];

		applyH(qnum, F, HC, eVars.C[q], n);

		HtildeCommF = eVars.Hsub[q]*F - F*eVars.Hsub[q];
		matrix Usqrtinv = invsqrt(Y[q]^O(Y[q]));
		ColumnBundle HCFUsqrtinv = HC*(F*Usqrtinv);
		gradq = HC*F*Usqrtinv;
		gradq -= O(eVars.C[q]*(eVars.C[q]^HCFUsqrtinv));
		gradq += 0.5*O(eVars.C[q]*HtildeCommF);
	    (*grad)[q]  = gradq*qnum.weight;
	}
}

void PerturbationSolver::step(const PerturbationGradient& dir, double alpha)
{	assert(dir.eInfo == &eInfo);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		axpy(alpha, dir.dY[q], e.vptInfo.dY[q]);
	}
	if (pInfo.incommensurate) {
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			axpy(alpha, dir.dY[q+eInfo.nStates], e.vptInfo.dY[q+eInfo.nStates]);
		}
	}
}

bool PerturbationSolver::report(int iter)
{
	/*
	debugBP("Iteration %d\n", iter);
	assert(eInfo.isMine(0));
	debugBP("%d\n", e.vptInfo.dY.size());
	complexScalarFieldTilde dat = e.vptInfo.dY[0].getColumn(0, 0);
	for (int i = 0; i < 10; i++) {
		debugBP("%f\n", (dat->data(false)[i]).x);
	}*/
	return false;
}

void PerturbationSolver::constrain(PerturbationGradient&)
{
	return;
}

double PerturbationSolver::minimize(const MinimizeParams& params)
{
	logPrintf("VPT solver is starting.\n");

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
		die("VPT not tested with spins.");

	for (auto sp: e.iInfo.species) {
		if (sp->isRelativistic())
			die("Nonrelativistic PPs only");
		if (sp->isUltrasoft() & pInfo.incommensurate)
			die("Ultrasoft for commensurate perturbations only.");
	}

	if (e.exCorr.needFiniteDifferencing())
	logPrintf("Using finite differencing for Vxc.\n");


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

	if (!eVars.wfnsFilename.empty()) {
		Energies ener;
		eVars.elecEnergyAndGrad(ener, nullptr, nullptr, true);
	}
	calcdGradTau();

	pInfo.testing = true;

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
	//TODO: Implement
}

void PerturbationSolver::applyH(QuantumNumber qnum, const diagMatrix& F, ColumnBundle& HC, ColumnBundle& C, ScalarFieldArray ntot)
{	assert(C); //make sure wavefunction is available for this state
	std::vector<matrix> HVdagCq(e.iInfo.species.size());
	std::vector<matrix> VdagCq(e.iInfo.species.size());
	HC.zero();

	//Propagate grad_n (Vscloc) to HCq (which is grad_Cq upto weights and fillings) if required
	ScalarFieldArray Vscloc = getVscloc(ntot);
	e.iInfo.augmentDensityGridGrad(Vscloc);

	//logPrintf("Vscloc\n");
	//printV(Vscloc[0]);

	HC += Idag_DiagV_I(C, Vscloc);

	//TODO track how augment/project depends on q num
	e.iInfo.project(C, VdagCq); //update the atomic projections
	e.iInfo.augmentDensitySphericalGrad(qnum, VdagCq, HVdagCq);

	HC += (-0.5) * L(C);

	//Nonlocal pseudopotentials:

	e.iInfo.EnlAndGrad(qnum, F, VdagCq, HVdagCq);
	e.iInfo.projectGrad(HVdagCq, C, HC);
}


void PerturbationSolver::applyH2(int q, const diagMatrix& Fq, ColumnBundle& HCq) //Used only for comparison against applyHamiltonian
{
	const QuantumNumber& qnum = e.eInfo.qnums[q];
	std::vector<matrix> HVdagCq(e.iInfo.species.size());

	HCq.zero();

	//Propagate grad_n (Vscloc) to HCq (which is grad_Cq upto weights and fillings) if required
	HCq += Idag_DiagV_I(eVars.C[q], eVars.Vscloc); //Accumulate Idag Diag(Vscloc) I C

	e.iInfo.augmentDensitySphericalGrad(qnum, eVars.VdagC[q], HVdagCq); //Contribution via pseudopotential density augmentation

	//Kinetic energy:
	ColumnBundle LCq = L(eVars.C[q]);
	HCq += (-0.5) * LCq;

	//Nonlocal pseudopotentials:
	e.iInfo.EnlAndGrad(qnum, Fq, eVars.VdagC[q], HVdagCq);
	e.iInfo.projectGrad(HVdagCq, eVars.C[q], HCq);
	//Compute subspace hamiltonian if needed:
	//Hsub[q] = C[q] ^ HCq;
}


void PerturbationSolver::dH_Vscloc(QuantumNumber qnum, ColumnBundle& HC, ColumnBundle& C, ScalarFieldArray dVscloc)
{	assert(C);

	std::vector<matrix> HVdagCq(e.iInfo.species.size());
	std::vector<matrix> VdagCq(e.iInfo.species.size());

	e.iInfo.augmentDensityGridGrad(dVscloc);

	//TODO fix derivative w.r.t. Vscloc
	//TODO: Use correct quantum number in incommensurate calculations

	HC.zero();

	HC += Idag_DiagV_I(C, dVscloc, &HC);

	e.iInfo.project(C, VdagCq); //update the atomic projections
	e.iInfo.augmentDensitySphericalGrad(qnum, VdagCq, HVdagCq);
	e.iInfo.projectGrad(HVdagCq, HC, HC);
}

//TODO complexify by ScalarFieldArray -> complexScalarFieldArray
void PerturbationSolver::dHpsi(QuantumNumber qnum, ColumnBundle& HC, ColumnBundle& C, ScalarFieldArray dn)
{
	ScalarFieldArray dVscloc(eVars.Vscloc.size());
	getdVsclocPsi(dn, dVscloc);

	//logPrintf("Test A %f\n", nrm2(HC));

	dH_Vscloc(qnum, HC, C, dVscloc);
}

void PerturbationSolver::dHpsi(QuantumNumber qnum, ColumnBundle& HC, ColumnBundle& C, complexScalarFieldArray dn, vector3<> q)
{
	complexScalarFieldArray dVscloc(eVars.Vscloc.size());

	getdVsclocPsi(dn, dVscloc, q);

	ColumnBundle HCr, HCi;
	HCr = HC.similar();
	HCi = HC.similar();
	HCr.zero(); HCi.zero();

	dH_Vscloc(qnum, HCr, C, Real(dVscloc));
	dH_Vscloc(qnum, HCi, C, Imag(dVscloc));

	//printCB(HCr);
	//printCB(HCi);

	HC = HCr;

	//logPrintf("Test A %f\n", nrm2(HCr));
	//logPrintf("Test B %f\n", nrm2(HCi));

	HC += complex(0,1)*HCi;

	//printCB(HC);
}

//TODO complexify
void PerturbationSolver::dHtau(QuantumNumber qnum, ColumnBundle& HC, ColumnBundle& C)
{
	assert(!pInfo.incommensurate);


	complexScalarFieldArray dVscloc_complex(eVars.Vscloc.size());
	ScalarFieldArray dVscloc(eVars.Vscloc.size());

	getdVsclocTau(dVscloc_complex);
	dVscloc = Real(dVscloc_complex);

	dH_Vscloc(qnum, HC, C, dVscloc);
}

void PerturbationSolver::dHtau(QuantumNumber qnum, ColumnBundle& HC, ColumnBundle& C, int qsign)
{	assert(C);

	complexScalarFieldArray dVscloc(eVars.Vscloc.size());

	getdVsclocTau(dVscloc, qsign == -1);

	ColumnBundle HCr, HCi;
	HCr = HC.similar();
	HCi = HC.similar();
	HCr.zero(); HCi.zero();

	dH_Vscloc(qnum, HCr, C, Real(dVscloc));
	dH_Vscloc(qnum, HCi, C, Imag(dVscloc));

	//logPrintf("AA");
	//printV(Real(dVscloc)[0]);
	//printV(Imag(dVscloc)[0]);

	HC = HCr;
	HC += complex(0,1)*HCi;
}

//simplified version
ScalarFieldArray PerturbationSolver::getVscloc(ScalarFieldArray ntot)
{	static StopWatch watch("EdensityAndVscloc"); watch.start();
	const IonInfo& iInfo = e.iInfo;

	ScalarFieldArray Vscloc(eVars.Vscloc.size());

	ScalarFieldTilde nTilde = J(ntot[0]); //TODO change

	// Local part of pseudopotential:
	ScalarFieldTilde VsclocTilde = clone(iInfo.Vlocps);

	// Hartree term:
	ScalarFieldTilde dH = (*e.coulomb)(nTilde); //Note: external charge and nuclear charge contribute to d_vac as well (see below)
	VsclocTilde += dH;

	// External charge:
	if(eVars.rhoExternal)
	{	ScalarFieldTilde phiExternal = (*e.coulomb)(eVars.rhoExternal);
		VsclocTilde += phiExternal;
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


void PerturbationSolver::getdVsclocPsi(ScalarFieldArray dn, ScalarFieldArray& dVscloc)
{	static StopWatch watch("EdensityAndVscloc"); watch.start();

	//logPrintf("getdVsc\n");
	//printV(dn[0]);


	ScalarFieldTilde dnTilde = J(dn[0]);
	//TODO implememtn get_nTot

	// Hartree term:
	ScalarFieldTilde dVsclocTilde = (*e.coulomb)(dnTilde); //Note: external charge and nuclear charge contribute to d_vac as well (see below)

	// Exchange and correlation, and store the real space Vscloc with the odd historic normalization factor of JdagOJ:
	//e.exCorr(eVars.get_nXC(), &eVars.Vxc, false, &eVars.tau, &eVars.Vtau);
	//Second derivative exc
	ScalarFieldArray dVxc(eVars.Vxc.size());
	ScalarField tmp;
	e.exCorr.getdVxc(eVars.get_nXC(), &dVxc, false, &tmp, &eVars.tau, &eVars.Vtau, dn);

	for(unsigned s=0; s<eVars.Vscloc.size(); s++)
	{	dVscloc[s] = JdagOJ(dVxc[s]);

		if(s<2) //Include all the spin-independent contributions along the diagonal alone
			dVscloc[s] += Jdag(O(dVsclocTilde), true);
	}

	//printV(dVscloc[0]);

	e.symm.symmetrize(dVscloc);
	watch.stop();
}


void PerturbationSolver::getdVsclocPsi(complexScalarFieldArray dn, complexScalarFieldArray& dVscloc, vector3<> k)
{	static StopWatch watch("EdensityAndVscloc"); watch.start();

	//logPrintf("getdVsc\n");
	//printV(dn[0]);


	complexScalarFieldTilde dnTilde = J(dn[0]);
	//TODO implememtn get_nTot

	// Hartree term:
	complexScalarFieldTilde dVsclocTilde;
	dVsclocTilde = -4*M_PI*Linv(O(dnTilde), &k);

	// Exchange and correlation, and store the real space Vscloc with the odd historic normalization factor of JdagOJ:
	//e.exCorr(eVars.get_nXC(), &eVars.Vxc, false, &eVars.tau, &eVars.Vtau);
	//Second derivative exc
	ScalarFieldArray dVxcRe(eVars.Vxc.size());
	ScalarFieldArray dVxcIm(eVars.Vxc.size());
	ScalarField tmpr;
	ScalarField tmpi;
	e.exCorr.getdVxc(eVars.get_nXC(), &dVxcRe, false, &tmpr, &eVars.tau, &eVars.Vtau, Real(dn));
	e.exCorr.getdVxc(eVars.get_nXC(), &dVxcIm, false, &tmpi, &eVars.tau, &eVars.Vtau, Imag(dn));

	for(unsigned s=0; s<eVars.Vscloc.size(); s++)
	{	dVscloc[s] = Complex(JdagOJ(dVxcRe[s]), JdagOJ(dVxcIm[s]));

		if(s<2) //Include all the spin-independent contributions along the diagonal alone
			dVscloc[s] += Jdag(O(dVsclocTilde), true);
	}

	//printV(dVscloc[0]);

	e.symm.symmetrize(dVscloc);
	watch.stop();
}



void PerturbationSolver::getdVsclocTau(complexScalarFieldArray& dVscloc, bool conjugate)
{	static StopWatch watch("EdensityAndVscloc"); watch.start();

	//debugBP("Vscloc size %d ", eVars.Vscloc.size());
	//debugBP("Dvext size %d ", pInfo.dVext.size());

	for(unsigned s=0; s<eVars.Vscloc.size(); s++)
	{
		if (!conjugate)
			dVscloc[s] = JdagOJ(pInfo.dVext[s]);
		else
			dVscloc[s] = JdagOJ(conj(pInfo.dVext[s]));
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


ScalarFieldArray PerturbationSolver::addnXC(ScalarFieldArray n)
{	if(e.iInfo.nCore)
	{	ScalarFieldArray nXC = clone(n);
		int nSpins = std::min(int(nXC.size()), 2); //1 for unpolarized and 2 for polarized
		for(int s=0; s<nSpins; s++) //note that off-diagonal components of spin-density matrix are excluded
			nXC[s] += (1./nSpins) * e.iInfo.nCore; //add core density
		return nXC;
	}
	else return n; //no cores
}
