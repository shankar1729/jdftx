/*-------------------------------------------------------------------
Copyright 2022 Brandon Li

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

#include <electronic/PerturbationSolver.h>
#include <electronic/TestPerturbation.h>
#include <electronic/Everything.h>
#include <electronic/ElecMinimizer.h>
#include <electronic/ColumnBundle.h>
#include <electronic/ExCorr.h>
#include <core/matrix.h>
#include <core/Units.h>
#include <core/Minimize.h>
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
	if (pInfo->commensurate)
		X.resize(eInfo->nStates);
	else
		X.resize(eInfo->nStates*2);
}

PerturbationGradient& PerturbationGradient::operator*=(double alpha)
{
	int nStates = eInfo->nStates;

	for(int q=eInfo->qStart; q<eInfo->qStop; q++)
	{
		if( X[q]) X[q] *= alpha;
		if (!pInfo->commensurate) if( X[q+nStates]) X[q+nStates] *= alpha;
	}

	return *this;
}

void axpy(double alpha, const PerturbationGradient& x, PerturbationGradient& y)
{	assert(x.eInfo == y.eInfo);

	int nStates = x.eInfo->nStates;

	for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++)
	{
		if(x.X[q]) { if(y.X[q]) axpy(alpha, x.X[q], y.X[q]); else y.X[q] = alpha*x.X[q]; }
		if (!x.pInfo->commensurate) if(x.X[q+nStates]) { if(y.X[q+nStates]) axpy(alpha, x.X[q+nStates], y.X[q+nStates]); else y.X[q+nStates] = alpha*x.X[q+nStates]; }
	}
}

double dot(const PerturbationGradient& x, const PerturbationGradient& y)
{	assert(x.eInfo == y.eInfo);

	int nStates = x.eInfo->nStates;
	double result = 0.0;

	for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++)
	{
		if(x.X[q] && y.X[q]) result += dotc(x.X[q], y.X[q]).real()*2.0;
		if (!x.pInfo->commensurate) if(x.X[q+nStates] && y.X[q+nStates]) result += dotc(x.X[q+nStates], y.X[q+nStates]).real()*2.0;
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
		if (!pInfo->commensurate) if(x[q+nStates] && y[q+nStates]) result += dotc(x[q+nStates], y[q+nStates]).real()*2.0;
	}

	mpiWorld->allReduce(result, MPIUtil::ReduceSum);
	return result;
}

PerturbationGradient clone(const PerturbationGradient& x) { return x; } //!< create a copy

void randomize(PerturbationGradient& x) { randomize(x.X, *x.eInfo); } //!< initialize with random numbers

ScalarField randomRealSpaceVector(ColumnBundle& ref, const GridInfo* gInfo) {
	ColumnBundle res = ref;
	res.randomize(0, 1);
	ScalarField vec = Real(I(res.getColumn(0, 0)));
	ScalarFieldTilde vtilde = J(vec);
	vtilde->setGzero(0);
	return changeGrid(I(vtilde), *gInfo);
}


PerturbationSolver::PerturbationSolver(Everything& e) : e(e), eVars(e.eVars), eInfo(e.eInfo), iInfo(e.iInfo), pInfo(e.vptInfo) {}


//TODO Test GPU code for ultrasoft derivs
//TODO Eliminate redundant zeros
// Determine changes to existing code
void PerturbationSolver::solvePerturbation() {
	if (!e.vptParams.nIterations)
		die("Error: Must specify nIterations in command solve-perturbation to enable the variational perturbation solver.\n")
	
	logPrintf("Variational perturbation solver is starting.\n");
	
	if (pInfo.datom) {
		AtomicMode m = pInfo.datom->mode;
		if (m.sp < 0 || m.sp >= iInfo.species.size() || m.at < 0 || m.at >= iInfo.species[m.sp]->atpos.size())
			die("Perturbed atom does not exist! Check species and atom number.\n");
		logPrintf("Perturbing species %d atom %d in direction %g %g %g.\n", m.sp, m.at, m.dirCartesian[0], m.dirCartesian[1], m.dirCartesian[2]);
	}
	
	if (pInfo.dVext)
		logPrintf("Perturbing Vexternal.\n");
	
	if(pInfo.commensurate && e.eVars.wfnsFilename.empty() && !e.eVars.skipWfnsInit) {
		logPrintf("No input wavefunctions given. Performing electronic structure calculation first.\n");
		elecFluidMinimize(e);
	}
	
	state.init(e);
	if (pInfo.commensurate) {
		init(state.X, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
			state.X[q].zero();
	} else {
		pInfo.initInc(state.X, 2*eInfo.nStates, eInfo.nBands, &eInfo);
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			state.X[q].zero();
			state.X[q+eInfo.nStates].zero();
		}
	}
	
	eVars.elecEnergyAndGrad(e.ener, 0, 0, true);

	updateExcorrCache();
	updateHC();
	updateNonlocalDerivs();
	
	pInfo.dGradTau.init(e);
	pInfo.dGradPsi.init(e);
	
	calcdGradTau();
	
	solve(pInfo.dGradTau, e.vptParams);
	
	hessian(pInfo.dGradPsi, state);
}

void PerturbationSolver::precondition(PerturbationGradient& v) {
	//Apply inverse kinetic preconditioner
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		precond_inv_kinetic(v.X[q], 1);
		if (!pInfo.commensurate) precond_inv_kinetic(v.X[q + eInfo.nStates], 1);
	}
}

void PerturbationSolver::Ksqrt(PerturbationGradient& v) {
	//Apply inverse kinetic preconditioner
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		precond_inv_kinetic(v.X[q], 1, true);
		if (!pInfo.commensurate) precond_inv_kinetic(v.X[q + eInfo.nStates], 1, true);
	}
}

void PerturbationSolver::hessian(PerturbationGradient& Av, const PerturbationGradient& v) {
	static StopWatch watch("dGradPsi"); watch.start();
	Av.init(e);
	
	if (pInfo.commensurate) {
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			//pInfo.dU[q] = (pInfo.dY[q]^O(eVars.C[q])) + (eVars.C[q]^O(pInfo.dY[q]));
			pInfo.dU[q] = 2*dagger_symmetrize(v.X[q]^pInfo.OC[q]);
			pInfo.dC[q] = v.X[q];
			pInfo.dC[q] -= 0.5*(eVars.C[q]*pInfo.dU[q]);
		}
	} else {
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			int kpq_k = q;
			int kmq_k = q+e.eInfo.nStates;
			int kpq = q ;
			int kmq = q + e.eInfo.nStates;

			pInfo.dC[kpq_k] = v.X[kpq_k] - O(pInfo.Cinc[kpq]*(pInfo.Cinc[kpq]^v.X[kpq_k]));
			pInfo.dC[kmq_k] = v.X[kmq_k] - O(pInfo.Cinc[kmq]*(pInfo.Cinc[kmq]^v.X[kmq_k]));
		}
	}
	if (pInfo.commensurate) {
		getdn(pInfo.dn, &pInfo.dC);
		getdVsclocPsi(pInfo.dn, pInfo.dVscloc);
		e.iInfo.augmentDensityGridGrad(pInfo.dVscloc);
		pInfo.E_nAug_dVsclocpsi = iInfo.getE_nAug();
	} else {
		getdnInc(&pInfo.dC, 0, pInfo.dnpq, pInfo.dnmq);
		getdVsclocPsi(pInfo.dnpq, pInfo.dVsclocpq, &pInfo.qvec);
		pInfo.dVsclocmq = conj(pInfo.dVsclocpq);
	}

	if (pInfo.commensurate) {
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			const QuantumNumber& qnum = e.eInfo.qnums[q];
			ColumnBundle dwGradEq, dHCq = eVars.C[q].similar(), HdCq = eVars.C[q].similar();

			diagMatrix F = eVars.F[q];
			applyH(e.eInfo.qnums[q], eVars.F[q], HdCq, pInfo.dC[q]);
			dHpsi(e.eInfo.qnums[q], dHCq, eVars.C[q], pInfo.dVscloc, &eVars.VdagC[q]);

			pInfo.CdagdHC[q] = eVars.C[q]^dHCq;
			pInfo.dHsub[q] = (pInfo.dC[q]^pInfo.HC[q]) + pInfo.CdagdHC[q] + (eVars.C[q]^HdCq);
			
			dwGradEq = dHCq;
			dwGradEq += HdCq;
			dwGradEq -= O(eVars.C[q]*pInfo.dHsub[q]+pInfo.dC[q]*eVars.Hsub[q]);
			
			dwGradEq = dwGradEq*eVars.F[q];
			
			dwGradEq += (-0.5)*pInfo.grad[q]*(eVars.F[q]*pInfo.dU[q]); //Second order correction term to be used if C does not minimize energy exactly

			Av.X[q] = dwGradEq*qnum.weight;
		}
	} else {
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {

			const QuantumNumber& qnum = e.eInfo.qnums[q];

			ColumnBundle dHC_kpq_k, dHC_kmq_k,
			HdC_kpq_k, HdC_kmq_k,
			dwGradEq_kpq_k, dwGradEq_kmq_k;

			diagMatrix F = eVars.F[q];
			int kpq_k = q;
			int kmq_k = q + e.eInfo.nStates;
			int kpq = q ;
			int kmq = q + e.eInfo.nStates;

			dHC_kpq_k = pInfo.Cinc[kpq].similar();
			dHC_kmq_k = pInfo.Cinc[kmq].similar();
			dHC_kpq_k.zero();
			dHC_kmq_k.zero();
			
			dHpsi(pInfo.kplusq_vectors[q], dHC_kpq_k, eVars.C[q], pInfo.dVsclocpq);
			applyH(pInfo.kplusq_vectors[q], F, HdC_kpq_k, pInfo.dC[kpq_k]);
			dwGradEq_kpq_k = dHC_kpq_k;
			dwGradEq_kpq_k += HdC_kpq_k;
			dwGradEq_kpq_k -= O(pInfo.Cinc[kpq]*(pInfo.Cinc[kpq]^dHC_kpq_k));
			dwGradEq_kpq_k -= O(pInfo.dC[kpq_k]*eVars.Hsub[q]);

			dHpsi(pInfo.kminusq_vectors[q], dHC_kmq_k, eVars.C[q], pInfo.dVsclocmq);
			applyH(pInfo.kminusq_vectors[q], F, HdC_kmq_k, pInfo.dC[kmq_k]);
			dwGradEq_kmq_k = dHC_kmq_k;
			dwGradEq_kmq_k += HdC_kmq_k;
			dwGradEq_kmq_k -= O(pInfo.Cinc[kmq]*(pInfo.Cinc[kmq]^dHC_kmq_k));
			dwGradEq_kmq_k -= O(pInfo.dC[kmq_k]*eVars.Hsub[q]);

			Av.X[kpq_k] = (dwGradEq_kpq_k - O(pInfo.Cinc[kpq]*(pInfo.Cinc[kpq]^dwGradEq_kpq_k)))*F*qnum.weight;
			Av.X[kmq_k] = (dwGradEq_kmq_k - O(pInfo.Cinc[kmq]*(pInfo.Cinc[kmq]^dwGradEq_kmq_k)))*F*qnum.weight;
		}
	}
	
	watch.stop();
}

void PerturbationSolver::calcdGradTau() {
	static StopWatch watch("dGradTau"); watch.start();
	if (pInfo.commensurate) {
		
		if (pInfo.datom && pInfo.datom->isUltrasoft(iInfo)) {
			getdnatom(pInfo.datom->dnatom);
			
			for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
				ColumnBundle dOC = eVars.C[q].similar();
				dOC.zero();
				
				iInfo.species[pInfo.datom->mode.sp]->augmentOverlapDeriv(eVars.C[q], dOC, pInfo.datom->Vatom_cached[q], pInfo.datom->dVatom_cached[q]);
				pInfo.dUmhalfatom[q] = (-0.5)*eVars.C[q]^dOC;
				
				pInfo.datom->dCatom[q] = eVars.C[q]*pInfo.dUmhalfatom[q];
			}
			
			ScalarFieldArray dnCatom(eVars.n.size());
			nullToZero(dnCatom, e.gInfo);
			getdn(dnCatom, &pInfo.datom->dCatom);
			pInfo.datom->dnatom += dnCatom;
		}
		
		nullToZero(pInfo.dVsclocTau, e.gInfo);
		
		if (pInfo.datom && pInfo.datom->isUltrasoft(iInfo))
			getdVsclocTau(pInfo.dVsclocTau, &pInfo.datom->dnatom);
		else
			getdVsclocTau(pInfo.dVsclocTau);
		
		iInfo.augmentDensityGridGrad(pInfo.dVsclocTau);
		pInfo.E_nAug_dVscloctau = iInfo.getE_nAug();
			
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			const QuantumNumber& qnum = e.eInfo.qnums[q];
			ColumnBundle dwGradEq, dHCq = eVars.C[q].similar(), HdCq = eVars.C[q].similar();
			diagMatrix F = eVars.F[q];

			dHtau(e.eInfo.qnums[q], dHCq, eVars.C[q], pInfo.dVsclocTau);
			pInfo.dHsubatom[q] = eVars.C[q]^dHCq;
			pInfo.CdagdHCatom[q] = pInfo.dHsubatom[q];
			
			if(pInfo.datom && pInfo.datom->isUltrasoft(iInfo)) {
				applyH(e.eInfo.qnums[q], eVars.F[q], HdCq, pInfo.datom->dCatom[q]);
				pInfo.dHsubatom[q] += (pInfo.datom->dCatom[q]^pInfo.HC[q]) + (eVars.C[q]^HdCq);
			}

			dwGradEq = dHCq;
			dwGradEq -= pInfo.OC[q]*pInfo.dHsubatom[q];
			
			if(pInfo.datom && pInfo.datom->isUltrasoft(iInfo)) {
				//Atom perturbation causes wfns to become unnormalized, so dC needs to be factored in the gradient
				dwGradEq += HdCq;
				dwGradEq -= O(pInfo.datom->dCatom[q]*eVars.Hsub[q]);
				
				//Derivative of overlap operator w.r.t atpos
				iInfo.species[pInfo.datom->mode.sp]->augmentOverlapDeriv((-1)*eVars.C[q]*eVars.Hsub[q], dwGradEq, pInfo.datom->Vatom_cached[q], pInfo.datom->dVatom_cached[q]);
			}
			
			dwGradEq = dwGradEq*F;
			
			if(pInfo.datom && pInfo.datom->isUltrasoft(iInfo))
				dwGradEq += pInfo.grad[q]*(pInfo.dUmhalfatom[q]*F); //Correction term 
			
			pInfo.dGradTau.X[q] = dwGradEq*qnum.weight;
		}
	} else {
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			const QuantumNumber& qnum = e.eInfo.qnums[q];
			ColumnBundle dHC_kpq_k, dHC_kmq_k,
			dwGradEq_kpq_k, dwGradEq_kmq_k;
			diagMatrix F = eVars.F[q];

			int kpq_k = q;
			int kmq_k = q + e.eInfo.nStates;
			int kpq = q ;
			int kmq = q + e.eInfo.nStates;

			dHC_kpq_k = pInfo.Cinc[kpq].similar();
			dHC_kmq_k = pInfo.Cinc[kmq].similar();
			dHC_kpq_k.zero();
			dHC_kmq_k.zero();

			dHtau(pInfo.kplusq_vectors[q], dHC_kpq_k, eVars.C[q], 1);
			dwGradEq_kpq_k = dHC_kpq_k - O(pInfo.Cinc[kpq]*(pInfo.Cinc[kpq]^dHC_kpq_k));

			dHtau(pInfo.kminusq_vectors[q], dHC_kmq_k, eVars.C[q], -1);
			dwGradEq_kmq_k = dHC_kmq_k - O(pInfo.Cinc[kmq]*(pInfo.Cinc[kmq]^dHC_kmq_k));

			pInfo.dGradTau.X[kpq_k] = dwGradEq_kpq_k*F*qnum.weight;
			pInfo.dGradTau.X[kmq_k] = dwGradEq_kmq_k*F*qnum.weight;
		}
	}
	
	watch.stop();
}

void PerturbationSolver::getGrad(std::vector<ColumnBundle> *grad, std::vector<ColumnBundle> Y) {
	assert(pInfo.testing);
	
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		const QuantumNumber& qnum = e.eInfo.qnums[q];
		ColumnBundle HC, gradq, HCFUmhalf;
		matrix dHtilde, HtildeCommF;

		HCFUmhalf.zero();
		gradq.zero();

		diagMatrix F = eVars.F[q];

		applyH(qnum, F, HC, eVars.C[q]);

		//HtildeCommF = eVars.Hsub[q]*F - F*eVars.Hsub[q];
		matrix Usqrtinv = invsqrt(Y[q]^O(Y[q]));
		//ColumnBundle HCFUsqrtinv = HC*(F*Usqrtinv);
		gradq = HC;
		gradq -= O(eVars.C[q]*eVars.Hsub[q]);
		//gradq += 0.5*O(eVars.C[q]*HtildeCommF);
		(*grad)[q]  = gradq*F*Usqrtinv*qnum.weight;
	}
}

void PerturbationSolver::computeIncommensurateWfns()
{
	die("Feature has not been implemented yet\n");
	//TODO: Implement
}

void PerturbationSolver::updateExcorrCache()
{
	ScalarFieldArray nXC = eVars.get_nXC();
	
	if (nXC.size() > 1)
		return;
	
	int iDirStart, iDirStop;
	TaskDivision(3, mpiWorld).myRange(iDirStart, iDirStop);
	{
		for(int i=iDirStart; i<iDirStop; i++) {
			pInfo.IDJn_cached[i] = I(D(J(nXC[0]),i));
		}

		nullToZero(pInfo.sigma_cached, e.gInfo);
		initZero(pInfo.sigma_cached);

		for(int i=iDirStart; i<iDirStop; i++)
			pInfo.sigma_cached += pInfo.IDJn_cached[i] * pInfo.IDJn_cached[i];
		pInfo.sigma_cached->allReduceData(mpiWorld, MPIUtil::ReduceSum);
	}
	
	e.exCorr.getSecondDerivatives(nXC[0], pInfo.e_nn_cached, pInfo.e_sigma_cached, pInfo.e_nsigma_cached, pInfo.e_sigmasigma_cached, 1e-9, &pInfo.sigma_cached);
}

void PerturbationSolver::updateHC() {
	if (!pInfo.commensurate) return;
	
	iInfo.augmentDensityGridGrad(eVars.Vscloc);
	pInfo.E_nAug_cached = iInfo.getE_nAug();
	
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		QuantumNumber qnum = eInfo.qnums[q];
		applyH(qnum, eVars.F[q], pInfo.HC[q], eVars.C[q]);
		pInfo.OC[q] = O(eVars.C[q]);
		
		pInfo.grad[q] = pInfo.HC[q];
		pInfo.grad[q] -= pInfo.OC[q]*(eVars.Hsub[q]);
	}
}

void PerturbationSolver::updateNonlocalDerivs()
{
	if (!pInfo.datom || !pInfo.commensurate) return;
	
	AtomicMode m = pInfo.datom->mode;
	auto sp = iInfo.species[m.sp];

	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		pInfo.datom->Vatom_cached[q] = *sp->getV(eVars.C[q], m.at);
		pInfo.datom->dVatom_cached[q] = -D(pInfo.datom->Vatom_cached[q], m.dirCartesian);

		pInfo.datom->VdagCatom_cached[q] = pInfo.datom->Vatom_cached[q]^eVars.C[q];
		pInfo.datom->dVdagCatom_cached[q] = pInfo.datom->dVatom_cached[q]^eVars.C[q];
	}

	sp->augmentDensityGridGradDeriv(eVars.Vscloc, m.at, &m.dirLattice);
	pInfo.datom->E_nAug_datom = sp->getE_nAug();
}



void PerturbationSolver::getdn(ScalarFieldArray& dn, const std::vector<ColumnBundle>* dC, const std::vector<ColumnBundle>* C)
{	static StopWatch watch("getdn"); watch.start();

	initZero(dn);
	//Runs over all states and accumulates density to the corresponding spin channel of the total density
	e.iInfo.augmentDensityInit();
	
	std::vector<matrix> VdagdCq;

	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{
		if (pInfo.densityAugRequired(e)) {
			VdagdCq.resize(e.iInfo.species.size());
			e.iInfo.project((*dC)[q], VdagdCq);
		}
		
		if (C) {
			std::vector<matrix> VdagCq(e.iInfo.species.size());
			e.iInfo.project((*C)[q], VdagCq); //update the atomic projections
			
			dn += 2.0*e.eInfo.qnums[q].weight * Real(diagouterI(eVars.F[q], (*C)[q], (*dC)[q], dn.size(), &e.gInfo));
			e.iInfo.augmentDensitySpherical(e.eInfo.qnums[q], eVars.F[q], VdagCq, &VdagdCq);
		} else {
			dn += 2.0*e.eInfo.qnums[q].weight * Real(diagouterI(eVars.F[q], eVars.C[q], (*dC)[q], dn.size(), &e.gInfo));
			e.iInfo.augmentDensitySpherical(e.eInfo.qnums[q], eVars.F[q], eVars.VdagC[q], &VdagdCq);
		}
	}
	e.iInfo.augmentDensityGrid(dn);
	for(ScalarField& ns: dn)
	{	nullToZero(ns, e.gInfo);
		ns->allReduceData(mpiWorld, MPIUtil::ReduceSum);
	}
	e.symm.symmetrize(dn);
	watch.stop();
}

void PerturbationSolver::getdnInc(const std::vector<ColumnBundle>* dC, const std::vector<ColumnBundle>* C, complexScalarFieldArray& dnpq, complexScalarFieldArray& dnmq)
{   static StopWatch watch("getdnInc"); watch.start();

	initZero(dnpq);
	dnpq[0] = 0*dnpq[0];
	initZero(dnmq);

	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{
		int kpq_k = q;
		int kmq_k = q+e.eInfo.nStates;

		if (C) {
			dnpq += e.eInfo.qnums[q].weight * diagouterI(eVars.F[q], (*C)[q], (*dC)[kmq_k], dnpq.size(), &e.gInfo);
			dnpq += e.eInfo.qnums[q].weight * diagouterI(eVars.F[q], (*dC)[kpq_k], (*C)[q], dnpq.size(), &e.gInfo);
		} else {
			dnpq += e.eInfo.qnums[q].weight * diagouterI(eVars.F[q], eVars.C[q], (*dC)[kmq_k], dnpq.size(), &e.gInfo);
			dnpq += e.eInfo.qnums[q].weight * diagouterI(eVars.F[q], (*dC)[kpq_k], eVars.C[q], dnpq.size(), &e.gInfo);
		}
	}

	for(complexScalarField& ns: dnpq)
	{	nullToZero(ns, e.gInfo);
		ns->allReduceData(mpiWorld, MPIUtil::ReduceSum);
	}

	e.symm.symmetrize(dnpq);
	
	dnmq = conj(dnpq);
	
	watch.stop();
}

void PerturbationSolver::getdnatom(ScalarFieldArray& dnatom)
{	static StopWatch watch("getdn"); watch.start();
	nullToZero(dnatom, e.gInfo, eVars.n.size());
	for (unsigned s = 0; s < dnatom.size(); s++)
		initZero(dnatom[s]);
	
	if (pInfo.datom && pInfo.datom->isUltrasoft(iInfo)) {
		AtomicMode m = pInfo.datom->mode;
		auto sp = iInfo.species[m.sp];
		//Runs over all states and accumulates density to the corresponding spin channel of the total density
		sp->augmentDensityInit(m.at);
		
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		{
			sp->augmentDensitySpherical(e.eInfo.qnums[q], eVars.F[q], pInfo.datom->VdagCatom_cached[q], &pInfo.datom->dVdagCatom_cached[q], 0, m.at);
		}
		sp->augmentDensityGrid(dnatom, m.at);
		
		
		sp->augmentDensityInit(m.at);
		
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		{
			sp->augmentDensitySpherical(e.eInfo.qnums[q], eVars.F[q], pInfo.datom->VdagCatom_cached[q], 0, 0, m.at);
		}
		sp->augmentDensityGrid(dnatom, m.at, &m.dirLattice);
		
		for(ScalarField& ns: dnatom)
		{	nullToZero(ns, e.gInfo);
			ns->allReduceData(mpiWorld, MPIUtil::ReduceSum);
		}
		e.symm.symmetrize(dnatom);
	}
	watch.stop();
}


void PerturbationSolver::applyH(const QuantumNumber& qnum, const diagMatrix& Fq, ColumnBundle& HCq, const ColumnBundle& Cq)
{	static StopWatch watch("applyH"); watch.start();
	assert(Cq); //make sure wavefunction is available for this state
	
	std::vector<matrix> HVdagCq(e.iInfo.species.size());
	std::vector<matrix> VdagCq(e.iInfo.species.size());
	
	HCq.zero();
	
	iInfo.setE_nAug(pInfo.E_nAug_cached);

	HCq += Idag_DiagV_I(Cq, eVars.Vscloc);

	e.iInfo.project(Cq, VdagCq); //update the atomic projections
	e.iInfo.augmentDensitySphericalGrad(qnum, VdagCq, HVdagCq);

	HCq += (-0.5) * L(Cq);

	//Nonlocal pseudopotentials:

	e.iInfo.EnlAndGrad(qnum, Fq, VdagCq, HVdagCq);
	e.iInfo.projectGrad(HVdagCq, Cq, HCq);
	
	watch.stop();
}

void PerturbationSolver::dH_Vscloc(const QuantumNumber& qnum, ColumnBundle& HCq, const ColumnBundle& Cq, const ScalarFieldArray& dVscloc, const std::vector<matrix>* VdagCq)
{	static StopWatch watch("dH_Vscloc"); watch.start();
	assert(Cq);

	std::vector<matrix> HVdagCq(e.iInfo.species.size());

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
	assert(pInfo.commensurate);

	iInfo.setE_nAug(pInfo.E_nAug_dVscloctau);
	
	dH_Vscloc(qnum, HCq, Cq, dVscloc);
	
	if (pInfo.datom) {
		//dAugatom(qnum, HC, C);
		
		matrix VatomdagCq, dVatomdagCq, HVatomdagCq, HdVatomdagCq;
		ColumnBundle dVHVdagCq, VHdVdagCq;
		
		AtomicMode m = pInfo.datom->mode;
		auto sp = iInfo.species[m.sp];
		auto Vatom = sp->getV(Cq, m.at);
		ColumnBundle dVatom = -D(*Vatom, m.dirCartesian);
		
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
		
		if (pInfo.datom->isUltrasoft(iInfo)) {
			matrix dHVatomdagCq;
			sp->setE_nAug(pInfo.datom->E_nAug_datom);
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

	getdVsclocTau(dVscloc, qsign);

	ColumnBundle HCr, HCi;
	HCr = HCq.similar();
	HCi = HCq.similar();
	HCr.zero(); HCi.zero();

	dH_Vscloc(qnum, HCr, Cq, Real(dVscloc));
	dH_Vscloc(qnum, HCi, Cq, Imag(dVscloc));
	
	assert(!pInfo.datom);

	HCq = HCr;
	HCq += complex(0,1)*HCi;
	watch.stop();
}

void PerturbationSolver::getdVsclocPsi(const ScalarFieldArray dn, ScalarFieldArray& dVscloc)
{	static StopWatch watch("getdVsclocPsi"); watch.start();
	assert(pInfo.commensurate);

	ScalarFieldTilde dnTilde = J(get_nTot(dn));

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

	assert(!pInfo.commensurate);

	complexScalarFieldTilde dnTilde = J(get_nTot(dn));

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


void PerturbationSolver::getdVsclocTau(ScalarFieldArray& dVscloc, ScalarFieldArray* dn)
{   static StopWatch watch("getdVsclocTau"); watch.start();
	
	nullToZero(dVscloc, e.gInfo);
	initZero(dVscloc);
	
	assert(!e.exCorr.orbitalDep);
		
	if (pInfo.datom && pInfo.datom->isUltrasoft(iInfo)) {
		assert(dn);
		getdVsclocPsi(*dn, dVscloc);
	}
	
	for(unsigned s=0; s<eVars.Vscloc.size(); s++) {
		if (pInfo.dVext) dVscloc[s] += JdagOJ(pInfo.dVext->dVext[s]);
		if (pInfo.drhoExt) dVscloc[s] += Jdag(O((*e.coulomb)(pInfo.drhoExt->drhoExt)), true);
		if (pInfo.dChargeBall) dVscloc[s] += Jdag(O((*e.coulomb)(pInfo.dChargeBall->drhoExt)), true);
		if (pInfo.dElectricField) {
			CoulombParams params = e.coulombParams;
			params.Efield = pInfo.dElectricField->Efield;
			pInfo.dElectricField->C = params.createCoulomb(e.gInfo, e.gInfoWfns, e.coulombWfns);
			dVscloc[s] += JdagOJ(pInfo.dElectricField->C->getEfieldPotential());
		}
	}
		
		
	if (pInfo.datom) {
		ScalarField dVlocps, drhoIon, dnCore;
		nullToZero(dVlocps, e.gInfo);
		nullToZero(drhoIon, e.gInfo);
		
		assert(!iInfo.nChargeball);
		
		ScalarFieldTilde rhoIon, nCoreTilde, tauCoreTilde, nChargeball;
		nullToZero(pInfo.datom->Vlocps, e.gInfo);
		initZero(pInfo.datom->Vlocps);
		nullToZero(rhoIon, e.gInfo);
		
		AtomicMode m = pInfo.datom->mode;
		auto sp = iInfo.species[m.sp];
		sp->updateLocal(pInfo.datom->Vlocps, rhoIon, nChargeball, nCoreTilde, tauCoreTilde, m.at);
		
		//Add long-range part to Vlocps and smoothen rhoIon:
		pInfo.datom->Vlocps += (*(e.coulomb))(rhoIon, Coulomb::PointChargeRight);
		assert(!(iInfo.computeStress and iInfo.ionWidth));
		//rhoIon = gaussConvolve(rhoIon, ionWidth);
		
		assert(!tauCoreTilde);
		
		//Process partial core density:
		if(nCoreTilde) {
			dnCore = -I(D(nCoreTilde, m.dirCartesian));
			ScalarFieldArray dVxc(eVars.Vxc.size());
			e.exCorr.getdVxc(eVars.get_nXC(), &dVxc, false, &eVars.tau, &eVars.Vtau, getdnXC(dnCore));
			
			
			for(unsigned s=0; s<eVars.Vscloc.size(); s++)
				dVscloc[s] += JdagOJ(dVxc[s]);
			
			pInfo.dnCoreA = dnCore;
			//dVxcA = dVxc;
		}
		
		for(unsigned s=0; s<eVars.Vscloc.size(); s++)
			dVscloc[s] += -Jdag(O(D(pInfo.datom->Vlocps, m.dirCartesian)), true);
		
		//nullToZero(dVlocpsA, e.gInfo);
		//initZero(dVlocpsA, e.gInfo);
		//dVlocpsA = -D(Vlocps, m.dirCartesian);
	}

	e.symm.symmetrize(dVscloc);
	watch.stop();
}

void PerturbationSolver::getdVsclocTau(complexScalarFieldArray& dVscloc, int qsign)
{	static StopWatch watch("getdVsclocTau"); watch.start();

	assert(!pInfo.commensurate && !pInfo.datom);
	
	nullToZero(dVscloc, e.gInfo);
	initZero(dVscloc);

	for(unsigned s=0; s<eVars.Vscloc.size(); s++)
	{
		vector3<> qvec = pInfo.qvec*qsign;
		
		if (pInfo.dVext) {
			if (qsign == 1) dVscloc[s] += JdagOJ(pInfo.dVext->dVextpq[s]);
			else if (qsign == -1) dVscloc[s] += JdagOJ(pInfo.dVext->dVextmq[s]);
		}
		
		if (pInfo.drhoExt) {
			if (qsign == 1) dVscloc[s] += Jdag(O(-4*M_PI*Linv(O(pInfo.drhoExt->drhoExtpq), &qvec)), true);
			else if (qsign == -1) dVscloc[s] += Jdag(O(-4*M_PI*Linv(O(pInfo.drhoExt->drhoExtmq), &qvec)), true);
		}
		
		if (pInfo.dChargeBall) {
			if (qsign == 1) dVscloc[s] += Jdag(O(-4*M_PI*Linv(O(pInfo.dChargeBall->drhoExtpq), &qvec)), true);
			else if (qsign == -1) dVscloc[s] += Jdag(O(-4*M_PI*Linv(O(pInfo.dChargeBall->drhoExtmq), &qvec)), true);
		}
	}

	e.symm.symmetrize(dVscloc);
	watch.stop();
}

ScalarFieldArray PerturbationSolver::getdnXC(const ScalarField dnCore) const
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

