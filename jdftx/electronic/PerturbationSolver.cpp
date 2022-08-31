/*
 * PerturbationSolver.cpp
 *
 *  Created on: Jul 22, 2022
 *      Author: brandon
 */

#include <electronic/PerturbationSolver.h>
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
	dY.resize(eInfo->nStates);
}

PerturbationGradient& PerturbationGradient::operator*=(double alpha)
{	for(int q=eInfo->qStart; q<eInfo->qStop; q++)
	{
		if(dY[q]) dY[q] *= alpha;
	}
	return *this;
}

void axpy(double alpha, const PerturbationGradient& x, PerturbationGradient& y)
{	assert(x.eInfo == y.eInfo);
	for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++)
	{	if(x.dY[q]) { if(y.dY[q]) axpy(alpha, x.dY[q], y.dY[q]); else y.dY[q] = alpha*x.dY[q]; }
	}
}

double dot(const PerturbationGradient& x, const PerturbationGradient& y)
{	assert(x.eInfo == y.eInfo);
	double result = 0.0;
	for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++)
	{	if(x.dY[q] && y.dY[q]) result += dotc(x.dY[q], y.dY[q]).real()*2.0;
	}
	mpiWorld->allReduce(result, MPIUtil::ReduceSum);
	return result;
}

double dot(const std::vector<ColumnBundle>& x, const std::vector<ColumnBundle>& y, int qStart, int qStop)
{
	double result = 0.0;
	for(int q=qStart; q<qStop; q++)
	{	if(x[q] && y[q]) result += dotc(x[q], y[q]).real()*2.0;
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

void printCB(ColumnBundle C) {
	double *dat = (double*)(C.getColumn(0, 0)->dataPref());
	logPrintf("ColumnBundle values %g %g %g %g %g %g %g %g %g %g\n", dat[0], dat[1], dat[2], dat[3], dat[4], dat[5], dat[6], dat[7], dat[8], dat[9]);
}

void printV(ScalarField V) {
	double *dat = V->dataPref();
	logPrintf("ColumnBundle values %g %g %g %g %g %g %g %g %g %g\n", dat[0], dat[1], dat[2], dat[3], dat[4], dat[5], dat[6], dat[7], dat[8], dat[9]);
}


ScalarFieldArray PerturbationSolver::getdn(std::vector<ColumnBundle>* dC, std::vector<ColumnBundle>* C)
{	ScalarFieldArray density(eVars.n.size());
	//Runs over all states and accumulates density to the corresponding spin channel of the total density
	//e.iInfo.augmentDensityInit();
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{
		if (C)
			density += 2.0*e.eInfo.qnums[q].weight * Real(diagouterI(eVars.F[q], (*C)[q], (*dC)[q], density.size(), &e.gInfo));
		else
			density += 2.0*e.eInfo.qnums[q].weight * Real(diagouterI(eVars.F[q], eVars.C[q], (*dC)[q], density.size(), &e.gInfo));
		//e.iInfo.augmentDensitySpherical(e.eInfo.qnums[q], eVars.F[q], eVars.VdagC[q]); //pseudopotential contribution
	}
	//e.iInfo.augmentDensityGrid(density);
	for(ScalarField& ns: density)
	{	nullToZero(ns, e.gInfo);
		ns->allReduceData(mpiWorld, MPIUtil::ReduceSum);
	}
	e.symm.symmetrize(density);
	return density;
}



ScalarFieldArray PerturbationSolver::getn(std::vector<ColumnBundle>& C)
{	ScalarFieldArray density(eVars.n.size());
	nullToZero(density, e.gInfo);

	debugBP("Density norm %g.\n", nrm2(density[0]));
	//Runs over all states and accumulates density to the corresponding spin channel of the total density
	//e->iInfo.augmentDensityInit();
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{	density += e.eInfo.qnums[q].weight * diagouterI(eVars.F[q], C[q], density.size(), &e.gInfo);
		//e->iInfo.augmentDensitySpherical(e->eInfo.qnums[q], F[q], VdagC[q]); //pseudopotential contribution
	}
	//e.iInfo.augmentDensityGrid(density);
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

double PerturbationSolver::compute(PerturbationGradient* grad, PerturbationGradient* Kgrad)
{
	debugBP("testC\n");
	if(grad) grad->init(e);
	if(Kgrad) Kgrad->init(e);

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

	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		pInfo.dU[q] = (pInfo.dY[q]^O(eVars.C[q])) + (eVars.C[q]^O(pInfo.dY[q]));
		pInfo.dC[q] = pInfo.dY[q];
		pInfo.dC[q] -= 0.5*(eVars.C[q]*pInfo.dU[q]);
	}

	ScalarFieldArray dn = getdn(&pInfo.dC);
	ScalarFieldArray n = eVars.n; //TODO change
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		debugBP("Compute %d", q);
		const QuantumNumber& qnum = e.eInfo.qnums[q];
		ColumnBundle dHC, HC, HdC, dwGradEq, d_HCFUmhalf;
		matrix dHtilde, HtildeCommF, dHtildeCommF;

		dwGradEq.zero(); d_HCFUmhalf.zero();

		diagMatrix F = eVars.F[q];

		debugBP("Compute A0.\n");

		dHpsi(dHC, eVars.C[q], dn);
		//dHtau(dHC, eVars.C[q]);

		applyH(q, F, HC, eVars.C[q], n);
		applyH(q, F, HdC, pInfo.dC[q], n);

		debugBP("Compute A1.\n");

		HtildeCommF = eVars.Hsub[q]*F - F*eVars.Hsub[q];
	    dHtilde = (pInfo.dC[q]^HC) + (eVars.C[q]^dHC) + (eVars.C[q]^HdC);
	    dHtildeCommF = dHtilde*F - F*dHtilde;

		//debugBP("Compute A2");
	    //TODO
	    //d_HCFUmhalf = (dHC*F) + (HdC*F) - (0.5*(HC*(F*pInfo.dU[q])));
	    d_HCFUmhalf = dHC*F;
	    d_HCFUmhalf += HdC*F;
	    d_HCFUmhalf -= 0.5*(HC*(F*pInfo.dU[q]));

		//debugBP("Compute A3");
	    //dwGradEq = -O(pInfo.dC[q]*(eVars.C[q]^HC*F) + eVars.C[q]*(pInfo.dC[q]^HC*F)) + d_HCFUmhalf - O(eVars.C[q]*(eVars.C[q]^d_HCFUmhalf))
	    //    + O(0.5*pInfo.dC[q]*HtildeCommF + 0.5*eVars.C[q]*dHtildeCommF - 0.125*eVars.C[q]*(HtildeCommF*pInfo.dU[q]+pInfo.dU[q]*HtildeCommF));
	    dwGradEq = -O(pInfo.dC[q]*(eVars.C[q]^HC*F));
	    dwGradEq -= O(eVars.C[q]*(pInfo.dC[q]^HC*F));
	    dwGradEq += d_HCFUmhalf;
		//debugBP("Compute A4");
	    dwGradEq -= O(eVars.C[q]*(eVars.C[q]^d_HCFUmhalf));
	    dwGradEq += O(0.5*pInfo.dC[q]*HtildeCommF);
	    dwGradEq += O(0.5*eVars.C[q]*dHtildeCommF);
		//debugBP("Compute A5");
	    dwGradEq -= O(0.125*eVars.C[q]*(HtildeCommF*pInfo.dU[q]+pInfo.dU[q]*HtildeCommF));
		debugBP("Compute A6 Norm %g %g.\n", nrm2(HtildeCommF), nrm2(dHtildeCommF));
	    pInfo.dGradPsi[q] = dwGradEq*qnum.weight;
	}

	debugBP("Compute A.\n");
	//Objective function is 1/2 x^T*A*x - b^T*x
	double ener = 0.5*dot(pInfo.dY, pInfo.dGradPsi, eInfo.qStart, eInfo.qStop) - dot(pInfo.dY, pInfo.dGradTau, eInfo.qStart, eInfo.qStop);
	//double ener = 0;

	debugBP("Compute B0");
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		if (grad) {
			grad->dY[q] = pInfo.dGradPsi[q] - pInfo.dGradTau[q];
			//grad->dY[q] = pInfo.dGradPsi[q];
			if (Kgrad) {
				//const QuantumNumber& qnum = eInfo.qnums[q];
				//double Nq = qnum.weight*trace(eVars.F[q]);
				//double KErollover = 2. * (Nq>1e-3 ? KEq/Nq : 1.);
				Kgrad->dY[q] = grad->dY[q];
				precond_inv_kinetic(Kgrad->dY[q], 1); //apply preconditioner
			}
		}
	}
	debugBP("Compute B");

	//if (grad)
	//	logPrintf("Norm Ax-b: %g", nrm2(grad->dY[0]));

	return ener;
}

void PerturbationSolver::calcdGradTau() {

	ScalarFieldArray n = eVars.n; //TODO change
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		const QuantumNumber& qnum = e.eInfo.qnums[q];
		ColumnBundle dHC, HC, dwGradEq, d_HCFUmhalf;
		matrix dHtilde, HtildeCommF, dHtildeCommF;
		diagMatrix F = eVars.F[q];

	    debugBP("dGradtau A.\n");
		//dHpsi(dHC, eVars.C[q], dn);
		dHtau(dHC, eVars.C[q]);
	    debugBP("dGradtau A1.\n");
		applyH(q, F, HC, eVars.C[q], n);
	    debugBP("dGradtau A2.\n");

		HtildeCommF = eVars.Hsub[q]*F - F*eVars.Hsub[q];
	    dHtilde = eVars.C[q]^dHC;
	    dHtildeCommF = dHtilde*F - F*dHtilde;

	    //TODO
	    //d_HCFUmhalf = (dHC*F) + (HdC*F) - (0.5*(HC*(F*pInfo.dU[q])));
	    d_HCFUmhalf = dHC*F;


	    debugBP("dGradtau B.\n");
	    //dwGradEq = -O(pInfo.dC[q]*(eVars.C[q]^HC*F) + eVars.C[q]*(pInfo.dC[q]^HC*F)) + d_HCFUmhalf - O(eVars.C[q]*(eVars.C[q]^d_HCFUmhalf))
	    //    + O(0.5*pInfo.dC[q]*HtildeCommF + 0.5*eVars.C[q]*dHtildeCommF - 0.125*eVars.C[q]*(HtildeCommF*pInfo.dU[q]+pInfo.dU[q]*HtildeCommF));
	    dwGradEq = d_HCFUmhalf;
	    dwGradEq -= O(eVars.C[q]*(eVars.C[q]^d_HCFUmhalf));
	    dwGradEq += O(0.5*eVars.C[q]*dHtildeCommF);
	    pInfo.dGradTau[q] = dwGradEq*qnum.weight;
	}

}

void PerturbationSolver::getGrad(std::vector<ColumnBundle> *grad, std::vector<ColumnBundle> Y) {
	debugBP("Getgrad A.\n");

	ScalarFieldArray n = eVars.n; //TODO change
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		const QuantumNumber& qnum = e.eInfo.qnums[q];
		ColumnBundle HC, gradq, HCFUmhalf;
		matrix dHtilde, HtildeCommF;

		HCFUmhalf.zero();
		gradq.zero();

		diagMatrix F = eVars.F[q];

		applyH(q, F, HC, eVars.C[q], n);

		HtildeCommF = eVars.Hsub[q]*F - F*eVars.Hsub[q];
		matrix Usqrtinv = invsqrt(Y[q]^O(Y[q]));
		ColumnBundle HCFUsqrtinv = HC*(F*Usqrtinv);
		gradq = HC*F*Usqrtinv;
		gradq -= O(eVars.C[q]*(eVars.C[q]^HCFUsqrtinv));
		gradq += 0.5*O(eVars.C[q]*HtildeCommF);
	    (*grad)[q]  = gradq*qnum.weight;
	}
	debugBP("Getgrad B.\n");
}

void PerturbationSolver::step(const PerturbationGradient& dir, double alpha)
{	assert(dir.eInfo == &eInfo);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		axpy(alpha, dir.dY[q], e.vptInfo.dY[q]);
	}
	debugBP("testD.\n");
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

ScalarField randomRealSpaceVector(ColumnBundle& ref) {
	ColumnBundle res = ref;
	res.randomize(0, 1);
	ScalarField vec = Real(I(res.getColumn(0, 0)));
	ScalarFieldTilde vtilde = J(vec);
	vtilde->setGzero(0);
	//res.getColumn(0, 0)->setGzero(0);
	return I(vtilde);
	//vec->dataPref()[0] = 0;
	//logPrintf("Norm %g", nrm2(vec));
	//return vec;
}

double PerturbationSolver::minimize(const MinimizeParams& params)
{
	std::cout <<"Starting";
	debugBP("testB.\n");

	bool testWhole = true;

	ScalarFieldArray Vext(eVars.n.size()), dVext(eVars.n.size());
	ScalarFieldArray dnNum, dnAnal;
	if (testWhole) {
		double h = 1e-2;
		nullToZero(Vext, e.gInfo);
		nullToZero(dVext, e.gInfo);
		Vext[0] = randomRealSpaceVector(eVars.C[0]);
		dVext[0] = randomRealSpaceVector(eVars.C[0]);
		//initRandomFlat(Vext[0]);
		//initRandomFlat(dVext[0]);
		complexScalarFieldTilde VextTilde, dVextTilde;
		nullToZero(VextTilde, e.gInfo);
		nullToZero(dVextTilde, e.gInfo);
		/*(VextTilde->dataPref())[69] = 1;
		(VextTilde->dataPref())[23] = 0.5;
		(VextTilde->dataPref())[44] = 0.6;
		(VextTilde->dataPref())[21] = -0.6;
		(dVextTilde->dataPref())[4] = 0.2;
		(dVextTilde->dataPref())[5] = -0.7;
		(dVextTilde->dataPref())[6] = 0.2;*/
		//Vext[0] = Real(I(VextTilde));
		//dVext[0] = Real(I(dVextTilde));
		pInfo.dVext = clone(dVext);

		//Vext[0] = Vext[0];
		logPrintf("Integral %g.\n", integral(Vext[0]));
		dVext[0] = dVext[0]*h;

		eVars.Vexternal = Vext+dVext;
		elecFluidMinimize(e);
		ScalarFieldArray n2 = clone(eVars.n);

		eVars.Vexternal = Vext;
		elecFluidMinimize(e);
		ScalarFieldArray n = clone(eVars.n);

		dnNum = (n2 - n)*(1/h);
	} else if(eVars.wfnsFilename.length() == 0) {
		elecFluidMinimize(e);
	}

	pInfo.dY.resize(eInfo.nStates);
	pInfo.dC.resize(eInfo.nStates);
	pInfo.dU.resize(eInfo.nStates);
	pInfo.dGradPsi.resize(eInfo.nStates);
	pInfo.dGradTau.resize(eInfo.nStates);

	init(pInfo.dY, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	init(pInfo.dC, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	init(pInfo.dGradPsi, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	init(pInfo.dGradTau, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	randomize(pInfo.dY, eInfo);
	//for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
	//	e.vptInfo.dY[q] *= 0;
	//}
	calcdGradTau();

	//testVPT();

	double E = Minimizable<PerturbationGradient>::minimize(params);

	if (testWhole) {
		dnAnal = (-1)*getdn(&pInfo.dC, &eVars.C);

		logPrintf("Difference %g.\n", nrm2(dnNum[0]-dnAnal[0])/nrm2(dnAnal[0]));
		printV(dnNum[0]);
		printV(dnAnal[0]);

		eVars.n = Vext;
		//die("Test complete.\n");
	}

	return E;
}


void PerturbationSolver::applyH(int q, const diagMatrix& F, ColumnBundle& HC, ColumnBundle& C, ScalarFieldArray ntot)
{	assert(C); //make sure wavefunction is available for this state
	const QuantumNumber& qnum = e.eInfo.qnums[q];
	std::vector<matrix> HVdagC(e.iInfo.species.size());
	std::vector<matrix> VdagC(e.iInfo.species.size());
	HC.zero();

	debugBP("applyH a.\n");
	//Propagate grad_n (Vscloc) to HCq (which is grad_Cq upto weights and fillings) if required
	ScalarFieldArray Vscloc = getVscloc(ntot);

	HC += Idag_DiagV_I(C, Vscloc);

	e.iInfo.project(C, VdagC); //update the atomic projections
	e.iInfo.augmentDensitySphericalGrad(e.eInfo.qnums[q], VdagC, HVdagC);

	//Kinetic energy:
	HC += (-0.5) * L(C);

	debugBP("applyH b.\n");
	//Nonlocal pseudopotentials:
	e.iInfo.EnlAndGrad(qnum, F, VdagC, HVdagC);
	e.iInfo.projectGrad(HVdagC, C, HC);
}


void PerturbationSolver::applyH2(int q, const diagMatrix& Fq, ColumnBundle& HCq) //Do not use this
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


void PerturbationSolver::dHpsi(ColumnBundle& HC, ColumnBundle& C, ScalarFieldArray dn)
{	assert(C);

	ScalarFieldArray dVscloc(eVars.Vscloc.size());

	getdVsclocPsi(dn, dVscloc);

	debugBP("dVscloc norm: %g.\n", nrm2(dVscloc[0]));

	HC.zero();
	HC += Idag_DiagV_I(C, dVscloc);
}

void PerturbationSolver::dHtau(ColumnBundle& HC, ColumnBundle& C)
{	assert(C);

	ScalarFieldArray dVscloc(eVars.Vscloc.size());

	getdVsclocTau(dVscloc);

	HC.zero();
	HC += Idag_DiagV_I(C, dVscloc);

	//Nonlocal pseudopotentials:
	//e.iInfo.project(C, VdagC); //update the atomic projections
	//e.iInfo.EnlAndGrad(qnum, F, VdagC, HVdagC);
	//e.iInfo.projectGrad(HVdagC, C, HC);
}

//simplified version
ScalarFieldArray PerturbationSolver::getVscloc(ScalarFieldArray ntot)
{	static StopWatch watch("EdensityAndVscloc"); watch.start();
	const IonInfo& iInfo = e.iInfo;

	ScalarFieldArray Vscloc(eVars.Vscloc.size());

	debugBP("N size %d.\n", (int)ntot.size());
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

	debugBP("VsclocPsi A.\n");
	for(unsigned s=0; s<eVars.Vscloc.size(); s++)
	{	dVscloc[s] = JdagOJ(dVxc[s]);
	debugBP("VsclocPsi B.\n");
		if(s<2) //Include all the spin-independent contributions along the diagonal alone
			dVscloc[s] += Jdag(O(dVsclocTilde), true);
	}

	debugBP("dVscloc norm: %g.\n", nrm2(dVscloc[0]));

	debugBP("VsclocPsi C.\n");
	e.symm.symmetrize(dVscloc);
	watch.stop();
}



void PerturbationSolver::getdVsclocTau(ScalarFieldArray& dVscloc)
{	static StopWatch watch("EdensityAndVscloc"); watch.start();

	//debugBP("Vscloc size %d ", eVars.Vscloc.size());
	//debugBP("Dvext size %d ", pInfo.dVext.size());

	for(unsigned s=0; s<eVars.Vscloc.size(); s++)
	{
		dVscloc[s] = JdagOJ(pInfo.dVext[s]);
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

void PerturbationSolver::testVPT() {
	int test = 10;
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

	std::vector<ColumnBundle> W;
	std::vector<ColumnBundle> dW;
	std::vector<ColumnBundle> W2;
	std::vector<ColumnBundle> Y;
	std::vector<ColumnBundle> dY;
	std::vector<ColumnBundle> Y2;
	std::vector<ColumnBundle> H1Y;
	std::vector<ColumnBundle> H2Y;
	std::vector<ColumnBundle> H3Y;
	std::vector<ColumnBundle> Grad1;
	std::vector<ColumnBundle> Grad2;
	PerturbationGradient pg;
	std::vector<ColumnBundle> dHY;
	pg.init(e);
	//std::vector<ColumnBundle> dGrad;
	init(W, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	init(dW, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	init(W2, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	init(Y, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	init(dY, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	init(Y2, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	init(H1Y, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	init(H2Y, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	init(H3Y, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	init(Grad1, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	init(Grad2, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	init(pg.dY, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	init(dHY, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	randomize(W, eInfo);
	randomize(dW, eInfo);
	double h = 1e-7;
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		H1Y[q].zero();
		H2Y[q].zero();
		H3Y[q].zero();
		dW[q] = dW[q]*h;
		orthonormalize(q, W[q]);
		Y[q] = W[q];
		W2[q] = W[q] + dW[q];
		//orthonormalize(q, Y[q]);
		Y2[q] = W2[q];
		orthonormalize(q, Y2[q]);
		dY[q] = Y2[q] - Y[q];
		dHY[q].zero();
	}
	logPrintf("Norm %g.\n", nrm2(H1Y[0]));
	logPrintf("Norm %g.\n", nrm2(H2Y[0]));

	if (test == 1) { //Test Hamtiltonian
		ColumnBundle tmpC;
		Energies ener;
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			//eVars.n = getn(Y);
			applyH(q, eVars.F[q], H1Y[q], Y[q], getn(Y));

			eVars.C[q] = Y[q];
			e.iInfo.project(eVars.C[q], eVars.VdagC[q]);
			eVars.n = eVars.calcDensity();
			eVars.EdensityAndVscloc(ener);
			//eVars.elecEnergyAndGrad(ener,0,0,true);
			e.iInfo.augmentDensityGridGrad(eVars.Vscloc);
			eVars.applyHamiltonian(q, eVars.F[q], H2Y[q], ener, true, false);
			applyH2(q, eVars.F[q], H3Y[q]);
		}
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			logPrintf("Difference %g.\n", nrm2(H1Y[q]-H2Y[q])/nrm2(H2Y[q]));
			logPrintf("Difference %g.\n", nrm2(H2Y[q]-H3Y[q])/nrm2(H2Y[q]));
			printCB(H1Y[q]);
			printCB(H2Y[q]);
			printCB(H3Y[q]);
		}
	} else if (test == 2) { //FD test d(Vexc)
		int ncount = eVars.n.size();
		ScalarFieldArray Vxc(ncount), Vxc2(ncount), n(ncount), n2(ncount), dn(ncount), dVxc(ncount), Dvxc(ncount);
		ScalarField v1,v2,dv;
		n = getn(Y);
		n2 =getn(Y2);
		dn = n2-n;
		e.exCorr.getVxcSimplified(addnXC(n), &Vxc, false, &v1);
		e.exCorr.getVxcSimplified(addnXC(n2), &Vxc2, false, &v2);
		e.exCorr.getdVxc(addnXC(n), &dVxc, false, &dv, 0, 0, dn);
		Dvxc = Vxc2-Vxc;
		logPrintf("Difference %g.\n", nrm2(Dvxc[0]-dVxc[0])/nrm2(dVxc[0]));
		printV(Dvxc[0]*(1/h));
		printV(dVxc[0]*(1/h));
		logPrintf("Difference %g.\n", nrm2((v2-v1)-dv)/nrm2(dv));
		printV((v2-v1)*(1/h));
		printV(dv*(1/h));
		logPrintf("Density");
		printV(n[0]);
		printV(dn[0]);
	} else if (test == 3) { //FD test dn
		ScalarFieldArray Vxc, n, n2, dnnum, dnanal;
		n = getn(Y);
		n2 =getn(Y2);
		dnnum = n2 - n;
		dnanal = getdn(&dY, &Y);
		logPrintf("Difference %g.\n", nrm2(dnnum[0]-dnanal[0])/nrm2(dnnum[0]));
		printV(dnnum[0]);
		printV(dnanal[0]);
		printCB(Y[0]);
		printCB(dY[0]);
	} else if (test == 4) { //FD test dE_psipsi
		Energies ener;
		eVars.C = Y;
		//for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			//e.iInfo.project(eVars.C[q], eVars.VdagC[q]);
		//}
		eVars.n = eVars.calcDensity();
		//eVars.n = eVars.calcDensity();
		//eVars.EdensityAndVscloc(ener);
		//eVars.elecEnergyAndGrad(ener);
		getGrad(&Grad1, Y);

		eVars.C = Y2;
		eVars.n = eVars.calcDensity();
		//for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			//e.iInfo.project(eVars.C[q], eVars.VdagC[q]);
		//}
		//eVars.n = eVars.calcDensity();
		//eVars.EdensityAndVscloc(ener);
		//eVars.elecEnergyAndGrad(ener);
		getGrad(&Grad2, W2);

		pInfo.dY = dW;
		eVars.C = Y;
		eVars.n = eVars.calcDensity();
		//for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			//e.iInfo.project(eVars.C[q], eVars.VdagC[q]);
		//}
		//eVars.n = eVars.calcDensity();
		//eVars.EdensityAndVscloc(ener);
		//eVars.elecEnergyAndGrad(ener);

		logPrintf("Test wut");
		compute(&pg, 0);


		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			logPrintf("Difference %g", nrm2((Grad2[q]-Grad1[q])-pg.dY[q])/nrm2(pg.dY[q]));
			printCB(Grad2[q]-Grad1[q]);
			printCB(pg.dY[q]);
		}
	} else if (test == 5) { //Test grad C = 0
		getGrad(&Y, Y);
		logPrintf("Norm %g.\n", nrm2(Y[0]));
		logPrintf("Norm %g.\n", nrm2(eVars.C[0]));
		logPrintf("Norm %g.\n", nrm2(Y[0])/nrm2(eVars.C[0]));
		printCB(Y[0]);
	} else if (test == 6) { //Compare Vxc
		int ncount = eVars.n.size();
		ScalarFieldArray Vxc(ncount), Vxc2(ncount), n(ncount);
		n = getn(Y);
		e.exCorr(addnXC(n), &Vxc, false);
		e.exCorr.getVxcSimplified(addnXC(n), &Vxc2, false);
		logPrintf("Difference %g.\n", nrm2(Vxc[0]-Vxc2[0])/nrm2(Vxc[0]));
		printV(Vxc[0]);
		printV(Vxc2[0]);
	} else if (test == 7) { //FD test hamiltonian

		/*for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			Energies ener;
			eVars.C = Y;
			eVars.n = getn(Y);
			applyH(q, eVars.F[q], H1Y[q], Y[q]);
			eVars.C = Y2;
			eVars.n = getn(Y2);
			applyH(q, eVars.F[q], H2Y[q], Y[q]);

			eVars.C = Y;
			eVars.n = getn(Y);
			eVars.EdensityAndVscloc(ener);
			dHpsi(dHY[q], Y[q], getdn(&dY, &Y));
			logPrintf("Norm %g.\n", nrm2(H2Y[q]-H1Y[q]));
			logPrintf("Difference %g.\n", nrm2(H2Y[q]-H1Y[q]-dHY[0])/nrm2(dHY[q]));
			printCB(H2Y[q]-H1Y[q]);
			printCB(dHY[q]);
		}*/
	} else if (test == 8) { //Compare getn
		ScalarFieldArray n1;
		ScalarFieldArray n2;
		eVars.C = Y;
		n1 = eVars.calcDensity();
		n2 = getn(Y);

		printV(n1[0]);
		printV(n2[0]);
	} else if (test == 9) { //FD test dC

		pInfo.dY = dW;
		eVars.C = Y;
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			pInfo.dU[q] = (pInfo.dY[q]^O(eVars.C[q])) + (eVars.C[q]^O(pInfo.dY[q]));
			pInfo.dC[q] = pInfo.dY[q];
			pInfo.dC[q] -= 0.5*(eVars.C[q]*pInfo.dU[q]);
			logPrintf("Difference %g.\n", nrm2(Y2[q]-Y[q]-pInfo.dC[q])/nrm2(pInfo.dC[q]));
			printCB(Y2[q]-Y[q]);
			printCB(pInfo.dC[q]);
		}
	} else if (test == 10) { //FD test dPsiTau
		Energies ener;
		eVars.C = Y;
		eVars.n = eVars.calcDensity();
		ScalarFieldArray Vext(eVars.n.size()), dVext(eVars.n.size());
		nullToZero(Vext, e.gInfo);
		nullToZero(dVext, e.gInfo);
		initRandomFlat(Vext[0]);
		initRandomFlat(dVext[0]);
		dVext[0] = dVext[0]*h;
		pInfo.dVext = dVext;

	    logPrintf("Test 10 A.\n");
		eVars.Vexternal = Vext;
		getGrad(&Grad1, Y);

		eVars.Vexternal = Vext+dVext;
		getGrad(&Grad2, Y);

		eVars.Vexternal = Vext;
		logPrintf("Test wut");
		calcdGradTau();

	    logPrintf("Test 10 B.\n");

		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			logPrintf("Difference %g", nrm2((Grad2[q]-Grad1[q])-pInfo.dGradTau[q])/nrm2(pInfo.dGradTau[q]));
			printCB(Grad2[q]-Grad1[q]);
			printCB(pInfo.dGradTau[q]);
		}
	}
	die("Test complete.\n");
}

//Compare applyH to applyHamiltonian
//Compare excorr functions
//Finite difference test dH
//Finite difference test dn
//Finite difference test dGrad psi/tau against Grad
//Test Grad(C) = 0 at converged wfns
