#include <electronic/TestPerturbation.h>
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


TestPerturbation::TestPerturbation(Everything& e) : e(e), eVars(e.eVars), eInfo(e.eInfo), pInfo(e.vptInfo) {}

bool TestPerturbation::compareHamiltonians() {
	logPrintf("\nTesting Hamiltonian function.\n");
	ColumnBundle tmpC;
	Energies ener;

	std::vector<ColumnBundle> H1C, H2C;

	setup(H1C); setup(H2C);

	setState(Cmin);

	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		eVars.n = ps->getn(C1);
		ps->applyH(eInfo.qnums[q], eVars.F[q], H1C[q], C1[q]);
	}

	setState(C1);

	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		eVars.applyHamiltonian(q, eVars.F[q], H2C[q], ener, true, false);
	}

	double delta = 0;

	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		delta += nrm2(H1C[q]-H2C[q])/nrm2(H2C[q]);
	}

	mpiWorld->allReduce(delta, MPIUtil::ReduceSum);

	delta /= eInfo.nStates;

	logPrintf("Difference %g.\n", delta);

	return delta < 1e-6;
}

bool TestPerturbation::compareVxc() {
	logPrintf("\nTesting simplified exchange-correlation.\n");
	int ncount = eVars.n.size();
	ScalarFieldArray Vxc(ncount), Vxc2(ncount), n(ncount);

	setState(Cmin);
	n = ps->getn(C1);
	e.exCorr(ps->addnXC(n), &Vxc, false);
	e.exCorr.getVxcSimplified(ps->addnXC(n), &Vxc2, false);

	double delta = nrm2(Vxc[0]-Vxc2[0])/nrm2(Vxc[0]);

	ps->printV(Vxc[0]);
	ps->printV(Vxc2[0]);
	logPrintf("Difference %g.\n", delta);
	return delta < 1e-6;
}

bool TestPerturbation::compareVscloc() {
	logPrintf("\nTesting self consistent potential.\n");
	int ncount = eVars.n.size();
	ScalarFieldArray Vscloc1(ncount), Vscloc2(ncount), n(ncount);

	setState(Cmin);
	Vscloc2 = ps->getVscloc(ps->getn(C1));

	setState(C1);
	Vscloc1 = eVars.Vscloc;


	double delta = nrm2(Vscloc1[0]-Vscloc2[0])/nrm2(Vscloc2[0]);

	ps->printV(Vscloc1[0]);
	ps->printV(Vscloc2[0]);
	logPrintf("Difference %g.\n", delta);
	return delta < 1e-6;
}

bool TestPerturbation::compareGetn() {
	logPrintf("\nTesting density.\n");
	ScalarFieldArray n1;
	ScalarFieldArray n2;

	setState(Cmin);
	n2 = ps->getn(C1);

	setState(C1);
	n1 = eVars.calcDensity();

	double delta = nrm2(n1[0]-n2[0])/nrm2(n2[0]);

	ps->printV(n1[0]);
	ps->printV(n2[0]);
	logPrintf("Difference %g.\n", delta);
	return delta < 1e-6;
}

bool TestPerturbation::testGradientIsZero() {
	logPrintf("\nTesting gradient of energy is zero at minimum.\n");
	std::vector<ColumnBundle> grad, HC;
	setup(grad); setup(HC);

	setState(Cmin);

	ps->getGrad(&grad, eVars.C);

	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		ps->applyH(eInfo.qnums[q], eVars.F[q], HC[q], Cmin[q]);
	}

	logPrintf("||\u2207E|| %g.\n", nrm2(grad[0]));
	logPrintf("||HC|| %g.\n", nrm2(HC[0]));

	double ratio = nrm2(grad[0])/nrm2(HC[0]);
	logPrintf("||\u2207E||/||HC|| %g.\n", ratio);
	//ps->printCB(grad[0]);
	return ratio < 1e-4;
}

bool TestPerturbation::FDTest_dVxc() {
	logPrintf("\nCommencing finite difference test of exchange-correlation.\n");
	int ncount = eVars.n.size();
	ScalarFieldArray Vxc1(ncount), Vxc2(ncount), n1(ncount), n2(ncount), dn(ncount), dVxc_anal(ncount), dVxc_num(ncount);
	ScalarField v1,v2,dv,dvnum;
	nullToZero(v1, e.gInfo);
	nullToZero(v2, e.gInfo);
	nullToZero(dv, e.gInfo);
	nullToZero(dvnum, e.gInfo);
	setState(Cmin);
	n1 = ps->getn(C1);
	n2 = ps->getn(C2);
	dn = n2-n1;

	e.exCorr.getVxcSimplified(ps->addnXC(n1), &Vxc1, false, &v1);
	e.exCorr.getVxcSimplified(ps->addnXC(n2), &Vxc2, false, &v2);
	e.exCorr.getdVxc(ps->addnXC(n1), &dVxc_anal, false, &dv, 0, 0, dn);
	dVxc_num = Vxc2-Vxc1;
	dvnum = v2-v1;
	double delta = nrm2(dVxc_num[0]-dVxc_anal[0])/nrm2(dVxc_anal[0]);
	//double deltav = nrm2(dvnum-dv)/nrm2(dv);

	logPrintf("Difference %g.\n", delta);
	//logPrintf("Difference test %g.\n", deltav);
	return delta < 1e-6;
}

bool TestPerturbation::FDTest_dn() {
	/*logPrintf("\nCommencing finite difference test of density.\n");
	ScalarFieldArray Vxc, n1, n2, dnnum, dnanal;
	n1 = ps->getn(C1);
	n2 = ps->getn(C2);
	dnnum = n2 - n1;
	dnanal = ps->getdn(&dC, &C1);
	double delta = nrm2(dnnum[0]-dnanal[0])/nrm2(dnnum[0]);
	logPrintf("Difference %g.\n", delta);
	ps->printV(dnnum[0]);
	ps->printV(dnanal[0]);
	ps->printCB(C1[0]);
	ps->printCB(dC[0]);
	return delta < 1e-6;*/


	logPrintf("\nCommencing finite difference test of density.\n");
	ScalarFieldArray Vxc, n1, n2, dnnum, dnanal;
	setState(C1);
	n1 = eVars.n;
	setState(C2);
	n2 = eVars.n;
	dnnum = n2 - n1;
	dnanal = ps->getdn(&dC, &C1);
	double delta = nrm2(dnnum[0]-dnanal[0])/nrm2(dnnum[0]);
	logPrintf("Difference %g.\n", delta);
	ps->printV(dnnum[0]);
	ps->printV(dnanal[0]);
	ps->printCB(C1[0]);
	ps->printCB(dC[0]);
	return delta < 1e-6;
}

bool TestPerturbation::FDTest_dVscloc() {
	logPrintf("\nCommencing finite difference test of Vscloc.\n");
	ScalarFieldArray Vscloc1, Vscloc2, dVscloc_num, dVscloc_anal;
	setState(C1);
	Vscloc1 = ps->getVscloc(eVars.n);
	setState(C2);
	Vscloc2 = ps->getVscloc(eVars.n);
	dVscloc_num = Vscloc2 - Vscloc1;

	setState(C1);
	dVscloc_anal.resize(Vscloc1.size());
	ps->getdVsclocPsi(ps->getdn(&dC, &C1), dVscloc_anal);
	double delta = nrm2(dVscloc_num[0]-dVscloc_anal[0])/nrm2(dVscloc_anal[0]);
	logPrintf("Difference %g.\n", delta);
	ps->printV(dVscloc_num[0]);
	ps->printV(dVscloc_anal[0]);
	return delta < 1e-6;
}

bool TestPerturbation::FDTest_dgradpsi() {
	logPrintf("\nCommencing finite difference test of gradient w.r.t. psi.\n");
	std::vector<ColumnBundle> Grad1, Grad2;
	PerturbationGradient dGrad_anal;
	std::vector<ColumnBundle> dGrad_num;
	setup(Grad1); setup(Grad2); setup(dGrad_num);

	Energies ener;
	setState(C1);
	ps->getGrad(&Grad1, C1);

	setState(C2);
	ps->getGrad(&Grad2, C2);

	setState(C1);
	pInfo.dY = dC;
	ps->compute(&dGrad_anal, 0);
	double delta = 0;

	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		dGrad_num[q] = Grad2[q]-Grad1[q];
		delta += nrm2(pInfo.dGradPsi[q] - dGrad_num[q])/nrm2(pInfo.dGradPsi[q]);
		ps->printCB(dGrad_num[q]);
		ps->printCB(pInfo.dGradPsi[q]);
	}


	mpiWorld->allReduce(delta, MPIUtil::ReduceSum);

	delta /= eInfo.nStates;
	logPrintf("Difference %g\n", delta);
	return delta < 1e-6;
}

bool TestPerturbation::FDTest_Hamiltonian() {
	logPrintf("\nCommencing finite difference test of Hamiltonian.\n");

	std::vector<ColumnBundle> H1Y, H2Y, dHY;
	setup(H1Y); setup(H2Y); setup(dHY);

	double delta = 0;
	setState(C1);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		ps->applyH(eInfo.qnums[q], eVars.F[q], H1Y[q], C1[q]);
	}
	setState(C2);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		ps->applyH(eInfo.qnums[q], eVars.F[q], H2Y[q], C1[q]);
	}

	setState(C1);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		ps->dHpsi(eInfo.qnums[q], dHY[q], C1[q], ps->getdn(&dC, &C1));

		delta += nrm2(H2Y[q]-H1Y[q]-dHY[q])/nrm2(dHY[q]);
		ps->printCB(H2Y[q]-H1Y[q]);
		ps->printCB(dHY[q]);
	}

	mpiWorld->allReduce(delta, MPIUtil::ReduceSum);
	delta /= eInfo.nStates;
	logPrintf("Difference %g.\n", delta);
	return delta < 1e-6;
}

bool TestPerturbation::FDTest_dC() {
	logPrintf("\nCommencing finite difference test of orthonormalized wfns.\n");
	pInfo.dY = dY;
	eVars.C = C1;
	double delta = 0;
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		pInfo.dU[q] = (pInfo.dY[q]^O(eVars.C[q])) + (eVars.C[q]^O(pInfo.dY[q]));
		pInfo.dC[q] = pInfo.dY[q];
		pInfo.dC[q] -= 0.5*(eVars.C[q]*pInfo.dU[q]);
		delta += nrm2(C2[q]-C1[q]-pInfo.dC[q])/nrm2(pInfo.dC[q]);
		ps->printCB(C2[q]-C1[q]);
		ps->printCB(pInfo.dC[q]);
	}
	mpiWorld->allReduce(delta, MPIUtil::ReduceSum);
	delta /= eInfo.nStates;
	logPrintf("Difference %g.\n", delta);
	return delta < 1e-6;
}

bool TestPerturbation::FDTest_dgradtau() {
	logPrintf("\nCommencing finite difference test of gradient w.r.t. tau.\n");
	std::vector<ColumnBundle> Grad1, Grad2;
	setup(Grad1); setup(Grad2);

	Energies ener;
	eVars.C = C1;
	eVars.n = eVars.calcDensity();
	ScalarFieldArray Vext(eVars.n.size()), dVext(eVars.n.size());
	nullToZero(Vext, e.gInfo);
	nullToZero(dVext, e.gInfo);
	initRandomFlat(Vext[0]);
	initRandomFlat(dVext[0]);
	dVext[0] = dVext[0]*h;
	pInfo.dVext = Complex(dVext);

	eVars.Vexternal = Vext;
	ps->getGrad(&Grad1, C1);

	eVars.Vexternal = Vext+dVext;
	ps->getGrad(&Grad2, C1);

	eVars.Vexternal = Vext;
	ps->calcdGradTau();


	double delta = 0;
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		delta += nrm2((Grad2[q]-Grad1[q])-pInfo.dGradTau[q])/nrm2(pInfo.dGradTau[q]);
		ps->printCB(Grad2[q]-Grad1[q]);
		ps->printCB(pInfo.dGradTau[q]);
	}

	mpiWorld->allReduce(delta, MPIUtil::ReduceSum);
	delta /= eInfo.nStates;
	logPrintf("Difference %g\n", delta);
	return delta < 1e-6;
}

void TestPerturbation::testVPT() {
	PerturbationGradient pg;
	pg.init(e);
	//std::vector<ColumnBundle> dGrad;
	Cmin = eVars.C;
	setup(Y1);
	setup(Y2);
	setup(dY);
	setup(C1);
	setup(C2);
	setup(dC);

	//logPrintf("%d", eInfo.qStart);
	//ps->printCB(Y1[0]);

	randomize(Y1, eInfo);
	randomize(dY, eInfo);

	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		dY[q] = dY[q]*h;
		Y2[q] = Y1[q] + dY[q];

		C1[q] = Y1[q];
		C2[q] = Y2[q];

		ps->orthonormalize(q, C1[q]);
		ps->orthonormalize(q, C2[q]);

		dC[q] = C2[q] - C1[q];
	}

	bool testPassed = true;

	/*Uncomment one of the following tests:*/

	testPassed &= compareHamiltonians();
	testPassed &= compareVxc();
	testPassed &= compareVscloc();
	testPassed &= compareGetn();
	testPassed &= testGradientIsZero();
	testPassed &= FDTest_dVxc();
	testPassed &= FDTest_dn();
	testPassed &= FDTest_dVscloc();
	testPassed &= FDTest_dgradpsi();
	testPassed &= FDTest_Hamiltonian();
	//testPassed &= FDTest_dC(); //Do not use this test
	testPassed &= FDTest_dgradtau();

	if (testPassed) {
		logPrintf("\nAll tests completed successfully.\n");
	}
	else {
		logPrintf("\nOne or more of the tests have failed.\n");
	}
}

void TestPerturbation::setup(std::vector<ColumnBundle> &Y) {
	init(Y, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		Y[q].zero();
	}
}


void TestPerturbation::setState(std::vector<ColumnBundle> &C) {
	eVars.C = C;
	Energies energ;
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		e.iInfo.project(eVars.C[q], eVars.VdagC[q]);
	}
	eVars.elecEnergyAndGrad(energ, nullptr, nullptr, true);
}
