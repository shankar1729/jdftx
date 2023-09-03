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


TestPerturbation::TestPerturbation(Everything& e, PerturbationSolver& ps) : e(e), eVars(e.eVars), eInfo(e.eInfo), pInfo(e.vptInfo), ps(ps), spring(e) {}

bool TestPerturbation::compareHamiltonians() {
	logPrintf("\nTesting Hamiltonian function.\n");
	ColumnBundle tmpC;
	Energies ener;

	std::vector<ColumnBundle> H1C, H2C;

	setup(H1C); setup(H2C);

	setState(C1);

	e.iInfo.augmentDensityGridGrad(eVars.Vscloc);
	pInfo.E_nAug_cached = e.iInfo.getE_nAug();
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		ps.applyH(eInfo.qnums[q], eVars.F[q], H1C[q], C1[q]);
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

bool TestPerturbation::testGradientIsZero() {
	logPrintf("\nTesting gradient of energy is zero at minimum.\n");
	std::vector<ColumnBundle> grad, HC;
	setup(grad); setup(HC);

	setState(Cmin);

	ps.getGrad(&grad, eVars.C);

	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		ps.applyH(eInfo.qnums[q], eVars.F[q], HC[q], Cmin[q]);
	}

	logPrintf("||\u2207E|| %g.\n", nrm2(grad[0]));
	logPrintf("||HC|| %g.\n", nrm2(HC[0]));

	double ratio = nrm2(grad[0])/nrm2(HC[0]);
	logPrintf("||\u2207E||/||HC|| %g.\n", ratio);
	//pInfo.sampleCB(grad[0]);
	return ratio < 1e-3;
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
	setState(C1);
	n1 = eVars.n;
	e.exCorr(eVars.get_nXC(), &Vxc1, false, &eVars.tau, &eVars.Vtau);
	
	setState(C2);
	n2 = eVars.n;
	e.exCorr(eVars.get_nXC(), &Vxc2, false, &eVars.tau, &eVars.Vtau);
	
	dn = n2-n1;
	
	setState(C1);
	e.exCorr.getdVxc(eVars.get_nXC(), &dVxc_anal, false, 0, 0, dn);
	dVxc_num = Vxc2-Vxc1;
	//dvnum = v2-v1;
	double delta = nrm2(dVxc_num[0]-dVxc_anal[0])/nrm2(dVxc_anal[0]);
	//double deltav = nrm2(dvnum-dv)/nrm2(dv);

	logPrintf("Difference %g.\n", delta);
	//logPrintf("Difference test %g.\n", deltav);
	return delta < 1e-6;
}

bool TestPerturbation::FDTest_dn() {
	/*logPrintf("\nCommencing finite difference test of density.\n");
	ScalarFieldArray Vxc, n1, n2, dnnum, dnanal;
	n1 = ps.getn(C1);
	n2 = ps.getn(C2);
	dnnum = n2 - n1;
	dnanal = ps.getdn(&dC, &C1);
	double delta = nrm2(dnnum[0]-dnanal[0])/nrm2(dnnum[0]);
	logPrintf("Difference %g.\n", delta);
	pInfo.sampleField(dnnum[0]);
	pInfo.sampleField(dnanal[0]);
	pInfo.sampleCB(C1[0]);
	pInfo.sampleCB(dC[0]);
	return delta < 1e-6;*/


	logPrintf("\nCommencing finite difference test of density.\n");
	ScalarFieldArray Vxc, n1, n2, dnnum, dnanal;
	
	setState(C1);
	n1 = clone(eVars.n);
	
	setState(C2);
	n2 = clone(eVars.n);
	
	dnnum = n2 - n1;
	
	dnanal.resize(eVars.n.size());
	nullToZero(dnanal, e.gInfo);
	ps.getdn(dnanal, &dC, &C1);
	double delta = nrm2(dnnum[0]-dnanal[0])/nrm2(dnnum[0]);
	logPrintf("Difference %g.\n", delta);
	pInfo.sampleField(dnnum[0], "dn_num");
	pInfo.sampleField(dnanal[0], "dn_anal");
	pInfo.sampleCB(C1[0], "C");
	pInfo.sampleCB(dC[0], "dC");
	return delta < 1e-6;
}

bool TestPerturbation::FDTest_dVscloc() {
	logPrintf("\nCommencing finite difference test of Vscloc.\n");
	ScalarFieldArray Vscloc1, Vscloc2, dn, dVscloc_num, dVscloc_anal;
	
	setState(C1);
	Vscloc1 = clone(eVars.Vscloc);
	
	setState(C2);
	Vscloc2 = clone(eVars.Vscloc);
	
	dVscloc_num = Vscloc2 - Vscloc1;

	setState(C1);
	dVscloc_anal.resize(Vscloc1.size());
	dn.resize(eVars.n.size());
	nullToZero(dn, e.gInfo);
	ps.getdn(dn, &dC, &C1);
	ps.getdVsclocPsi(dn, dVscloc_anal);
	double delta = nrm2(dVscloc_num[0]-dVscloc_anal[0])/nrm2(dVscloc_anal[0]);
	logPrintf("Difference %g.\n", delta);
	pInfo.sampleField(dVscloc_num[0], "dVscloc_num");
	pInfo.sampleField(dVscloc_anal[0], "dVscloc_anal");
	return delta < 1e-6;
}


bool TestPerturbation::FDTest_dnatom() {
	if (!mode->isUltrasoft(e.iInfo))
		return true;
	
	logPrintf("\nCommencing finite difference test of density w.r.t. atomic positions.\n");
	ScalarFieldArray dn1, dn2, dn_num, dn_anal;
	pInfo.datom = mode;
	
	setAtpos1();
	dn1 = clone(eVars.n);
	
	setAtpos2();
	dn2 = clone(eVars.n);
	dn_num = dn2-dn1;

	setAtpos1();
	ps.calcdGradTau();
	dn_anal = mode->dnatom;
	//dn_anal = ps.getdnatom();
	pInfo.datom = 0;
	
	double delta = nrm2(dn_num[0]-dn_anal[0])/nrm2(dn_anal[0]);
	logPrintf("Difference %g.\n", delta);
	pInfo.sampleField(dn_num[0], "dn_num");
	pInfo.sampleField(dn_anal[0], "dn_anal");
	return delta < 1e-6;
}

bool TestPerturbation::FDTest_dVsclocatom() {
	logPrintf("\nCommencing finite difference test of Vscloc w.r.t. atomic positions.\n");
	ScalarFieldArray Vscloc1, Vscloc2, dVscloc_num, dVscloc_anal;
	pInfo.datom = mode;
	
	setAtpos1();
	Vscloc1 = clone(eVars.Vscloc);
	
	setAtpos2();
	Vscloc2 = clone(eVars.Vscloc);
	
	dVscloc_num = Vscloc2 - Vscloc1;

	setAtpos1();
	ps.calcdGradTau();
	dVscloc_anal = clone(pInfo.dVsclocTau);
	//dVscloc_anal.resize(Vscloc1.size());
	//ps.getdVsclocTau(dVscloc_anal);
	pInfo.datom = 0;
	
	double delta = nrm2(dVscloc_num[0]-dVscloc_anal[0])/nrm2(dVscloc_anal[0]);
	logPrintf("Difference %g.\n", delta);
	pInfo.sampleField(dVscloc_num[0], "dVscloc_num");
	pInfo.sampleField(dVscloc_anal[0], "dVscloc_anal");
	return delta < 1e-6;
}

bool TestPerturbation::FDTest_dVlocpsatom() {
	
	/*logPrintf("\nCommencing finite difference test of Vlocps w.r.t. atomic positions.\n");
	ScalarFieldTilde Vlocps1, Vlocps2, dVlocps_num, dVlocps_anal;
	ScalarFieldArray Vxc1, Vxc2, dVxc_num, dVxc_anal;
	ScalarField nCore1, nCore2, dnCore_num, dnCore_anal;
	matrix E_naug1, E_naug2, dE_naug_num,  dE_naug_anal;
	ScalarFieldArray tmp, Vscloc;
	
	auto sp = e.iInfo.species[mode->atomdisplacement.sp];
	pInfo.datom = mode;
	
	setAtpos1();
	
	Vlocps1 = clone(e.iInfo.Vlocps);
	Vxc1 = clone(eVars.Vxc);
	if (ps.ultrasoftDensityAugRequired) nCore1 = clone(e.iInfo.nCore);
	Vscloc = clone(eVars.Vscloc);
	sp->augmentDensityGridGrad(Vscloc);
	if (ps.ultrasoftDensityAugRequired) E_naug1 = sp->getE_nAug();
	
	setAtpos2();
	
	if (ps.ultrasoftDensityAugRequired) nCore2 = clone(e.iInfo.nCore);
	sp->augmentDensityGridGrad(Vscloc);
	if (ps.ultrasoftDensityAugRequired) E_naug2 = sp->getE_nAug();
	
	Vlocps2 = clone(e.iInfo.Vlocps);
	Vxc2 = clone(eVars.Vxc);
	dVlocps_num = Vlocps2 - Vlocps1;
	dVxc_num = Vxc2 - Vxc1;
	if (ps.ultrasoftDensityAugRequired) dnCore_num = nCore2 - nCore1;
	if (ps.ultrasoftDensityAugRequired) dE_naug_num = E_naug2 - E_naug1;
	
	setAtpos1();
	
	tmp.resize(eVars.Vscloc.size());
	ps.getdVsclocTau(tmp);
	dVlocps_anal = ps.dVlocpsA;
	dVxc_anal = ps.dVxcA;
	if (ps.ultrasoftDensityAugRequired) dnCore_anal = ps.dnCoreA;
	sp->augmentDensityGridGradDeriv(Vscloc, mode->atomdisplacement.at, &m.dirLattice);
	if (ps.ultrasoftDensityAugRequired) dE_naug_anal = sp->getE_nAug();
	pInfo.datom = 0;
	
	double delta = nrm2(dVlocps_num-dVlocps_anal)/nrm2(dVlocps_anal);
	//double delta2 = nrm2(dVxc_num[0]-dVxc_anal[0])/nrm2(dVxc_anal[0]);
	double delta3;
	double delta4;
	if (ps.ultrasoftDensityAugRequired) delta3 = nrm2(dnCore_num-dnCore_anal)/nrm2(dnCore_anal);
	if (ps.ultrasoftDensityAugRequired)  delta4 = nrm2(dE_naug_num-dE_naug_anal)/nrm2(dE_naug_anal);
	logPrintf("Difference %g.\n", delta);
	//logPrintf("Difference %g.\n", delta2);
	if (ps.ultrasoftDensityAugRequired) logPrintf("Difference %g.\n", delta3);
	if (ps.ultrasoftDensityAugRequired) logPrintf("Difference %g.\n", delta4);
	pInfo.sampleField(dVlocps_num);
	pInfo.sampleField(dVlocps_anal);
	//pInfo.sampleField(dVxc_num[0]);
	//pInfo.sampleField(dVxc_anal[0]);
	if (ps.ultrasoftDensityAugRequired) pInfo.sampleField(dnCore_num);
	if (ps.ultrasoftDensityAugRequired) pInfo.sampleField(dnCore_anal);
	if (ps.ultrasoftDensityAugRequired) pInfo.sampleMat(dE_naug_num);
	if (ps.ultrasoftDensityAugRequired) pInfo.sampleMat(dE_naug_anal);
	return delta < 1e-6;*/
	
	return true;
}

bool TestPerturbation::FDTest_dgradpsi() {
	logPrintf("\nCommencing finite difference test of gradient w.r.t. psi.\n");
	std::vector<ColumnBundle> Grad1, Grad2;
	PerturbationGradient dGrad_anal, v;
	std::vector<ColumnBundle> dGrad_num;
	setup(Grad1); setup(Grad2); setup(dGrad_num);

	Energies ener;
	setState(C1);
	ps.getGrad(&Grad1, C1);

	setState(Y2);
	ps.getGrad(&Grad2, Y2);

	setState(C1);
	v.init(e);
	v.X = dY;
	ps.hessian(dGrad_anal, v);
	double delta = 0;

	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		dGrad_num[q] = Grad2[q]-Grad1[q];
		delta += nrm2(dGrad_anal.X[q] - dGrad_num[q])/nrm2(dGrad_anal.X[q]);
		pInfo.sampleCB(dGrad_num[q], "dGrad_num");
		pInfo.sampleCB(dGrad_anal.X[q], "dGrad_anal");
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
		ps.applyH(eInfo.qnums[q], eVars.F[q], H1Y[q], C1[q]);
	}
	setState(C2);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		ps.applyH(eInfo.qnums[q], eVars.F[q], H2Y[q], C1[q]);
	}

	setState(C1);
	ScalarFieldArray dVscloc(eVars.Vscloc.size()), dn;
	
	dn.resize(eVars.n.size());
	nullToZero(dn, e.gInfo);
	ps.getdn(dn, &dC, &C1);
	ps.getdVsclocPsi(dn, dVscloc);
	e.iInfo.augmentDensityGridGrad(dVscloc);
	pInfo.E_nAug_dVsclocpsi = e.iInfo.getE_nAug();
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		ps.dHpsi(eInfo.qnums[q], dHY[q], C1[q], dVscloc);

		delta += nrm2(H2Y[q]-H1Y[q]-dHY[q])/nrm2(dHY[q]);
		pInfo.sampleCB(H2Y[q]-H1Y[q], "dHC_num");
		pInfo.sampleCB(dHY[q], "dHC_anal");
	}

	mpiWorld->allReduce(delta, MPIUtil::ReduceSum);
	delta /= eInfo.nStates;
	logPrintf("Difference %g.\n", delta);
	return delta < 1e-6;
}


bool TestPerturbation::FDTest_Hamiltoniandatom() {
	logPrintf("\nCommencing finite difference test of Hamiltonian w.r.t atomic positions.\n");

	std::vector<ColumnBundle> H1Y, H2Y, dHY;
	setup(H1Y); setup(H2Y); setup(dHY);

	double delta = 0;
	pInfo.datom = mode;
	
	setAtpos1();
	
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		ps.applyH(eInfo.qnums[q], eVars.F[q], H1Y[q], C1[q]);
	}
	
	setAtpos2();
	
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		ps.applyH(eInfo.qnums[q], eVars.F[q], H2Y[q], C1[q]);
	}

	setAtpos1();
	
	ps.calcdGradTau();
	
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		ps.dHtau(eInfo.qnums[q], dHY[q], C1[q], pInfo.dVsclocTau);
		delta += nrm2(H2Y[q]-H1Y[q]-dHY[q])/nrm2(dHY[q]);
		pInfo.sampleCB(H2Y[q]-H1Y[q], "dHC_num");
		pInfo.sampleCB(dHY[q], "dHC_anal");
	}
	
	pInfo.datom = 0;

	mpiWorld->allReduce(delta, MPIUtil::ReduceSum);
	delta /= eInfo.nStates;
	logPrintf("Difference %g.\n", delta);
	return delta < 1e-6;
}



bool TestPerturbation::FDTest_Overlapatom() {
	if (!mode->isUltrasoft(e.iInfo))
		return true;
	
	logPrintf("\nCommencing finite difference test of Overlap w.r.t atomic positions.\n");

	std::vector<ColumnBundle> O1Y, O2Y, dOY_anal, dOY_num;
	setup(O1Y); setup(O2Y); setup(dOY_anal); setup(dOY_num);

	double delta = 0;
	pInfo.datom = mode;
	
	setAtpos1();
	
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		O1Y[q] = O(C1[q]);
	}
	
	setAtpos2();
	
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		O2Y[q] = O(C1[q]);
		dOY_num[q] = O2Y[q]-O1Y[q];
	}
	
	setAtpos1();
	
	ps.calcdGradTau();
	
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		dOY_anal[q] = C1[q].similar();
		dOY_anal[q].zero();
		e.iInfo.species[pInfo.datom->mode.sp]->augmentOverlapDeriv(C1[q], dOY_anal[q], pInfo.datom->Vatom_cached[q], pInfo.datom->dVatom_cached[q]);
		//ps.dHtau(eInfo.qnums[q], dOY_anal[q], C1[q], pInfo.dVsclocatom);
		delta += nrm2(dOY_num[q]-dOY_anal[q])/nrm2(dOY_anal[q]);
		pInfo.sampleCB(dOY_num[q], "dOY_num");
		pInfo.sampleCB(dOY_anal[q], "dOY_anal");
	}
	
	pInfo.datom = 0;

	mpiWorld->allReduce(delta, MPIUtil::ReduceSum);
	delta /= eInfo.nStates;
	logPrintf("Difference %g.\n", delta);
	return delta < 1e-6;
}


bool TestPerturbation::FDTest_dV() {
	logPrintf("\nCommencing finite difference test of nonlocal projectors.\n");

	std::vector<matrix> VdagC1(eInfo.nStates), VdagC2(eInfo.nStates), dVdagC(eInfo.nStates);
	
	double delta = 0;
	pInfo.datom = mode;
	
	auto sp = e.iInfo.species[mode->mode.sp];
	setAtpos1();
	
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		auto Vatom = sp->getV(C1[q], mode->mode.at);
		VdagC1[q] = (*Vatom)^C1[q];
		pInfo.sampleCB(*Vatom, "Vatom");
	}
	
	setAtpos2();
	
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		auto Vatom = sp->getV(C1[q], mode->mode.at);
		VdagC2[q] = (*Vatom)^C1[q];
		pInfo.sampleCB(*Vatom, "Vatom perturbed");
	}

	setAtpos1();
	
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		//auto dVatom = spe->getV(C1[q], mode->atomdisplacement.at, &m.dir);
		auto dVatom = sp->getV(C1[q], mode->mode.at);
		dVdagC[q] = -D(*dVatom, mode->mode.dirCartesian)^C1[q];
		
		
		delta += nrm2(VdagC2[q]-VdagC1[q]-dVdagC[q])/nrm2(dVdagC[q]);
		pInfo.sampleMat(VdagC2[q]-VdagC1[q], "dVdagC_num");
		pInfo.sampleMat(dVdagC[q], "dVdagC_anal");
	}
	pInfo.datom = 0;

	mpiWorld->allReduce(delta, MPIUtil::ReduceSum);
	delta /= eInfo.nStates;
	logPrintf("Difference %g.\n", delta);
	return delta < 1e-6;
}

bool TestPerturbation::FDTest_dC() {
	logPrintf("\nCommencing finite difference test of orthonormalized wfns.\n");
	eVars.C = C1;
	double delta = 0;
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		pInfo.dU[q] = 2*dagger_symmetrize(dY[q]^O(eVars.C[q]));
		pInfo.dC[q] = dY[q];
		pInfo.dC[q] -= 0.5*(eVars.C[q]*pInfo.dU[q]);
		delta += nrm2(dC[q]-pInfo.dC[q])/nrm2(pInfo.dC[q]);
		pInfo.sampleCB(dC[q], "dC_num");
		pInfo.sampleCB(pInfo.dC[q], "dC_anal");
	}
	mpiWorld->allReduce(delta, MPIUtil::ReduceSum);
	delta /= eInfo.nStates;
	logPrintf("Difference %g.\n", delta);
	return delta < 1e-5;
}

bool TestPerturbation::FDTest_dCatom() {
	if (!mode->isUltrasoft(e.iInfo))
		return true;
	
	logPrintf("\nCommencing finite difference test of dC w.r.t. atomic positions.\n");
	
	pInfo.datom = mode;
	setAtpos1();
	ps.calcdGradTau();
	pInfo.datom = 0;

	double delta = 0;
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		delta += nrm2((C2atom[q]-C1[q])-mode->dCatom[q])/nrm2(mode->dCatom[q]);
		pInfo.sampleCB(C2atom[q]-C1[q], "dC_num");
		pInfo.sampleCB(mode->dCatom[q], "dC_anal");
	}

	mpiWorld->allReduce(delta, MPIUtil::ReduceSum);
	delta /= eInfo.nStates;
	logPrintf("Difference %g\n", delta);
	return delta < 1e-5;
}

bool TestPerturbation::FDTest_dgradtau() {
	logPrintf("\nCommencing finite difference test of gradient w.r.t. tau.\n");
	std::vector<ColumnBundle> Grad1, Grad2;
	setup(Grad1); setup(Grad2);

	ScalarFieldArray Vext(eVars.n.size());
	nullToZero(Vext, e.gInfo);
	pInfo.dVext = VextMode;

	eVars.Vexternal = Vext;
	setState(C1);
	ps.getGrad(&Grad1, C1);

	eVars.Vexternal = Vext+VextMode->dVext;
	setState(C1);
	ps.getGrad(&Grad2, C1);

	eVars.Vexternal = Vext;
	setState(C1);
	ps.calcdGradTau();
	pInfo.dVext = 0;


	double delta = 0;
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		delta += nrm2((Grad2[q]-Grad1[q])-pInfo.dGradTau.X[q])/nrm2(pInfo.dGradTau.X[q]);
		pInfo.sampleCB(Grad2[q]-Grad1[q], "dGrad_num");
		pInfo.sampleCB(pInfo.dGradTau.X[q], "dGrad_anal");
	}

	mpiWorld->allReduce(delta, MPIUtil::ReduceSum);
	delta /= eInfo.nStates;
	logPrintf("Difference %g\n", delta);
	return delta < 1e-6;
}

bool TestPerturbation::FDTest_dgradtauatom() {
	logPrintf("\nCommencing finite difference test of gradient w.r.t. atomic positions.\n");
	std::vector<ColumnBundle> Grad1, Grad2;
	setup(Grad1); setup(Grad2);
	pInfo.datom = mode;

	setAtpos1();
	ps.getGrad(&Grad1, C1);
	
	setAtpos2();
	ps.getGrad(&Grad2, C1);
	
	setAtpos1();
	ps.calcdGradTau();
	pInfo.datom = 0;

	double delta = 0;
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		delta += nrm2((Grad2[q]-Grad1[q])-pInfo.dGradTau.X[q])/nrm2(pInfo.dGradTau.X[q]);
		pInfo.sampleCB(Grad2[q]-Grad1[q], "dGrad_num");
		pInfo.sampleCB(pInfo.dGradTau.X[q], "dGrad_anal");
	}

	mpiWorld->allReduce(delta, MPIUtil::ReduceSum);
	delta /= eInfo.nStates;
	logPrintf("Difference %g\n", delta);
	return delta < 1e-6;
}

bool TestPerturbation::FDTest_dsqV() {
	logPrintf("\nTesting second derivative of nonlocal projectors.\n");
	std::vector<ColumnBundle> V(eInfo.nStates), V1(eInfo.nStates), V2(eInfo.nStates), V12(eInfo.nStates), dsqV_num(eInfo.nStates), dsqV_anal(eInfo.nStates);
	
	auto sp = e.iInfo.species[modeA->mode.sp];
	
	setdsqposnn();
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		V[q] = *sp->getV(eVars.C[q], modeA->mode.at);
	
	setdsqposnp();
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		V1[q] = *sp->getV(eVars.C[q], modeA->mode.at);
	
	setdsqpospn();
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		V2[q] = *sp->getV(eVars.C[q], modeA->mode.at);
	
	setdsqpospp();
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		V12[q] = *sp->getV(eVars.C[q], modeA->mode.at);
	
	setdsqpos0();
	double delta = 0;
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		dsqV_anal[q] = D(D(V[q], modeA->mode.dirCartesian), modeB->mode.dirCartesian);
		dsqV_num[q] = V12[q]-V1[q]-V2[q]+V[q];
		
		delta += nrm2(dsqV_num[q] - dsqV_anal[q])/nrm2(dsqV_anal[q]);
		pInfo.sampleCB(dsqV_num[q], "dsqV_num");
		pInfo.sampleCB(dsqV_anal[q], "dsqV_anal");
	}
	
	mpiWorld->allReduce(delta, MPIUtil::ReduceSum);
	delta /= eInfo.nStates;
	logPrintf("Difference %g\n", delta);
	return delta < 1e-4;
}

bool TestPerturbation::FDTest_dsqEnl() {
	logPrintf("\nTesting second derivative of Enl.\n");
	double Enn, Enp, Epn, Epp, dsqEanal, dsqEnum;
	
	auto sp = e.iInfo.species[modeA->mode.sp];
	
	setdsqposnn();
	Enn = e.ener.E["Enl"];
	
	setdsqposnp();
	Enp = e.ener.E["Enl"];
	
	setdsqpospn();
	Epn = e.ener.E["Enl"];
	
	setdsqpospp();
	Epp = e.ener.E["Enl"];
	
	setdsqpos0();
	dsqEnum = Epp - Enp - Epn + Enn;
	dsqEanal = spring.dsqEnl(modeA, modeB);
	
	double delta = (dsqEnum-dsqEanal)/dsqEanal;
	
	logPrintf("Eanal: %g, Enum: %g\n", dsqEanal, dsqEnum);
	logPrintf("Difference %g\n", delta);
	return delta < 1e-4;
}

bool TestPerturbation::FDTest_dsqEloc() {
	logPrintf("\nTesting second derivative of Eloc.\n");
	double Enn, Enp, Epn, Epp, dsqEanal, dsqEnum;
	
	auto sp = e.iInfo.species[modeA->mode.sp];
	
	setdsqposnn();
	Enn = e.ener.E["Eloc"];
	
	setdsqposnp();
	Enp = e.ener.E["Eloc"];
	
	setdsqpospn();
	Epn = e.ener.E["Eloc"];
	
	setdsqpospp();
	Epp = e.ener.E["Eloc"];
	
	setdsqpos0();
	dsqEnum = Epp - Enp - Epn + Enn;
	pInfo.datom = modeA;
	ps.calcdGradTau();
	dsqEanal = spring.dsqEloc(modeA, modeB);
	pInfo.datom = 0;
	
	double delta = (dsqEnum-dsqEanal)/dsqEanal;
	
	logPrintf("Eanal: %g, Enum: %g\n", dsqEanal, dsqEnum);
	logPrintf("Difference %g\n", delta);
	return delta < 1e-4;
}

bool TestPerturbation::FDTest_dsqEH() {
	return true;
}

bool TestPerturbation::FDTest_dsqExc() {
	return true;
}

bool TestPerturbation::FDTest_dsqExccore() {
	return true;
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
	setup(C2atom);
	
	pInfo.dGradTau.init(e);
	
	int sp = 0;
	int atom = 1;
	pos1 = e.iInfo.species[sp]->atpos[atom];
	vector3<> deltapos = vector3<double>(1.0,0.0,0.0)*h;
	mode = std::make_shared<AtomPerturbation>(sp, atom, deltapos, e);
 	pos2 = pos1+mode->mode.dirLattice;
	
	vector3<> deltaposA = vector3<double>(1e-3,0,0);
	vector3<> deltaposB = vector3<double>(0,1e-3,0);
	modeA = std::make_shared<AtomPerturbation>(sp, atom, deltaposA, e);
	modeB = std::make_shared<AtomPerturbation>(sp, atom, deltaposB, e);

	//logPrintf("%d", eInfo.qStart);
	//pInfo.sampleCB(Y1[0]);

	randomize(Y1, eInfo);
	randomize(dY, eInfo);
	
	VextMode = std::make_shared<VextPerturbation>(e);
	ScalarFieldArray dVext(eVars.n.size());
	nullToZero(dVext, e.gInfo);
	initRandomFlat(dVext[0]);
	dVext[0] = dVext[0]*h;
	VextMode->dVext = dVext;

	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		matrix Usqrtinv;
		
		Usqrtinv = invsqrt(Y1[q]^O(Y1[q]));
		//Y1[q] = Cmin[q];
		Y1[q] = Y1[q]*Usqrtinv;
		C1[q] = Y1[q];
		
		
		dY[q] = dY[q]*h;
		Y2[q] = Y1[q] + dY[q];

		C2[q] = Y2[q];
		Usqrtinv = invsqrt(C2[q]^O(C2[q]));
		C2[q] = C2[q]*Usqrtinv;
		
		dC[q] = C2[q] - C1[q];
		
		C2atom[q] = C1[q];
	}

	setAtpos2();
	
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		matrix Usqrtinv = invsqrt(C2atom[q]^O(C2atom[q]));
		C2atom[q] = C2atom[q]*Usqrtinv;
	}
	
	setAtpos1();
	
	bool testPassed = true;

	testPassed &= compareHamiltonians();
	testPassed &= testGradientIsZero();
	testPassed &= FDTest_dVxc();
	testPassed &= FDTest_dn();
	testPassed &= FDTest_dVscloc();
	testPassed &= FDTest_dnatom();
	testPassed &= FDTest_dVsclocatom();
	//testPassed &= FDTest_dVlocpsatom();
	testPassed &= FDTest_dgradpsi();
	testPassed &= FDTest_Hamiltonian();
	testPassed &= FDTest_Hamiltoniandatom();
	testPassed &= FDTest_Overlapatom();
	testPassed &= FDTest_dV();
	testPassed &= FDTest_dC();
	testPassed &= FDTest_dCatom();
	testPassed &= FDTest_dgradtau();
	testPassed &= FDTest_dgradtauatom();
	//testPassed &= FDTest_dsqV();
	//testPassed &= FDTest_dsqEnl();
	//testPassed &= FDTest_dsqEloc();
	//testPassed &= FDTest_dsqEH();
	//testPassed &= FDTest_dsqExc();
	//testPassed &= FDTest_dsqExccore();

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
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		eVars.C[q] = C[q];
		eVars.orthonormalize(q, 0, true);
	}
	e.iInfo.update(e.ener);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		//eVars.VdagC[q].clear();
		//e.iInfo.project(eVars.C[q], eVars.VdagC[q]);
	}
	eVars.elecEnergyAndGrad(e.ener, 0, 0, true);
	
	ps.updateExcorrCache();
	ps.updateHC();
	ps.updateNonlocalDerivs();
}

void TestPerturbation::setAtpos1() {
	auto sp = e.iInfo.species[mode->mode.sp];
	sp->atpos[mode->mode.at] = pos1;
	mpiWorld->bcastData(sp->atpos);
	sp->sync_atpos();
	setState(C1);
}

void TestPerturbation::setAtpos2() {
	auto sp = e.iInfo.species[mode->mode.sp];
	sp->atpos[mode->mode.at] = pos2;
	mpiWorld->bcastData(sp->atpos);
	sp->sync_atpos();
	setState(C2atom);
}


void TestPerturbation::setdsqpos0() {
	auto sp = e.iInfo.species[modeA->mode.sp];
	sp->atpos[modeA->mode.at] = pos1;
	mpiWorld->bcastData(sp->atpos);
	sp->sync_atpos();
	setState(C1);
}

void TestPerturbation::setdsqposnn() {
	auto sp = e.iInfo.species[modeA->mode.sp];
	sp->atpos[modeA->mode.at] = pos1 - 0.5*modeA->mode.dirLattice - 0.5*modeB->mode.dirLattice;
	mpiWorld->bcastData(sp->atpos);
	sp->sync_atpos();
	setState(C1);
}

void TestPerturbation::setdsqposnp() {
	auto sp = e.iInfo.species[modeA->mode.sp];
	sp->atpos[modeA->mode.at] = pos1 - 0.5*modeA->mode.dirLattice + 0.5*modeB->mode.dirLattice;
	mpiWorld->bcastData(sp->atpos);
	sp->sync_atpos();
	setState(C1);
}

void TestPerturbation::setdsqpospn() {
	auto sp = e.iInfo.species[modeB->mode.sp];
	sp->atpos[modeB->mode.at] = pos1 + 0.5*modeA->mode.dirLattice - 0.5*modeB->mode.dirLattice;
	mpiWorld->bcastData(sp->atpos);
	sp->sync_atpos();
	setState(C1);
}

void TestPerturbation::setdsqpospp() {
	auto sp = e.iInfo.species[modeA->mode.sp];
	sp->atpos[modeA->mode.at] = pos1 + 0.5*modeA->mode.dirLattice + 0.5*modeB->mode.dirLattice;
	mpiWorld->bcastData(sp->atpos);
	sp->sync_atpos();
	setState(C1);
}
