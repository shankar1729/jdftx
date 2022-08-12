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

ScalarFieldArray PerturbationSolver::getdn()
{	ScalarFieldArray density(eVars.n.size());
	//Runs over all states and accumulates density to the corresponding spin channel of the total density
	//e.iInfo.augmentDensityInit();
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{	density += e.eInfo.qnums[q].weight * diagouterI(eVars.F[q], eVars.C[q], density.size(), &e.gInfo);
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



PerturbationSolver::PerturbationSolver(Everything& e) : e(e), eVars(e.eVars), eInfo(e.eInfo), pInfo(e.vptInfo) {}

double PerturbationSolver::compute(PerturbationGradient* grad, PerturbationGradient* Kgrad)
{
	logPrintf("testC\n");
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

	ScalarFieldArray dn = getdn();

	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		const QuantumNumber& qnum = e.eInfo.qnums[q];
		ColumnBundle dHC, HC, HdC, dwGradEq, d_HCFUmhalf;
		matrix dHtilde, HtildeCommF, dHtildeCommF;
		diagMatrix F = eVars.F[q];

		dHpsi(dHC, eVars.C[q], dn);
		dHtau(dHC, eVars.C[q]);
		applyH(q, F, HC, eVars.C[q]);
		applyH(q, F, HdC, pInfo.dC[q]);

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
	    dwGradEq += O(eVars.C[q]*(pInfo.dC[q]^HC*F));
	    dwGradEq += d_HCFUmhalf;
	    dwGradEq -= O(eVars.C[q]*(eVars.C[q]^d_HCFUmhalf));
	    dwGradEq += O(0.5*pInfo.dC[q]*HtildeCommF);
	    dwGradEq += O(0.5*eVars.C[q]*dHtildeCommF);
	    dwGradEq -= O(0.125*eVars.C[q]*(HtildeCommF*pInfo.dU[q]+pInfo.dU[q]*HtildeCommF));
	    dwGradEq = dwGradEq*qnum.weight;
	}


	/*for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		if (grad) {
			grad->dY[q] = e.vptInfo.dY[q] - eVars.C[q];
			if (Kgrad) {
				Kgrad->dY[q] = grad->dY[q];
			}
		}
	}*/

	//Objective function is 1/2 x^T*A*x - b^T*x
	double ener = 0.5*dot(e.vptInfo.dY, e.vptInfo.dY, eInfo.qStart, eInfo.qStop) - dot(eVars.C, e.vptInfo.dY, eInfo.qStart, eInfo.qStop);
	return ener;
}

void PerturbationSolver::step(const PerturbationGradient& dir, double alpha)
{	assert(dir.eInfo == &eInfo);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		axpy(alpha, dir.dY[q], e.vptInfo.dY[q]);
	}
	logPrintf("testD\n");
}

bool PerturbationSolver::report(int iter)
{
	/*
	logPrintf("Iteration %d\n", iter);
	assert(eInfo.isMine(0));
	logPrintf("%d\n", e.vptInfo.dY.size());
	complexScalarFieldTilde dat = e.vptInfo.dY[0].getColumn(0, 0);
	for (int i = 0; i < 10; i++) {
		logPrintf("%f\n", (dat->data(false)[i]).x);
	}*/
	return false;
}

void PerturbationSolver::constrain(PerturbationGradient&)
{
	return;
}

double PerturbationSolver::minimize(const MinimizeParams& params)
{
	std::cout <<"Starting";
	logPrintf("testB\n");
	if(eVars.wfnsFilename.length() == 0) {
		elecFluidMinimize(e);
	}
	pInfo.dY.resize(eInfo.nStates);
	pInfo.dC.resize(eInfo.nStates);
	pInfo.dU.resize(eInfo.nStates);

	init(e.vptInfo.dY, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	init(e.vptInfo.dC, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	randomize(e.vptInfo.dY, eInfo);
	//for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
	//	e.vptInfo.dY[q] *= 0;
	//}
	return Minimizable<PerturbationGradient>::minimize(params);
}


void PerturbationSolver::applyH(int q, const diagMatrix& F, ColumnBundle& HC, ColumnBundle& C)
{	assert(C); //make sure wavefunction is available for this state
	const QuantumNumber& qnum = e.eInfo.qnums[q];
	std::vector<matrix> HVdagC(e.iInfo.species.size());
	std::vector<matrix> VdagC(e.iInfo.species.size());

	//Propagate grad_n (Vscloc) to HCq (which is grad_Cq upto weights and fillings) if required
	HC += Idag_DiagV_I(C, eVars.Vscloc);

	//Kinetic energy:
	HC += (-0.5) * L(C);

	//Nonlocal pseudopotentials:
	e.iInfo.project(C, VdagC); //update the atomic projections
	e.iInfo.EnlAndGrad(qnum, F, VdagC, HVdagC);
	e.iInfo.projectGrad(HVdagC, C, HC);
}

void PerturbationSolver::dHpsi(ColumnBundle& HC, ColumnBundle& C, ScalarFieldArray dn)
{	assert(C);

	ScalarFieldArray dVscloc;

	getdVsclocPsi(dn, dVscloc);

	HC += Idag_DiagV_I(C, dVscloc);
}

void PerturbationSolver::dHtau(ColumnBundle& HC, ColumnBundle& C)
{	assert(C);

	ScalarFieldArray dVscloc;

	getdVsclocTau(dVscloc);

	HC -= Idag_DiagV_I(C, dVscloc);

	//Nonlocal pseudopotentials:
	//e.iInfo.project(C, VdagC); //update the atomic projections
	//e.iInfo.EnlAndGrad(qnum, F, VdagC, HVdagC);
	//e.iInfo.projectGrad(HVdagC, C, HC);
}

//simplified version
void PerturbationSolver::getVscloc()
{	static StopWatch watch("EdensityAndVscloc"); watch.start();
	const IonInfo& iInfo = e.iInfo;

	ScalarFieldTilde nTilde = J(eVars.get_nTot());

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
	e.exCorr(eVars.get_nXC(), &eVars.Vxc, false, &eVars.tau, &eVars.Vtau);

	for(unsigned s=0; s<eVars.Vscloc.size(); s++)
	{	eVars.Vscloc[s] = JdagOJ(eVars.Vxc[s]);
		if(s<2) //Include all the spin-independent contributions along the diagonal alone
			eVars.Vscloc[s] += Jdag(O(VsclocTilde), true);

		//External potential contributions:
		if(eVars.Vexternal.size())
			eVars.Vscloc[s] += JdagOJ(eVars.Vexternal[s]);
	}
	e.symm.symmetrize(eVars.Vscloc);
	if(eVars.Vtau[0]) e.symm.symmetrize(eVars.Vtau);
	watch.stop();
}


void PerturbationSolver::getdVsclocPsi(ScalarFieldArray dn, ScalarFieldArray& dVscloc)
{	static StopWatch watch("EdensityAndVscloc"); watch.start();

	ScalarFieldTilde dnTilde = J(eVars.get_nTot());

	// Hartree term:
	ScalarFieldTilde dVsclocTilde = (*e.coulomb)(dnTilde); //Note: external charge and nuclear charge contribute to d_vac as well (see below)

	// Exchange and correlation, and store the real space Vscloc with the odd historic normalization factor of JdagOJ:
	//e.exCorr(eVars.get_nXC(), &eVars.Vxc, false, &eVars.tau, &eVars.Vtau);
	//Second derivative exc
	ScalarFieldArray dVxc;
	e.exCorr.getdVxc(eVars.get_nXC(), &dVxc, false, &eVars.tau, &eVars.Vtau, dn);

	for(unsigned s=0; s<eVars.Vscloc.size(); s++)
	{	eVars.Vscloc[s] = JdagOJ(dVxc[s]);
		if(s<2) //Include all the spin-independent contributions along the diagonal alone
			dVscloc[s] += Jdag(O(dVsclocTilde), true);
	}
	e.symm.symmetrize(dVscloc);
	watch.stop();
}



void PerturbationSolver::getdVsclocTau(ScalarFieldArray& dVscloc)
{	static StopWatch watch("EdensityAndVscloc"); watch.start();

	for(unsigned s=0; s<eVars.Vscloc.size(); s++)
	{
		dVscloc[s] = pInfo.dVext[s];
	}

	e.symm.symmetrize(dVscloc);
	watch.stop();
}

//Compare applyH to applyHamiltonian
//Finite difference test dH
//Finite difference test dn
//Finite difference test dGrad psi/tau
