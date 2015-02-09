/*-------------------------------------------------------------------
Copyright 2015 Ravishankar Sundararaman

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

#include <electronic/ElectronScattering.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/ColumnBundleTransform.h>
#include <core/LatticeUtils.h>

ElectronScattering::ElectronScattering()
: eta(0.), Ecut(0.), fCut(1e-6), omegaMax(0.)
{
}

void ElectronScattering::dump(const Everything& everything)
{	Everything& e = (Everything&)everything; //may modify everything to save memory / optimize
	
	logPrintf("\n----- Electron-electron scattering Im(Sigma) -----\n"); logFlush();

	//Update default parameters:
	if(!eta)
	{	eta = e.eInfo.kT;
		if(!eta) die("eta must be specified explicitly since electronic temperature is zero.\n");
	}
	if(!Ecut) Ecut = e.cntrl.Ecut;
	double oMin = DBL_MAX, oMax = -DBL_MAX; //occupied energy range
	double uMin = DBL_MAX, uMax = -DBL_MAX; //unoccupied energy range
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		for(int b=0; b<e.eInfo.nBands; b++)
		{	double E = e.eVars.Hsub_eigs[q][b];
			double f = e.eVars.F[q][b];
			if(f > fCut) //sufficiently occupied
			{	oMin = std::min(oMin, E);
				oMax = std::max(oMax, E);
			}
			if(f < 1.-fCut) //sufficiently unoccupied
			{	uMin = std::min(uMin, E);
				uMax = std::max(uMax, E);
			}
		}
	mpiUtil->allReduce(oMin, MPIUtil::ReduceMin);
	mpiUtil->allReduce(oMax, MPIUtil::ReduceMax);
	mpiUtil->allReduce(uMin, MPIUtil::ReduceMin);
	mpiUtil->allReduce(uMax, MPIUtil::ReduceMax);
	if(!omegaMax) omegaMax = std::max(uMax-uMin, oMax-oMin);
	double Emin = uMin - omegaMax;
	double Emax = oMax + omegaMax;
	//--- print selected values after fixing defaults:
	logPrintf("Frequency resolution:    %lg\n", eta);
	logPrintf("Dielectric matrix Ecut:  %lg\n", Ecut);
	logPrintf("Maximum energy transfer: %lg\n", omegaMax);
	
	//Make necessary quantities available on all processes:
	std::vector<ColumnBundle> C(e.eInfo.nStates);
	std::vector<diagMatrix> E(e.eInfo.nStates);
	std::vector<diagMatrix> F(e.eInfo.nStates);
	for(int q=0; q<e.eInfo.nStates; q++)
	{	int procSrc = e.eInfo.whose(q);
		if(procSrc == mpiUtil->iProcess())
		{	std::swap(C[q], e.eVars.C[q]);
			std::swap(E[q], e.eVars.Hsub_eigs[q]);
			std::swap(F[q], e.eVars.F[q]);
		}
		else
		{	C[q].init(e.eInfo.nBands, e.basis[q].nbasis * e.eInfo.spinorLength(), &e.basis[q], &e.eInfo.qnums[q]);
			E[q].resize(e.eInfo.nBands);
			F[q].resize(e.eInfo.nBands);
		}
		C[q].bcast(procSrc);
		E[q].bcast(procSrc);
		F[q].bcast(procSrc);
	}
	
	//Report maximum nearest-neighbour eigenvalue change (to guide choice of eta)
	const Supercell& supercell = *(e.coulombParams.supercell);
	matrix3<> kBasisT = inv(supercell.Rsuper) * e.gInfo.R;
	vector3<> kBasis[3]; for(int j=0; j<3; j++) kBasis[j] = kBasisT.row(j);
	PeriodicLookup<vector3<>> plook(supercell.kmesh, e.gInfo.GGT);
	size_t ikStart, ikStop;
	TaskDivision(supercell.kmesh.size(), mpiUtil).myRange(ikStart, ikStop);
	double dEmax = 0.;
	for(size_t ik=ikStart; ik<ikStop; ik++)
	{	const diagMatrix& Ei = E[supercell.kmeshTransform[ik].iReduced];
		for(int j=0; j<3; j++)
		{	size_t jk = plook.find(supercell.kmesh[ik] + kBasis[j]);
			assert(jk != string::npos);
			const diagMatrix& Ej = E[supercell.kmeshTransform[jk].iReduced];
			for(int b=0; b<e.eInfo.nBands; b++)
				if(Emin <= Ei[b] && Ei[b] <= Emax)
					dEmax = std::max(dEmax, fabs(Ej[b]-Ei[b]));
		}
	}
	mpiUtil->allReduce(dEmax, MPIUtil::ReduceMax);
	logPrintf("Maximum k-neighbour dE: %lg (guide for selecting eta)\n", dEmax);
	
	//Initialize reduced q-Mesh:
	//--- q-mesh is a k-point dfference mesh, which could differ from k-mesh for off-Gamma meshes
	std::vector<QuantumNumber> qmesh(supercell.kmesh.size());
	for(size_t iq=0; iq<qmesh.size(); iq++)
	{	qmesh[iq].k = supercell.kmesh[iq] - supercell.kmesh[0]; //k-difference
		qmesh[iq].weight = 1./qmesh.size(); //uniform mesh
		qmesh[iq].spin = 0;
	}
	logPrintf("Symmetries reduced momentum transfers (q-mesh) from %d to ", int(qmesh.size()));
	qmesh = e.symm.reduceKmesh(qmesh);
	logPrintf("%d entries\n", int(qmesh.size())); logFlush();
	TaskDivision qDivision(qmesh.size(), mpiUtil);
	
	//Initialize polarizability/dielectric bases corresponding to qmesh:
	logPrintf("Setting up reduced polarizability bases at Ecut = %lg: ", Ecut); logFlush();
	std::vector<Basis> basisChi(qmesh.size());
	double avg_nbasis = 0.;
	const GridInfo& gInfoBasis = e.gInfoWfns ? *e.gInfoWfns : e.gInfo;
	logSuspend();
	for(size_t q=qDivision.start(); q<qDivision.stop(); q++)
	{	basisChi[q].setup(gInfoBasis, e.iInfo, Ecut, qmesh[q].k);
		avg_nbasis += qmesh[q].weight * basisChi[q].nbasis;
	}
	logResume();
	mpiUtil->allReduce(avg_nbasis, MPIUtil::ReduceSum);
	logPrintf("nbasis = %.2lf average, %.2lf ideal\n", avg_nbasis, pow(sqrt(2*Ecut),3)*(e.gInfo.detR/(6*M_PI*M_PI)));
	logFlush();

	//Initialize common wavefunction basis and ColumnBundle transforms for full k-mesh:
	logPrintf("Setting up k-mesh wavefunction transforms ... "); logFlush();
	double kMaxSq = 0;
	for(const vector3<>& k: supercell.kmesh)
		kMaxSq = std::max(kMaxSq, e.gInfo.GGT.metric_length_squared(k));
	double GmaxEff = sqrt(2.*e.cntrl.Ecut) + sqrt(kMaxSq);
	double EcutEff = 0.5*GmaxEff*GmaxEff * (1.+symmThreshold); //add some margin for round-off error safety
	logSuspend();
	Basis basis; basis.setup(e.gInfo, e.iInfo, EcutEff, vector3<>());
	logResume();
	ColumnBundleTransform::BasisWrapper basisWrapper(basis);
	std::vector<matrix3<int>> sym = e.symm.getMatrices();
	std::vector<std::shared_ptr<ColumnBundleTransform>> transform(supercell.kmesh.size());
	for(size_t ik=0; ik<supercell.kmesh.size(); ik++)
	{	const Supercell::KmeshTransform& kTransform = supercell.kmeshTransform[ik];
		const Basis& basisC = e.basis[kTransform.iReduced];
		const vector3<>& kC = e.eInfo.qnums[kTransform.iReduced].k;
		transform[ik] = std::make_shared<ColumnBundleTransform>(kC, basisC, supercell.kmesh[ik], basisWrapper,
			e.eInfo.spinorLength(), sym[kTransform.iSym], kTransform.invert);
	}
	logPrintf("done.\n"); logFlush();

	
	die("Not yet implemented.\n");
	logPrintf("\n"); logFlush();
}
