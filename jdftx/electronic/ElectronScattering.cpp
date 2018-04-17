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
#include <electronic/SpeciesInfo_internal.h>
#include <core/SphericalHarmonics.h>
#include <core/LatticeUtils.h>
#include <core/Random.h>

//Matrix imaginary part (in the operator sense; not elementwise):
matrix Im(const matrix& m)
{	return complex(0,-0.5)*(m - dagger(m)); //(m - dagger(m))/(2i)
}

ElectronScattering::ElectronScattering()
: eta(0.), Ecut(0.), fCut(1e-6), omegaMax(0.), slabResponse(false), EcutTransverse(0.)
{
}

void ElectronScattering::dump(const Everything& everything)
{	Everything& e = (Everything&)everything; //may modify everything to save memory / optimize
	this->e = &everything;
	nBands = e.eInfo.nBands;
	nSpinor = e.eInfo.spinorLength();
	nSpins = e.eInfo.nSpins();
	qCount = e.eInfo.nStates / nSpins; //reduce k-mesh size (without spin)
	
	if(slabResponse)
		logPrintf("\n----- Slab dielectric matrix calculation -----\n");
	else
		logPrintf("\n----- Electron-electron scattering Im(Sigma) -----\n");
	logFlush();
	
	//Update default parameters:
	if(!eta)
	{	eta = e.eInfo.smearingWidth;
		if(!eta) die("eta must be specified explicitly since electronic temperature is zero.\n");
	}
	if(!Ecut) Ecut = e.cntrl.Ecut;
	if(!EcutTransverse) EcutTransverse = e.cntrl.Ecut;
	double oMin = DBL_MAX, oMax = -DBL_MAX; //occupied energy range
	double uMin = DBL_MAX, uMax = -DBL_MAX; //unoccupied energy range
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		for(int b=0; b<nBands; b++)
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
	mpiWorld->allReduce(oMin, MPIUtil::ReduceMin);
	mpiWorld->allReduce(oMax, MPIUtil::ReduceMax);
	mpiWorld->allReduce(uMin, MPIUtil::ReduceMin);
	mpiWorld->allReduce(uMax, MPIUtil::ReduceMax);
	if(!omegaMax) omegaMax = std::max(uMax-uMin, oMax-oMin);
	Emin = uMin - omegaMax;
	Emax = oMax + omegaMax;
	//--- print selected values after fixing defaults:
	logPrintf("Frequency resolution:    %lg\n", eta);
	logPrintf("Dielectric matrix Ecut:  %lg\n", Ecut);
	logPrintf("Maximum energy transfer: %lg\n", omegaMax);
	
	//Initialize frequency grid:
	diagMatrix omegaGrid, wOmega;
	omegaGrid.push_back(0.);
	wOmega.push_back(0.5*eta); //integration weight (halved at endpoint)
	while(omegaGrid.back()<omegaMax + 10*eta) //add margin for covering enough of the Lorentzians
	{	omegaGrid.push_back(omegaGrid.back() + eta);
		wOmega.push_back(eta);
	}
	int iOmegaStart, iOmegaStop; //split dielectric computation over frequency grid
	TaskDivision omegaDiv(omegaGrid.size(), mpiWorld);
	omegaDiv.myRange(iOmegaStart, iOmegaStop);
	logPrintf("Initialized frequency grid with resolution %lg and %d points.\n", eta, omegaGrid.nRows());

	//Handle simpler slab response case, which only requires zero momentum transfer, separately:
	if(slabResponse)
	{	std::swap(C, e.eVars.C);
		std::swap(E, e.eVars.Hsub_eigs);
		std::swap(F, e.eVars.F);
		std::swap(VdagC, e.eVars.VdagC);
		dumpSlabResponse(e, omegaGrid);
		return;
	}
	
	//Make necessary quantities available on all processes:
	C.resize(e.eInfo.nStates);
	E.resize(e.eInfo.nStates);
	F.resize(e.eInfo.nStates);
	VdagC.resize(e.eInfo.nStates);
	for(int q=0; q<e.eInfo.nStates; q++)
	{	int procSrc = e.eInfo.whose(q);
		if(e.eInfo.isMine(q))
		{	std::swap(C[q], e.eVars.C[q]);
			std::swap(E[q], e.eVars.Hsub_eigs[q]);
			std::swap(F[q], e.eVars.F[q]);
			std::swap(VdagC[q], e.eVars.VdagC[q]);
		}
		else
		{	C[q].init(nBands, e.basis[q].nbasis * nSpinor, &e.basis[q], &e.eInfo.qnums[q]);
			E[q].resize(nBands);
			F[q].resize(nBands);
			VdagC[q].resize(e.iInfo.species.size());
		}
		C[q].bcast(procSrc);
		E[q].bcast(procSrc);
		F[q].bcast(procSrc);
		for(unsigned iSp=0; iSp<e.iInfo.species.size(); iSp++)
			if(e.iInfo.species[iSp]->isUltrasoft())
			{	if(!e.eInfo.isMine(q))
					VdagC[q][iSp].init(e.iInfo.species[iSp]->nProjectors(), nBands);
				VdagC[q][iSp].bcast(procSrc);
			}
	}
	
	//Randomize supercell to improve load balancing on k-mesh:
	{	std::vector< vector3<> >& kmesh = e.coulombParams.supercell->kmesh;
		std::vector<Supercell::KmeshTransform>& kmeshTransform = e.coulombParams.supercell->kmeshTransform;
		for(size_t ik=0; ik<kmesh.size()-1; ik++)
		{	size_t jk = ik + floor(Random::uniform(kmesh.size()-ik));
			mpiWorld->bcast(jk);
			if(jk !=ik && jk < kmesh.size())
			{	std::swap(kmesh[ik], kmesh[jk]);
				std::swap(kmeshTransform[ik], kmeshTransform[jk]);
			}
		}
	}
	
	//Report maximum nearest-neighbour eigenvalue change (to guide choice of eta)
	supercell = e.coulombParams.supercell;
	matrix3<> kBasisT = inv(supercell->Rsuper) * e.gInfo.R;
	vector3<> kBasis[3]; for(int j=0; j<3; j++) kBasis[j] = kBasisT.row(j);
	plook = std::make_shared< PeriodicLookup< vector3<> > >(supercell->kmesh, e.gInfo.GGT);
	size_t ikStart, ikStop;
	TaskDivision(supercell->kmesh.size(), mpiWorld).myRange(ikStart, ikStop);
	double dEmax = 0.;
	for(size_t ik=ikStart; ik<ikStop; ik++)
	for(int iSpin=0; iSpin<nSpins; iSpin++)
	{	const diagMatrix& Ei = E[supercell->kmeshTransform[ik].iReduced+iSpin*qCount];
		for(int j=0; j<3; j++)
		{	size_t jk = plook->find(supercell->kmesh[ik] + kBasis[j]);
			assert(jk != string::npos);
			const diagMatrix& Ej = E[supercell->kmeshTransform[jk].iReduced+iSpin*qCount];
			for(int b=0; b<nBands; b++)
				if(Emin <= Ei[b] && Ei[b] <= Emax)
					dEmax = std::max(dEmax, fabs(Ej[b]-Ei[b]));
		}
	}
	mpiWorld->allReduce(dEmax, MPIUtil::ReduceMax);
	logPrintf("Maximum k-neighbour dE: %lg (guide for selecting eta)\n", dEmax);
	
	//Initialize reduced q-Mesh:
	//--- q-mesh is a k-point dfference mesh, which could differ from k-mesh for off-Gamma meshes
	qmesh.resize(supercell->kmesh.size());
	for(size_t iq=0; iq<qmesh.size(); iq++)
	{	qmesh[iq].k = supercell->kmesh[iq] - supercell->kmesh[0]; //k-difference
		qmesh[iq].weight = 1./qmesh.size(); //uniform mesh
		qmesh[iq].spin = 0;
	}
	logPrintf("Symmetries reduced momentum transfers (q-mesh) from %d to ", int(qmesh.size()));
	qmesh = e.symm.reduceKmesh(qmesh);
	logPrintf("%d entries\n", int(qmesh.size())); logFlush();
	
	//Initialize polarizability/dielectric bases corresponding to qmesh:
	logPrintf("Setting up reduced polarizability bases at Ecut = %lg: ", Ecut); logFlush();
	basisChi.resize(qmesh.size());
	double avg_nbasis = 0.;
	const GridInfo& gInfoBasis = e.gInfoWfns ? *e.gInfoWfns : e.gInfo;
	logSuspend();
	for(size_t iq=0; iq<qmesh.size(); iq++)
	{	basisChi[iq].setup(gInfoBasis, e.iInfo, Ecut, qmesh[iq].k);
		avg_nbasis += qmesh[iq].weight * basisChi[iq].nbasis;
	}
	logResume();
	logPrintf("nbasis = %.2lf average, %.2lf ideal\n", avg_nbasis, pow(sqrt(2*Ecut),3)*(e.gInfo.detR/(6*M_PI*M_PI)));
	logFlush();

	//Initialize common wavefunction basis and ColumnBundle transforms for full k-mesh:
	logPrintf("Setting up k-mesh wavefunction transforms ... "); logFlush();
	double kMaxSq = 0;
	for(const vector3<>& k: supercell->kmesh)
	{	kMaxSq = std::max(kMaxSq, e.gInfo.GGT.metric_length_squared(k));
		for(const QuantumNumber& qnum: qmesh)
			kMaxSq = std::max(kMaxSq, e.gInfo.GGT.metric_length_squared(k + qnum.k));
	}
	double kWeight = double(e.eInfo.spinWeight) / supercell->kmesh.size();
	double GmaxEff = sqrt(2.*e.cntrl.Ecut) + sqrt(kMaxSq);
	double EcutEff = 0.5*GmaxEff*GmaxEff * (1.+symmThreshold); //add some margin for round-off error safety
	logSuspend();
	basis.setup(gInfoBasis, e.iInfo, EcutEff, vector3<>());
	logResume();
	ColumnBundleTransform::BasisWrapper basisWrapper(basis);
	std::vector<SpaceGroupOp> sym = e.symm.getMatrices();
	for(size_t ik=ikStart; ik<ikStop; ik++)
	{	const vector3<>& k = supercell->kmesh[ik];
		for(const QuantumNumber& qnum: qmesh)
		{	vector3<> k2 = k + qnum.k; double roundErr;
			vector3<int> k2sup = round((k2 - supercell->kmesh[0]) * supercell->super, &roundErr);
			assert(roundErr < symmThreshold);
			auto iter = transform.find(k2sup);
			if(iter == transform.end())
			{	size_t ik2 = plook->find(k2); assert(ik2 != string::npos);
				const Supercell::KmeshTransform& kTransform = supercell->kmeshTransform[ik2];
				const Basis& basisC = e.basis[kTransform.iReduced];
				const vector3<>& kC = e.eInfo.qnums[kTransform.iReduced].k;
				transform[k2sup] = std::make_shared<ColumnBundleTransform>(kC, basisC, k2, basisWrapper,
					nSpinor, sym[kTransform.iSym], kTransform.invert);
				//Initialize corresponding quantum number:
				QuantumNumber qnum;
				qnum.k = k2;
				qnum.spin = 0;
				qnum.weight = kWeight;
				qnumMesh[k2sup] = qnum;
			}
		}
	}
	logPrintf("done.\n"); logFlush();

	//Main loop over momentum transfers:
	diagMatrix ImKscrHead(omegaGrid.size(), 0.);
	std::vector<diagMatrix> ImSigma(e.eInfo.nStates, diagMatrix(nBands,0.));
	for(size_t iq=0; iq<qmesh.size(); iq++)
	{	logPrintf("\nMomentum transfer %d of %d: q = ", int(iq+1), int(qmesh.size()));
		qmesh[iq].k.print(globalLog, " %+.5lf ");
		int nbasis = basisChi[iq].nbasis;
		nAugRhoAtomInit(iq);
		
		//Construct Coulomb operator (regularizes G=0 using the tricks developed for EXX):
		matrix invKq = inv(coulombMatrix(iq));
		
		//Calculate chi_KS:
		std::vector<matrix> chiKS(omegaGrid.nRows());
		logPrintf("\tComputing chi_KS ...  "); logFlush(); 
		size_t nkMine = ikStop-ikStart;
		int ikInterval = std::max(1, int(round(nkMine/20.))); //interval for reporting progress
		for(size_t ik=ikStart; ik<ikStop; ik++)
		for(int iSpin=0; iSpin<nSpins; iSpin++)
		{	//Report progress:
			size_t ikDone = ik-ikStart+1;
			if(ikDone % ikInterval == 0 && !iSpin)
			{	logPrintf("%d%% ", int(round(ikDone*100./nkMine)));
				logFlush();
			}
			//Get events:
			size_t jk; matrix nij;
			std::vector<Event> events = getEvents(true, iSpin, ik, iq, jk, nij);
			if(!events.size()) continue;
			//Collect contributions for each frequency:
			for(int iOmega=0; iOmega<omegaGrid.nRows(); iOmega++)
			{	double omega = omegaGrid[iOmega];
				complex omegaTilde(omega, 2*eta);
				complex one(1,0);
				std::vector<complex> Xks; Xks.reserve(events.size());
				for(const Event& event: events)
					Xks.push_back(-e.gInfo.detR * kWeight * event.fWeight
						* (one/(event.Eji - omegaTilde) + one/(event.Eji + omegaTilde)) );
				chiKS[iOmega] += (nij * Xks) * dagger(nij);
			}
		}
		for(int iOmega=0; iOmega<omegaGrid.nRows(); iOmega++)
		{	chiKS[iOmega].allReduce(MPIUtil::ReduceSum);
			if(!omegaDiv.isMine(iOmega)) chiKS[iOmega] = 0; //no longer needed on this process
		}
		logPrintf("done.\n"); logFlush();
		
		//Figure out head entry index:
		int iHead = 0;
		for(const vector3<int>& iG: basisChi[iq].iGarr)
		{	if(!iG.length_squared()) break;
			iHead++;
		}
		assert(iHead < nbasis);
		
		//Calculate Im(screened Coulomb operator):
		logPrintf("\tComputing Im(Kscreened) ... "); logFlush();
		std::vector<matrix> ImKscr(omegaGrid.nRows());
		for(int iOmega=iOmegaStart; iOmega<iOmegaStop; iOmega++)
		{	ImKscr[iOmega] = Im(inv(invKq - chiKS[iOmega]));
			chiKS[iOmega] = 0; //free to save memory
			ImKscrHead[iOmega] += qmesh[iq].weight * ImKscr[iOmega](iHead,iHead).real(); //accumulate head of ImKscr
		}
		chiKS.clear(); //free memory; no longer needed
		for(int iOmega=0; iOmega<omegaGrid.nRows(); iOmega++)
		{	if(!omegaDiv.isMine(iOmega)) ImKscr[iOmega] = zeroes(nbasis,nbasis);
			ImKscr[iOmega].bcast(omegaDiv.whose(iOmega));
		}
		logPrintf("done.\n"); logFlush();
		
		//Calculate ImSigma contributions:
		logPrintf("\tComputing ImSigma ... "); logFlush(); 
		for(size_t ik=ikStart; ik<ikStop; ik++)
		for(int iSpin=0; iSpin<nSpins; iSpin++)
		{	//Report progress:
			size_t ikDone = ik-ikStart+1;
			if(ikDone % ikInterval == 0 && !iSpin)
			{	logPrintf("%d%% ", int(round(ikDone*100./nkMine)));
				logFlush();
			}
			//Get events:
			size_t jk; matrix nij;
			std::vector<Event> events = getEvents(false, iSpin, ik, iq, jk, nij);
			if(!events.size()) continue;
			//Integrate over frequency for event contributions to linewidth:
			diagMatrix eventContrib(events.size(), 0);
			for(int iOmega=0; iOmega<omegaGrid.nRows(); iOmega++)
			{	//Construct energy conserving delta-function:
				double omega = omegaGrid[iOmega];
				complex omegaTilde(omega, 2*eta);
				diagMatrix delta; delta.reserve(events.size());
				for(const Event& event: events)
					delta.push_back(e.gInfo.detR * event.fWeight //overlap and sign for electron / hole
						* (2*eta/M_PI) * ( 1./(event.Eji - omegaTilde).norm() - 1./(event.Eji + omegaTilde).norm()) ); //Normalized Lorentzians
				eventContrib += wOmega[iOmega] * delta * diag(dagger(nij) * ImKscr[iOmega] * nij);
			}
			//Accumulate contributions to linewidth:
			int iReduced = supercell->kmeshTransform[ik].iReduced; //directly collect to reduced k-point
			double symFactor = e.eInfo.spinWeight / (supercell->kmesh.size() * e.eInfo.qnums[iReduced].weight); //symmetrization factor = 1 / |orbit of iReduced|
			double qWeight = qmesh[iq].weight;
			for(size_t iEvent=0; iEvent<events.size(); iEvent++)
			{	const Event& event = events[iEvent];
				ImSigma[iReduced+iSpin*qCount][event.i] += symFactor * qWeight * eventContrib[iEvent];
			}
		}
		logPrintf("done.\n"); logFlush();
	}
	logPrintf("\n");
	
	ImKscrHead.allReduce(MPIUtil::ReduceSum);
	for(diagMatrix& IS: ImSigma)
		IS.allReduce(MPIUtil::ReduceSum);
	for(int q=0; q<e.eInfo.nStates; q++)
		for(int b=0; b<nBands; b++)
		{	double Eqb = E[q][b];
			if(Eqb<Emin || Eqb>Emax)
				ImSigma[q][b] = NAN; //clearly mark as invalid
		}
	
	string fname = e.dump.getFilename("ImSigma_ee");
	logPrintf("Dumping %s ... ", fname.c_str()); logFlush();
	e.eInfo.write(ImSigma, fname.c_str());
	logPrintf("done.\n");

	fname = e.dump.getFilename("ImKscrHead");
	logPrintf("Dumping %s ... ", fname.c_str()); logFlush();
	if(mpiWorld->isHead())
	{	FILE* fp = fopen(fname.c_str(), "w");
		for(int iOmega=0; iOmega<omegaGrid.nRows(); iOmega++)
			fprintf(fp, "%lf %le\n", omegaGrid[iOmega], ImKscrHead[iOmega]);
		fclose(fp);
	}
	logPrintf("done.\n");

	logPrintf("\n"); logFlush();
}

//Calculate diag(A*dagger(B)) without constructing large intermediate matrix
diagMatrix diagouter(const matrix& A, const matrix& B)
{	assert(A.nRows()==B.nRows());
	assert(A.nCols()==B.nCols());
	//Elementwise multiply A and conj(B);
	matrix ABconj = conj(B);
	callPref(eblas_zmul)(ABconj.nData(), A.dataPref(),1, ABconj.dataPref(),1);
	//Add columns of ABconj:
	for(int col=1; col<ABconj.nCols(); col++)
		callPref(eblas_zaxpy)(ABconj.nRows(), 1., ABconj.dataPref()+ABconj.index(0,col),1, ABconj.dataPref(),1);
	//Return real part as diagMatrix:
	diagMatrix result(ABconj.nRows(), 0.);
	eblas_daxpy(result.nRows(), 1., (double*)ABconj.data(),2, result.data(),1);
	return result;
}

std::vector<ElectronScattering::Event> ElectronScattering::getEvents(bool chiMode, int iSpin, size_t ik, size_t iq, size_t& jk, matrix& nij) const
{	static StopWatch watchI("ElectronScattering::getEventsI"), watchJ("ElectronScattering::getEventsJ"), watchAug("ElectronScattering::nAug");
	//Find target k-point:
	const vector3<>& ki = slabResponse ? e->eInfo.qnums[ik].k : supercell->kmesh[ik];
	const vector3<> kj = ki + qmesh[iq].k;
	if(slabResponse)
	{	assert(iq == 0);
		jk = ik;
	}
	else	
	{	jk = plook->find(kj);
		assert(jk != string::npos);
	}
	
	//Compile list of events:
	int iReduced = slabResponse ? ik: supercell->kmeshTransform[ik].iReduced+iSpin*qCount;
	int jReduced = slabResponse ? jk: supercell->kmeshTransform[jk].iReduced+iSpin*qCount;
	const diagMatrix &Ei = E[iReduced], &Fi = F[iReduced];
	const diagMatrix &Ej = E[jReduced], &Fj = F[jReduced];
	std::vector<Event> events; events.reserve((nBands*nBands)/2);
	std::vector<bool> iUsed(nBands,false), jUsed(nBands,false); //sets of i and j actually referenced
	Event event;
	for(event.i=0; event.i<nBands; event.i++)
	for(event.j=0; event.j<nBands; event.j++)
	{	event.fWeight = chiMode ? 0.5*(Fi[event.i] - Fj[event.j]) : (1. - Fi[event.i] - Fj[event.j]);
		double Eii = Ei[event.i];
		double Ejj = Ej[event.j];
		event.Eji = Ejj - Eii;
		if(!chiMode)
		{	if(Eii<Emin || Eii>Emax) event.fWeight = 0.; //state out of relevant range
			if(event.fWeight * (Eii-Ejj) <= 0) event.fWeight = 0; //wrong sign for energy transfer
		}
		bool needEvent = (fabs(event.fWeight) > fCut);
		if(needEvent)
		{	events.push_back(event);
			iUsed[event.i] = true;
			jUsed[event.j] = true;
		}
	}
	if(!events.size()) return events;
	
	//Get wavefunctions in real space:
	ColumnBundle Ci, Cj;
	std::vector<matrix> VdagCi, VdagCj;
	if(!slabResponse)
	{	Ci = getWfns(ik, iSpin, ki, &VdagCi);
		Cj = getWfns(jk, iSpin, kj, &VdagCj);
	}
	else
	{	VdagCi = VdagC[ik];
		VdagCj = VdagC[jk];
	}
	std::vector< std::vector<complexScalarField> > conjICi(nBands), ICj(nBands);
	watchI.start();
	for(int i=0; i<nBands; i++) if(iUsed[i])
	{	conjICi[i].resize(nSpinor);
		for(int s=0; s<nSpinor; s++)
			conjICi[i][s] = conj(I((slabResponse ? C[ik] : Ci).getColumn(i,s))); 
	}
	for(int j=0; j<nBands; j++) if(jUsed[j])
	{	ICj[j].resize(nSpinor);
		for(int s=0; s<nSpinor; s++)
			ICj[j][s] = I((slabResponse ? C[jk] : Cj).getColumn(j,s));
	}
	watchI.stop();
	
	//Initialize pair densities:
	watchJ.start();
	const Basis& basis_q = basisChi[iq];
	int nbasis = basis_q.nbasis;
	nij = zeroes(nbasis, events.size());
	complex* nijData = nij.dataPref();
	for(const Event& event: events)
	{	complexScalarField Inij;
		for(int s=0; s<nSpinor; s++)
			Inij += conjICi[event.i][s] * ICj[event.j][s];
		callPref(eblas_gather_zdaxpy)(nbasis, 1., basis_q.index.dataPref(), J(Inij)->dataPref(), nijData);
		nijData += nbasis;
	}
	watchJ.stop();
	
	//Augmentation:
	watchAug.start();
	for(unsigned iSp=0; iSp<e->iInfo.species.size(); iSp++)
	{	const SpeciesInfo& sp = *(e->iInfo.species[iSp]);
		if(sp.isUltrasoft())
		{	int nProjSp = sp.MnlAll.nRows(); //number of projectors per atom including spinor indices
			int nProj = nProjSp / nSpinor;
			//Collect atomic density matrices for all events:
			matrix RhoAll(nProj*nProj*sp.atpos.size(), events.size());
			complex* RhoAllData = RhoAll.dataPref();
			for(const Event& event: events)
			{	for(unsigned atom=0; atom<sp.atpos.size(); atom++)
				{	//Atomic density matrix:
					matrix atomVdagCi = VdagCi[iSp](atom*nProjSp,(atom+1)*nProjSp, event.i,event.i+1);
					matrix atomVdagCj = VdagCj[iSp](atom*nProjSp,(atom+1)*nProjSp, event.j,event.j+1);
					matrix Rho = atomVdagCj * dagger(atomVdagCi); //density matrix in projector basis on this atom (Includes dagger(Vj) at inner index j, Vi at outer index i)
					if(sp.isRelativistic()) Rho = sp.fljAll * Rho * sp.fljAll; //transformation for relativistic pseudopotential
					if(nSpinor==2) Rho = Rho(0,2,nProjSp, 0,2,nProjSp) + Rho(1,2,nProjSp, 1,2,nProjSp); //contract spinorial index
					callPref(eblas_copy)(RhoAllData, Rho.dataPref(), Rho.nData());
					RhoAllData += Rho.nData();
				}
			}
			//Augment:
			nij += nAugRhoAtom[iSp] * RhoAll;
		}
	}
	watchAug.stop();
	
	return events;
}

ColumnBundle ElectronScattering::getWfns(size_t ik, int iSpin, const vector3<>& k, std::vector<matrix>* VdagCi) const
{	static StopWatch watch("ElectronScattering::getWfns"); watch.start();
	double roundErr;
	vector3<int> kSup = round((k - supercell->kmesh[0]) * supercell->super, &roundErr);
	assert(roundErr < symmThreshold);
	ColumnBundle result(nBands, basis.nbasis * nSpinor, &basis, &qnumMesh.find(kSup)->second, isGpuEnabled());
	result.zero();
	const ColumnBundleTransform& cbt = *(transform.find(kSup)->second);
	const Supercell::KmeshTransform& kmt = supercell->kmeshTransform[ik];
	cbt.scatterAxpy(1., C[kmt.iReduced+iSpin*qCount], result,0,1);
	if(VdagCi) *VdagCi = cbt.transformVdagC(VdagC[kmt.iReduced+iSpin*qCount], kmt.iSym);
	watch.stop();
	return result;
}

matrix ElectronScattering::coulombMatrix(size_t iq) const
{	//Use function implemented in Polarizability:
	matrix coulombMatrix(const ColumnBundle& V, const Everything& e, vector3<> dk);
	const Basis& basis_q = basisChi[iq];
	ColumnBundle V(basis_q.nbasis, basis_q.nbasis, &basis_q);
	V.zero();
	complex* Vdata = V.data();
	double normFac = 1./sqrt(e->gInfo.detR);
	for(size_t b=0; b<basis_q.nbasis; b++)
		Vdata[V.index(b,b)] = normFac;
	return coulombMatrix(V, *e, qmesh[iq].k);
}

void ElectronScattering::nAugRhoAtomInit(size_t iq)
{	static StopWatch watch("ElectronScattering::nAugInit");
	nAugRhoAtom.clear();
	nAugRhoAtom.resize(e->iInfo.species.size());
	const Basis& basis_q = basisChi[iq];
	//For each ultrasoft species:
	for(unsigned iSp=0; iSp<e->iInfo.species.size(); iSp++)
	{	const SpeciesInfo& sp = *(e->iInfo.species[iSp]);
		if(sp.isUltrasoft())
		{	watch.start();
			//Count  number of Qij,lm entries:
			int nQijlm = 0;
			for(const auto& entry: sp.Qradial)
				nQijlm += (2*entry.first.l+1);
			//Initialize matrix with Qij,lm functions in basis:
			matrix& n = nAugRhoAtom[iSp];
			n = zeroes(basis_q.nbasis, nQijlm*sp.atpos.size());
			std::map<std::pair<SpeciesInfo::QijIndex,int>, int> ijlmMap;
			complex* nData = n.dataPref();
			int atomStride = basis_q.nbasis * nQijlm;
			int ijlmIndex = 0;
			for(const auto& entry: sp.Qradial)
			{	const int& l = entry.first.l;
				for(int m=-l; m<=l; m++)
				{	callPref(Vnl)(basis_q.nbasis, atomStride, sp.atpos.size(), l, m, qmesh[iq].k, basis_q.iGarr.dataPref(),
						e->gInfo.G, sp.atposManaged.dataPref(), entry.second, nData);
					ijlmMap[std::make_pair(entry.first,m)] = ijlmIndex;
					nData += basis_q.nbasis;
					ijlmIndex++;
				}
			}
			assert(ijlmIndex == nQijlm);
			//Convert to augmentation functions per density-matrix entry:
			int nProj = sp.MnlAll.nRows()/nSpinor;
			const double invDetR = 1./e->gInfo.detR;
			matrix ijlmTOVV = zeroes(nQijlm, nProj*nProj);
			//--- Triple loop over first projector:
			int i1 = 0;
			for(int l1=0; l1<int(sp.VnlRadial.size()); l1++)
			for(int p1=0; p1<int(sp.VnlRadial[l1].size()); p1++)
			for(int m1=-l1; m1<=l1; m1++)
			{	//--- Triple loop over second projector:
				int i2 = 0;
				for(int l2=0; l2<int(sp.VnlRadial.size()); l2++)
				for(int p2=0; p2<int(sp.VnlRadial[l2].size()); p2++)
				for(int m2=-l2; m2<=l2; m2++)
				{	std::vector<YlmProdTerm> terms = expandYlmProd(l1,m1, l2,m2);
					for(const YlmProdTerm& term: terms)
					{	//Find index into ijlm function for current term (if any):
						SpeciesInfo::QijIndex qIndex = { l1, p1, l2, p2, term.l };
						auto ijlmIter = ijlmMap.find(std::make_pair(qIndex,term.m));
						if(ijlmIter == ijlmMap.end()) continue;  //no entry at this l
						const int& ijlmIndex = ijlmIter->second;
						ijlmTOVV.set(ijlmIndex, i1*nProj+i2, (term.coeff*invDetR)*cis(0.5*M_PI*(l2-l1-term.l))); //phase=(-i)^(l+l1-l2)
					}
					i2++;
				}
				assert(i2==nProj);
				i1++;
			}
			assert(i1==nProj);
			n = n * tiledBlockMatrix(ijlmTOVV, sp.atpos.size());
			watch.stop();
		}
	}
}


//-----  Slab response related; eventually merge these with Polarizability instead -----

//Setup an ellipsoidal basis with different Ecuts along iDir and remaining tranverse directions.
//The plane waves with no transverse component are listed first, and their count is returned.
int setupEllipsoidalBasis(Basis& basis, const GridInfo& gInfo, const IonInfo& iInfo, double Ecut, double EcutTransverse, int iDir, const vector3<> k)
{
	//Bounding box:
	vector3<int> iGbox;
	for(int i=0; i<3; i++)
	{	const double& Ecut_i = i==iDir ? Ecut : EcutTransverse;
		iGbox[i] = 1 + int(sqrt(2*Ecut_i) * gInfo.R.column(i).length() / (2*M_PI));
	}
	//Add the indices with normal components alone first:
	std::vector<int> indexVec;
	vector3<int> iG;
	for(iG[iDir]=-iGbox[iDir]; iG[iDir]<=iGbox[iDir]; iG[iDir]++)
		if(0.5*gInfo.GGT.metric_length_squared(iG+k) <= Ecut)
			indexVec.push_back(gInfo.fullGindex(iG));
	int nSlab = indexVec.size();
	//Add remaining:
	double iDirScale = sqrt(EcutTransverse/Ecut); //scale factor on iDir to make ellipsoid spherical
	for(iG[0]=-iGbox[0]; iG[0]<=iGbox[0]; iG[0]++)
	for(iG[1]=-iGbox[1]; iG[1]<=iGbox[1]; iG[1]++)
	for(iG[2]=-iGbox[2]; iG[2]<=iGbox[2]; iG[2]++)
	{	if(iG.length_squared()==iG[iDir]*iG[iDir]) continue; //normal-component-only already included above
		vector3<> iGkScaled = iG+k;
		iGkScaled[iDir] *= iDirScale;
		if(0.5*gInfo.GGT.metric_length_squared(iGkScaled) <= EcutTransverse)
			indexVec.push_back(gInfo.fullGindex(iG));
	}
	basis.setup(gInfo, iInfo, indexVec);
	return nSlab; //number of normal-component only bases
}


void ElectronScattering::dumpSlabResponse(Everything& e, const diagMatrix& omegaGrid)
{
	//Zero momentum transfer quantum number:
	QuantumNumber qnum;
	qnum.k = vector3<>(); //only zero momentum transfer for slab response
	qnum.weight = 1.;
	qnum.spin = 0;
	qmesh.assign(1, qnum);
	
	//Basis:
	const GridInfo& gInfoBasis = e.gInfoWfns ? *e.gInfoWfns : e.gInfo;
	logPrintf("Setting up reduced polarizability bases at Ecut = %lg and EcutTransverse = %lg: ", Ecut, EcutTransverse); logFlush();
	assert(e.coulombParams.geometry == CoulombParams::Slab);
	basisChi.resize(qmesh.size());
	int nBasisSlab = setupEllipsoidalBasis(basisChi[0], gInfoBasis, e.iInfo, Ecut, EcutTransverse, e.coulombParams.iDir, qmesh[0].k);
	logPrintf("nbasis = %lu, %.2lf ideal\n", basisChi[0].nbasis, sqrt(2*Ecut)*(2*EcutTransverse)*(e.gInfo.detR/(6*M_PI*M_PI)));
	nAugRhoAtomInit(0);
	
	//Construct Coulomb operator (regularizes G=0 using the tricks developed for EXX):
	matrix invKq = inv(coulombMatrix(0));
	
	//Calculate chi_KS:
	std::vector<matrix> chiKS(omegaGrid.nRows());
	logPrintf("\tComputing chi_KS ...  "); logFlush(); 
	int nqMine = e.eInfo.qStop - e.eInfo.qStart;
	int iqInterval = std::max(1, int(round(nqMine/20.))); //interval for reporting progress
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{	//Report progress:
		size_t iqDone = q-e.eInfo.qStart+1;
		if(iqDone % iqInterval == 0)
		{	logPrintf("%d%% ", int(round(iqDone*100./nqMine)));
			logFlush();
		}
		//Get events:
		size_t jk; matrix nij;
		std::vector<Event> events = getEvents(true, 0, q, 0, jk, nij);
		if(!events.size()) continue;
		//Collect contributions for each frequency:
		for(int iOmega=0; iOmega<omegaGrid.nRows(); iOmega++)
		{	double omega = omegaGrid[iOmega];
			complex omegaTilde(omega, 2*eta);
			complex one(1,0);
			std::vector<complex> Xks; Xks.reserve(events.size());
			for(const Event& event: events)
				Xks.push_back(-e.gInfo.detR * e.eInfo.qnums[q].weight * event.fWeight
					* (one/(event.Eji - omegaTilde) + one/(event.Eji + omegaTilde)) );
			chiKS[iOmega] += (nij * Xks) * dagger(nij);
		}
	}
	int iOmegaStart, iOmegaStop; //split remaining computation over frequency grid
	TaskDivision omegaDiv(omegaGrid.size(), mpiWorld);
	omegaDiv.myRange(iOmegaStart, iOmegaStop);
	for(int iOmega=0; iOmega<omegaGrid.nRows(); iOmega++)
	{	chiKS[iOmega].allReduce(MPIUtil::ReduceSum);
		if(!omegaDiv.isMine(iOmega)) chiKS[iOmega] = 0; //no longer needed on this process
	}
	logPrintf("done.\n"); logFlush();

	//Output result:
	string fname = e.dump.getFilename("slabResponse");
	logPrintf("Dumping %s ... ", fname.c_str()); logFlush();
	MPIUtil::File fp;
	mpiWorld->fopenWrite(fp, fname.c_str());
	mpiWorld->fseek(fp, iOmegaStart*nBasisSlab*nBasisSlab*sizeof(complex), SEEK_SET);
	for(int iOmega=iOmegaStart; iOmega<iOmegaStop; iOmega++)
	{	matrix chiExt = (invKq*inv(invKq - chiKS[iOmega])*invKq - invKq)(0,nBasisSlab, 0,nBasisSlab); //external field susceptibility
		mpiWorld->fwrite(chiExt.data(), sizeof(complex), chiExt.nData(), fp);
	}
	mpiWorld->fclose(fp); logPrintf("done.\n");
	
	if(mpiWorld->isHead())
	{
		//Output frequency list:
		string fname = e.dump.getFilename("slabResponseOmega");
		logPrintf("Dumping %s ... ", fname.c_str()); logFlush();
		FILE* fp = fopen(fname.c_str(), "w");
		if(!fp) die("Failed to open '%s' for writing.\n", fname.c_str());
		omegaGrid.print(fp);
		fclose(fp); logPrintf("done.\n");
		
		//Output basis:
		fname = e.dump.getFilename("slabResponseBasis");
		logPrintf("Dumping %s ... ", fname.c_str()); logFlush();
		fp = fopen(fname.c_str(), "w");
		if(!fp) die("Failed to open '%s' for writing.\n", fname.c_str());
		for(const vector3<int>& iGn: basisChi[0].iGarr)
			fprintf(fp, "%d %d %d\n", iGn[0], iGn[1], iGn[2]);
		fclose(fp); logPrintf("done.\n");
	}
	logPrintf("\n");
}
