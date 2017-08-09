/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver, Deniz Gunceler
Copyright 1996-2003 Sohrab Ismail-Beigi

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

#include <electronic/Everything.h>
#include <electronic/ElecMinimizer.h>
#include <electronic/ColumnBundle.h>
#include <electronic/ExCorr.h>
#include <electronic/ExactExchange.h>
#include <fluid/FluidSolver.h>
#include <core/matrix.h>
#include <core/Units.h>
#include <core/ScalarFieldIO.h>
#include <cstdio>
#include <cmath>
#include <limits.h>

ElecVars::ElecVars()
: isRandom(true), initLCAO(true), skipWfnsInit(false), HauxInitialized(false), lcaoIter(-1), lcaoTol(1e-6)
{
}

// Kernel for generating box shaped external potentials
void applyBoxPot(int i, vector3<> r, matrix3<>& R, const ElecVars::BoxPotential* bP, double* Vbox)
{	// Map lattice intervals [0,1) to [-0.5,0.5)
	r = inv(R)*r;
	for(int j=0; j<3; j++) r[j] -= floor(0.5+r[j]);
	r = R*r;
	//Set potential
	for(int j=0; j<3; j++)
	{	if((r[j] < bP->min[j]) or (r[j] > bP->max[j]))
		{	Vbox[i] = bP->Vout;
			return;
		}
	}
	Vbox[i] = bP->Vin;
}

void ElecVars::setup(const Everything &everything)
{	
	this->e = &everything;
	logPrintf("\n---------- Allocating electronic variables ----------\n"); logFlush();

	const ElecInfo &eInfo = e->eInfo;
	const IonInfo& iInfo = e->iInfo;
	const std::vector<Basis>& basis = e->basis;
	const GridInfo& gInfo = e->gInfo;

	n.resize(eInfo.nDensities);
	if(e->ionDynamicsParams.tMax)
	  nAccumulated.resize(eInfo.nDensities);
	Vscloc.resize(n.size());
	if(eInfo.hasU)
	{	iInfo.rhoAtom_initZero(rhoAtom);
		iInfo.rhoAtom_initZero(U_rhoAtom);
	}
	
	if(VexternalFilename.size())
	{	Vexternal.resize(VexternalFilename.size());
		for(unsigned s=0; s<Vexternal.size(); s++)
		{	Vexternal[s] = ScalarFieldData::alloc(gInfo);
			logPrintf("Reading external potential from '%s'\n", VexternalFilename[s].c_str());
			loadRawBinary(Vexternal[s], VexternalFilename[s].c_str());
		}
		if(Vexternal.size()==1 && n.size()==2) //Replicate potential for second spin:
			Vexternal.push_back(Vexternal[0]->clone());
	}

	//Vexternal contributions from boxPot's.
	for(size_t j = 0; j<boxPot.size(); j++)
	{	//Create potential
		ScalarField temp; nullToZero(temp, e->gInfo);
		applyFunc_r(e->gInfo, applyBoxPot, gInfo.R, &boxPot[j], temp->data());
		temp = I(gaussConvolve(J(temp), boxPot[j].convolve_radius));
		//Add to Vexternal
		if(!Vexternal.size()) Vexternal.resize(n.size());
		for(unsigned s=0; s<n.size(); s++) Vexternal[s] += temp;
	}

	//Vexternal contributions due to external field
	if(e->coulombParams.Efield.length_squared())
	{	ScalarField temp = e->coulomb->getEfieldPotential();
		//Add to Vexternal
		if(!Vexternal.size()) Vexternal.resize(n.size());
		for(unsigned s=0; s<n.size(); s++) Vexternal[s] += temp;
	}

	if(rhoExternalFilename.length())
	{	logPrintf("Reading external charge from '%s'\n", rhoExternalFilename.c_str());
		ScalarField temp(ScalarFieldData::alloc(gInfo));
		loadRawBinary(temp, rhoExternalFilename.c_str());
		rhoExternal = J(temp);
	}
	
	
	//Initialize matrix arrays if required:
	Hsub.resize(eInfo.nStates);
	Hsub_evecs.resize(eInfo.nStates);
	Hsub_eigs.resize(eInfo.nStates);
	if(eInfo.fillingsUpdate==ElecInfo::FillingsHsub)
		Haux_eigs.resize(eInfo.nStates);
	if(eigsFilename.length())
	{	eInfo.read(Hsub_eigs, eigsFilename.c_str());
		if(eInfo.fillingsUpdate==ElecInfo::FillingsHsub)
		{	Haux_eigs = Hsub_eigs;
			HauxInitialized = true;
		}
	}
	VdagC.resize(eInfo.nStates, std::vector<matrix>(e->iInfo.species.size()));
	
	//Read in electron (spin) density if needed
	if(e->cntrl.fixed_H)
	{	string fnamePattern = nFilenamePattern.length() ? nFilenamePattern : VFilenamePattern; //Command ensures that the pattern has a "$VAR" in it
		#define READchannel(var, suffix) \
		{	string fname = fnamePattern; \
			size_t pos = fname.find("$VAR"); \
			assert(pos != string::npos); \
			fname.replace(pos,4, suffix); \
			logPrintf("Reading " #suffix " from file '%s' ... ", fname.c_str()); logFlush(); \
			nullToZero(var, e->gInfo); \
			loadRawBinary(var, fname.c_str()); \
			logPrintf("done\n"); logFlush(); \
		}
		#define READ(var) \
		{	var.resize(eInfo.nDensities); \
			if(var.size()==1) { READchannel(var[0], #var) } \
			else \
			{	READchannel(var[0], #var "_up") \
				READchannel(var[1], #var "_dn") \
				if(var.size()==4) \
				{	READchannel(var[2], #var "_re") \
					READchannel(var[3], #var "_im") \
				} \
			} \
		}
		#define READrhoAtom(var) \
		{	string fname = fnamePattern; \
			size_t pos = fname.find("$VAR"); \
			assert(pos != string::npos); \
			fname.replace(pos,4, #var); \
			logPrintf("Reading " #var " from file '%s' ... ", fname.c_str()); logFlush(); \
			FILE* fp = fopen(fname.c_str(), "r"); \
			for(matrix& m: var) m.read(fp); \
			fclose(fp); \
			logPrintf("done\n"); logFlush(); \
		}
		if(nFilenamePattern.length())
		{	READ(n)
			if(e->exCorr.needsKEdensity()) READ(tau)
			if(eInfo.hasU) READrhoAtom(rhoAtom)
		}
		else
		{	READ(Vscloc)
			if(e->exCorr.needsKEdensity()) READ(Vtau)
			if(eInfo.hasU) READrhoAtom(U_rhoAtom)
		}
		#undef READ
		#undef READchannel
		//Handle hybrid functionals
		if(e->exCorr.exxFactor())
			die("Band-structure (fixed electron density) calculations temporarily\n"
				"unsupported with hybrid functionals: use v1.1.2 or earlier.\n");
	}
	
	//Wavefunction initialiation (bypass in dry runs and phonon supercell calculations)
	if(skipWfnsInit)
	{	C.resize(eInfo.nStates); //skip memory allocation, but initialize array
		logPrintf("Skipped wave function initialization.\n");
	}
	else
	{
		// Initialize ColumnBundle arrays for the electronic wave-functions:
		logPrintf("Initializing wave functions:  ");
		init(C, eInfo.nStates, eInfo.nBands, &basis[0], &eInfo);

		//Initial wave functions
		int nBandsInited = 0;
		if(wfnsFilename.length())
		{	logPrintf("reading from '%s'\n", wfnsFilename.c_str()); logFlush();
			if(readConversion) readConversion->Ecut = e->cntrl.Ecut;
			read(C, wfnsFilename.c_str(), eInfo, readConversion.get());
			nBandsInited = (readConversion && readConversion->nBandsOld) ? readConversion->nBandsOld : eInfo.nBands;
			isRandom = (nBandsInited<eInfo.nBands);
		}
		else if(initLCAO)
		{	nBandsInited = LCAO();
		}
		
		//Randomize and orthogonalize any uninitialized bands:
		if(nBandsInited < eInfo.nBands)
		{	if(nBandsInited) logPrintf("Setting upper %d bands to ", eInfo.nBands-nBandsInited);
			logPrintf("bandwidth-limited random numbers\n"); logFlush();
			for(int q=eInfo.qStart; q<eInfo.qStop; q++)
			{	//randomize uninitialized bands:
				C[q].randomize(nBandsInited, eInfo.nBands);
				//don't mix (during orthogonalization) the random columns with the init'd ones:
				if(nBandsInited)
				{	//Ensure that the initalized bands are orthonormal:
					ColumnBundle Cfixed = C[q].getSub(0, nBandsInited);
					ColumnBundle OCfixed = O(Cfixed);
					matrix ortho = invsqrt(Cfixed^OCfixed);
					Cfixed = Cfixed * ortho;
					OCfixed = OCfixed * ortho;
					//Project out initalized band directions from the rest:
					C[q] -= Cfixed * (OCfixed^C[q]);
					C[q].setSub(0, Cfixed);
				}
			}
		}
		
		//Orthogonalize initial wavefunctions:
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{	C[q] = C[q] * invsqrt(C[q]^O(C[q]));
			iInfo.project(C[q], VdagC[q]);
		}
	}
	
	//Fluid setup:
	if(fluidParams.fluidType != FluidNone)
	{	logPrintf("----- createFluidSolver() ----- (Fluid-side solver setup)\n");
		fluidSolver = std::shared_ptr<FluidSolver>(createFluidSolver(*e, fluidParams));
		if(!fluidSolver) die("Failed to create fluid solver.\n");
		if(fluidInitialStateFilename.length())
		{	logPrintf("Reading fluid state from '%s'\n", fluidInitialStateFilename.c_str()); logFlush();
			fluidSolver->loadState(fluidInitialStateFilename.c_str());
		}
	}
	
	//Citations:
	if(!e->cntrl.scf)
	{	if(eInfo.fillingsUpdate==ElecInfo::FillingsHsub)
			Citations::add("Total energy minimization with Auxiliary Hamiltonian",
				"C. Freysoldt, S. Boeck, and J. Neugebauer, Phys. Rev. B 79, 241103(R) (2009)");
		else
			Citations::add("Total energy minimization",
				"T.A. Arias, M.C. Payne and J.D. Joannopoulos, Phys. Rev. Lett. 69, 1077 (1992)");
	}
	if(!std::isnan(eInfo.mu))
		Citations::add("Grand-canonical (fixed-potential) DFT",
			"R. Sundararaman, W. A. Goddard III and T. A. Arias, J. Chem. Phys. 146, 114104 (2017)");
}

ScalarFieldArray ElecVars::get_nXC() const
{	if(e->iInfo.nCore)
	{	ScalarFieldArray nXC = clone(n);
		int nSpins = std::min(int(nXC.size()), 2); //1 for unpolarized and 2 for polarized
		for(int s=0; s<nSpins; s++) //note that off-diagonal components of spin-density matrix are excluded
			nXC[s] += (1./nSpins) * e->iInfo.nCore; //add core density
		return nXC;
	}
	else return n; //no cores
}

//Electronic density functional and gradient
void ElecVars::EdensityAndVscloc(Energies& ener, const ExCorr* alternateExCorr)
{	static StopWatch watch("EdensityAndVscloc"); watch.start();
	const ElecInfo& eInfo = e->eInfo;
	const IonInfo& iInfo = e->iInfo;
	
	ScalarFieldTilde nTilde = J(get_nTot());
	
	// Local part of pseudopotential:
	ener.E["Eloc"] = dot(nTilde, O(iInfo.Vlocps));
	ScalarFieldTilde VsclocTilde = clone(iInfo.Vlocps);
	
	// Hartree term:
	ScalarFieldTilde dH = (*e->coulomb)(nTilde); //Note: external charge and nuclear charge contribute to d_vac as well (see below)
	ener.E["EH"] = 0.5*dot(nTilde, O(dH));
	VsclocTilde += dH;

	// External charge:
	ener.E["Eexternal"] = 0.;
	if(rhoExternal)
	{	ScalarFieldTilde phiExternal = (*e->coulomb)(rhoExternal);
		ener.E["Eexternal"] += dot(nTilde + iInfo.rhoIon, O(phiExternal));
		if(rhoExternalSelfEnergy)
			ener.E["Eexternal"] += 0.5 * dot(rhoExternal, O(phiExternal));
		VsclocTilde += phiExternal;
	}
	
	//Fluid contributions
	ScalarFieldTilde VtauTilde;
	if(fluidParams.fluidType != FluidNone)
	{	//Compute n considered for cavity formation (i.e-> include chargeball and partial cores)
		ScalarFieldTilde nCavityTilde = clone(nTilde);
		if(iInfo.nCore) nCavityTilde += J(iInfo.nCore);
		if(iInfo.nChargeball) nCavityTilde += iInfo.nChargeball;

		//Net electric charge:
		ScalarFieldTilde rhoExplicitTilde = nTilde + iInfo.rhoIon + rhoExternal;
		if(!fluidSolver->k2factor) rhoExplicitTilde->setGzero(0.); //No screening => apply neutralizing background charge
		fluidSolver->set(rhoExplicitTilde, nCavityTilde);
		// If the fluid doesn't have a gummel loop, minimize it each time:
		if(!fluidSolver->useGummel()) fluidSolver->minimizeFluid();
		
		// Compute the energy and accumulate gradients:
		ener.E["A_diel"] = fluidSolver->get_Adiel_and_grad(&d_fluid, &V_cavity);
		VsclocTilde += d_fluid;
		VsclocTilde += V_cavity;

		//Chemical-potential correction due to potential of electron in bulk fluid
		double bulkPotential = fluidSolver->bulkPotential();
		//Chemical-potential correction due to finite nuclear width in fluid interaction:
		double muCorrection = fluidSolver->ionWidthMuCorrection();
		ener.E["MuShift"] = (eInfo.nElectrons - iInfo.getZtot()) * (muCorrection - bulkPotential);
		VsclocTilde->setGzero(muCorrection - bulkPotential + VsclocTilde->getGzero());
	}

	//Atomic density-matrix contributions: (DFT+U)
	if(eInfo.hasU)
		ener.E["U"] = iInfo.rhoAtom_computeU(rhoAtom, U_rhoAtom);
	
	// Exchange and correlation, and store the real space Vscloc with the odd historic normalization factor of JdagOJ:
	ScalarFieldArray Vxc;
	const ExCorr& exCorr = alternateExCorr ? *alternateExCorr : e->exCorr;
	ener.E["Exc"] = exCorr(get_nXC(), &Vxc, false, &tau, &Vtau);
	if(!exCorr.hasEnergy() && !e->cntrl.scf)
		die("Potential functionals do not support total-energy minimization; use SCF instead.\n")
	if(exCorr.orbitalDep)
	{	if(!e->cntrl.scf)
		{	if(e->cntrl.fixed_H) die("Orbital-dependent potential functionals do not support fix-density; use fix-potential instead.\n")
			else die("Orbital-dependent potential functionals do not support total-energy minimization; use SCF instead.\n")
		}
		Vxc += exCorr.orbitalDep->getPotential();
	}
	if(VtauTilde) Vtau.resize(n.size());
	for(unsigned s=0; s<Vscloc.size(); s++)
	{	Vscloc[s] = JdagOJ(Vxc[s]);
		if(s<2) //Include all the spin-independent contributions along the diagonal alone
			Vscloc[s] += Jdag(O(VsclocTilde), true);
		//External potential contributions:
		if(Vexternal.size())
		{	ener.E["Eexternal"] += e->gInfo.dV * dot(n[s], Vexternal[s]);
			Vscloc[s] += JdagOJ(Vexternal[s]);
		}
		e->symm.symmetrize(Vscloc[s]);
		if(VtauTilde) Vtau[s] += I(VtauTilde);
		if(Vtau[s]) e->symm.symmetrize(Vtau[s]);
	}
	watch.stop();
}


//-----  Electronic energy and (preconditioned) gradient calculation ----------

double ElecVars::elecEnergyAndGrad(Energies& ener, ElecGradient* grad, ElecGradient* Kgrad, bool calc_Hsub)
{	
	const ElecInfo& eInfo = e->eInfo;
	
	//Cleanup old gradients:
	if(grad) { grad->C.assign(eInfo.nStates, ColumnBundle()); }
	if(Kgrad) { Kgrad->C.assign(eInfo.nStates, ColumnBundle()); }
	
	//Determine whether Hsub and hence HC needs to be calculated:
	bool need_Hsub = calc_Hsub || grad;
	double mu = 0., Bz = 0.;
	
	if(eInfo.fillingsUpdate==ElecInfo::FillingsHsub)
	{	//Update nElectrons from mu, or mu from nElectrons as appropriate:
		if(std::isnan(eInfo.mu)) mu = eInfo.findMu(Haux_eigs, eInfo.nElectrons, Bz);
		else { mu = eInfo.mu; ((ElecInfo&)eInfo).nElectrons = eInfo.nElectronsCalc(mu, Haux_eigs, Bz); }
		//Compute fillings from aux hamiltonian eigenvalues:
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
			F[q] = eInfo.smear(eInfo.muEff(mu,Bz,q), Haux_eigs[q]);
		//Update TS and muN:
		eInfo.updateFillingsEnergies(Haux_eigs, ener);
		//Report for SCF (ElecMinimizer handles for minimize):
		if(e->cntrl.scf && n[0]) eInfo.smearReport();
	}
	
	//Update the density and density-dependent pieces if required:
	n = calcDensity();
	if(e->exCorr.needsKEdensity()) tau = KEdensity();
	if(eInfo.hasU) e->iInfo.rhoAtom_calc(F, C, rhoAtom); //Atomic density matrix contributions for DFT+U
	EdensityAndVscloc(ener); //Calculate density functional and its gradient
	if(need_Hsub) e->iInfo.augmentDensityGridGrad(Vscloc); //Update Vscloc projected onto spherical functions for ultrasoft psps
	
	//--------- Wavefunction dependent parts -----------
	
	std::vector<ColumnBundle> HC(eInfo.nStates); //gradient w.r.t C (upto weights and fillings)
	
	//Exact exchange if required:
	ener.E["EXX"] = 0.;
	if(e->exCorr.exxFactor())
	{	double aXX = e->exCorr.exxFactor();
		double omega = e->exCorr.exxRange();
		assert(e->exx);
		ener.E["EXX"] = (*e->exx)(aXX, omega, F, C, need_Hsub ? &HC : 0);
	}
	
	//Do the single-particle contributions one state at a time to save memory (and for better cache warmth):
	ener.E["KE"] = 0.;
	ener.E["Enl"] = 0.;
	for(int q=eInfo.qStart; q<e->eInfo.qStop; q++)
	{	double KEq = applyHamiltonian(q, F[q], HC[q], ener, need_Hsub);
		if(grad) //Calculate wavefunction gradients:
		{	const QuantumNumber& qnum = eInfo.qnums[q];
			HC[q] -= O(C[q]) * Hsub[q]; //Include orthonormality contribution
			grad->C[q] = HC[q] * (F[q]*qnum.weight);
			if(Kgrad)
			{	double Nq = qnum.weight*trace(F[q]);
				double KErollover = 2. * (Nq>1e-3 ? KEq/Nq : 1.);
				precond_inv_kinetic(HC[q], KErollover); //apply preconditioner
				std::swap(Kgrad->C[q], HC[q]); //this frees HC[q]
			}
		}
	}
	mpiUtil->allReduce(ener.E["KE"], MPIUtil::ReduceSum);
	mpiUtil->allReduce(ener.E["Enl"], MPIUtil::ReduceSum);
	
	double dmuContrib = 0., dBzContrib = 0.;
	if(grad and eInfo.fillingsUpdate==ElecInfo::FillingsHsub and (std::isnan(eInfo.mu) or eInfo.Mconstrain)) //contribution due to N/M constraint via the mu/Bz gradient 
	{	double dmuNum[2] = {0.,0.}, dmuDen[2] = {0.,0.}; //numerator and denominator of dmuContrib resolved by spin channels (if any)
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{	diagMatrix fprime = eInfo.smearPrime(eInfo.muEff(mu,Bz,q), Haux_eigs[q]);
			double w = eInfo.qnums[q].weight;
			int sIndex = eInfo.qnums[q].index();
			dmuNum[sIndex] += w * trace(fprime * (diag(Hsub[q])-Haux_eigs[q]));
			dmuDen[sIndex] += w * trace(fprime);
		}
		mpiUtil->allReduce(dmuNum, 2, MPIUtil::ReduceSum);
		mpiUtil->allReduce(dmuDen, 2, MPIUtil::ReduceSum);
		if(std::isnan(eInfo.mu) and eInfo.Mconstrain)
		{	//Fixed N and M (effectively independent constraints on Nup and Ndn)
			double dmuContribUp = dmuNum[0]/dmuDen[0];
			double dmuContribDn = dmuNum[1]/dmuDen[1];
			dmuContrib = 0.5*(dmuContribUp + dmuContribDn);
			dBzContrib = 0.5*(dmuContribUp - dmuContribDn);
		}
		else if(eInfo.Mconstrain)
		{	//Fixed M only
			dmuContrib = 0.;
			dBzContrib = (dmuNum[0]-dmuNum[1])/(dmuDen[0]-dmuDen[1]);
		}
		else
		{	//Fixed N only
			dmuContrib = (dmuNum[0]+dmuNum[1])/(dmuDen[0]+dmuDen[1]);
			dBzContrib = 0.;
		}
	}
	
	//Auxiliary hamiltonian gradient:
	if(grad && eInfo.fillingsUpdate==ElecInfo::FillingsHsub)
	{	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{	const QuantumNumber& qnum = eInfo.qnums[q];
			matrix gradF0 = Hsub[q]-Haux_eigs[q]; //gradient w.r.t fillings except for constraint contributions
			matrix gradF = gradF0 - eye(eInfo.nBands)*eInfo.muEff(dmuContrib,dBzContrib,q); //gradient w.r.t fillings
			grad->Haux[q] = qnum.weight * dagger_symmetrize(eInfo.smearGrad(eInfo.muEff(mu,Bz,q), Haux_eigs[q], gradF));
			if(Kgrad) //Drop the fermiPrime factors in preconditioned gradient:
				Kgrad->Haux[q] = (-e->cntrl.subspaceRotationFactor) * gradF0;
		}
	}
	
	return relevantFreeEnergy(*e);
}

//Make phase (and degenerate-subspace rotations) of wavefunctions reproducible 
void fixPhase(matrix& evecs, const diagMatrix& eigs, const ColumnBundle& C)
{	const double tol = 1e-4;
	//Pick out the head elements:
	const std::vector<int>& head = C.basis->head;
	int nSpinor = C.spinorLength();
	matrix Chead = zeroes(head.size()*nSpinor, C.nCols());
	for(size_t n=0; n<head.size(); n++)
		for(int s=0; s<nSpinor; s++)
			callPref(eblas_zaxpy)(C.nCols(), 1.,
				C.dataPref() + n + s*C.basis->nbasis, C.colLength(), 
				Chead.dataPref() + nSpinor*n + s, Chead.nRows() );
	Chead = Chead * evecs;
	//Hamiltonian for degeneracy breaking (some definite (and unlikely to be symmetric) diagonal in the low G head)
	diagMatrix headH(Chead.nRows());
	for(int n=0; n<Chead.nRows(); n++)
		headH[n] = pow(n, M_PI);
	//Resolve degeneracies
	matrix degFix = eye(C.nCols());
	bool degFound = false;
	for(int bStart=0; bStart<C.nCols()-1;)
	{	int bStop=bStart;
		while(bStop<eigs.nCols() && eigs[bStop]<eigs[bStart]+tol) bStop++;
		if(bStop-bStart > 1) //degeneracy
		{	degFound = true;
			matrix CheadSub = Chead(0,Chead.nRows(), bStart,bStop);
			matrix degEvecs; diagMatrix degEigs;
			(dagger(CheadSub) * headH * CheadSub).diagonalize(degEvecs, degEigs);
			degFix.set(bStart,bStop, bStart,bStop, degEvecs);
		}
		bStart = bStop;
	}
	if(degFound)
	{	evecs = evecs*degFix;
		Chead = Chead*degFix;
	}
	//Fix phase:
	matrix phaseFix = eye(C.nCols());
	for(int b=0; b<C.nCols(); b++)
	{	complex phase;
		double normPrev = 0;
		for(int n=0; n<Chead.nRows(); n++)
		{	const complex c = Chead(n,b);
			if(c.norm() > normPrev)
			{	phase = c.conj()/c.abs();
				normPrev = c.norm();
			}
		}
		phaseFix.set(b,b, phase);
	}
	evecs = evecs*phaseFix;
}

void ElecVars::setEigenvectors()
{	const ElecInfo& eInfo = e->eInfo;
	logPrintf("Setting wave functions to eigenvectors of Hamiltonian\n"); logFlush();
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
	{
		fixPhase(Hsub_evecs[q], Hsub_eigs[q], C[q]);
		C[q] = C[q] * Hsub_evecs[q];
		for(matrix& VdagCq_sp: VdagC[q])
			if(VdagCq_sp) VdagCq_sp = VdagCq_sp * Hsub_evecs[q];
		
		if(eInfo.fillingsUpdate==ElecInfo::FillingsHsub && !e->cntrl.scf)
			Haux_eigs[q] = Hsub_eigs[q];
		
		//Apply corresponding changes to Hsub:
		Hsub[q] = Hsub_eigs[q]; //now diagonal
		Hsub_evecs[q] = eye(eInfo.nBands);
	}
}

ScalarFieldArray ElecVars::KEdensity() const
{	ScalarFieldArray tau(n.size());
	//Compute KE density from valence electrons:
	for(int q=e->eInfo.qStart; q<e->eInfo.qStop; q++)
		for(int iDir=0; iDir<3; iDir++)
			tau += (0.5*C[q].qnum->weight) * diagouterI(F[q], D(C[q],iDir), tau.size(), &e->gInfo);
	for(ScalarField& tau_s: tau)
	{	nullToZero(tau_s, e->gInfo);
		e->symm.symmetrize(tau_s); //Symmetrize
		tau_s->allReduce(MPIUtil::ReduceSum);
	}
	//Add core KE density model:
	if(e->iInfo.tauCore)
	{	for(unsigned s=0; s<tau.size(); s++)
			tau[s] += (1.0/tau.size()) * e->iInfo.tauCore; //add core KE density
	}
	return tau;
}

ScalarFieldArray ElecVars::calcDensity() const
{	ScalarFieldArray density(n.size());
	//Runs over all states and accumulates density to the corresponding spin channel of the total density
	e->iInfo.augmentDensityInit();
	for(int q=e->eInfo.qStart; q<e->eInfo.qStop; q++)
	{	density += e->eInfo.qnums[q].weight * diagouterI(F[q], C[q], density.size(), &e->gInfo);
		e->iInfo.augmentDensitySpherical(e->eInfo.qnums[q], F[q], VdagC[q]); //pseudopotential contribution
	}
	e->iInfo.augmentDensityGrid(density);
	
	for(ScalarField& ns: density)
	{	nullToZero(ns, e->gInfo);
		e->symm.symmetrize(ns); //Symmetrize
		ns->allReduce(MPIUtil::ReduceSum);
	}
	return density;
}

void ElecVars::orthonormalize(int q, matrix* extraRotation)
{	assert(e->eInfo.isMine(q));
	VdagC[q].clear();
	matrix rot = invsqrt(C[q]^O(C[q], &VdagC[q])); //Compute U:
	if(extraRotation) *extraRotation = (rot = rot * (*extraRotation)); //set rot and extraRotation to the net transformation
	C[q] = C[q] * rot;
	e->iInfo.project(C[q], VdagC[q], &rot); //update the atomic projections
}

double ElecVars::applyHamiltonian(int q, const diagMatrix& Fq, ColumnBundle& HCq, Energies& ener, bool need_Hsub)
{	assert(C[q]); //make sure wavefunction is available for this states
	const QuantumNumber& qnum = e->eInfo.qnums[q];
	std::vector<matrix> HVdagCq(e->iInfo.species.size());
	
	//Propagate grad_n (Vscloc) to HCq (which is grad_Cq upto weights and fillings) if required
	if(need_Hsub)
	{	HCq += Idag_DiagV_I(C[q], Vscloc); //Accumulate Idag Diag(Vscloc) I C
		e->iInfo.augmentDensitySphericalGrad(qnum, Fq, VdagC[q], HVdagCq); //Contribution via pseudopotential density augmentation
		if(e->exCorr.needsKEdensity() && Vtau[qnum.index()]) //Contribution via orbital KE:
		{	for(int iDir=0; iDir<3; iDir++)
				HCq -= (0.5*e->gInfo.dV) * D(Idag_DiagV_I(D(C[q],iDir), Vtau), iDir);
		}
		if(e->eInfo.hasU) //Contribution via atomic density matrix projections (DFT+U)
			e->iInfo.rhoAtom_grad(C[q], U_rhoAtom, HCq);
	}

	//Kinetic energy:
	double KEq;
	{	ColumnBundle LCq = L(C[q]);
		if(HCq) HCq += (-0.5) * LCq;
		KEq = qnum.weight * (-0.5) * traceinner(Fq, C[q], LCq).real();
		ener.E["KE"] += KEq;
	}
	
	//Nonlocal pseudopotentials:
	ener.E["Enl"] += qnum.weight * e->iInfo.EnlAndGrad(qnum, Fq, VdagC[q], HVdagCq);
	if(HCq) e->iInfo.projectGrad(HVdagCq, C[q], HCq);
	
	//Compute subspace hamiltonian if needed:
	if(need_Hsub)
	{	Hsub[q] = C[q] ^ HCq;
		Hsub[q].diagonalize(Hsub_evecs[q], Hsub_eigs[q]);
	}
	return KEq;
}
