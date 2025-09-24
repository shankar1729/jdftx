/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver

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

#include <electronic/Dump.h>
#include <electronic/Dump_internal.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/SpeciesInfo.h>
#include <electronic/ExactExchange.h>
#include <electronic/DOS.h>
#include <electronic/Polarizability.h>
#include <electronic/ElectronScattering.h>
#include <electronic/LatticeMinimizer.h>
#include <fluid/FluidSolver.h>
#include <core/VectorField.h>
#include <core/ScalarFieldIO.h>
#include <ctime>

Dump::Dump()
: potentialSubtraction(true), bandProjectionOrtho(false), bandProjectionNorm(true), Munfold(1,1,1), curIter(0)
{
}

void Dump::setup(const Everything& everything)
{	e = &everything;
	if(dos) dos->setup(everything);
	
	//Add some citations here so that they are included in a dry run:
	for(auto dumpPair: *this)
		switch(dumpPair.second)
		{	case DumpQMC:
			{	Citations::add("Quantum Monte Carlo solvation",
					"K.A. Schwarz, R. Sundararaman, K. Letchworth-Weaver, T.A. Arias and R. Hennig, Phys. Rev. B 85, 201102(R) (2012)");
				break;
			}
			case DumpDvac:
			case DumpDtot:
			case DumpBulkEpsilon:
			case DumpSlabEpsilon:
			case DumpChargedDefect:
			{	if(potentialSubtraction)
					Citations::add("Smooth electrostatic potentials by atom-potential subtraction",
						"R. Sundararaman and Y. Ping, J. Chem. Phys. 146, 104109 (2017)");
				break;
			}
			default:; //No action necessary (needed only to suppress compiler warnings)
		}
}


void Dump::operator()(DumpFrequency freq, int iter)
{
	if(!checkInterval(freq, iter)) return; // => don't dump this time
	curIter = iter; curFreq = freq; //used by getFilename()
	
	bool foundVars = false; //whether any variables are to be dumped at this frequency
	for(auto entry: *this)
		if(entry.first==freq && entry.second!=DumpNone)
		{	foundVars = true;
			break;
		}
	if(!foundVars) return;
	logPrintf("\n");
	
	const ElecInfo &eInfo = e->eInfo;
	const ElecVars &eVars = e->eVars;
	const IonInfo &iInfo = e->iInfo;

	bool isCevec = (freq==DumpFreq_Gummel) || (freq==DumpFreq_Ionic) || (freq==DumpFreq_End); //whether C has been set to eigenfunctions
	bool hasFluid = (eVars.fluidParams.fluidType != FluidNone);
	
	//Macro to determine which variables should be dumped:
	//(Drop the "Dump" from the DumpVariables enum name to use as var)
	#define ShouldDump(var) count(std::make_pair(freq,Dump##var))

	#define StartDump(prefix) \
		string fname = getFilename(prefix); \
		logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();

	#define EndDump \
		logPrintf("done\n"); logFlush();

	#define DUMP_nocheck(object, prefix) \
		{	StartDump(prefix) \
			if(mpiWorld->isHead()) saveRawBinary(object, fname.c_str()); \
			EndDump \
		}
	
	#define DUMP_spinCollection(object, prefix) \
		{	if(object.size()==1) DUMP_nocheck(object[0], prefix) \
			else if(object.size()!=0) \
			{	DUMP_nocheck(object[0], (prefix+string("_up")).c_str()) \
				DUMP_nocheck(object[1], (prefix+string("_dn")).c_str()) \
				if(object.size()==4) \
				{	DUMP_nocheck(object[2], (prefix+string("_re")).c_str()) \
					DUMP_nocheck(object[3], (prefix+string("_im")).c_str()) \
				} \
			} \
		}

	#define DUMP(object, prefix, varname) \
		if(ShouldDump(varname) && object) \
		{	DUMP_nocheck(object, prefix) \
		}
	
	// Set up date/time stamp
	time_t timenow = time(0);
	tm *mytm = localtime(&timenow);
	ostringstream stampStream;
	stampStream << (mytm->tm_mon+1) << '.' << mytm->tm_mday << '.'
		<< mytm->tm_hour << '-' << mytm->tm_min << '-' << mytm->tm_sec;
	stamp = stampStream.str();
	
	if((ShouldDump(State) and eInfo.fillingsUpdate==ElecInfo::FillingsHsub) or ShouldDump(Fillings))
	{	//Dump fillings
		double wInv = eInfo.spinType==SpinNone ? 0.5 : 1.0; //normalization factor from external to internal fillings
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) ((ElecVars&)eVars).F[q] *= (1./wInv);
		StartDump("fillings")
		eInfo.write(eVars.F, fname.c_str());
		EndDump
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) ((ElecVars&)eVars).F[q] *= wInv;
	}
	
	if(ShouldDump(State))
	{
		//Dump wave functions
		StartDump("wfns")
		eInfo.write(eVars.C, fname.c_str());
		EndDump
		
		if(hasFluid)
		{	//Dump state of fluid:
			StartDump("fluidState")
			if(mpiWorld->isHead()) eVars.fluidSolver->saveState(fname.c_str());
			EndDump
		}
	}

	if(ShouldDump(IonicPositions)
		|| (ShouldDump(State)
			&& (e->ionicMinParams.nIterations>0
				|| e->latticeMinParams.nIterations>0
				|| e->ionicDynParams.nSteps>0 ) ) )
	{	StartDump("ionpos")
		FILE* fp = mpiWorld->isHead() ? fopen(fname.c_str(), "w") : nullLog;
		if(!fp) die("Error opening %s for writing.\n", fname.c_str());
		iInfo.printPositions(fp);  //needs to be called from all processes (for magnetic moment computation)
		if(mpiWorld->isHead())fclose(fp);
		EndDump
		if(iInfo.metadynamicsBond) iInfo.metadynamicsBond->save();
	}
	if(ShouldDump(Forces))
	{	StartDump("force")
		if(mpiWorld->isHead())
		{	FILE* fp = fopen(fname.c_str(), "w");
			if(!fp) die("Error opening %s for writing.\n", fname.c_str());
			iInfo.forces.print(*e, fp);
			fclose(fp);
		}
		EndDump
	}
	if(ShouldDump(Lattice) || (ShouldDump(State) && e->latticeMinParams.nIterations>0))
	{	StartDump("lattice")
		if(mpiWorld->isHead()) 
		{	FILE* fp = fopen(fname.c_str(), "w");
			if(!fp) die("Error opening %s for writing.\n", fname.c_str());
			fprintf(fp, "lattice");
			for(int j=0; j<3; j++)
			{	fprintf(fp, " \\\n\t");
				for(int k=0; k<3; k++)
					fprintf(fp, "%20.15lf ", e->gInfo.R(j,k));
			}
			fprintf(fp, "#Note: latt-scale has been absorbed into these lattice vectors.\n");
			fclose(fp);
		}
		EndDump
	}
	DUMP(I(iInfo.rhoIon), "Nion", IonicDensity)
	
	if(ShouldDump(ElecDensity))
		DUMP_spinCollection(eVars.n, "n")
	if(ShouldDump(ElecDensityAccum))
		DUMP_spinCollection(eVars.nAccum, "nAccum")
	if(iInfo.nCore) DUMP(iInfo.nCore, "nCore", CoreDensity)
	
	if((ShouldDump(KEdensity) or (e->exCorr.needsKEdensity() and ShouldDump(ElecDensity))))
	{	const auto& tau = (e->exCorr.needsKEdensity() ? e->eVars.tau : e->eVars.KEdensity());
		 DUMP_spinCollection(tau, "tau")
	}

	//Electrostatic and fluid potentials:
	ScalarFieldTilde d_vac; ScalarField d_tot;
	bool needDtot = ShouldDump(Dtot)
		|| ShouldDump(SlabEpsilon)
		|| ShouldDump(BulkEpsilon)
		|| ShouldDump(ChargedDefect);
	if(ShouldDump(Dvac) || needDtot)
	{	d_vac = iInfo.Vlocps + (*e->coulomb)(J(eVars.get_nTot())); //local pseudopotential + Hartree term
		//Subtract neutral-atom reference potential (gives smoother result):
		if(potentialSubtraction)
		{	ScalarFieldTilde dAtomic;
			for(auto sp: e->iInfo.species) if(sp->atpos.size()) sp->accumulateAtomicPotential(dAtomic);
			d_vac -= dAtomic;
		}
		if(eVars.rhoExternal) d_vac += (*e->coulomb)(eVars.rhoExternal); //potential due to external charge (if any)
		if(e->coulombParams.Efield.length_squared()) d_vac += J(e->coulomb->getEfieldPotential());
	}
	DUMP(I(d_vac), "d_vac", Dvac);
	if(hasFluid)
	{	if(ShouldDump(Dfluid) || needDtot)
		{	double GzeroCorrection = eVars.fluidSolver->ionWidthMuCorrection() - eVars.fluidSolver->bulkPotential();
			ScalarFieldTilde d_fluid = clone(eVars.d_fluid);
			if(eVars.fluidSolver->A_rhoNonES)
				d_fluid -= eVars.fluidSolver->A_rhoNonES;
			DUMP(I(d_fluid), "d_fluid", Dfluid);
			if(needDtot) d_tot = I(d_vac + d_fluid) + GzeroCorrection;
			DUMP(d_tot, "d_tot", Dtot);
		}
		DUMP(I(eVars.V_cavity), "V_cavity", Vcavity);
		DUMP(I(eVars.V_cavity + eVars.d_fluid), "V_fluidTot", VfluidTot);
	}
	else
	{	if(needDtot)
		{	d_tot = I(d_vac);
			DUMP(d_tot, "d_tot", Dtot);
		}
	}
	if(ShouldDump(SlabEpsilon))
		if(slabEpsilon)
			slabEpsilon->dump(*e, d_tot);
	if(ShouldDump(BulkEpsilon))
		if(bulkEpsilon)
			bulkEpsilon->dump(*e, d_tot);
	if(ShouldDump(ChargedDefect))
		if(chargedDefect)
			chargedDefect->dump(*e, d_tot);
	d_tot = 0;

	DUMP(I(iInfo.Vlocps), "Vlocps", Vlocps)
	if(ShouldDump(Vscloc))
		DUMP_spinCollection(eVars.Vscloc, "Vscloc")
	if(ShouldDump(Vscloc) and e->exCorr.needsKEdensity())
		DUMP_spinCollection(eVars.Vtau, "Vtau")
	
	if(ShouldDump(BandEigs) ||
		(ShouldDump(State) &&
			( (eInfo.fillingsUpdate == ElecInfo::FillingsHsub)
			|| (e->exCorr.orbitalDep && isCevec) ) ) )
	{	StartDump("eigenvals")
		eInfo.write(eVars.Hsub_eigs, fname.c_str());
		EndDump
	}
	
	if(ShouldDump(BandProjections) && isCevec)
	{	StartDump("bandProjections")
		dumpProjections(*e, fname.c_str(), bandProjectionOrtho, bandProjectionNorm);
		EndDump
	}
	
	if(ShouldDump(EigStats))
	{	StartDump("eigStats")
		double Emin = +INFINITY; int qEmin = 0;
		double Emax = -INFINITY; int qEmax = 0;
		double HOMO = -INFINITY; int qHOMO = 0;
		double LUMO = +INFINITY; int qLUMO = 0;
		double gap = +INFINITY; int qGap = 0;
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{	double curEmin = +INFINITY, curEmax = -INFINITY;
			double curHOMO = -INFINITY, curLUMO = +INFINITY;
			for(int b=0; b<eInfo.nBands; b++)
			{	double F = eVars.F[q][b], E = eVars.Hsub_eigs[q][b];
				curEmin = std::min(curEmin, E);
				curEmax = std::max(curEmax, E);
				if(F>=0.5) curHOMO = std::max(curHOMO, E);
				if(F<=0.5) curLUMO = std::min(curLUMO, E);
			}
			double curGap = curLUMO - curHOMO; //optical gap at current q
			if(curEmin < Emin) { Emin=curEmin; qEmin=q; }
			if(curEmax > Emax) { Emax=curEmax; qEmax=q; }
			if(curHOMO > HOMO) { HOMO=curHOMO; qHOMO=q; }
			if(curLUMO < LUMO) { LUMO=curLUMO; qLUMO=q; }
			if(curGap < gap) { gap=curGap; qGap=q; }
		}
		mpiWorld->allReduce(Emin, qEmin, MPIUtil::ReduceMin);
		mpiWorld->allReduce(Emax, qEmax, MPIUtil::ReduceMax);
		mpiWorld->allReduce(HOMO, qHOMO, MPIUtil::ReduceMax);
		mpiWorld->allReduce(LUMO, qLUMO, MPIUtil::ReduceMin);
		mpiWorld->allReduce(gap, qGap, MPIUtil::ReduceMin);
		double gapIndirect = LUMO - HOMO;
		double mu = NAN, Bz;
		if(std::isfinite(HOMO) && std::isfinite(LUMO))
			mu = (!std::isnan(eInfo.mu)) ? eInfo.mu : eInfo.findMu(e->eVars.Hsub_eigs, eInfo.nElectrons, Bz);
		//Print results:
		FILE* fp = 0;
		if(mpiWorld->isHead()) fp = fopen(fname.c_str(), "w");
		logPrintf("\n");
		#define teePrintf(...) \
			{	fprintf(globalLog, "\t" __VA_ARGS__); \
				if(mpiWorld->isHead()) fprintf(fp, __VA_ARGS__); \
			}
		#define printQuantity(name, value, q) \
			if(std::isfinite(value)) \
			{	const QuantumNumber& qnum = eInfo.qnums[q]; \
				teePrintf(name ": %+.6lf at state %d ( [ %+.6lf %+.6lf %+.6lf ] spin %d )\n", \
					value, q, qnum.k[0], qnum.k[1], qnum.k[2], qnum.spin); \
			} \
			else teePrintf(name ": unavailable\n");
		printQuantity("eMin", Emin, qEmin)
		printQuantity("HOMO", HOMO, qHOMO)
		if(std::isfinite(mu)) teePrintf("mu  : %+.6lf\n", mu) else teePrintf("mu  : unavailable\n")
		printQuantity("LUMO", LUMO, qLUMO)
		printQuantity("eMax", Emax, qEmax)
		if(std::isfinite(gapIndirect)) teePrintf("HOMO-LUMO gap: %+.6lf\n", gapIndirect) else teePrintf("HOMO-LUMO gap: unavailable\n")
		printQuantity("Optical gap  ", gap, qGap)
		#undef printQuantity
		#undef teePrintf
		if(mpiWorld->isHead()) fclose(fp);
		logFlush();
	}
	
	if(eInfo.hasU && (ShouldDump(RhoAtom) || ShouldDump(ElecDensity)))
	{	StartDump("rhoAtom")
		if(mpiWorld->isHead())
		{	FILE* fp = fopen(fname.c_str(), "w");
			for(const matrix& m: eVars.rhoAtom) m.write(fp);
			fclose(fp);
		}
		EndDump
	}
	
	if(eInfo.hasU && ShouldDump(Vscloc))
	{	StartDump("U_rhoAtom")
		if(mpiWorld->isHead())
		{	FILE* fp = fopen(fname.c_str(), "w");
			for(const matrix& m: eVars.U_rhoAtom) m.write(fp);
			fclose(fp);
		}
		EndDump
	}
	
	if(ShouldDump(Ecomponents))
	{	StartDump("Ecomponents")
		if(mpiWorld->isHead())
		{	FILE* fp = fopen(fname.c_str(), "w");
			if(!fp) die("Error opening %s for writing.\n", fname.c_str());	
			e->ener.print(fp);
			fclose(fp);
		}
		EndDump
	}
	
	if(ShouldDump(BandUnfold))
		dumpUnfold();
	
	if((ShouldDump(BoundCharge) || ShouldDump(SolvationRadii)) && hasFluid)
	{	ScalarFieldTilde nboundTilde = (-1.0/(4*M_PI*e->gInfo.detR)) * L(eVars.d_fluid);
		if(iInfo.ionWidth) nboundTilde = gaussConvolve(nboundTilde, iInfo.ionWidth);
		double rhoTot_Gzero = J(eVars.get_nTot())->getGzero() + iInfo.rhoIon->getGzero();
		if(eVars.rhoExternal) rhoTot_Gzero += eVars.rhoExternal->getGzero();
		nboundTilde->setGzero(-rhoTot_Gzero); //total bound charge will neutralize system
		if(ShouldDump(SolvationRadii))
		{	StartDump("Rsol")
			if(mpiWorld->isHead()) dumpRsol(I(nboundTilde), fname);
			EndDump
		}
		DUMP(I(nboundTilde), "nbound", BoundCharge)
	}
	
	if(ShouldDump(FluidDensity) && hasFluid)
		eVars.fluidSolver->dumpDensities(getFilename("fluid%s").c_str());
	
	if(ShouldDump(QMC) && isCevec) dumpQMC();
	if(ShouldDump(Ocean) && isCevec) dumpOcean();
	if(ShouldDump(BGW) && isCevec) dumpBGW();
	
	if(ShouldDump(Dipole))
	{	StartDump("Moments")
		if(mpiWorld->isHead()) dumpMoment(*e, fname.c_str());
		EndDump
	}
	
	if(ShouldDump(SIC))
	{	DumpSelfInteractionCorrection selfInteractionCorrection(*e);
		selfInteractionCorrection.dump(getFilename("SIC").c_str());
	}

	if(ShouldDump(Symmetries))
	{	StartDump("sym")
		if(mpiWorld->isHead())
		{	FILE* fp = fopen(fname.c_str(), "w");
			const std::vector<SpaceGroupOp>& sym = e->symm.getMatrices();
			for(const SpaceGroupOp& op: sym)
			{	op.rot.print(fp, "%2d ", false);
				fprintf(fp, "%lg %lg %lg\n", op.a[0], op.a[1], op.a[2]);
				fprintf(fp, "\n");
			}
			fclose(fp);
		}
		EndDump
	}
	
	if(ShouldDump(Kpoints))
	{	StartDump("kPts")
		if(mpiWorld->isHead())
		{	FILE* fp = fopen(fname.c_str(), "w");
			eInfo.kpointsPrint(fp, true);
			fclose(fp);
		}
		EndDump
		if(e->symm.mode != SymmetriesNone)
		{	StartDump("kMap")
			if(mpiWorld->isHead())
			{	FILE* fp = fopen(fname.c_str(), "w");
				e->symm.printKmap(fp);
				fclose(fp);
			}
			EndDump
		}
	}
	
	if(ShouldDump(Gvectors))
	{	StartDump("Gvectors")
		if(mpiWorld->isHead())
		{	FILE* fp = fopen(fname.c_str(), "w");
			for(int q=0; q<eInfo.nStates; q++)
			{	//Header:
				fprintf(fp, "# ");
				eInfo.kpointPrint(fp, q, true);
				fprintf(fp, "\n");
				//G-vectors
				const Basis& basis = e->basis[q];
				for(const vector3<int>& iG: basis.iGarr)
					fprintf(fp, "%d %d %d\n", iG[0], iG[1], iG[2]);
				fprintf(fp, "\n");
			}
			fclose(fp);
		}
		EndDump
	}
	
	if(ShouldDump(OrbitalDep) && e->exCorr.orbitalDep && isCevec)
		e->exCorr.orbitalDep->dump();
	
	if(ShouldDump(Excitations))
	{	if(e->eInfo.isNoncollinear()) logPrintf("WARNING: Excitations dump not supported with noncollinear spins (Skipping)\n");
		else
		{	StartDump("Excitations")
			dumpExcitations(*e, fname.c_str());
			EndDump
		}
	}
	
	if(ShouldDump(FCI))
	{	StartDump("FCI")
		dumpFCI(*e, fname.c_str());
		EndDump
	}
	
	if(ShouldDump(Spin) and eInfo.isNoncollinear())
	{	StartDump("S")
		std::vector<matrix> S(eInfo.nStates);
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) //kpoint/spin
		{	//Calculate spin overlap:
			vector3<matrix> Sq = spinOverlap(eVars.C[q]);
			//Reshape to a N x 3N matrix:
			S[q] = zeroes(eInfo.nBands, eInfo.nBands*3);
			for(int k=0; k<3; k++) //cartesian direction
				S[q].set(0,eInfo.nBands, eInfo.nBands*k,eInfo.nBands*(k+1), Sq[k]);
		}
		eInfo.write(S, fname.c_str(), eInfo.nBands, eInfo.nBands*3);
		EndDump
	}
	
	if(ShouldDump(Momenta) or ShouldDump(Velocities) or ShouldDump(FermiVelocity))
	{
		//Common code: compute matrix elements:
		std::vector<matrix> momenta(eInfo.nStates); //full matrices
		std::vector<std::vector<vector3<>>> v(eInfo.nStates, std::vector<vector3<>>(eInfo.nBands)); //diagonal parts
		bool needFull = ShouldDump(Momenta); //whether full matrices are needed
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) //kpoint/spin
		{	if(needFull) momenta[q] = zeroes(eInfo.nBands, eInfo.nBands*3);
			for(int iDir=0; iDir<3; iDir++) //cartesian direction
			{	matrix Pqk = complex(0,-1) * iInfo.rHcommutator(eVars.C[q], iDir, eVars.Hsub_eigs[q]);
				if(needFull)
					momenta[q].set(0,eInfo.nBands, eInfo.nBands*iDir,eInfo.nBands*(iDir+1), Pqk);
				for(int b=0; b<eInfo.nBands; b++)
					v[q][b][iDir] = Pqk(b,b).real();
			}
		}
		
		//Output full momentum matrix elements:
		if(ShouldDump(Momenta))
		{	StartDump("momenta")
			eInfo.write(momenta, fname.c_str(), eInfo.nBands, eInfo.nBands*3);
			EndDump
		}
		
		//Output diagonal band-velocity matrix elements:
		if(ShouldDump(Velocities))
		{	StartDump("velocities")
			FILE* fp;
			if(mpiWorld->isHead())
			{	fp = fopen(fname.c_str(), "w");
				if(!fp) die_alone("Error opening %s for writing.\n", fname.c_str());
				for(int q=0; q<eInfo.nStates; q++)
				{	if(not eInfo.isMine(q)) mpiWorld->recvData(v[q], eInfo.whose(q), q);
					fwrite(v[q].data(), sizeof(vector3<>), eInfo.nBands, fp);
				}
				fclose(fp);
			}
			else
			{	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
					mpiWorld->sendData(v[q], 0, q);
			} 
			EndDump
		}
		
		//Print Fermi velocity and related statistics:
		if(ShouldDump(FermiVelocity))
		{	StartDump("FermiVelocity")
			if(eInfo.fillingsUpdate == ElecInfo::FillingsHsub)
			{	FILE* fp;
				const int nMu = 7;
				const double dmu[nMu] = {0., -eInfo.smearingWidth, +eInfo.smearingWidth,
					-2.*eInfo.smearingWidth, +2.*eInfo.smearingWidth,
					-3.*eInfo.smearingWidth, +3.*eInfo.smearingWidth };
				std::vector<double> gEf(nMu); //density of states per unit volume at Fermi level with each dmu
				std::vector<vector3<>> vFabsSum(nMu); //g(Ef)<|v_i|> for each Cartesian direction i at Fermi level with each dmu
				std::vector<matrix3<>> vFsum(nMu), vFsqSum(nMu); //g(Ef)<(1/v)v.vT> and g(Ef)<v.vT> at Fermi level with each dmu
				double BzUnused = 0.;
				double mu = (!std::isnan(eInfo.mu))
					? eInfo.mu
					: eInfo.findMu(e->eVars.Hsub_eigs, eInfo.nElectrons, BzUnused);
				//Get Cartesian symmetry matrices:
				std::vector<matrix3<>> symCart;
				for(const SpaceGroupOp& op: e->symm.getMatrices())
					symCart.push_back(e->gInfo.R * op.rot * e->gInfo.invR);
				//Compute sum:
				for(int q=eInfo.qStart; q<eInfo.qStop; q++)
				{	double dosPrefac = eInfo.qnums[q].weight/e->gInfo.detR;
					for(int b=0; b<eInfo.nBands; b++)
					{	const double& Eqb = eVars.Hsub_eigs[q][b];
						const vector3<>& vqb = v[q][b];
						double vSq = vqb.length_squared();
						matrix3<> vv = outer(vqb, vqb);
						for(int iMu=0; iMu<nMu; iMu++)
						{	double w = dosPrefac * (-eInfo.smearPrime(mu+dmu[iMu], Eqb));
							gEf[iMu] += w;
							//Landauer version based on |v| along each Cartesian axis:
							for(const matrix3<>& rot: symCart)
							{	vector3<> vRot = rot * vqb;
								for(int iDir=0; iDir<3; iDir++)
									vFabsSum[iMu][iDir] += w * fabs(vRot[iDir]);
							}
							//Tensorial versions:
							if(vSq > symmThresholdSq)
							{	vFsum[iMu] += (w/sqrt(vSq)) * vv;
								vFsqSum[iMu] += w * vv;
							}
						}
					}
				}
				mpiWorld->allReduceData(gEf, MPIUtil::ReduceSum);
				mpiWorld->allReduceData(vFabsSum, MPIUtil::ReduceSum);
				mpiWorld->allReduceData(vFsum, MPIUtil::ReduceSum);
				mpiWorld->allReduceData(vFsqSum, MPIUtil::ReduceSum);
				for(vector3<>& m: vFabsSum) m *= (0.5/symCart.size()); //all rotated versions accumulated above (0.5 from Landauer formula)
				for(matrix3<>& m: vFsum) e->symm.symmetrize(m);
				for(matrix3<>& m: vFsqSum) e->symm.symmetrize(m);
				//Write from head:
				if(mpiWorld->isHead())
				{	fp = fopen(fname.c_str(), "w");
					if(!fp) die_alone("Error opening %s for writing.\n", fname.c_str());
					for(int iMu=0; iMu<nMu; iMu++)
					{	fprintf(fp, "Results for dmu = %7lg Eh:\n", dmu[iMu]);
						fprintf(fp, "-----------------------------\n");
						fprintf(fp, "\nvF [Eh-a0]: %12lg\n", sqrt(trace(vFsqSum[iMu])/gEf[iMu]));
						fprintf(fp, "\ng(Ef) [1/(Eh-a0^3)]: %12lg\n", gEf[iMu]);
						fprintf(fp, "\ng|v|(Ef) [1/a0^2]: "); vFabsSum[iMu].print(fp, "%12lg ");
						fprintf(fp, "\ngv(Ef) [1/a0^2]: %12lg\n", trace(vFsum[iMu])/3); vFsum[iMu].print(fp, "%12lg ", true, 1e-12);
						fprintf(fp, "\ngvv(Ef) [Eh/a0]: %12lg\n", trace(vFsqSum[iMu])/3); vFsqSum[iMu].print(fp, "%12lg ", true, 1e-12);
						fprintf(fp, "\n\n");
					}
					fprintf(fp, "\nUnit conversions:\n");
					fprintf(fp, "1/a0^2 = %lg 1/(Ohm nm^2)\n", Ohm*std::pow(10.*Angstrom,2));
					fprintf(fp, "Eh/a0 = %lg 1/(Ohm nm fs)\n", Ohm*(10.*Angstrom)*fs);
					fclose(fp);
				}
				EndDump
			}
			else logPrintf("skipping (no smearing)\n");
		}
	}
	
	if(ShouldDump(RealSpaceWfns))
	{	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{	int nSpinor = eVars.C[q].spinorLength();
			for(int b=0; b<eInfo.nBands; b++) for(int s=0; s<nSpinor; s++)
			{	ostringstream prefixStream;
				prefixStream << "wfns_" << q << '_' << (b*nSpinor+s) << ".rs";
				StartDump(prefixStream.str())
				saveRawBinary(I(eVars.C[q].getColumn(b,s)), fname.c_str());
				EndDump
			}
		}
	}

	if(ShouldDump(ExcCompare))
	{	StartDump("ExcCompare") logPrintf("\n");
		FILE* fp = fopen(fname.c_str(), "w");
		if(!fp) die("Error opening %s for writing.\n", fname.c_str());
		//Compute exact exchange if required:
		std::map<double,double> unscaledEXX; //exact exchange energies for each range parameter
		if(e->exCorr.exxFactor()) //use the main functional Exx when available
			unscaledEXX[e->exCorr.exxRange()] = e->ener.E["EXX"] / e->exCorr.exxFactor();
		for(auto exc: e->exCorrDiff)
			if(exc->exxFactor() && unscaledEXX.find(exc->exxRange())==unscaledEXX.end())
			{	double omega = exc->exxRange();
				//Exact exchange not yet computed for this range parameter:
				logPrintf("\tComputing exact exchange with range-parameter %lg  ... ", omega);
				unscaledEXX[omega] = (*e->exx)(1., omega, eVars.F, eVars.C);
				logPrintf("EXX = %.16lf (unscaled)\n", unscaledEXX[omega]); logFlush();
			}
		//KE density for meta-GGAs, if required
		ScalarFieldArray tau;
		bool needTau = false;
		for(auto exc: e->exCorrDiff)
			needTau |= exc->needsKEdensity();
		if(needTau)
		{	if(e->exCorr.needsKEdensity()) //Used by main functional i.e. already computed in ElecVars
				tau = eVars.tau;
			else //Not used by main functional, need to compute:
				tau = eVars.KEdensity();
		}
		//Print main functional energy:
		fprintf(fp, "%s(%s) = %.16lf\n", relevantFreeEnergyName(*e),
				e->exCorr.getName().c_str(), relevantFreeEnergy(*e));
		//Compute all the functionals in the list:
		for(auto exc: e->exCorrDiff)
		{	double Enew = relevantFreeEnergy(*e) - (e->ener.E["Exc"] + e->ener.E["EXX"] + e->ener.E["Exc_core"])
				+ (*exc)(eVars.get_nXC(), 0, false, &tau)
				+ (iInfo.nCore ? -(*exc)(iInfo.nCore, 0, false, &iInfo.tauCore) : 0)
				+ exc->exxFactor() * unscaledEXX[exc->exxRange()];
			logPrintf("\t%s = %.16lf for '%s'\n", relevantFreeEnergyName(*e), Enew, exc->getName().c_str());
			fprintf(fp, "%s(%s) = %.16lf\n", relevantFreeEnergyName(*e), exc->getName().c_str(), Enew);
		}
		fclose(fp);
		logPrintf("\t"); EndDump
	}
	
	if(ShouldDump(FluidDebug) && hasFluid)
		eVars.fluidSolver->dumpDebug(getFilename("fluid%s").c_str());

	if(ShouldDump(DOS))
	{	dos->dump();
	}
	
	if(ShouldDump(Stress))
	{	StartDump("stress")
		if(mpiWorld->isHead())
		{	FILE* fp = fopen(fname.c_str(), "w");
			if(!fp) die("Error opening %s for writing.\n", fname.c_str());
			fprintf(fp, "# Stress tensor [Eh/a0^3]:\n");
			e->iInfo.stress.print(fp, "%12lg ", false);
			fclose(fp);
		}
		EndDump
	}
	
	if(ShouldDump(XCanalysis))
	{	ScalarFieldArray tauW = XC_Analysis::tauWeizsacker(*e); DUMP_spinCollection(tauW, "tauW");
		ScalarFieldArray spness = XC_Analysis::spness(*e); DUMP_spinCollection(spness, "spness");
		ScalarFieldArray sVh = XC_Analysis::sHartree(*e); DUMP_spinCollection(sVh, "sHartree");
	}

	if(ShouldDump(EresolvedDensity))
	{	std::vector<diagMatrix>& F = ((ElecVars&)e->eVars).F;
		std::vector<diagMatrix> Forig = F; //backup fillings
		int iRange=0;
		for(const auto& Erange: densityErange)
		{	//Set fillings based on current energy range:
			for(int q=eInfo.qStart; q<eInfo.qStop; q++)
				for(int b=0; b<eInfo.nBands; b++)
				{	double e_qb = eVars.Hsub_eigs[q][b];
					F[q][b] = (e_qb >= Erange.first && e_qb <= Erange.second) ? 1. : 0.;
				}
			//Calculate and dump density:
			ScalarFieldArray density = eVars.calcDensity();
			ostringstream oss; oss << "EresolvedDensity." << iRange; iRange++;
			DUMP_spinCollection(density, oss.str())
		}
		F = Forig; //restore fillings
	}
	
	if(ShouldDump(FermiDensity)) 
	{	std::vector<diagMatrix>& F = ((ElecVars&)e->eVars).F;
		std::vector<diagMatrix> Forig = F; //backup fillings
		int iRange=0;
		for(const auto& muLevel: fermiDensityLevels)
		{	//Set fillings based on derivative of smearing function evaluated at muLevel
			double muF, Bz;
			//calculate mu if not set
			muF = ( !std::isnan(muLevel) ? muLevel : eInfo.findMu(eVars.Hsub_eigs,eInfo.nElectrons,Bz) );
			for(int q=eInfo.qStart; q<eInfo.qStop; q++)
				for(int b=0; b<eInfo.nBands; b++)
					F[q][b] = -eInfo.smearPrime(muF,eVars.Hsub_eigs[q][b]);
			//Calculate and dump density
			ScalarFieldArray density = eVars.calcDensity();
			ostringstream oss; oss << "FermiDensity." << iRange; iRange++;
			DUMP_spinCollection(density, oss.str())
		}
		F = Forig; //restore fillings
	}
	
	if(freq==DumpFreq_End and (ShouldDump(R) or ShouldDump(L) or ShouldDump(Q) or ShouldDump(Berry)))
	{	if(!dumpCprime) dumpCprime = std::make_shared<DumpCprime>(); //default parameters
		dumpCprime->dump((Everything&)*e);
	}
	
	if(ShouldDump(DWfns))
	{	//Dump wave functions
		if (!e->pertInfo.commensurate)
			logPrintf("WARNING! Incommensurate wfns are not dumped correctly.");
		StartDump("dwfns")
		eInfo.write(e->pertInfo.dC, fname.c_str());
		EndDump
	}

	if(ShouldDump(Dn))
	{	if (e->pertInfo.commensurate)
		{	ScalarFieldArray dn = e->pertInfo.dn;
			if(e->pertInfo.datom) dn = dn + e->pertInfo.datom->dnatom;
			DUMP_spinCollection(dn, "dn")
		}
		else
		{	DUMP_spinCollection(e->pertInfo.dnmq, "dn-q")
			DUMP_spinCollection(e->pertInfo.dnpq, "dn+q")
		}
	}

	if (ShouldDump(DVext) && e->pertInfo.dVext)
	{	if (e->pertInfo.commensurate)
		{	DUMP_spinCollection(e->pertInfo.dVext->dVext, "dVext")
		}
		else
		{	DUMP_spinCollection(e->pertInfo.dVext->dVextmq, "dVext-q")
			DUMP_spinCollection(e->pertInfo.dVext->dVextpq, "dVext+q")
		}
	}
	
	if(ShouldDump(DVscloc))
	{	if (e->pertInfo.commensurate)
		{	ScalarFieldArray dVscloc = e->pertInfo.dVsclocTau + e->pertInfo.dVscloc;
			DUMP_spinCollection(dVscloc, "dVscloc")
		}
		else
		{	DUMP_spinCollection(e->pertInfo.dVsclocmq, "dVscloc-q")
			DUMP_spinCollection(e->pertInfo.dVsclocpq, "dVscloc+q")
		}
	}

	//----------------------------------------------------------------------
	//The following compute-intensive things are free to clear wavefunctions
	//to conserve memory etc. and should therefore happen at the very end
	
	if(freq==DumpFreq_End && ShouldDump(Polarizability))
	{	polarizability->dump(*e);
	}

	if(freq==DumpFreq_End && ShouldDump(ElectronScattering))
	{	electronScattering->dump(*e);
	}
}

bool Dump::checkInterval(DumpFrequency freq, int iter) const
{	//Check if dumping should occur at this iteration
	auto intervalFreq = interval.find(freq);
	return (intervalFreq==interval.end() //cound not find interval
		|| ((iter+1) % intervalFreq->second == 0)); //or iteration number is divisible by it
}

string Dump::getFilename(string varName) const
{	//Create a map of substitutions:
	std::map<string,string> subMap;
	ostringstream ossIter; ossIter << curIter;
	subMap["$VAR"] = varName;
	subMap["$ITER"] = ossIter.str();
	subMap["$STAMP"] = stamp;
	subMap["$INPUT"] = inputBasename;
	//Find the relevant pattern:
	string fname = format;
	auto iter = formatFreq.find(curFreq);
	if(iter != formatFreq.end())
		fname = iter->second; //override with frequency-dependent format
	//Apply the substitutions:
	while(true)
	{	bool found = false;
		for(auto sub: subMap)
		{	size_t pos = fname.find(sub.first);
			if(pos != string::npos)
			{	fname.replace(pos, sub.first.length(), sub.second);
				found = true;
			}
		}
		if(!found) break; //no matches in this pass
	}
	return fname;
}
