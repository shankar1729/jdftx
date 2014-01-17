/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver
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

#include <cmath>
#include <cstdio>
#include <ctime>

#include <electronic/Dump.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/SpeciesInfo.h>
#include <electronic/operators.h>
#include <electronic/ExactExchange.h>
#include <electronic/SelfInteractionCorrection.h>
#include <electronic/DOS.h>
#include <electronic/Polarizability.h>
#include <electronic/LatticeMinimizer.h>
#include <fluid/FluidSolver.h>
#include <core/DataMultiplet.h>
#include <core/WignerSeitz.h>
#include <core/DataIO.h>
#include <core/Units.h>
#include <sstream>

void dumpExcitations(const Everything& e, const char* filename);
namespace Moments{void dumpMoment(const Everything& e, const char* filename, int n, vector3<> origin);}

namespace XC_Analysis{
	DataRptrCollection tauWeizsacker(const Everything& e);
	DataRptrCollection spness(const Everything& e);
	DataRptrCollection sHartree(const Everything& e);
}

void Dump::setup(const Everything& everything)
{	e = &everything;
	if(wannier.group.size())
		wannier.setup(everything);
	
	if(dos) dos->setup(everything);
	
	//Add citation for QMC coupling if required (done here so that it works in dry run):
	for(auto dumpPair: *this)
		if(dumpPair.second == DumpQMC)
			Citations::add("Quantum Monte Carlo solvation",
				"K.A. Schwarz, R. Sundararaman, K. Letchworth-Weaver, T.A. Arias and R. Hennig, Phys. Rev. B 85, 201102(R) (2012)");
}


void Dump::operator()(DumpFrequency freq, int iter)
{
	//Check if dumping should occur at this iteration
	auto intervalFreq = interval.find(freq);
	if(intervalFreq!=interval.end() //interval found
		&& ((iter+1) % intervalFreq->second != 0)) //and iteration number is not divisible by it
		return; // => don't dump this time

	const ElecInfo &eInfo = e->eInfo;
	const ElecVars &eVars = e->eVars;
	const IonInfo &iInfo = e->iInfo;

	bool isCevec = (freq==DumpFreq_Gummel) || (freq==DumpFreq_Ionic) || (freq==DumpFreq_End); //whether C has been set to eigenfunctions
	
	//Macro to determine which variables should be dumped:
	//(Drop the "Dump" from the DumpVariables enum name to use as var)
	#define ShouldDump(var) (  ShouldDumpNoAll(var) || count(std::make_pair(freq,DumpAll)) )

	//Determine whether to dump, but don't include in the 'All' collection
	#define ShouldDumpNoAll(var) count(std::make_pair(freq,Dump##var))

	#define StartDump(prefix) \
		string fname = getFilename(prefix); \
		logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();

	#define EndDump \
		logPrintf("done\n"); logFlush();
	

	#define DUMP_nocheck(object, prefix) \
		{	StartDump(prefix) \
			if(mpiUtil->isHead()) saveRawBinary(object, fname.c_str()); \
			EndDump \
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
	
	if((ShouldDump(State) and eInfo.fillingsUpdate!=ElecInfo::ConstantFillings) or ShouldDump(Fillings))
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
		write(eVars.C, fname.c_str(), eInfo);
		EndDump
		
		if(eVars.fluidSolver)
		{	//Dump state of fluid:
			StartDump("fluidState")
			if(mpiUtil->isHead()) eVars.fluidSolver->saveState(fname.c_str());
			EndDump
		}
		
		if(eInfo.fillingsUpdate!=ElecInfo::ConstantFillings and eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux)
		{	//Dump auxilliary hamiltonian
			StartDump("Haux")
			eInfo.write(eVars.B, fname.c_str());
			EndDump
		}
	}

	if(ShouldDump(IonicPositions) || (ShouldDump(State) && e->ionicMinParams.nIterations>0))
	{	StartDump("ionpos")
		FILE* fp = mpiUtil->isHead() ? fopen(fname.c_str(), "w") : nullLog;
		if(!fp) die("Error opening %s for writing.\n", fname.c_str());
		iInfo.printPositions(fp);  //needs to be called from all processes (for magnetic moment computation)
		if(mpiUtil->isHead())fclose(fp);
		EndDump
	}
	if(ShouldDump(Forces))
	{	StartDump("force")
		if(mpiUtil->isHead()) 
		{	FILE* fp = fopen(fname.c_str(), "w");
			if(!fp) die("Error opening %s for writing.\n", fname.c_str());
			iInfo.forces.print(*e, fp);
			fclose(fp);
		}
		EndDump
	}
	if(ShouldDump(Lattice) || (ShouldDump(State) && e->latticeMinParams.nIterations>0))
	{	StartDump("lattice")
		if(mpiUtil->isHead()) 
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
	
	if(eInfo.spinType == SpinZ)
	{	DUMP(eVars.n[0], "n_up", ElecDensity)
		DUMP(eVars.n[1], "n_dn", ElecDensity)
	}
	else DUMP(eVars.n[0], "n", ElecDensity)
	if(iInfo.nCore) DUMP(iInfo.nCore, "nCore", CoreDensity)
	
	DUMP(I(eVars.d_vac), "d_vac", Dvac);
	if(eVars.fluidParams.fluidType != FluidNone)
	{	DUMP(I(eVars.d_fluid), "d_fluid", Dfluid);
		DUMP(I(eVars.d_vac + eVars.d_fluid), "d_tot", Dtot);
		DUMP(I(eVars.V_cavity), "V_cavity", Vcavity);
		DUMP(I(eVars.V_cavity + eVars.d_fluid), "V_fluidTot", VfluidTot);
	}

	DUMP(I(iInfo.Vlocps), "Vlocps", Vlocps)
	if(eInfo.spinType == SpinZ)
	{	DUMP(eVars.Vscloc[0], "Vscloc_up", Vscloc)
		DUMP(eVars.Vscloc[1], "Vscloc_dn", Vscloc)
	}
	else DUMP(eVars.Vscloc[0], "Vscloc", Vscloc)
	
	if(ShouldDump(HsubEvecs))
	{	StartDump("Hsub_evecs")
		eInfo.write(eVars.Hsub_evecs, fname.c_str());
		EndDump
	}
	
	// Dumps tau (positive Kinetic Energy density)
	if((ShouldDump(KEdensity) or (e->exCorr.needsKEdensity() and ShouldDump(ElecDensity))))
	{	const auto& tau = (e->exCorr.needsKEdensity() ? e->eVars.tau : e->eVars.KEdensity());
		if(eInfo.spinType == SpinZ)
		{	DUMP_nocheck(tau[0], "KEdensity_up");
			DUMP_nocheck(tau[1], "KEdensity_dn");
			
			DataRptrCollection tauW = XC_Analysis::tauWeizsacker(*e);
			DUMP_nocheck(tauW[0], "tauW_up");
			DUMP_nocheck(tauW[1], "tauW_dn");
			
			DataRptrCollection spness = XC_Analysis::spness(*e);
			DUMP_nocheck(spness[0], "spness_up");
			DUMP_nocheck(spness[1], "spness_dn");
			
			DataRptrCollection sVh = XC_Analysis::sHartree(*e);
			DUMP_nocheck(sVh[0], "sHartree_up");
			DUMP_nocheck(sVh[1], "sHartree_dn");
			
		}
		else DUMP_nocheck(tau[0], "KEdensity");
	}
	
	if(ShouldDump(BandEigs))
	{	StartDump("eigenvals")
		eInfo.write(eVars.Hsub_eigs, fname.c_str());
		EndDump
	}

	if(eInfo.hasU && ShouldDump(RhoAtom))
	{	StartDump("rhoAtom")
		FILE* fp = mpiUtil->isHead() ? fopen(fname.c_str(), "w") : nullLog;
		iInfo.computeU(eVars.F, eVars.C, 0, 0, fp); //needs to be called from all processes
		if(mpiUtil->isHead()) fclose(fp);
		EndDump
	}
	
	if(ShouldDump(Ecomponents))
	{	StartDump("Ecomponents")
		if(mpiUtil->isHead())
		{	FILE* fp = fopen(fname.c_str(), "w");
			if(!fp) die("Error opening %s for writing.\n", fname.c_str());	
			e->ener.print(fp);
			fclose(fp);
		}
		EndDump
	}
	
	if((ShouldDump(BoundCharge) || ShouldDump(SolvationRadii)) && eVars.fluidParams.fluidType!=FluidNone)
	{	DataGptr nboundTilde = (-1.0/(4*M_PI*e->gInfo.detR)) * L(eVars.d_fluid);
		if(iInfo.ionWidth) nboundTilde = gaussConvolve(nboundTilde, iInfo.ionWidth);
		nboundTilde->setGzero(-(J(eVars.get_nTot())+iInfo.rhoIon)->getGzero()); //total bound charge will neutralize system
		if(ShouldDump(SolvationRadii))
		{	StartDump("Rsol")
			if(mpiUtil->isHead()) dumpRsol(I(nboundTilde), fname);
			EndDump
		}
		DUMP(I(nboundTilde), "nbound", BoundCharge)
	}
	
	if(ShouldDump(FluidDensity))
		if(eVars.fluidSolver)
			eVars.fluidSolver->dumpDensities(getFilename("fluid%s").c_str());
	
	if(ShouldDump(QMC) && isCevec) dumpQMC();
	
	if(ShouldDump(Dipole))
	{	StartDump("Moments")
		if(mpiUtil->isHead()) Moments::dumpMoment(*e, fname.c_str(), 1, vector3<>(0.,0.,0.));
		EndDump
	}
	
	if(ShouldDump(SIC))
	{	SelfInteractionCorrection selfInteractionCorrection(*e);
		selfInteractionCorrection.dump(getFilename("SIC").c_str());
	}

	if(ShouldDump(Symmetries))
	{	StartDump("sym")
		if(mpiUtil->isHead())
		{	FILE* fp = fopen(fname.c_str(), "w");
			const std::vector< matrix3<int> >& sym = e->symm.getMatrices();
			for(const matrix3<int>& m: sym)
			{	m.print(fp, "%2d ", false);
				fprintf(fp, "\n");
			}
			fclose(fp);
		}
		EndDump
	}
	
	if(ShouldDump(Kpoints))
	{	StartDump("kPts")
		if(mpiUtil->isHead())
		{	FILE* fp = fopen(fname.c_str(), "w");
			eInfo.kpointsPrint(fp, true);
			fclose(fp);
		}
		EndDump
		if(e->symm.mode != SymmetriesNone)
		{	StartDump("kMap")
			if(mpiUtil->isHead())
			{	FILE* fp = fopen(fname.c_str(), "w");
				e->symm.printKmap(fp);
				fclose(fp);
			}
			EndDump
		}
	}
	
	if(ShouldDump(OrbitalDep) && e->exCorr.orbitalDep && isCevec)
		e->exCorr.orbitalDep->dump();
	
	//----------------- The following are not included in 'All' -----------------------

	if(ShouldDumpNoAll(Excitations))
	{	StartDump("Excitations")
		dumpExcitations(*e, fname.c_str());
		EndDump
	}
	
	if(ShouldDumpNoAll(Momenta))
	{	StartDump("momenta")
		std::vector<matrix> momenta(eInfo.nStates);
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) //kpoint/spin
		{	momenta[q] = zeroes(eInfo.nBands, eInfo.nBands*3);
			for(int k=0; k<3; k++) //cartesian direction
				momenta[q].set(0,eInfo.nBands, eInfo.nBands*k,eInfo.nBands*(k+1),
					complex(0,e->gInfo.detR) * (eVars.C[q] ^ D(eVars.C[q], k)) );
		}
		eInfo.write(momenta, fname.c_str(), eInfo.nBands, eInfo.nBands*3);
		EndDump
	}
	
	if(ShouldDumpNoAll(RealSpaceWfns))
	{	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
			for(int b=0; b<eInfo.nBands; b++)
			{	ostringstream prefixStream;
				prefixStream << "wfns_" << q << '_' << b << ".rs";
				StartDump(prefixStream.str())
				saveRawBinary(I(eVars.C[q].getColumn(b)), fname.c_str());
				EndDump
			}
	}

	if(ShouldDumpNoAll(ExcCompare))
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
		DataRptrCollection tau;
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
	
	if(ShouldDumpNoAll(FluidDebug))
		if(eVars.fluidSolver)
			eVars.fluidSolver->dumpDebug(getFilename("fluid%s").c_str());

	if(ShouldDumpNoAll(Wannier))
	{	wannier.saveMLWF();
	}
	
	if(ShouldDumpNoAll(OptVext))
	{	if(eInfo.spinType == SpinZ)
		{	DUMP(eVars.Vexternal[0], "optVextUp", OptVext)
			DUMP(eVars.Vexternal[1], "optVextDn", OptVext)
		}
		else DUMP(eVars.Vexternal[0], "optVext", OptVext)
	}
	
	if(ShouldDumpNoAll(DOS))
	{	dos->dump();
	}
	
	if(ShouldDumpNoAll(Stress))
	{	
		StartDump("stress")
		
		//This part needs to happen on all processes (since it calls ElecVars::energyAndGrad())
		logSuspend();
		LatticeMinimizer lattMin(*((Everything*) e));
		auto stress = lattMin.calculateStress();
		matrix3<> stressTensor;
		for(size_t i=0; i<lattMin.strainBasis.size(); i++)
			stressTensor += stress[i]*lattMin.strainBasis[i];
		lattMin.restore();
		logResume();
		
		if(mpiUtil->isHead())
		{	FILE* fp = fopen(fname.c_str(), "w");
			if(!fp) die("Error opening %s for writing.\n", fname.c_str());
			
			// Dump stress in strain basis units
			fprintf(fp, "%zu strain basis elements\n", lattMin.strainBasis.size());
			for(const matrix3<>& s: lattMin.strainBasis)
			{	s.print(fp, " %lg ");
				fprintf(fp, "\n");
			}
			fprintf(fp, "\n\n");
			
			fprintf(fp, "stress (in strain units, magnitudes along directions above)\n");
			for(size_t j=0; j<lattMin.strainBasis.size(); j++)
			{	fprintf(fp, "%.5e \t ", stress[j]);
			}
			fprintf(fp, "\n\n");
			
			// Dump stress in lattice units
			fprintf(fp, "stress (in strain units)");
			for(int j=0; j<3; j++)
			{	fprintf(fp, " \\\n\t");
				for(int k=0; k<3; k++)
					fprintf(fp, "%20.15lf ", stressTensor(j,k));
			}

			fclose(fp);
		}
		EndDump
	}

	//Polarizability dump deletes wavefunctions to free memory and should therefore happen at the very end
	if(freq==DumpFreq_End && ShouldDumpNoAll(Polarizability))
	{	polarizability->dump(*e);
	}
}

string Dump::getFilename(string varName) const
{	//Create a map of substitutions:
	std::map<string,string> subMap;
	subMap["$VAR"] = varName;
	subMap["$STAMP"] = stamp;
	//Apply the substitutions:
	string fname = format;
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

inline void set_rInv(size_t iStart, size_t iStop, const vector3<int>& S, const matrix3<>& RTR, const WignerSeitz* ws, double* rInv)
{	vector3<> invS; for(int k=0; k<3; k++) invS[k] = 1./S[k];
	matrix3<> meshMetric = Diag(invS) * RTR * Diag(invS);
	const double rWidth = 1.;
	THREAD_rLoop( 
		double r = sqrt(meshMetric.metric_length_squared(ws->restrict(iv, S, invS)));
		if(r < rWidth) //Polynomial with f',f" zero at origin and f,f',f" matched at rWidth
		{	double t = r/rWidth;
			rInv[i] = (1./rWidth) * (1. + 0.5*(1.+t*t*t*(-2.+t)) + (1./12)*(1.+t*t*t*(-4.+t*3.))*2. );
		}
		else rInv[i] = 1./r;
	)
}

void Dump::dumpRsol(DataRptr nbound, string fname)
{	
	//Compute normalization factor for the partition:
	int nAtomsTot = 0; for(const auto& sp: e->iInfo.species) nAtomsTot += sp->atpos.size();
	const double nFloor = 1e-5/nAtomsTot; //lower cap on densities to prevent Nyquist noise in low density regions
	DataRptr nAtomicTot;
	for(const auto& sp: e->iInfo.species)
	{	RadialFunctionG nRadial;
		logSuspend(); sp->getAtom_nRadial(0,0, nRadial); logResume();
		for(unsigned atom=0; atom<sp->atpos.size(); atom++)
		{	DataRptr nAtomic = radialFunction(e->gInfo, nRadial, sp->atpos[atom]);
			double nMin, nMax; callPref(eblas_capMinMax)(e->gInfo.nr, nAtomic->dataPref(), nMin, nMax, nFloor);
			nAtomicTot += nAtomic;
		}
	}
	DataRptr nboundByAtomic = (nbound*nbound) * inv(nAtomicTot);

	DataRptr rInv0(DataR::alloc(e->gInfo));
	{	logSuspend(); WignerSeitz ws(e->gInfo.R); logResume();
		threadLaunch(set_rInv, e->gInfo.nr, e->gInfo.S, e->gInfo.RTR, &ws, rInv0->data());
	}
	
	//Compute bound charge 1/r and 1/r^2 expectation values weighted by atom-density partition:
	FILE* fp = fopen(fname.c_str(), "w");
	fprintf(fp, "#Species   rMean +/- rSigma [bohrs]   (rMean +/- rSigma [Angstrom])   sqrt(Int|nbound^2|) in partition\n");
	for(const auto& sp: e->iInfo.species)
	{	RadialFunctionG nRadial;
		logSuspend(); sp->getAtom_nRadial(0,0, nRadial); logResume();
		for(unsigned atom=0; atom<sp->atpos.size(); atom++)
		{	DataRptr w = radialFunction(e->gInfo, nRadial, sp->atpos[atom]) * nboundByAtomic;
			//Get r centered at current atom:
			DataGptr trans; nullToZero(trans, e->gInfo); initTranslation(trans, e->gInfo.R*sp->atpos[atom]);
			DataRptr rInv = I(trans * J(rInv0), true);
			//Compute moments:
			double wNorm = integral(w);
			double rInvMean = integral(w * rInv) / wNorm;
			double rInvSqMean = integral(w * rInv * rInv) / wNorm;
			double rInvSigma = sqrt(rInvSqMean - rInvMean*rInvMean);
			double rMean = 1./rInvMean;
			double rSigma = rInvSigma / (rInvMean*rInvMean);
			//Print stats:
			fprintf(fp, "Rsol %s    %.2lf +/- %.2lf    ( %.2lf +/- %.2lf A )   Qrms: %.1le\n", sp->name.c_str(),
				rMean, rSigma, rMean/Angstrom, rSigma/Angstrom, sqrt(wNorm));
		}
	}
	fclose(fp);
}



namespace Moments{
	inline double map_to_interval(double x)
	{	double temp =  modf(x, new double);
		if(temp>0.5)
			return temp-1.;
		else if(temp<-0.5)
			return temp+1.;
		else
			return temp;
	}
	inline vector3<> map_to_interval(vector3<> x)
	{	for(int j=0;j<3;j++)
			x[j] = map_to_interval(x[j]);
		return x;
	}

	void rn_pow_x(int i, vector3<> r, int dir, matrix3<> R, double moment, vector3<> r0, double* rx)
	{	vector3<> lx = map_to_interval(inv(R)*(r-r0));
		rx[i] = pow((R*lx)[dir], moment);
	}
	
	void dumpMoment(const Everything& e, const char* filename, int moment, vector3<> origin)
	{	const GridInfo& g = e.gInfo;
		
		FILE* fp = fopen(filename, "w");
		if(!fp) die("Error opening %s for writing.\n", filename);	
		
		// Calculates the electronic moment about origin
		DataRptr r0, r1, r2;
		nullToZero(r0, g); 	nullToZero(r1, g); 	nullToZero(r2, g);
		applyFunc_r(g, rn_pow_x, 0, g.R, moment, origin, r0->data());
		applyFunc_r(g, rn_pow_x, 1, g.R, moment, origin, r1->data());
		applyFunc_r(g, rn_pow_x, 2, g.R, moment, origin, r2->data());
		vector3<> elecMoment;
		elecMoment[0] = integral(e.eVars.n[0]*r0);
		elecMoment[1] = integral(e.eVars.n[0]*r1);
		elecMoment[2] = integral(e.eVars.n[0]*r2);
		fprintf(fp, "Electron moment of order %i: %f\t%f\t%f", moment, elecMoment[0], elecMoment[1], elecMoment[2]);
		
		// Calculates the ionic moment about the origin
		vector3<> ionMoment(0., 0., 0.);
		for(auto sp: e.iInfo.species)
			for(unsigned n=0; n < sp->atpos.size(); n++)
			{	vector3<> cartesianCoord = g.R * map_to_interval(sp->atpos[n]);
				for(int j=0; j<3; j++)
					ionMoment[j] += -sp->Z * pow(cartesianCoord[j], moment); 
			}
		fprintf(fp, "\nIon moment of order %i: %f\t%f\t%f", moment, ionMoment[0], ionMoment[1], ionMoment[2]);

		
		// Calculates the total (elec+ion) dipole moment
		fprintf(fp, "\nTotal moment of order %i: %f\t%f\t%f", moment, ionMoment[0]+elecMoment[0], ionMoment[1]+elecMoment[1], ionMoment[2]+elecMoment[2]);		
	}
}

void dumpExcitations(const Everything& e, const char* filename)
{
	const GridInfo& g = e.gInfo;

	struct excitation
	{	int q,o,u;
		double dE;
		double dreal, dimag, dnorm;
		
		excitation(int q, int o, int u, double dE, double dreal, double dimag, double dnorm): q(q), o(o), u(u), dE(dE), dreal(dreal), dimag(dimag), dnorm(dnorm){};
		
		inline bool operator<(const excitation& other) const {return dE<other.dE;}
		void print(FILE* fp) const { fprintf(fp, "%5i %3i %3i %12.5e %12.5e %12.5e %12.5e\n", q, o, u, dE, dreal, dimag, dnorm); }
	};
	std::vector<excitation> excitations;

	double maxHOMO=-DBL_MAX, minLUMO=DBL_MAX; // maximum (minimum) of all HOMOs (LUMOs) in all qnums
	int maxHOMOq=0, minLUMOq=0, maxHOMOn=0, minLUMOn=0; //Indices and energies for the indirect gap
	
	//Select relevant eigenvals:
	std::vector<diagMatrix> eigsQP;
	if(e.exCorr.orbitalDep && e.dump.count(std::make_pair(DumpFreq_End, DumpOrbitalDep)))
	{	//Search for an eigenvalsQP file:
		string fname = e.dump.getFilename("eigenvalsQP");
		FILE* fp = fopen(fname.c_str(), "r");
		if(fp)
		{	fclose(fp);
			eigsQP.resize(e.eInfo.nStates);
			e.eInfo.read(eigsQP, fname.c_str());
		}
	}
	const std::vector<diagMatrix>& eigs = eigsQP.size() ? eigsQP : e.eVars.Hsub_eigs;
	
	// Integral kernel's for Fermi's golden rule
	DataRptr r0, r1, r2;
	nullToZero(r0, g); 	nullToZero(r1, g); 	nullToZero(r2, g);
	applyFunc_r(g, Moments::rn_pow_x, 0, g.R, 1, vector3<>(0.,0.,0.), r0->data());
	applyFunc_r(g, Moments::rn_pow_x, 1, g.R, 1, vector3<>(0.,0.,0.), r1->data());
	applyFunc_r(g, Moments::rn_pow_x, 2, g.R, 1, vector3<>(0.,0.,0.), r2->data());
	
	//Find and cache all excitations in system (between same qnums)
	bool insufficientBands = false;
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{	//Find local HOMO and check band sufficiency:
		int HOMO = e.eInfo.findHOMO(q);
		if(HOMO+1>=e.eInfo.nBands) { insufficientBands=true; break; }
		
		//Update global HOMO and LUMO of current process:
		if(eigs[q][HOMO]   > maxHOMO) { maxHOMOq = q; maxHOMOn = HOMO;   maxHOMO = eigs[q][HOMO];   }
		if(eigs[q][HOMO+1] < minLUMO) { minLUMOq = q; minLUMOn = HOMO+1; minLUMO = eigs[q][HOMO+1]; }
		
		for(int o=HOMO; o>=0; o--)
		{	for(int u=(HOMO+1); u<e.eInfo.nBands; u++)
			{	complex x = integral(I(e.eVars.C[q].getColumn(u))*r0*I(e.eVars.C[q].getColumn(o)));
				complex y = integral(I(e.eVars.C[q].getColumn(u))*r1*I(e.eVars.C[q].getColumn(o)));
				complex z = integral(I(e.eVars.C[q].getColumn(u))*r2*I(e.eVars.C[q].getColumn(o)));
				vector3<> dreal(x.real(), y.real(),z.real());
				vector3<> dimag(x.imag(), y.imag(),z.imag());
				vector3<> dnorm(sqrt(x.norm()), sqrt(y.norm()),sqrt(z.norm()));
				double dE = eigs[q][u]-eigs[q][o]; //Excitation energy
				excitations.push_back(excitation(q, o, u, dE, dreal.length_squared(), dimag.length_squared(), dnorm.length_squared()));
			}
		}
	}
	mpiUtil->allReduce(insufficientBands, MPIUtil::ReduceLOr);
	if(insufficientBands)
	{	logPrintf("Insufficient bands to calculate excited states!\n");
		logPrintf("Increase the number of bands (elec-n-bands) and try again!\n");
		return;
	}
	
	//Transmit results to head process:
	if(mpiUtil->isHead())
	{	excitations.reserve(excitations.size() * mpiUtil->nProcesses());
		for(int jProcess=1; jProcess<mpiUtil->nProcesses(); jProcess++)
		{	//Receive data:
			size_t nExcitations; mpiUtil->recv(nExcitations, jProcess, 0);
			std::vector<int> msgInt(4 + nExcitations*3); 
			std::vector<double> msgDbl(2 + nExcitations*4);
			mpiUtil->recv(msgInt.data(), msgInt.size(), jProcess, 1);
			mpiUtil->recv(msgDbl.data(), msgDbl.size(), jProcess, 2);
			//Unpack:
			std::vector<int>::const_iterator intPtr = msgInt.begin();
			std::vector<double>::const_iterator dblPtr = msgDbl.begin();
			//--- globals:
			int j_maxHOMOq = *(intPtr++); int j_maxHOMOn = *(intPtr++); double j_maxHOMO = *(dblPtr++);
			int j_minLUMOq = *(intPtr++); int j_minLUMOn = *(intPtr++); double j_minLUMO = *(dblPtr++);
			if(j_maxHOMO > maxHOMO) { maxHOMOq=j_maxHOMOq; maxHOMOn=j_maxHOMOn; maxHOMO=j_maxHOMO; }
			if(j_minLUMO < minLUMO) { minLUMOq=j_minLUMOq; minLUMOn=j_minLUMOn; minLUMO=j_minLUMO; }
			//--- excitation array:
			for(size_t iExcitation=0; iExcitation<nExcitations; iExcitation++)
			{	int q = *(intPtr++); int o = *(intPtr++); int u = *(intPtr++);
				double dE = *(dblPtr++);
				double dreal = *(dblPtr++); double dimag = *(dblPtr++); double dnorm = *(dblPtr++);
				excitations.push_back(excitation(q, o, u, dE, dreal, dimag, dnorm));
			}
		}
	}
	else
	{	//Pack data:
		std::vector<int> msgInt; std::vector<double> msgDbl;
		size_t nExcitations = excitations.size();
		msgInt.reserve(4 + nExcitations*3);
		msgDbl.reserve(2 + nExcitations*4);
		msgInt.push_back(maxHOMOq); msgInt.push_back(maxHOMOn); msgDbl.push_back(maxHOMO);
		msgInt.push_back(minLUMOq); msgInt.push_back(minLUMOn); msgDbl.push_back(minLUMO);
		for(const excitation& e: excitations)
		{	msgInt.push_back(e.q); msgInt.push_back(e.o); msgInt.push_back(e.u);
			msgDbl.push_back(e.dE);
			msgDbl.push_back(e.dreal); msgDbl.push_back(e.dimag); msgDbl.push_back(e.dnorm);
		}
		//Send data:
		mpiUtil->send(nExcitations, 0, 0);
		mpiUtil->send(msgInt.data(), msgInt.size(), 0, 1);
		mpiUtil->send(msgDbl.data(), msgDbl.size(), 0, 2);
	}

	//Process and print excitations:
	if(!mpiUtil->isHead()) return;
	
	FILE* fp = fopen(filename, "w");
	if(!fp) die("Error opening %s for writing.\n", filename);
	
	std::sort(excitations.begin(), excitations.end());
	const excitation& opt = excitations.front();
	fprintf(fp, "Using %s eigenvalues.\n", eigsQP.size() ? "discontinuity-corrected QP" : "KS");
	fprintf(fp, "Optical (direct) gap: %.5e (from n = %i to %i in qnum = %i)\n", opt.dE, opt.o, opt.u, opt.q);
	fprintf(fp, "Indirect gap: %.5e (from (%i, %i) to (%i, %i))\n\n", minLUMO-maxHOMO, maxHOMOq, maxHOMOn, minLUMOq, minLUMOn);
	
	fprintf(fp, "Optical excitation energies and corresponding electric dipole transition strengths\n");
	fprintf(fp, "qnum   i   f      dE        |<psi1|r|psi2>|^2 (real, imag, norm)\n");
	for(const excitation& e: excitations) e.print(fp);
	fclose(fp);
}


namespace XC_Analysis{

	#define cutoff 1e-8
	
	void regularize(int i, vector3<> r, double* tau)
	{	tau[i] = (tau[i]<cutoff ? 0. : tau[i]);
	}
	void spness_kernel(int i, vector3<> r, double* tauW, double* tau, double* spness)
	{	spness[i] = (tau[i]>cutoff ? tauW[i]/tau[i] : 1.);
	}
	
	DataRptrCollection tauWeizsacker(const Everything& e)
	{	DataRptrCollection tauW(e.eVars.n.size());
		for(size_t j=0; j<e.eVars.n.size(); j++)
			tauW[j] = (1./8.)*lengthSquared(gradient(e.eVars.n[j]))*pow(e.eVars.n[j], -1);
		return tauW;
	}
	
	DataRptrCollection spness(const Everything& e)
	{
			const auto& tau = (e.exCorr.needsKEdensity() ? e.eVars.tau : e.eVars.KEdensity());
		
			DataRptrCollection tauW = tauWeizsacker(e);
			DataRptrCollection spness(e.eVars.n.size()); nullToZero(spness, e.gInfo);
			for(size_t j=0; j<e.eVars.n.size(); j++)
			{	applyFunc_r(e.gInfo, regularize, tau[j]->data());
				applyFunc_r(e.gInfo, regularize, tauW[j]->data());
				//spness[j] = tauW[j]*pow(tau[j], -1);
				applyFunc_r(e.gInfo, spness_kernel, tauW[j]->data(), tau[j]->data(), spness[j]->data());
			}

			return spness;
	}
	
	DataRptrCollection sHartree(const Everything& e)
	{
		DataRptrCollection sVh(e.eVars.n.size());
		DataGptrCollection nTilde = J(e.eVars.n);
		for(size_t s=0; s<e.eVars.n.size(); s++)
			sVh[s] = I((*(e.coulomb))(nTilde[s]));
		
		return sVh;
	}
}