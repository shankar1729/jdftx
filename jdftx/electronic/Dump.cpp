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
#include <electronic/matrix.h>
#include <electronic/SpeciesInfo.h>
#include <electronic/operators.h>
#include <electronic/ExactExchange.h>
#include <core/DataMultiplet.h>
#include <sstream>

void dumpSpinorbit(const Everything& e);

void Dump::setup(const Everything& everything)
{	e = &everything;
	if(wannier.group.size())
		wannier.setup(everything);
}


void Dump::operator()(DumpFrequency freq)
{
	const ElecInfo &eInfo = e->eInfo;
	const ElecVars &eVars = e->eVars;
	const IonInfo &iInfo = e->iInfo;

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
	

	#define DUMP(object, prefix, varname) \
		if(ShouldDump(varname)) \
		{	StartDump(prefix) \
			saveRawBinary(object, fname.c_str()); \
			EndDump \
		}

	// Set up date/time stamp
	time_t timenow = time(0);
	tm *mytm = localtime(&timenow);
	ostringstream stampStream;
	stampStream << (mytm->tm_mon+1) << '.' << mytm->tm_mday << '.'
		<< mytm->tm_hour << '-' << mytm->tm_min << '-' << mytm->tm_sec;
	stamp = stampStream.str();
	
	if(ShouldDump(State))
	{
		//Dump wave functions
		StartDump("wfns")
		write(eVars.C, fname.c_str());
		EndDump
		
		if(eVars.fluidSolver)
		{	//Dump state of fluid:
			StartDump("fluidState")
			eVars.fluidSolver->saveState(fname.c_str());
			EndDump
		}
		
		if(eInfo.fillingsUpdate!=ElecInfo::ConstantFillings)
		{	//Dump fillings
			StartDump("fillings")
			FILE* fp = fopen(fname.c_str(), "w");
			if(!fp) die("Error opening %s for writing.\n", fname.c_str());
			eInfo.printFillings(fp);
			fclose(fp);
			EndDump
			
			if(eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux)
			{	//Dump auxilliary hamiltonian
				StartDump("Haux")
				FILE* fp = fopen(fname.c_str(), "w");
				if(!fp) die("Error opening %s for writing.\n", fname.c_str());
				write(eVars.B, fname.c_str());
				fclose(fp);
				EndDump
			}
		}
	}

	if(ShouldDump(IonicPositions) || (ShouldDump(State) && e->ionicMinParams.nIterations>0))
	{	StartDump("ionpos")
		FILE* fp = fopen(fname.c_str(), "w");
		if(!fp) die("Error opening %s for writing.\n", fname.c_str());
		iInfo.printPositions(fp);
		fclose(fp);
		EndDump
	}
	if(ShouldDump(Forces))
	{	StartDump("force")
		FILE* fp = fopen(fname.c_str(), "w");
		if(!fp) die("Error opening %s for writing.\n", fname.c_str());
		iInfo.forces.print(*e, fp);
		fclose(fp);
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
	if(eVars.fluidType != FluidNone)
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
		write(eVars.Hsub_evecs, fname.c_str());
		EndDump
	}

	if(ShouldDump(BandEigs))
	{	StartDump("eigenvals")
		FILE* fp = fopen(fname.c_str(), "w");
		if(!fp) die("Error opening %s for writing.\n", fname.c_str());
		for(int q=0; q < eInfo.nStates; q++)
		{	for (int b=0; b < eInfo.nBands; b++)
				fprintf(fp, "%18.10le ", eVars.Hsub_eigs[q][b]);
			fprintf(fp, "\n");
		}
		fclose(fp);
		EndDump
	}

	if(ShouldDump(Ecomponents))
	{	StartDump("Ecomponents")
		FILE* fp = fopen(fname.c_str(), "w");
		if(!fp) die("Error opening %s for writing.\n", fname.c_str());
		e->ener.print(fp);
		fclose(fp);
		EndDump
	}
	
	if(ShouldDump(BoundCharge) && eVars.fluidType!=FluidNone)
	{	StartDump("nbound")
		DataGptr nboundTilde = (-1.0/(4*M_PI*e->gInfo.detR)) * L(eVars.d_fluid);
		nboundTilde->data()[0] = -(J(eVars.get_nTot())+iInfo.rhoIon)->data()[0]; //total bound charge will neutralize system
		saveRawBinary(I(nboundTilde), fname.c_str());
		EndDump
	}
	
	if(ShouldDump(FluidDensity))
		if(eVars.fluidSolver)
			eVars.fluidSolver->dumpDensities(getFilename("fluid%s").c_str());
	
	if(ShouldDump(QMC)) dumpQMC();

	//----------------- The following are not included in 'All' -----------------------
	
	if(ShouldDumpNoAll(RealSpaceWfns))
	{	for(int q=0; q<eInfo.nStates; q++)
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
			unscaledEXX[e->exCorr.exxRange()] = e->ener.EXX / e->exCorr.exxFactor();
		for(auto exc: e->exCorrDiff)
			if(exc->exxFactor() && unscaledEXX.find(exc->exxRange())==unscaledEXX.end())
			{	double omega = exc->exxRange();
				//Exact exchange not yet computed for this range parameter:
				logPrintf("\tComputing exact exchange with range-parameter %lg  ... ", omega);
				logFlush();
				auto relevantExx = e->exx.find(omega);
				unscaledEXX[omega] = (*relevantExx->second)(eVars.F, eVars.C);
				logPrintf("EXX = %.16lf (unscaled)\n", unscaledEXX[omega]);
			}
		//KE density for meta-GGAs, if required
		DataRptrCollection tau, Vxc, Vtau;
		bool needTau = false;
		for(auto exc: e->exCorrDiff)
			needTau |= exc->needsKEdensity();
		if(needTau)
		{	if(e->exCorr.needsKEdensity()) //Used by main functional i.e. already computed in ElecVars
				tau = eVars.tau;
			else //Not used by main functional, need to compute:
			{	tau = eVars.KEdensity();
				if(iInfo.tauCore)
					for(unsigned s=0; s<tau.size(); s++)
						tau[s] += (1.0/tau.size()) * iInfo.tauCore;
			}
		}
		//Print main functional energy:
		fprintf(fp, "%s(%s) = %.16lf\n", relevantFreeEnergyName(*e),
				e->exCorr.getName().c_str(), relevantFreeEnergy(*e));
		//Compute all the functionals in the list:
		for(auto exc: e->exCorrDiff)
		{	double Enew = relevantFreeEnergy(*e) - (e->ener.Exc + e->ener.EXX + e->ener.Exc_core)
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

	if(ShouldDumpNoAll(SpinOrbit)) dumpSpinorbit(*e);
	
	if(ShouldDumpNoAll(Projectors))
	{	StartDump("proj");
		FILE *fp = fopen(fname.c_str(), "w");
		if(!fp) die("Error opening %s for writing.\n", fname.c_str());
		for(int q=0; q<eInfo.nStates; q++)
			for(auto sp: iInfo.species)
				sp->writeProjectors(q, fp);
		fclose(fp);
		EndDump
	}
	
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
			{	fname.replace(pos, pos + sub.first.length(), sub.second);
				found = true;
			}
		}
		if(!found) break; //no matches in this pass
	}
	return fname;
}

bool Dump::shouldDump(DumpFrequency freq, int iter) const
{	auto intervalFreq = interval.find(freq);
	if(intervalFreq==interval.end()) return true; //no interval => every iteration
	else return ((iter+1) % intervalFreq->second == 0); //every interval[freq] iterations
}


//Dumps 4 matrices, corresponding to the 4 spin subblocks of Hso(s,k,s',k)
void dumpSpinorbit(const Everything& e)
{	const ElecInfo &eInfo = e.eInfo;
	const ElecVars &eVars = e.eVars;
	if(eInfo.spinType != SpinZ) die("Spin-orbit calculation requires spinType==SpinZ.\n");
	DataRptrVec dVsc_up = gradient(eVars.Vscloc[0]);
	DataRptrVec dVsc_dn = gradient(eVars.Vscloc[1]);
	
	matrix Hso(eInfo.nBands, eInfo.nBands, isGpuEnabled());
	string fname = e.dump.getFilename("Hso_%d.%s");
	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
	char filename[1024];
	for(int q_up = 0; q_up<eInfo.nStates/2; q_up++)
	{	int q_dn = q_up + eInfo.nStates/2;
		complex i(0,1);
		// up up
		Hso.zero();
		Hso -= (eVars.C[q_up] ^ Idag_DiagV_I(D(eVars.C[q_up], 1), dVsc_up[0]));
		Hso += (eVars.C[q_up] ^ Idag_DiagV_I(D(eVars.C[q_up], 0), dVsc_up[1]));
		sprintf(filename, fname.c_str(), q_up, "uu");
		Hso.write(filename);
		// up dn
		Hso.zero();
		Hso -= (eVars.C[q_up] ^ Idag_DiagV_I(D(eVars.C[q_dn], 2), dVsc_dn[1]));
		Hso += (eVars.C[q_up] ^ Idag_DiagV_I(D(eVars.C[q_dn], 1), dVsc_dn[2]));
		Hso -= i*(eVars.C[q_up] ^ Idag_DiagV_I(D(eVars.C[q_dn], 2), dVsc_dn[0]));
		Hso += i*(eVars.C[q_up] ^ Idag_DiagV_I(D(eVars.C[q_dn], 0), dVsc_dn[2]));
		sprintf(filename, fname.c_str(), q_up, "ud");
		Hso.write(filename);
		// dn up
		Hso.zero();
		Hso -= (eVars.C[q_dn]^Idag_DiagV_I(D(eVars.C[q_up], 2), dVsc_up[1]));
		Hso += (eVars.C[q_dn]^Idag_DiagV_I(D(eVars.C[q_up], 1), dVsc_up[2]));
		Hso += i*(eVars.C[q_dn]^Idag_DiagV_I(D(eVars.C[q_up], 2), dVsc_up[0]));
		Hso -= i*(eVars.C[q_dn]^Idag_DiagV_I(D(eVars.C[q_up], 0), dVsc_up[2]));
		sprintf(filename, fname.c_str(), q_up, "du");
		Hso.write(filename);
		// dn dn
		Hso.zero();
		Hso += (eVars.C[q_dn] ^ Idag_DiagV_I(D(eVars.C[q_dn], 1), dVsc_dn[0]));
		Hso -= (eVars.C[q_dn] ^ Idag_DiagV_I(D(eVars.C[q_dn], 0), dVsc_dn[1]));
		sprintf(filename, fname.c_str(), q_up, "dd");
		Hso.write(filename);
	}
	EndDump
}
