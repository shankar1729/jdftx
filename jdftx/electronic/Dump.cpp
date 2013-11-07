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
#include <core/DataIO.h>
#include <sstream>

void dumpExcitations(const Everything& e, const char* filename);
namespace Moments{void dumpMoment(const Everything& e, const char* filename, int n, vector3<> origin);}

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
	
	if((ShouldDump(State) and eInfo.fillingsUpdate!=ElecInfo::ConstantFillings) or ShouldDump(Fillings))
	{	//Dump fillings
		StartDump("fillings")
		FILE* fp = fopen(fname.c_str(), "w");
		if(!fp) die("Error opening %s for writing.\n", fname.c_str());
		eInfo.printFillings(fp);
		fclose(fp);
		EndDump
	}
	
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
		
		if(eInfo.fillingsUpdate!=ElecInfo::ConstantFillings and eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux)
		{	//Dump auxilliary hamiltonian
			StartDump("Haux")
			FILE* fp = fopen(fname.c_str(), "w");
			if(!fp) die("Error opening %s for writing.\n", fname.c_str());
			write(eVars.B, fname.c_str());
			fclose(fp);
			EndDump
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
	if(ShouldDump(Lattice) || (ShouldDump(State) && e->latticeMinParams.nIterations>0))
	{	StartDump("lattice")
		FILE* fp = fopen(fname.c_str(), "w");
		if(!fp) die("Error opening %s for writing.\n", fname.c_str());
		fprintf(fp, "lattice");
		for(int j=0; j<3; j++)
		{	fprintf(fp, " \\\n\t");
			for(int k=0; k<3; k++)
				fprintf(fp, "%20.15lf ", e->gInfo.R(j,k));
		}
		fprintf(fp, "#Note: latt-scale has been absorbed into these lattice vectors.\n");
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
		write(eVars.Hsub_evecs, fname.c_str());
		EndDump
	}
	
	// Dumps tau (positive Kinetic Energy density)
	if(ShouldDump(KEdensity) or (e->exCorr.needsKEdensity() and ShouldDump(ElecDensity)))
	{	const auto& tau = (e->exCorr.needsKEdensity() ? e->eVars.tau : e->eVars.KEdensity());
		if(eInfo.spinType == SpinZ)
		{	{	StartDump("KEdensity_up");
				saveRawBinary(tau[0], fname.c_str());
				EndDump;
			}
			{	StartDump("KEdensity_dn");
				saveRawBinary(tau[1], fname.c_str());
				EndDump;
			}
		}
		else
		{	StartDump("KEdensity");
			saveRawBinary(tau[0], fname.c_str());
			EndDump;
		}
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

	if(eInfo.hasU && ShouldDump(RhoAtom))
	{	StartDump("rhoAtom")
		FILE* fp = fopen(fname.c_str(), "w");
		iInfo.computeU(eVars.F, eVars.C, 0, 0, fp);
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
	
	if(ShouldDump(BoundCharge) && eVars.fluidParams.fluidType!=FluidNone)
	{	StartDump("nbound")
		DataGptr nboundTilde = (-1.0/(4*M_PI*e->gInfo.detR)) * L(eVars.d_fluid);
		nboundTilde->data()[0] = -(J(eVars.get_nTot())+iInfo.rhoIon)->data()[0]; //total bound charge will neutralize system
		saveRawBinary(I(nboundTilde), fname.c_str());
		EndDump
	}
	
	if(ShouldDump(FluidDensity))
		if(eVars.fluidSolver)
			eVars.fluidSolver->dumpDensities(getFilename("fluid%s").c_str());
	
	if(ShouldDump(QMC) && isCevec) dumpQMC();
	
	if(ShouldDump(Dipole))
	{	StartDump("Moments")
		Moments::dumpMoment(*e, fname.c_str(), 1, vector3<>(0.,0.,0.));
		EndDump
	}
	
	if(ShouldDump(SIC))
	{	SelfInteractionCorrection selfInteractionCorrection(*e);
		selfInteractionCorrection.dump(getFilename("SIC").c_str());
	}

	//----------------- The following are not included in 'All' -----------------------

	if(ShouldDumpNoAll(Excitations))
	{	StartDump("Excitations")
		dumpExcitations(*e, fname.c_str());
		EndDump
	}
	
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
		// Cache the type of calculation (fixed-n)
		bool fixed_n = e->cntrl.fixed_n;
		((Everything*) e)->cntrl.fixed_n = false;
		
		logSuspend();
		LatticeMinimizer lattMin(*((Everything*) e));
		auto stress = lattMin.calculateStress();
		matrix3<> stressTensor;
		for(size_t i=0; i<lattMin.strainBasis.size(); i++)
			stressTensor += stress[i]*lattMin.strainBasis[i];
		lattMin.restore();
		logResume();
		
		StartDump("stress")
		FILE* fp = fopen(fname.c_str(), "w");
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
		
		// Reset fixed_n variable
		((Everything*) e)->cntrl.fixed_n = fixed_n;
		
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
	
		FILE* fp = fopen(filename, "w");
		if(!fp) die("Error opening %s for writing.\n", filename);
	
		struct excitation
		{	int q,o,u;
			double dE;
			double dreal, dimag, dnorm;
			
			excitation(int q, int o, int u, double dE, double dreal, double dimag, double dnorm): q(q), o(o), u(u), dE(dE), dreal(dreal), dimag(dimag), dnorm(dnorm){};
			
			inline bool operator<(const excitation& other) const {return dE<other.dE;}
		};
		std::vector<excitation> excitations;

		double maxHOMO, minLUMO; // maximum (minimum) of all HOMOs (LUMOs) in all qnums
		double opticalGap; // optical gap of the system (minimum energy excitation for all qnums

		// Indices and energies for the indirect gap
		int maxHOMOq = 0;
		int minLUMOq = 0;
		int maxHOMOn = e.eInfo.findHOMO(maxHOMOq);
		if(maxHOMOn >= (e.eInfo.nBands-1))
		{	fprintf(fp, "Insufficient bands to calculate excited states!\n");
			fprintf(fp, "Increase the number of bands (elec-n-bands) and try again!\n");
			return;
		}
		int minLUMOn = maxHOMOn+1;
		maxHOMO = e.eVars.Hsub_eigs[maxHOMOq][maxHOMOn];
		minLUMO = e.eVars.Hsub_eigs[minLUMOq][minLUMOn];
		
		// Indices and energies for the direct (optical) gap
		int opticalq = maxHOMOq, opticaln = maxHOMOn;
		opticalGap = minLUMO - maxHOMO;
		
		// Integral kernel's for Fermi's golden rule
		DataRptr r0, r1, r2;
		nullToZero(r0, g); 	nullToZero(r1, g); 	nullToZero(r2, g);
		applyFunc_r(g, Moments::rn_pow_x, 0, g.R, 1, vector3<>(0.,0.,0.), r0->data());
		applyFunc_r(g, Moments::rn_pow_x, 1, g.R, 1, vector3<>(0.,0.,0.), r1->data());
		applyFunc_r(g, Moments::rn_pow_x, 2, g.R, 1, vector3<>(0.,0.,0.), r2->data());
		
		// Find and cache all excitations in system (between same qnums)
		for(int q=0; q<e.eInfo.nStates; q++)
		{	
			int HOMO = e.eInfo.findHOMO(q);
					
			// Update globa HOMO and LUMO
			if(e.eVars.Hsub_eigs[q][HOMO] > maxHOMO)
			{	maxHOMOq = q;
				maxHOMOn = HOMO;
				maxHOMO = e.eVars.Hsub_eigs[maxHOMOq][maxHOMOn];
			}
			if(e.eVars.Hsub_eigs[q][HOMO+1] < minLUMO)	
			{	minLUMOq = q;
				minLUMOn = HOMO+1;
				minLUMO = e.eVars.Hsub_eigs[minLUMOq][minLUMOn];
			}
			
			for(int o=HOMO; o>=0; o--)
			{	for(int u=(HOMO+1); u<e.eInfo.nBands; u++)
				{	complex x = integral(I(e.eVars.C[q].getColumn(u))*r0*I(e.eVars.C[q].getColumn(o)));
					complex y = integral(I(e.eVars.C[q].getColumn(u))*r1*I(e.eVars.C[q].getColumn(o)));
					complex z = integral(I(e.eVars.C[q].getColumn(u))*r2*I(e.eVars.C[q].getColumn(o)));
					vector3<> dreal(x.real(), y.real(),z.real());
					vector3<> dimag(x.imag(), y.imag(),z.imag());
					vector3<> dnorm(sqrt(x.norm()), sqrt(y.norm()),sqrt(z.norm()));
					
					// Calculate the excitation energy, update optical gap if necessary
					double dE = e.eVars.Hsub_eigs[q][u]-e.eVars.Hsub_eigs[q][o];
					if(dE < opticalGap) 
					{	opticalq = q; opticaln = HOMO;
						opticalGap = dE;
					}
					
					excitations.push_back(excitation(q, o, u, dE,dreal.length_squared(), dimag.length_squared(), dnorm.length_squared()));
				}
			}
		}
		std::sort(excitations.begin(), excitations.end());

		fprintf(fp, "Optical (direct) gap: %.5e (from n = %i to %i in qnum = %i)\n", opticalGap, opticaln, opticaln+1, opticalq);
		fprintf(fp, "Indirect gap: %.5e (from (%i, %i) to (%i, %i))\n\n", minLUMO-maxHOMO, maxHOMOq, maxHOMOn, minLUMOq, minLUMOn);
		
		fprintf(fp, "Optical excitation energies and corresponding electric dipole transition strengths\n");
		fprintf(fp, "qnum,\tinitial,\tfinal,\tdE,\t|<psi1|r|psi2>|^2 (real, imag, norm)\n");
		for(size_t i=0; i<excitations.size(); i++)
			fprintf(fp, "%i \t %i \t %i \t %.5e \t %.5e \t %.5e \t %.5e\n", excitations[i].q, excitations[i].o, excitations[i].u, excitations[i].dE, excitations[i].dreal, excitations[i].dimag, excitations[i].dnorm);
}
