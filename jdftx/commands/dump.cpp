/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

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

#include <commands/command.h>
#include <electronic/Everything.h>
#include <electronic/Polarizability.h>
#include <electronic/ElectronScattering.h>
#include <electronic/Vibrations.h>
#include <electronic/Dump_internal.h>
#include <core/Units.h>

struct CommandDumpOnly : public Command
{
	CommandDumpOnly() : Command("dump-only")
	{
		comments = 
			"Bypass all minimization, perform a single energy evaluation at\n"
			"the initial state, and process dump commands at the end state.\n"
			"Useful for post-processing a converged calculation.";
		
		forbid("fix-electron-potential");
		forbid("fix-electron-density");
	}

	void process(ParamList& pl, Everything& e)
	{	e.cntrl.dumpOnly = true;
	}

	void printStatus(Everything& e, int iRep)
	{
	}
}
commandDumpOnly;


EnumStringMap<DumpFrequency> freqMap
(	DumpFreq_End, "End",
	DumpFreq_Lattice, "Lattice",
	DumpFreq_Ionic, "Ionic",
	DumpFreq_Gummel, "Gummel",
	DumpFreq_Fluid, "Fluid",
	DumpFreq_Electronic, "Electronic"
);
EnumStringMap<DumpFrequency> freqDescMap
(	DumpFreq_End, "Dump specified vars at the end of the calculation",
	DumpFreq_Lattice, "Dump specified vars every (few) lattice minimization step(s)",
	DumpFreq_Ionic, "Dump specified vars every (few) ionic step(s)",
	DumpFreq_Gummel, "Dump specified vars every (few) fluid+electron minimize of the gummel loop",
	DumpFreq_Fluid, "Dump specified vars every (few) fluid step(s)",
	DumpFreq_Electronic, "Dump specified vars every (few) electronic step(s)"
);


EnumStringMap<DumpVariable> varMap
(	DumpAll, "All",
	DumpNone, "None",
	DumpState, "State",
	DumpIonicPositions, "IonicPositions",
	DumpForces, "Forces",
	DumpLattice, "Lattice",
	DumpIonicDensity, "IonicDensity",
	DumpElecDensity, "ElecDensity",
	DumpCoreDensity, "CoreDensity",
	DumpKEdensity, "KEdensity",
	DumpFluidDensity, "FluidDensity",
	DumpDvac, "Dvac",
	DumpDfluid, "Dfluid",
	DumpDtot, "Dtot",
	DumpVcavity, "Vcavity",
	DumpVfluidTot, "VfluidTot",
	DumpVlocps, "Vlocps",
	DumpVscloc, "Vscloc",
	DumpBandEigs, "BandEigs",
	DumpEigStats, "EigStats",
	DumpFillings, "Fillings",
	DumpRhoAtom, "RhoAtom",
	DumpEcomponents, "Ecomponents",
	DumpExcCompare, "ExcCompare",
	DumpBoundCharge, "BoundCharge",
	DumpSolvationRadii, "SolvationRadii",
	DumpQMC, "QMC",
       	DumpOcean, "Ocean",
	DumpRealSpaceWfns, "RealSpaceWfns",
	DumpFluidDebug, "FluidDebug",
	DumpSlabEpsilon, "SlabEpsilon",
	DumpChargedDefect, "ChargedDefect",
	DumpOptVext, "optVext",
	DumpDOS, "DOS",
	DumpSIC, "SelfInteractionCorrection",
	DumpDipole, "Dipole",
	DumpStress, "Stress",
	DumpExcitations, "Excitations",
	DumpMomenta, "Momenta",
	DumpSymmetries, "Symmetries",
	DumpKpoints, "Kpoints",
	DumpGvectors, "Gvectors",
	DumpOrbitalDep, "OrbitalDep",
	DumpXCanalysis, "XCanalysis",
	DumpEresolvedDensity, "EresolvedDensity",
	DumpFermiDensity, "FermiDensity" 
);
EnumStringMap<DumpVariable> varDescMap
(	DumpAll,            "Dump most things (except those marked not in All)",
	DumpNone,           "Dump nothing",
	DumpState,          "All variables needed to restart calculation: wavefunction and fluid state/fillings if any",
	DumpIonicPositions, "Ionic positions in the same format (and coordinate system) as the input file",
	DumpForces,         "Forces on the ions in the coordinate system selected by command forces-output-coords",
	DumpLattice,        "Lattice vectors in the same format as the input file",
	DumpIonicDensity,   "Nuclear charge density (with gaussians)",
	DumpElecDensity,    "Electronic densities (n or nup,ndn)",
	DumpCoreDensity,    "Total core electron density (from partial core corrections)",
	DumpKEdensity,      "Kinetic energy density of the valence electrons",
	DumpFluidDensity,   "Fluid densities (NO,NH,nWater for explicit fluids, cavity function for PCMs)",
	DumpDvac,           "Electrostatic potential due to explicit system alone",
	DumpDfluid,         "Electrostatic potential due to fluid alone",
	DumpDtot,           "Total electrostatic potential",
	DumpVcavity,        "Fluid cavitation potential on the electron density that determines the cavity",
	DumpVfluidTot,      "Total contribution of fluid to the electron potential",
	DumpVlocps,         "Local part of pseudopotentials",
	DumpVscloc,         "Self-consistent potential",
	DumpBandEigs,       "Band Eigenvalues",
	DumpEigStats,       "Band eigenvalue statistics: HOMO, LUMO, min, max and Fermi level",
	DumpFillings,       "Fillings",
	DumpRhoAtom,        "Atomic-orbital projected density matrices (only for species with +U enabled)",
	DumpEcomponents,    "Components of the energy",
	DumpBoundCharge,    "Bound charge in the fluid",
	DumpSolvationRadii, "Effective solvation radii based on fluid bound charge distribution",
	DumpQMC,            "Blip'd orbitals and potential for CASINO",
	DumpOcean,           "Wave functions for Ocean code",
	DumpRealSpaceWfns,  "Real-space wavefunctions (one column per file) [not in All]",
	DumpExcCompare,     "Energies for other exchange-correlation functionals (see command elec-ex-corr-compare) [not in All]",
	DumpFluidDebug,     "Fluid specific debug output if any  [not in All]",
	DumpSlabEpsilon,    "Local dielectric function of a slab (see command slab-epsilon)  [not in All]",
	DumpChargedDefect,  "Calculate energy correction for charged defect (see command charged-defect)  [not in All]",
	DumpOptVext,        "Optimized external potentials (see command invertKohnSham) [not in All]",
	DumpDOS,            "Density of States (see command density-of-states) [not in All]",
	DumpSIC,            "Calculates Perdew-Zunger self-interaction corrected Kohn-Sham eigenvalues",
	DumpDipole,         "Dipole moment of explicit charges (ionic and electronic)",
	DumpStress,         "Dumps dE/dR_ij where R_ij is the i'th component of the j'th lattice vector",
	DumpExcitations,    "Dumps dipole moments and transition strength (electric-dipole) of excitations",
	DumpMomenta,        "Momentum matrix elements in a binary file (indices outer to inner: state, cartesian direction, band1, band2)",
	DumpSymmetries,     "List of symmetry matrices (in covariant lattice coordinates)",
	DumpKpoints,        "List of reduced k-points in calculation, and mapping to the unreduced k-point mesh",
	DumpGvectors,       "List of G vectors in reciprocal lattice basis, for each k-point",
	DumpOrbitalDep,     "Custom output from orbital-dependent functionals (eg. quasi-particle energies, discontinuity potential)",
	DumpXCanalysis,     "Debug VW KE density, single-particle-ness and spin-polarzied Hartree potential",
	DumpEresolvedDensity, "Electron density from bands within specified energy ranges",
	DumpFermiDensity,	"Electron density from fermi-derivative at specified energy" 
);

struct CommandDump : public Command
{
	CommandDump() : Command("dump")
	{
		format = "<freq> <var> <var> ...";
		comments =
			"<freq> is one of:"
			+ addDescriptions(freqMap.optionList(), linkDescription(freqMap, freqDescMap))
			+ "\n\nand each <var> is one of:"
			+ addDescriptions(varMap.optionList(), linkDescription(varMap, varDescMap))
			+ "\n\nList of dumped variables from multiple instances will be accumulated for each <freq>."
			"\nUse command dump-interval to dump at regular intervals instead of every iteration.";
		allowMultiple = true;
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	DumpFrequency freq;
		pl.get(freq, DumpFreq_Delim, freqMap, "freq");
		//Handle the default call
		if(freq==DumpFreq_Delim)
		{	e.dump.insert(std::make_pair(DumpFreq_End,DumpState));
			return;
		}
		//For any real dump frequency:
		while(true)
		{	DumpVariable var;
			pl.get(var, DumpDelim, varMap, "var");
			if(var==DumpDelim) break; //will happen at end of command line
			e.dump.insert(std::make_pair(freq,var));
		}
	}

	void printStatus(Everything& e, int iRep)
	{	//Coealesce dump outputs for each frequency into a single command
		std::multimap<DumpFrequency,DumpVariable> dumpMap(e.dump.begin(), e.dump.end());
		typedef std::multimap<DumpFrequency,DumpVariable>::iterator Iter;
		Iter i = dumpMap.begin();
		int iDump=0;
		while(i!=dumpMap.end())
		{	if(iDump==iRep) logPrintf("%s", freqMap.getString(i->first));
			std::pair<Iter,Iter> range = dumpMap.equal_range(i->first);
			for(i=range.first; i!=range.second; i++)
				if(iDump==iRep) logPrintf(" %s", varMap.getString(i->second));
			iDump++;
		}
	}
}
commandDump;


struct CommandDumpInterval : public Command
{
	CommandDumpInterval() : Command("dump-interval")
	{
		format = "<freq> <interval>";
		comments = 
			"Dump every <interval> iterations of type <freq>=Ionic|Electronic|Fluid|Gummel.\n"
			"Without this command, the behavior defaults to <interval>=1 for each <freq>.\n"
			"(See command dump for more details)";
		allowMultiple = true;
	}

	void process(ParamList& pl, Everything& e)
	{	//get the frequency:
		DumpFrequency freq;
		pl.get(freq, DumpFreq_Delim, freqMap, "freq", true);
		if(freq==DumpFreq_End)
			throw string("<freq> must be one of Ionic|Electronic|Fluid|Gummel");
		if(e.dump.interval.find(freq) != e.dump.interval.end())
			throw string("dump-interval has been specified multiple times for <freq>=") + freqMap.getString(freq);
		//get the interval:
		int interval;
		pl.get(interval, 1, "interval", true);
		if(interval<1)
			throw string("<interval> must be a positive integer");
		//Set the interval
		e.dump.interval[freq] = interval;
	}

	void printStatus(Everything& e, int iRep)
	{	std::map<DumpFrequency,int>::const_iterator iter = e.dump.interval.begin();
		for(int i=0; i<iRep; i++) iter++; //access the iRep'th entry of interval
		logPrintf("%s %d", freqMap.getString(iter->first), iter->second);
	}
}
commandDumpInterval;


struct CommandDumpName : public Command
{
	CommandDumpName() : Command("dump-name")
	{
		format = "<format>";
		comments = 
			"Control the filename pattern for dump output, where <format> is an\n"
			"arbitrary format string that will be substituted according to:\n"
			"+ $VAR   -> name of the variable being dumped (this must be present)\n"
			"+ $STAMP -> time-stamp at the start of dump";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.dump.format, string("$STAMP.$VAR"), "format");
		if(e.dump.format.find("$VAR")==string::npos)
			throw "<format> = " + e.dump.format + " doesn't contain the pattern $VAR";
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s", e.dump.format.c_str());
	}
}
commandDumpName;


EnumStringMap<Polarizability::EigenBasis> polarizabilityMap
(	Polarizability::NonInteracting, "NonInteracting",
	Polarizability::External, "External",
	Polarizability::Total, "Total"
);

struct CommandPolarizability : public Command
{
    CommandPolarizability() : Command("polarizability")
	{
		format = "<eigenBasis>=" + polarizabilityMap.optionList() + " [<Ecut>=0] [<nEigs>=0]";
		comments = "Output polarizability matrix in specified eigenBasis";
		
		forbid("electron-scattering"); //both are major operations that are given permission to destroy Everything if necessary
	}
	
	void process(ParamList& pl, Everything& e)
	{	e.dump.polarizability = std::make_shared<Polarizability>();
		pl.get(e.dump.polarizability->eigenBasis, Polarizability::NonInteracting, polarizabilityMap, "eigenBasis");
		pl.get(e.dump.polarizability->Ecut, 0., "Ecut");
		pl.get(e.dump.polarizability->nEigs, 0, "nEigs");
		e.dump.insert(std::make_pair(DumpFreq_End, DumpPolarizability));
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s %lg %d", polarizabilityMap.getString(e.dump.polarizability->eigenBasis),
			e.dump.polarizability->Ecut, e.dump.polarizability->nEigs);
	}
}
commandPolarizability;


enum ElectronScatteringMember
{	ESM_eta,
	ESM_Ecut,
	ESM_fCut,
	ESM_omegaMax,
	ESM_delim
};
EnumStringMap<ElectronScatteringMember> esmMap
(	ESM_eta, "eta",
	ESM_Ecut, "Ecut",
	ESM_fCut, "fCut",
	ESM_omegaMax, "omegaMax"
);

struct CommandElectronScattering : public Command
{
    CommandElectronScattering() : Command("electron-scattering")
	{
		format = "<key1> <value1> ...";
		comments = "Calculate electron-electron scattering rates (expensive!)\n"
			"and output contribution to imaginary part of electron self-energy\n"
			"(calculated effectively using full-frequency G0W0).\n"
			"\n"
			"The following key-value pairs can appear in any order:\n"
			"\n+ eta <eta>\n\n"
			"   <eta> in Eh specifies frequency resolution and half the imaginary part\n"
			"   ascribed to probe frequency (if zero, the electron temperature is used.)\n"
			"\n+ Ecut <Ecut>\n\n"
			"   <Ecut> in Eh specifies energy cut-off for dielectric matrices.\n"
			"   (If zero, the wavefunction cutoff from elec-cutoff is used instead.)\n"
			"\n+ fCut <fCut>\n\n"
			"   <fCut> specifies threshold for considering states fully occupied or\n"
			"   unoccupied in optimizing sums over states (default: 1e-6)\n"
			"\n+ omegaMax <omegaMax>\n\n"
			"   <omegaMax> in Eh is the maximum energy transfer to account for\n"
			"   and hence the maximum frequency in dielectric function frequency grid.\n"
			"   (if zero, autodetermine from available eigenvalues)";
		
		forbid("polarizability"); //both are major operations that are given permission to destroy Everything if necessary
	}
	
	void process(ParamList& pl, Everything& e)
	{	e.dump.electronScattering = std::make_shared<ElectronScattering>();
		e.dump.insert(std::make_pair(DumpFreq_End, DumpElectronScattering));
		ElectronScattering& es = *(e.dump.electronScattering);
		while(true)
		{	ElectronScatteringMember key;
			pl.get(key, ESM_delim, esmMap, "key");
			switch(key)
			{	case ESM_eta: pl.get(es.eta, 0., "eta", true); break;
				case ESM_Ecut: pl.get(es.Ecut, 0., "Ecut", true); break;
				case ESM_fCut: pl.get(es.fCut, 0., "fCut", true); break;
				case ESM_omegaMax: pl.get(es.omegaMax, 0., "omegaMax", true); break;
				case ESM_delim: return; //end of input
			}
		}
	}

	void printStatus(Everything& e, int iRep)
	{	const ElectronScattering& es = *(e.dump.electronScattering);
		logPrintf(" \\\n\teta      %lg", es.eta);
		logPrintf(" \\\n\tEcut     %lg", es.Ecut);
		logPrintf(" \\\n\tfCut     %lg", es.fCut);
		logPrintf(" \\\n\tomegaMax %lg", es.omegaMax);
	}
}
commandElectronScattering;


struct CommandPolarizabilityKdiff : public Command
{
    CommandPolarizabilityKdiff() : Command("polarizability-kdiff")
	{
		format = "<dk0> <dk1> <dk2> [<dkFilenamePattern>]";
		comments = "Select k-point difference (in reciprocal lattice coords) for polarizability output.\n"
			"\n"
			"<dkFilenamePattern> may be specified to read offset band structure calcualations when <dk>\n"
			"does not belong to the k-point mesh. This string should be a filename pattern containing\n"
			"$VAR (to be replaced by eigenvals and wfns) and $q (to be replaced by state index)";
		
		require("polarizability");
	}
	
	void process(ParamList& pl, Everything& e)
	{	pl.get(e.dump.polarizability->dk[0], 0., "dk0", true);
		pl.get(e.dump.polarizability->dk[1], 0., "dk1", true);
		pl.get(e.dump.polarizability->dk[2], 0., "dk2", true);
		//Optional filename for offset states:
		string& dkFilenamePattern = e.dump.polarizability->dkFilenamePattern;
		pl.get(dkFilenamePattern, string(), "dkFilenamePattern");
		if(dkFilenamePattern.length())
		{	if(dkFilenamePattern.find("$VAR")==string::npos) throw "<dkFilenamePattern> = " + dkFilenamePattern + " doesn't contain the pattern $VAR";
			if(dkFilenamePattern.find("$q")==string::npos) throw "<dkFilenamePattern> = " + dkFilenamePattern + " doesn't contain the pattern $q";
		}
	}

	void printStatus(Everything& e, int iRep)
	{	for(int i=0; i<3; i++) logPrintf("%lg ", e.dump.polarizability->dk[i]);
		logPrintf("%s", e.dump.polarizability->dkFilenamePattern.c_str());
	}
}
commandPolarizabilityKdiff;


struct CommandDumpEresolvedDensity : public Command
{
    CommandDumpEresolvedDensity() : Command("dump-Eresolved-density")
	{
		format = "<Emin> <Emax>";
		comments =
			"Output electron density from bands within a specified energy range\n"
			"[Emin,Emax] (in Eh).\n"
			"\n"
			"When issued multiple times, the outputs will be\n"
			"numbered sequenetially EresolvedDensity.0 etc.\n"
			"This automatically invokes dump at End; dumping at\n"
			"other frequencies may be requested using the dump command.";
		
		allowMultiple = true;
	}
	
	void process(ParamList& pl, Everything& e)
	{	double Emin, Emax;
		pl.get(Emin, 0., "Emin", true);
		pl.get(Emax, 0., "Emax", true);
		if(Emin >= Emax) throw string("Emin must be < Emax");
		e.dump.densityErange.push_back(std::make_pair(Emin, Emax));
		e.dump.insert(std::make_pair(DumpFreq_End, DumpEresolvedDensity));
	}
	
	void printStatus(Everything& e, int iRep)
	{	const auto& Erange = e.dump.densityErange[iRep];
		logPrintf("%lg %lg", Erange.first, Erange.second);
	}
}
commandDumpEresolvedDensity;

struct CommandDumpFermiDensity : public Command 
{
    CommandDumpFermiDensity() : Command("dump-fermi-density")
	{
		format = "[<muLevel>]";
		comments =
			"Output electron density calculated from derivative of Fermi-Dirac\n"
			"function evaluated at desired chemical potential <muLevel>.\n"
			"If unspecified, calculated chemical potential will be used.\n"
			"\n"
			"When issued multiple times, the outputs will be\n"
			"numbered sequenetially FermiDensity.0 etc.\n"
			"This automatically invokes dump at End; dumping at\n"
			"other frequencies may be requested using the dump command.";
		
		allowMultiple = true;
		require("elec-fermi-fillings"); //require a temperature
	}
	
	void process(ParamList& pl, Everything& e)
	{	double muLevel;
		pl.get(muLevel, std::numeric_limits<double>::quiet_NaN(), "muLevel", false);
		e.dump.fermiDensityLevels.push_back(muLevel);
		e.dump.insert(std::make_pair(DumpFreq_End, DumpFermiDensity));
	}
	
	void printStatus(Everything& e, int iRep)
	{	const auto& muLevel = e.dump.fermiDensityLevels[iRep];
		logPrintf("%lg", muLevel);
	}
}
commandDumpFermiDensity;


//---------- Vibrations ------------------------

//An enum corresponding to members of class Vibrations
enum VibrationsMember
{	VM_dr,
	VM_centralDiff,
	VM_useConstraints,
	VM_translationSym,
	VM_rotationSym,
	VM_omegaMin,
	VM_T,
	VM_omegaResolution,
	VM_Delim
};

EnumStringMap<VibrationsMember> vibMap
(	VM_dr, "dr",
	VM_centralDiff, "centralDiff",
	VM_useConstraints, "useConstraints",
	VM_translationSym, "translationSym",
	VM_rotationSym, "rotationSym",
	VM_omegaMin, "omegaMin",
	VM_T, "T",
	VM_omegaResolution, "omegaResolution"
);

struct CommandVibrations : public Command
{
    CommandVibrations() : Command("vibrations")
	{
		format = "<key1> <args1> ...";
		comments =
			"Calculate vibrational modes of the system using a finite difference method.\n"
			"Note that this command should typically be issued in a run with converged ionic\n"
			"positions; ionic (and lattice) minimization are bypassed by the vibrations module.\n"
			"Any number of the following subcommands and their arguments may follow:\n"
			"+ dr <dr>: perturbation amplitude in bohrs for force matrix calculation (default: 0.01).\n"
			"+ centralDiff yes|no: use a central difference formula for the second derivative\n"
			"   to achieve higher accuracy at twice the cost (default: no)\n"
			"+ useConstraints yes|no: restrict modes of motion as specified by move flags\n"
			"   and constraints in the ion command (default: no)\n"
			"+ traslationSym yes|no: whether to assume overall translation symmetry (default yes).\n"
			"   Can be turned off to get vibrational levels in an external potential.\n"
			"+ rotationSym yes|no: project out rotational modes (default no). Improves reliability for\n"
			"   molecular calculations. Valid only for geometries with an unambiguous center of mass.\n"
			"+ omegaMin <omegaMin>: frequency cutoff (in Eh) for free energy calculation (default: 2e-4)\n"
			"+ T <T>: temperature (in Kelvin) for free energy calculation (default: 298)\n"
			"+ omegaResolution <omegaResolution>: resolution for detecting and reporting degeneracies\n"
			"   in modes (default: 1e-4). Does not affect free energies and all modes are still printed.";
		forbid("fix-electron-density");
		forbid("fix-electron-potential");
	}
	
	void process(ParamList& pl, Everything& e)
	{	e.vibrations = std::make_shared<Vibrations>();
		while(true)
		{	VibrationsMember key;
			pl.get(key, VM_Delim, vibMap, "key");
			switch(key)
			{	case VM_dr: pl.get(e.vibrations->dr, 0.01, "dr", true); break;
				case VM_centralDiff: pl.get(e.vibrations->centralDiff, false, boolMap, "centralDiff", true); break;
				case VM_useConstraints: pl.get(e.vibrations->useConstraints, false, boolMap, "useConstraints", true); break;
				case VM_translationSym: pl.get(e.vibrations->translationSym, true, boolMap, "translationSym", true); break;
				case VM_rotationSym: pl.get(e.vibrations->rotationSym, false, boolMap, "rotationSym", true); break;
				case VM_omegaMin: pl.get(e.vibrations->omegaMin, 2e-4, "omegaMin", true); break;
				case VM_T: pl.get(e.vibrations->T, 298., "T", true); e.vibrations->T *= Kelvin; break;
				case VM_omegaResolution: pl.get(e.vibrations->omegaResolution, 1e-4, "omegaResolution", true); break;
				case VM_Delim: return; //end of input
			}
		}
		
	}
	
	void printStatus(Everything& e, int iRep)
	{	logPrintf("\\\n\tdr %g", e.vibrations->dr);
		logPrintf("\\\n\tcentralDiff %s", boolMap.getString(e.vibrations->centralDiff));
		logPrintf("\\\n\tuseConstraints %s", boolMap.getString(e.vibrations->useConstraints));
		logPrintf("\\\n\ttranslationSym %s", boolMap.getString(e.vibrations->translationSym));
		logPrintf("\\\n\trotationSym %s", boolMap.getString(e.vibrations->rotationSym));
		logPrintf("\\\n\tomegaMin %g", e.vibrations->omegaMin);
		logPrintf("\\\n\tT %g", e.vibrations->T/Kelvin);
		logPrintf("\\\n\tomegaResolution %g", e.vibrations->omegaResolution);
	}
}
commandVibrations;

//-------------------------------------------------------------------------------------------------

struct CommandSlabEpsilon : public Command
{
	CommandSlabEpsilon() : Command("slab-epsilon")
	{
		format = "<DtotFile> <sigma> [<Ex>=0] [<Ey>=0] [<Ez>=0]";
		comments = 
			"Calculate dielectric function of a slab given the electrostatic potential\n"
			"output from another calculation on same system with a different electric field.\n"
			"+ <DtotFile> contains the electrostatic potential from the other calculation\n"
			"+ <sigma> in bohrs specifies gaussian smoothing in output\n"
			"+ optional <Ex>,<Ey>,Ez> specify the electric-field applied\n"
			"  in the calculation that generated <DtotFile>.";
		
		require("coulomb-truncation-embed");
	}

	void process(ParamList& pl, Everything& e)
	{	if(e.coulombParams.geometry != CoulombParams::Slab)
			throw string("coulomb-interaction must be in Slab mode");
		e.dump.slabEpsilon = std::make_shared<SlabEpsilon>();
		SlabEpsilon& se = *(e.dump.slabEpsilon);
		pl.get(se.dtotFname, string(), "DtotFile", true);
		pl.get(se.sigma, 0., "sigma", true);
		pl.get(se.Efield[0], 0., "Ex");
		pl.get(se.Efield[1], 0., "Ey");
		pl.get(se.Efield[2], 0., "Ez");
		e.dump.insert(std::make_pair(DumpFreq_End, DumpSlabEpsilon)); //dump at end by default
	}

	void printStatus(Everything& e, int iRep)
	{	const SlabEpsilon& se = *(e.dump.slabEpsilon);
		logPrintf("%s %lg %lg %lg %lg", se.dtotFname.c_str(), se.sigma, se.Efield[0], se.Efield[1], se.Efield[2]);
	}
}
commandSlabEpsilon;

//-------------------------------------------------------------------------------------------------

struct CommandChargedDefectCorrection : public Command
{
	CommandChargedDefectCorrection() : Command("charged-defect-correction")
	{
		format = "<DtotFile> <bulkEps>|<slabEpsFile> <rMin> <rSigma>";
		comments = 
			"Calculate energy correction for bulk or surface charged defects.\n"
			"The correction is calculated assuming the defects to be model\n"
			"charges specified using command charged-defect.\n"
			"\n"
			"<DtotFile> contains the electrostatic potential from a reference\n"
			"neutral calculation with similar geometry (lattice vectors and grid\n"
			"must match exactly).\n"
			"\n"
			"For bulk defect calculations (coulomb-interaction Periodic),\n"
			"<bulkEps> is the bulk dielectric constant to use in the correction.\n"
			"\n"
			"For surface defect calculations (coulomb-interaction Slab ...)\n"
			"<slabEpsFile> specifies a dielectric profile calculated using command\n"
			"slab-epsilon in a similar geometry (the number of points along the slab\n"
			"normal direction must match exactly). Note that coulomb-truncation-embed\n"
			"must be specified for charged-defect correction in Slab geometry.\n"
			"\n"
			"<rMin> specifies the distance away from the defect center to use in\n"
			"the determination of the alignment potential, with rSigma specifying an\n"
			"error function turn-on distance. The code wil generate a text file with\n"
			"the spherically averaged model and DFT electrostatic potentials, which\n"
			"can be used to check the calculated alignment and refine rMin and rSigma.";
		
		require("latt-scale");
		require("coords-type");
		require("coulomb-interaction");
	}

	void process(ParamList& pl, Everything& e)
	{	e.dump.chargedDefect = std::make_shared<ChargedDefect>();
		ChargedDefect& cd = *(e.dump.chargedDefect);
		//Ref potential:
		pl.get(cd.dtotFname, string(), "DtotFile", true);
		//Dielectric function:
		switch(e.coulombParams.geometry)
		{	case CoulombParams::Periodic: pl.get(cd.bulkEps, 1., "bulkEps", true); break;
			case CoulombParams::Slab: pl.get(cd.slabEpsFname, string(), "slabEpsFile", true); break;
			default: throw string("coulomb-interaction must be either Slab or Periodic");
		}
		//Alignment potential ranges:
		pl.get(cd.rMin, 0., "rMin", true);
		pl.get(cd.rSigma, 0., "rSigma", true);
		e.dump.insert(std::make_pair(DumpFreq_End, DumpChargedDefect)); //dump at end by default
	}

	void printStatus(Everything& e, int iRep)
	{	const ChargedDefect& cd = *(e.dump.chargedDefect);
		logPrintf("%s ", cd.dtotFname.c_str());
		switch(e.coulombParams.geometry)
		{	case CoulombParams::Periodic: logPrintf("%lg", cd.bulkEps); break;
			case CoulombParams::Slab: logPrintf("%s", cd.slabEpsFname.c_str()); break;
			default:; //should never be encountered
		}
		logPrintf(" %lg %lg", cd.rMin, cd.rSigma);
	}
}
CommandChargedDefectCorrection;

struct CommandChargedDefect : public Command
{
	CommandChargedDefect() : Command("charged-defect")
	{
		format = "<x0> <x1> <x2> <q> <sigma>";
		comments = 
			"Specify model charge distribution(s) for charged-defect-correction.\n"
			"Issue the command once for each defect in the calculation cell.\n"
			"The defect charge density will be modeled by a Gaussian of width\n"
			"<sigma> with net electron count <q> located at <x0>,<x1>,<x2>\n"
			"(in the coordinate system specified by coords-type). The Gaussian must\n"
			"be resolvable on the plane-wave grid; recommend rSigma >= 1 bohr.";
		
		require("charged-defect-correction");
		allowMultiple = true;
	}

	void process(ParamList& pl, Everything& e)
	{	ChargedDefect::Center cdc;
		//Position:
		vector3<> pos;
		pl.get(pos[0], 0., "x0", true);
		pl.get(pos[1], 0., "x1", true);
		pl.get(pos[2], 0., "x2", true);
		cdc.pos = (e.iInfo.coordsType==CoordsLattice) ? pos : inv(e.gInfo.R)*pos;
		//Charge and width:
		pl.get(cdc.q, 0., "q", true);
		pl.get(cdc.sigma, 0., "sigma", true);
		e.dump.chargedDefect->center.push_back(cdc);
	}

	void printStatus(Everything& e, int iRep)
	{	const ChargedDefect::Center& cdc = e.dump.chargedDefect->center[iRep];
		vector3<> pos = (e.iInfo.coordsType==CoordsLattice) ? cdc.pos : e.gInfo.R*cdc.pos;
		logPrintf("%lg %lg %lg  %+lg %lg", pos[0], pos[1], pos[2], cdc.q, cdc.sigma);
	}
}
commandChargedDefect;
