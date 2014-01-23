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
#include <electronic/Vibrations.h>
#include <core/Units.h>

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
	DumpFillings, "Fillings",
	DumpRhoAtom, "RhoAtom",
	DumpEcomponents, "Ecomponents",
	DumpExcCompare, "ExcCompare",
	DumpBoundCharge, "BoundCharge",
	DumpSolvationRadii, "SolvationRadii",
	DumpQMC, "QMC",
	DumpRealSpaceWfns, "RealSpaceWfns",
	DumpFluidDebug, "FluidDebug",
	DumpOptVext, "optVext",
	DumpDOS, "DOS",
	DumpSIC, "SelfInteractionCorrection",
	DumpDipole, "Dipole",
	DumpStress, "Stress",
	DumpExcitations, "Excitations",
	DumpMomenta, "Momenta",
	DumpSymmetries, "Symmetries",
	DumpKpoints, "Kpoints",
	DumpOrbitalDep, "OrbitalDep"
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
	DumpFillings,       "Fillings",
	DumpRhoAtom,        "Atomic-orbital projected density matrices (only for species with +U enabled)",
	DumpEcomponents,    "Components of the energy",
	DumpBoundCharge,    "Bound charge in the fluid",
	DumpSolvationRadii, "Effective solvation radii based on fluid bound charge distribution",
	DumpQMC,            "Blip'd orbitals and potential for CASINO",
	DumpRealSpaceWfns,  "Real-space wavefunctions (one column per file) [not in All]",
	DumpExcCompare,     "Energies for other exchange-correlation functionals (see elec-ex-corr-compare) [not in All]",
	DumpFluidDebug,     "Fluid specific debug output if any  [not in All]",
	DumpOptVext,        "Optimized external potentials (see invertKohnSham) [not in All]",
	DumpDOS,            "Density of States (see density-of-states) [not in All]",
	DumpSIC,            "Calculates Perdew-Zunger self-interaction corrected Kohn-Sham eigenvalues",
	DumpDipole,         "Dipole moment of explicit charges (ionic and electronic)",
	DumpStress,         "Dumps dE/dR_ij where R_ij is the i'th component of the j'th lattice vector",
	DumpExcitations,    "Dumps dipole moments and transition strength (electric-dipole) of excitations",
	DumpMomenta,        "Momentum matrix elements in a binary file (indices outer to inner: state, cartesian direction, band1, band2)",
	DumpSymmetries,     "List of symmetry matrices (in covariant lattice coordinates)",
	DumpKpoints,        "List of reduced k-points in calculation, and mapping to the unreduced k-point mesh",
	DumpOrbitalDep,     "Custom output from orbital-dependent functionals (eg. quasi-particle energies, discontinuity potential)"
);

struct CommandDump : public Command
{
	CommandDump() : Command("dump")
	{
		format = "<freq> <var> <var> ...";
		comments =
			"<freq> is one of:"
			+ addDescriptions(freqMap.optionList(), linkDescription(freqMap, freqDescMap))
			+ "\nand each <var*> is one of:"
			+ addDescriptions(varMap.optionList(), linkDescription(varMap, varDescMap))
			+ "\nList of dumped variables from multiple instances will be accumulated for each <freq>."
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
			"Dump every <interval> iterations of type <freq>=Ionic|Electronic|Fluid|Gummel\n"
			"Without this command, the behavior defaults to <interval>=1 for each <freq>.";
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
			"  Control the filename pattern for dump output:\n"
			"    <format> is an arbitrary format string that will be substituted according to:\n"
			"       $VAR -> name of the variable being dumped (this must be present somewhere in the string)\n"
			"       $STAMP -> time-stamp at the start of dump";
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
		comments = "Output polarizability matrix in specified eigeneigenBasis";
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

struct CommandPolarizabilityKdiff : public Command
{
    CommandPolarizabilityKdiff() : Command("polarizability-kdiff")
	{
		format = "<dk0> <dk1> <dk2> [<dkFilenamePattern>]";
		comments = "Select k-point difference (in reciprocal lattice coords) for polarizability output.\n"
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
			"  dr <dr>: perturbation amplitude in bohrs for force matrix calculation (default: 0.01).\n"
			"  centralDiff yes|no: use a central difference formula for the second derivative\n"
			"     to achieve higher accuracy at twice the cost (default: no)\n"
			"  useConstraints yes|no: restrict modes of motion as specified by move flags\n"
			"     and constraints in the ion command (default: no)\n"
			"  traslationSym yes|no: whether to assume overall translation symmetry (default yes).\n"
			"     Can be turned off to get vibrational levels in an external potential.\n"
			"  rotationSym yes|no: project out rotational modes (default no). Improves reliability for\n"
			"     molecular calculations. Valid only for geometries with an unambiguous center of mass.\n"
			"  omegaMin <omegaMin>: frequency cutoff (in Eh) for free energy calculation (default: 2e-4)\n"
			"  T <T>: temperature (in Kelvin) for free energy calculation (default: 298)\n"
			"  omegaResolution <omegaResolution>: resolution for detecting and reporting degeneracies\n"
			"     in modes (default: 1e-4). Does not affect free energies and all modes are still printed.";
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
