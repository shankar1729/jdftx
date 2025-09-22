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
#include <electronic/ColumnBundle.h>
#include <electronic/Polarizability.h>
#include <electronic/ElectronScattering.h>
#include <electronic/Vibrations.h>
#include <electronic/Dump_internal.h>
#include <electronic/DumpBGW_internal.h>
#include <core/Units.h>

struct CommandDumpOnly : public Command
{
	CommandDumpOnly() : Command("dump-only", "jdftx/Output")
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

//Variables allowed at dump frequency Init
EnumStringMap<DumpVariable> varInitMap
(	DumpNone, "None",
	DumpIonicPositions, "IonicPositions",
	DumpLattice, "Lattice",
	DumpCoreDensity, "CoreDensity",
	DumpVlocps, "Vlocps",
	DumpSymmetries, "Symmetries",
	DumpKpoints, "Kpoints",
	DumpGvectors, "Gvectors"
);

EnumStringMap<DumpFrequency> freqMap
(	DumpFreq_End, "End",
	DumpFreq_Init, "Init",
	DumpFreq_Ionic, "Ionic",
	DumpFreq_Gummel, "Gummel",
	DumpFreq_Fluid, "Fluid",
	DumpFreq_Electronic, "Electronic"
);
EnumStringMap<DumpFrequency> freqDescMap
(	DumpFreq_End, "Dump specified vars at the end of the calculation",
	DumpFreq_Init, "Dump specified vars from " + varInitMap.optionList() + " after initialization (even in dry run)",
	DumpFreq_Ionic, "Dump specified vars every (few) ionic / lattice step(s)",
	DumpFreq_Gummel, "Dump specified vars every (few) fluid+electron minimize of the gummel loop",
	DumpFreq_Fluid, "Dump specified vars every (few) fluid step(s)",
	DumpFreq_Electronic, "Dump specified vars every (few) electronic step(s)"
);


EnumStringMap<DumpVariable> varMap
(	DumpNone, "None",
	DumpState, "State",
	DumpIonicPositions, "IonicPositions",
	DumpForces, "Forces",
	DumpLattice, "Lattice",
	DumpIonicDensity, "IonicDensity",
	DumpElecDensity, "ElecDensity",
	DumpElecDensityAccum, "ElecDensityAccum",
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
	DumpBandProjections, "BandProjections",
	DumpEigStats, "EigStats",
	DumpFillings, "Fillings",
	DumpRhoAtom, "RhoAtom",
	DumpBandUnfold, "BandUnfold",
	DumpEcomponents, "Ecomponents",
	DumpExcCompare, "ExcCompare",
	DumpBoundCharge, "BoundCharge",
	DumpSolvationRadii, "SolvationRadii",
	DumpQMC, "QMC",
	DumpOcean, "Ocean",
	DumpBGW, "BGW",
	DumpRealSpaceWfns, "RealSpaceWfns",
	DumpFluidDebug, "FluidDebug",
	DumpSlabEpsilon, "SlabEpsilon",
	DumpBulkEpsilon, "BulkEpsilon",
	DumpChargedDefect, "ChargedDefect",
	DumpDOS, "DOS",
	DumpSIC, "SelfInteractionCorrection",
	DumpDipole, "Dipole",
	DumpStress, "Stress",
	DumpExcitations, "Excitations",
	DumpFCI, "FCI",
	DumpSpin, "Spin",
	DumpMomenta, "Momenta",
	DumpVelocities, "Velocities",
	DumpFermiVelocity, "FermiVelocity",
	DumpR, "R",
	DumpL, "L",
	DumpQ, "Q",
	DumpBerry, "Berry",
	DumpSymmetries, "Symmetries",
	DumpKpoints, "Kpoints",
	DumpGvectors, "Gvectors",
	DumpOrbitalDep, "OrbitalDep",
	DumpXCanalysis, "XCanalysis",
	DumpEresolvedDensity, "EresolvedDensity",
	DumpFermiDensity, "FermiDensity",
	DumpDWfns, "DWfns",
	DumpDn, "Dn",
	DumpDVext, "DVext",
	DumpDVscloc, "DVscloc",
	DumpHC, "HC"
);
EnumStringMap<DumpVariable> varDescMap
(	DumpNone,           "Dump nothing",
	DumpState,          "All variables needed to restart calculation: wavefunction and fluid state/fillings if any",
	DumpIonicPositions, "Ionic positions in the same format (and coordinate system) as the input file",
	DumpForces,         "Forces on the ions in the coordinate system selected by command forces-output-coords",
	DumpLattice,        "Lattice vectors in the same format as the input file",
	DumpIonicDensity,   "Nuclear charge density (with gaussians)",
	DumpElecDensity,    "Electronic densities (n or nup,ndn)",
	DumpElecDensityAccum, "Electronic densities (n or nup,ndn) accumulated over MD trajectory",
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
	DumpBandProjections,"Projections of each band state against each atomic orbital",
	DumpEigStats,       "Band eigenvalue statistics: HOMO, LUMO, min, max and Fermi level",
	DumpFillings,       "Fillings",
	DumpRhoAtom,        "Atomic-orbital projected density matrices (only for species with +U enabled)",
	DumpBandUnfold,     "Unfold band structure from supercell to unit cell (see command band-unfold)",
	DumpEcomponents,    "Components of the energy",
	DumpBoundCharge,    "Bound charge in the fluid",
	DumpSolvationRadii, "Effective solvation radii based on fluid bound charge distribution",
	DumpQMC,            "Blip'd orbitals and potential for CASINO \\cite Katie-QMC",
	DumpOcean,          "Wave functions for Ocean code",
	DumpBGW,            "G-space wavefunctions, density and potential for Berkeley GW (requires HDF5 support)",
	DumpRealSpaceWfns,  "Real-space wavefunctions (one column per file)",
	DumpExcCompare,     "Energies for other exchange-correlation functionals (see command elec-ex-corr-compare)",
	DumpFluidDebug,     "Fluid specific debug output if any ",
	DumpSlabEpsilon,    "Local dielectric function of a slab (see command slab-epsilon) ",
	DumpBulkEpsilon,    "Dielectric constant of a periodic solid (see command bulk-epsilon) ",
	DumpChargedDefect,  "Calculate energy correction for charged defect (see command charged-defect) ",
	DumpDOS,            "Density of States (see command density-of-states)",
	DumpSIC,            "Calculates Perdew-Zunger self-interaction corrected Kohn-Sham eigenvalues",
	DumpDipole,         "Dipole moment of explicit charges (ionic and electronic)",
	DumpStress,         "Dumps dE/dR_ij where R_ij is the i'th component of the j'th lattice vector",
	DumpExcitations,    "Dumps dipole moments and transition strength (electric-dipole) of excitations",
	DumpFCI,            "Output Coulomb matrix elements in FCIDUMP format",
	DumpSpin,           "Spin matrix elements from non-collinear calculations in a binary file (indices outer to inner: state, cartesian direction, band1, band2)",
	DumpMomenta,        "Momentum matrix elements in a binary file (indices outer to inner: state, cartesian direction, band1, band2)",
	DumpVelocities,     "Diagonal momentum/velocity matrix elements in a binary file  (indices outer to inner: state, band, cartesian direction)",
	DumpFermiVelocity,  "Fermi velocity, density of states at Fermi level and related quantities",
	DumpR,              "Position operator matrix elements, only allowed at End (see command Cprime-params)",
	DumpL,              "Angular momentum matrix elements, only allowed at End (see command Cprime-params)",
	DumpQ,              "Quadrupole r*p matrix elements, only allowed at End (see command Cprime-params)",
	DumpBerry,          "Berry curvature i <dC/dk| X |dC/dk>, only allowed at End (see command Cprime-params)",
	DumpSymmetries,     "List of symmetry matrices (in covariant lattice coordinates)",
	DumpKpoints,        "List of reduced k-points in calculation, and mapping to the unreduced k-point mesh",
	DumpGvectors,       "List of G vectors in reciprocal lattice basis, for each k-point",
	DumpOrbitalDep,     "Custom output from orbital-dependent functionals (eg. quasi-particle energies, discontinuity potential)",
	DumpXCanalysis,     "Debug VW KE density, single-particle-ness and spin-polarzied Hartree potential",
	DumpEresolvedDensity, "Electron density from bands within specified energy ranges",
	DumpFermiDensity,	"Electron density from fermi-derivative at specified energy",
	DumpDWfns, "Perturbation Wavefunctions",
	DumpDn, 			"First order change in electronic density",
	DumpDVext, 			"External perturbation",
	DumpDVscloc, 		"First order change in local self-consistent potential",
	DumpHC,				"Hamiltonian applied on wavefunctions for non self-consistent forces"
);

struct CommandDump : public Command
{
	CommandDump() : Command("dump", "jdftx/Output")
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
			pl.get(var, DumpDelim, freq==DumpFreq_Init ? varInitMap : varMap, "var");
			if(var==DumpDelim) break; //will happen at end of command line
			e.dump.insert(std::make_pair(freq,var));
			//Check for unsupported features:
			#ifndef HDF5_ENABLED
			if(var==DumpBGW) throw string("BerkeleyGW interface requires HDF5 support (CMake option EnableHDF5)\n");
			#endif
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
	CommandDumpInterval() : Command("dump-interval", "jdftx/Output")
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
		if(freq==DumpFreq_End || freq==DumpFreq_Init)
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
	CommandDumpName() : Command("dump-name", "jdftx/Output")
	{
		format = "<format> [<freq1> <format1>] [<freq2> <format2>] ...";
		comments = 
			"Control the filename pattern for dump output, where <format> is an\n"
			"arbitrary format string that will be substituted according to:\n"
			"+ $VAR   -> name of the variable being dumped (this must be present)\n"
			"+ $ITER  -> iteration number of relevant dump frequency\n"
			"+ $INPUT -> base name of input file, or 'stdin'\n"
			"+ $STAMP -> time-stamp at the start of dump\n"
			"\n"
			"Optionally, a different <format> could be specified for some dump frequencies.";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.dump.format, string("$INPUT.$VAR"), "format");
		if(e.dump.format.find("$VAR")==string::npos)
			throw "<format> = " + e.dump.format + " doesn't contain the pattern $VAR";
		//Check for additional specifictaions:
		while(true)
		{	DumpFrequency freq; pl.get(freq, DumpFreq_Delim, freqMap, "<freqN>");
			if(freq==DumpFreq_Delim) break; //no more freq's
			string format; pl.get(format, string(), "<formatN>", true);
			if(format.find("$VAR")==string::npos)
				throw "<format> = " + format + " doesn't contain the pattern $VAR";
			e.dump.formatFreq[freq] = format;
		}
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s", e.dump.format.c_str());
		for(auto entry: e.dump.formatFreq)
			logPrintf(" \\\n\t%s %s", freqMap.getString(entry.first), entry.second.c_str());
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
    CommandPolarizability() : Command("polarizability", "jdftx/Output")
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
	ESM_RPA,
	ESM_dumpEpsilon,
	ESM_slabResponse,
	ESM_EcutTransverse,
	ESM_computeRange,
	ESM_delim
};
EnumStringMap<ElectronScatteringMember> esmMap
(	ESM_eta, "eta",
	ESM_Ecut, "Ecut",
	ESM_fCut, "fCut",
	ESM_omegaMax, "omegaMax",
	ESM_dumpEpsilon, "dumpEpsilon",
	ESM_RPA, "RPA",
	ESM_slabResponse, "slabResponse",
	ESM_EcutTransverse, "EcutTransverse",
	ESM_computeRange, "computeRange"
);

struct CommandElectronScattering : public Command
{
    CommandElectronScattering() : Command("electron-scattering", "jdftx/Output")
	{
		format = "<key1> <value1> ...";
		comments = "Calculate electron-electron scattering rates (expensive!)\n"
			"and output contribution to imaginary part of electron self-energy\n"
			"(calculated effectively using full-frequency G0W0).\n"
			"\n"
			"The following key-value pairs can appear in any order:\n"
			"\n+ eta <eta>\n\n"
			"   <eta> in Eh specifies frequency grid resolution (required)\n"
			"\n+ Ecut <Ecut>\n\n"
			"   <Ecut> in Eh specifies energy cut-off for dielectric matrices.\n"
			"   (If zero, the wavefunction cutoff from elec-cutoff is used instead.)\n"
			"\n+ fCut <fCut>\n\n"
			"   <fCut> specifies threshold for considering states fully occupied or\n"
			"   unoccupied in optimizing sums over states (default: 1e-6)\n"
			"\n+ omegaMax <omegaMax>\n\n"
			"   <omegaMax> in Eh is the maximum energy transfer to account for\n"
			"   and hence the maximum frequency in dielectric function frequency grid.\n"
			"   (if zero, autodetermine from available eigenvalues)\n"
			"\n+ RPA yes|no\n\n"
			"   If yes, use RPA response that ignores XC contribution. (default: no).\n"
			"\n+ dumpEpsilon yes|no\n\n"
			"   If yes, dump dielectric function in GG' basis. (default: no).\n"
			"\n+ slabResponse yes|no\n\n"
			"   Whether to output slab-normal-direction susceptibility instead.\n"
			"   This needs slab geometry in coulomb-interaction, and will bypass the\n"
			"   actual electron-electron scattering calculation and output.\n"
			"\n+ EcutTransverse <EcutTransverse>\n\n"
			"   <EcutTransverse> in Eh specifies energy cut-off for dielectric matrix in.\n"
			"   directions trasverse to the slab normal; only valid when slabResponse = yes.\n"
			"   (If zero, use the same value as Ecut above.)\n\n"
			"\n+ computeRange <iqStart> <iqStop>\n\n"
			"   If specified, only calculate momentum transfers in range [iqStart , iqStop] in\n"
			"   the current run, in order to split the overall calculation into smaller jobs.\n"
			"   Note that the indices are 1-based, and the range includes both end-points.\n"
			"   To combine the final results, perform a final run without computeRange specified.";
			
		require("coulomb-interaction");
		forbid("polarizability"); //both are major operations that are given permission to destroy Everything if necessary
	}
	
	void process(ParamList& pl, Everything& e)
	{	e.dump.electronScattering = std::make_shared<ElectronScattering>();
		e.dump.insert(std::make_pair(DumpFreq_End, DumpElectronScattering));
		ElectronScattering& es = *(e.dump.electronScattering);
		while(true)
		{	ElectronScatteringMember key;
			pl.get(key, ESM_delim, esmMap, "key");
			if(key == ESM_delim) break; //end of input
			switch(key)
			{	case ESM_eta: pl.get(es.eta, 0., "eta", true); break;
				case ESM_Ecut: pl.get(es.Ecut, 0., "Ecut", true); break;
				case ESM_fCut: pl.get(es.fCut, 0., "fCut", true); break;
				case ESM_omegaMax: pl.get(es.omegaMax, 0., "omegaMax", true); break;
				case ESM_RPA: pl.get(es.RPA, false, boolMap, "RPA", true); break;
				case ESM_dumpEpsilon: pl.get(es.dumpEpsilon, false, boolMap, "dumpEpsilon", true); break;
				case ESM_slabResponse: pl.get(es.slabResponse, false, boolMap, "slabResponse", true); break;
				case ESM_EcutTransverse: pl.get(es.EcutTransverse, 0., "EcutTransverse", true); break;
				case ESM_computeRange:
				{	es.computeRange = true; //Note that inputs are [start,stop] 1-based, while internally we have [start,stop) 0-based
					pl.get(es.iqStart, size_t(0), "iqStart", true); if(es.iqStart < 1) throw string("Must have iqStart >= 1");
					pl.get(es.iqStop, size_t(0), "iqStop", true); if(es.iqStop < es.iqStart) throw string("Must have iqStop >= iqStart");
					es.iqStart -= 1; //convert to 0-based index. Note that iqStop becomes a non-included 0-based index without change
					break;
				}
				case ESM_delim: break; //never encountered; to suppress compiler warning
			}
		}
		if(es.slabResponse)
		{	if(e.coulombParams.geometry != CoulombParams::Slab)
				throw string("slabResponse = yes requires slab geometry in coulomb-interaction");
		}
		else
		{	if(es.EcutTransverse) throw string("Cannot specify EcutTransverse when slabResponse = no");
		}
		if(es.eta <= 0.) throw string("Must specify frequency grid resolution eta > 0.");
	}

	void printStatus(Everything& e, int iRep)
	{	const ElectronScattering& es = *(e.dump.electronScattering);
		logPrintf(" \\\n\teta      %lg", es.eta);
		logPrintf(" \\\n\tEcut     %lg", es.Ecut);
		logPrintf(" \\\n\tfCut     %lg", es.fCut);
		logPrintf(" \\\n\tomegaMax %lg", es.omegaMax);
		logPrintf(" \\\n\tRPA      %s", boolMap.getString(es.RPA));
		logPrintf(" \\\n\tdumpEpsilon      %s", boolMap.getString(es.dumpEpsilon));
		logPrintf(" \\\n\tslabResponse %s", boolMap.getString(es.slabResponse));
		if(es.slabResponse) logPrintf(" \\\n\tEcutTransverse %lg", es.EcutTransverse);
		if(es.computeRange)  logPrintf(" \\\n\tcomputeRange %lu %lu", es.iqStart+1, es.iqStop);
	}
}
commandElectronScattering;


struct CommandPolarizabilityKdiff : public Command
{
    CommandPolarizabilityKdiff() : Command("polarizability-kdiff", "jdftx/Output")
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
    CommandDumpEresolvedDensity() : Command("dump-Eresolved-density", "jdftx/Output")
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
    CommandDumpFermiDensity() : Command("dump-fermi-density", "jdftx/Output")
	{
		format = "[<muLevel>]";
		comments =
			"Output electron density calculated from derivative of smearing\n"
			"function evaluated at desired chemical potential <muLevel>.\n"
			"If unspecified, calculated chemical potential will be used.\n"
			"\n"
			"When issued multiple times, the outputs will be\n"
			"numbered sequenetially FermiDensity.0 etc.\n"
			"This automatically invokes dump at End; dumping at\n"
			"other frequencies may be requested using the dump command.";
		
		allowMultiple = true;
		require("elec-smearing"); //require smearing
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
    CommandVibrations() : Command("vibrations", "jdftx/Output")
	{
		format = "<key1> <args1> ...";
		comments =
			"Calculate vibrational modes of the system using a finite difference method.\n"
			"Note that this command should typically be issued in a run with converged ionic\n"
			"positions; ionic (and lattice) minimization are bypassed by the vibrations module.\n"
			"\n"
			"Any number of the following subcommands and their arguments may follow:\n"
			"+ dr <dr>: perturbation amplitude in bohrs for force matrix calculation (default: 0.01).\n"
			"+ centralDiff yes|no: use a central difference formula for the second derivative\n"
			"   to achieve higher accuracy at twice the cost (default: no)\n"
			"+ useConstraints yes|no: restrict modes of motion as specified by move flags\n"
			"   and constraints in the ion command (default: no)\n"
			"+ translationSym yes|no: whether to assume overall translation symmetry (default yes).\n"
			"   Can be turned off to get vibrational levels in an external potential.\n"
			"+ rotationSym yes|no: project out rotational modes (default no). Improves reliability for\n"
			"   molecular calculations. Valid only for geometries with an unambiguous center of mass.\n"
			"+ omegaMin <omegaMin>: frequency cutoff (in Eh) for free energy calculation (default: 2e-4)\n"
			"+ T <T>: temperature (in Kelvin) for free energy calculation (default: 298)\n"
			"+ omegaResolution <omegaResolution>: resolution for detecting and reporting degeneracies\n"
			"   in modes (default: 1e-4). Does not affect free energies and all modes are still printed.\n"
			"\n"
			"Note that for a periodic system with k-points, wave functions may be incompatible\n"
			"with and without the vibrations command due to symmetry-breaking by the perturbations.\n"
			"To avoid this, perform electronic optimization (without initial-state) in the\n"
			"vibration calculation itself, or consider using the phonon code instead."
			;
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
	CommandSlabEpsilon() : Command("slab-epsilon", "jdftx/Output")
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
		if(!(se.Efield-e.coulombParams.Efield).length_squared())
			throw(string("Applied electric fields in reference and present calculations are equal"));
		e.dump.insert(std::make_pair(DumpFreq_End, DumpSlabEpsilon)); //dump at end by default
	}

	void printStatus(Everything& e, int iRep)
	{	const SlabEpsilon& se = *(e.dump.slabEpsilon);
		logPrintf("%s %lg %lg %lg %lg", se.dtotFname.c_str(), se.sigma, se.Efield[0], se.Efield[1], se.Efield[2]);
	}
}
commandSlabEpsilon;

//-------------------------------------------------------------------------------------------------

struct CommandBulkEpsilon : public Command
{
	CommandBulkEpsilon() : Command("bulk-epsilon", "jdftx/Output")
	{
		format = "<DtotFile> [<Ex>=0] [<Ey>=0] [<Ez>=0]";
		comments = 
			"Calculate dielectric constant of a bulk material given the electrostatic potential\n"
			"output from another calculation on same system with a different electric field.\n"
			"+ <DtotFile> contains the electrostatic potential from the other calculation\n"
			"+ optional <Ex>,<Ey>,Ez> specify the electric-field applied\n"
			"  in the calculation that generated <DtotFile>.\n"
			"It is recommended to apply field only to one reciprocal lattice direction,\n"
			"and use a supercell of the material along that direction.";
		
		require("electric-field");
	}

	void process(ParamList& pl, Everything& e)
	{	if(e.coulombParams.geometry != CoulombParams::Periodic)
			throw string("coulomb-interaction must be in Periodic mode");
		e.dump.bulkEpsilon = std::make_shared<BulkEpsilon>();
		BulkEpsilon& be = *(e.dump.bulkEpsilon);
		pl.get(be.dtotFname, string(), "DtotFile", true);
		pl.get(be.Efield[0], 0., "Ex");
		pl.get(be.Efield[1], 0., "Ey");
		pl.get(be.Efield[2], 0., "Ez");
		if(!(be.Efield-e.coulombParams.Efield).length_squared())
			throw(string("Applied electric fields in reference and present calculations are equal"));
		e.dump.insert(std::make_pair(DumpFreq_End, DumpBulkEpsilon)); //dump at end by default
	}

	void printStatus(Everything& e, int iRep)
	{	const BulkEpsilon& be = *(e.dump.bulkEpsilon);
		logPrintf("%s %lg %lg %lg", be.dtotFname.c_str(), be.Efield[0], be.Efield[1], be.Efield[2]);
	}
}
commandBulkEpsilon;

//-------------------------------------------------------------------------------------------------

extern EnumStringMap<int> truncationDirMap;

struct CommandChargedDefectCorrection : public Command
{
	CommandChargedDefectCorrection() : Command("charged-defect-correction", "jdftx/Output")
	{
		format = "[Slab <dir>=100|010|001] <DtotFile> <bulkEps>|<slabEpsFile> <rMin> <rSigma>";
		comments = 
			"Calculate energy correction for bulk or surface charged defects \\cite ElectrostaticPotential\n"
			"The correction is calculated assuming the defects to be model\n"
			"charges specified using command charged-defect.\n"
			"\n"
			"By default, the defect is assumed bulk for coulomb-interaction Periodic\n"
			"and surface for coulomb-interaction Slab.  However, for the Periodic case,\n"
			"the optional [Slab <dir>] overrides this to calculate surface defects\n"
			"without truncated Coulomb potentials.  Note that coulomb-truncation-embed\n"
			"must be specified when using truncated coulomb potentials in Slab mode.\n"
			"In Periodic mode, the correction assumes a slab centered at the origin\n"
			"(i.e. analogous to using xCenter 0 0 0 in the truncated mode).\n"
			"\n"
			"<DtotFile> contains the electrostatic potential from a reference\n"
			"neutral calculation with similar geometry (lattice vectors and grid\n"
			"must match exactly).\n"
			"\n"
			"For bulk defect calculations, <bulkEps> is the bulk dielectric constant.\n"
			"\n"
			"For surface defect calculations, <slabEpsFile> specifies a dielectric\n"
			"profile calculated using command slab-epsilon in a similar geometry\n"
			"(the number of points along the slab normal direction must match exactly).\n"
			"Optionally, the <slabEpsFile> may contain an additional column for the\n"
			"in-plane response (which is not computed by command slab-epsilon),\n"
			"in which case an anisotropic dielectric model is used.\n"
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
		//Check for geometry override:
		cd.geometry = e.coulombParams.geometry;
		cd.iDir = e.coulombParams.iDir;
		string slabSpec;
		pl.get(slabSpec, string(), "Slab|DtotFile", true);
		if(slabSpec == "Slab")
		{	if(cd.geometry != CoulombParams::Periodic)
				throw string("Slab geometry override should only be specified for coulomb-interaction Periodic");
			cd.geometry = CoulombParams::Slab;
			pl.get(cd.iDir, 0, truncationDirMap, "dir", true);
		}
		else pl.rewind();
		//Ref potential:
		pl.get(cd.dtotFname, string(), "DtotFile", true);
		//Dielectric function (and add citation depending on geometry):
		switch(cd.geometry)
		{	case CoulombParams::Periodic:
			{	pl.get(cd.bulkEps, 1., "bulkEps", true);
				Citations::add("Correction scheme for charged bulk defects",
					"C Freysoldt, J Neugebauer and C. van de Walle, Phys. Rev. Lett. 102, 016402 (2009)");
				break;
			}
			case CoulombParams::Slab:
			{	pl.get(cd.slabEpsFname, string(), "slabEpsFile", true);
				Citations::add("Correction scheme for charged surface defects",
					"H Komsa and A Pasquarello, Alfredo, Phys. Rev. Lett. 110, 095505 (2013)");
				break;
			}
			default: throw string("coulomb-interaction must be either Slab or Periodic");
		}
		//Alignment potential ranges:
		pl.get(cd.rMin, 0., "rMin", true);
		pl.get(cd.rSigma, 0., "rSigma", true);
		e.dump.insert(std::make_pair(DumpFreq_End, DumpChargedDefect)); //dump at end by default
	}

	void printStatus(Everything& e, int iRep)
	{	const ChargedDefect& cd = *(e.dump.chargedDefect);
		if(cd.geometry != e.coulombParams.geometry)
			logPrintf("Slab %s ", truncationDirMap.getString(cd.iDir));
		logPrintf("%s ", cd.dtotFname.c_str());
		switch(cd.geometry)
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
	CommandChargedDefect() : Command("charged-defect", "jdftx/Output")
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

struct CommandPotentialSubtraction : public Command
{
	CommandPotentialSubtraction() : Command("potential-subtraction", "jdftx/Output")
	{	format = "<subtract>=yes|no";
		comments = 
			"Whether to subtract neutral atom potentials in dumped potentials (Dtot and Dvac).\n"
			"This subtraction produces much smoother potentials and is enabled by default \\cite ElectrostaticPotential.";
	}
	
	void process(ParamList& pl, Everything& e)
	{	pl.get(e.dump.potentialSubtraction, true, boolMap, "subtract");
	}
	
	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s", boolMap.getString(e.dump.potentialSubtraction));
	}
}
commandPotentialSubtraction;


struct CommandBandUnfold : public Command
{
	CommandBandUnfold() : Command("band-unfold", "jdftx/Output")
	{
		format = " \\\n\t<M00> <M01> <M02> \\\n\t<M10> <M11> <M12> \\\n\t<M20> <M21> <M22>";
		comments =
			"Unfold band structure from a supercell calculation to a unit cell\n"
			"with lattice vectors Runit, defined by the integer matrix M such\n"
			"that current lattice vectors R = Runit * M.";
	}

	void process(ParamList& pl, Everything& e)
	{	matrix3<int>& M = e.dump.Munfold;
		for(int j=0; j<3; j++) for(int k=0; k<3; k++)
		{	ostringstream oss; oss << "s" << j << k;
			pl.get(M(j,k), 0, oss.str(), true);
		}
		e.dump.insert(std::make_pair(DumpFreq_End, DumpBandUnfold));
	}

	void printStatus(Everything& e, int iRep)
	{	for (int j=0; j < 3; j++)
		{	logPrintf(" \\\n\t");
			for(int k=0; k<3; k++) logPrintf("%d ",e.dump.Munfold(j,k));
		}
	}
}
commandBandUnfold;


enum BGWparamsMember
{	BGWpm_nBandsDense,
	BGWpm_blockSize,
	BGWpm_clusterSize,
	BGWpm_nBandsV,
	BGWpm_saveVxc,
	BGWpm_saveVxx,
	BGWpm_rpaExx,
	BGWpm_offDiagV,
	BGWpm_EcutChiFluid,
	BGWpm_elecOnly,
	BGWpm_q0,
	BGWpm_freqReMax_eV,
	BGWpm_freqReStep_eV,
	BGWpm_freqBroaden_eV,
	BGWpm_freqNimag,
	BGWpm_freqPlasma,
	BGWpm_Ecut_rALDA,
	BGWpm_kFcut_rALDA,
	BGWpm_kernelSym_rALDA,
	BGWpm_Delim
};
EnumStringMap<BGWparamsMember> bgwpmMap
(	BGWpm_nBandsDense, "nBandsDense",
	BGWpm_blockSize, "blockSize",
	BGWpm_clusterSize, "clusterSize",
	BGWpm_nBandsV, "nBandsV",
	BGWpm_saveVxc, "saveVxc",
	BGWpm_saveVxx, "saveVxx",
	BGWpm_rpaExx, "rpaExx",
	BGWpm_offDiagV, "offDiagV",
	BGWpm_EcutChiFluid, "EcutChiFluid",
	BGWpm_elecOnly, "elecOnly",
	BGWpm_q0, "q0",
	BGWpm_freqReMax_eV, "freqReMax_eV",
	BGWpm_freqReStep_eV, "freqReStep_eV",
	BGWpm_freqBroaden_eV, "freqBroaden_eV",
	BGWpm_freqNimag, "freqNimag",
	BGWpm_freqPlasma, "freqPlasma",
	BGWpm_Ecut_rALDA, "Ecut_rALDA",
	BGWpm_kFcut_rALDA, "kFcut_rALDA",
	BGWpm_kernelSym_rALDA, "kernelSym_rALDA"
);
EnumStringMap<BGWparamsMember> bgwpmDescMap
(	BGWpm_nBandsDense, "If non-zero, use a dense ScaLAPACK solver to calculate more bands",
	BGWpm_blockSize, "Block size for ScaLAPACK diagonalization (default: 32)",
	BGWpm_clusterSize, "Maximum eigenvalue cluster size to allocate extra ScaLAPACK workspace for (default: 10)",
	BGWpm_nBandsV, "If non-zero, number of bands for Vxc and Vxx output",
	BGWpm_saveVxc, "Whether to write exchange-correlation matrix elements (default: yes)",
	BGWpm_saveVxx, "Whether to write exact-exchange matrix elements (default: no)",
	BGWpm_rpaExx, "Whether to compute RPA-consistent exact-exchange energy (default: no)",
	BGWpm_offDiagV, "Whether to write off-diagonal matrix elements of Vxc and/or Vxx (default: no)",
	BGWpm_EcutChiFluid, "KE cutoff in hartrees for fluid polarizability output (default: 0; set non-zero to enable)",
	BGWpm_elecOnly, "Whether fluid polarizability output should only include electronic response (default: true)",
	BGWpm_q0, "Zero wavevector replacement to be used for polarizability output (default: (0,0,0))",
	BGWpm_freqReMax_eV, "Maximum real frequency in eV (default: 30.)",
	BGWpm_freqReStep_eV, "Real frequency grid spacing in eV (default: 1.)",
	BGWpm_freqBroaden_eV, "Broadening (imaginary part) of real frequency grid in eV (default: 0.1)",
	BGWpm_freqNimag, "Number of imaginary frequencies (default: 25)",
	BGWpm_freqPlasma, "Plasma frequency in Hartrees used in GW imaginary frequency grid (default: 1.), set to zero for RPA frequency grid",
	BGWpm_Ecut_rALDA, "KE cutoff in hartrees for rALDA polarizability output (default: 0; set non-zero to enable)",
	BGWpm_kFcut_rALDA, "kF cutoff (in 1/a0) for rALDA regularization (enabled if non-zero)",
	BGWpm_kernelSym_rALDA, "Kernel symmetrization for rALDA if true, and wavevector symmetrization if false (default)"
);

struct CommandBGWparams : public Command
{
	CommandBGWparams() : Command("bgw-params", "jdftx/Output")
	{	
		format = "<key1> <value1> <key2> <value2> ...";
		comments = "Control BGW output. Possible keys and value types are:"
			+ addDescriptions(bgwpmMap.optionList(), linkDescription(bgwpmMap, bgwpmDescMap))
			+ "\n\nAny number of these key-value pairs may be specified in any order.";
	}

	void process(ParamList& pl, Everything& e)
	{	e.dump.bgwParams = std::make_shared<BGWparams>();
		BGWparams& bgwp = *(e.dump.bgwParams);
		while(true)
		{	BGWparamsMember key;
			pl.get(key, BGWpm_Delim, bgwpmMap, "key");
			#define READ_AND_CHECK(param, op, val) \
				case BGWpm_##param: \
					pl.get(bgwp.param, val, #param, true); \
					if(!(bgwp.param op val)) throw string(#param " must be " #op " " #val); \
					break;
			#define READ_BOOL(param) \
				case BGWpm_##param: \
					pl.get(bgwp.param, false, boolMap, #param, true); \
					break;
			switch(key)
			{	READ_AND_CHECK(nBandsDense, >=, 0)
				READ_AND_CHECK(blockSize, >, 0)
				READ_AND_CHECK(clusterSize, >, 0)
				READ_AND_CHECK(nBandsV, >=, 0)
				READ_BOOL(saveVxc)
				READ_BOOL(saveVxx)
				READ_BOOL(rpaExx)
				READ_BOOL(offDiagV)
				READ_AND_CHECK(EcutChiFluid, >=, 0.)
				READ_BOOL(elecOnly)
				case BGWpm_q0:
					for(int dir=0; dir<3; dir++)
						pl.get(bgwp.q0[dir], 0., "q0", true);
					break;
				READ_AND_CHECK(freqReMax_eV, >, 0.)
				READ_AND_CHECK(freqReStep_eV, >, 0.)
				READ_AND_CHECK(freqBroaden_eV, >, 0.)
				READ_AND_CHECK(freqNimag, >, 0)
				READ_AND_CHECK(freqPlasma, >=, 0.)
				READ_AND_CHECK(Ecut_rALDA, >=, 0.)
				READ_AND_CHECK(kFcut_rALDA, >=, 0.)
				READ_BOOL(kernelSym_rALDA)
				case BGWpm_Delim: return; //end of input
			}
			#undef READ_AND_CHECK
			#undef READ_BOOL
		}
	}

	void printStatus(Everything& e, int iRep)
	{	assert(e.dump.bgwParams);
		const BGWparams& bgwp = *(e.dump.bgwParams);
		#define PRINT(param, format) logPrintf(" \\\n\t" #param " " format, bgwp.param);
		#define PRINT_BOOL(param) logPrintf(" \\\n\t" #param " %s", boolMap.getString(bgwp.param));
		PRINT(nBandsDense, "%d")
		PRINT(blockSize, "%d")
		PRINT(clusterSize, "%d")
		PRINT(nBandsV, "%d")
		PRINT_BOOL(saveVxc)
		PRINT_BOOL(saveVxx)
		PRINT_BOOL(rpaExx)
		PRINT_BOOL(offDiagV)
		PRINT(EcutChiFluid, "%lg")
		PRINT_BOOL(elecOnly)
		logPrintf(" \\\n\tq0 %lg %lg %lg", bgwp.q0[0], bgwp.q0[1], bgwp.q0[2]);
		PRINT(freqReMax_eV, "%lg")
		PRINT(freqReStep_eV, "%lg")
		PRINT(freqBroaden_eV, "%lg")
		PRINT(freqNimag, "%d")
		PRINT(freqPlasma, "%lg")
		PRINT(Ecut_rALDA, "%lg")
		PRINT(kFcut_rALDA, "%lg")
		PRINT_BOOL(kernelSym_rALDA)
		#undef PRINT
		#undef PRINT_BOOL
	}
}
commandBGWparams;


struct CommandBandProjectionParams : public Command
{
	CommandBandProjectionParams() : Command("band-projection-params", "jdftx/Output")
	{	
		format = "<ortho>=yes|no <norm>=yes|no";
		comments = "Control band-projections output:\n"
			"\t<ortho>: whether to use ortho-orbitals.\n"
			"\t<norm>: whether to output only norm or complex amplitude.";
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.dump.bandProjectionOrtho, false, boolMap, "ortho", true);
		pl.get(e.dump.bandProjectionNorm, true, boolMap, "norm", true);
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s %s",
			boolMap.getString(e.dump.bandProjectionOrtho),
			boolMap.getString(e.dump.bandProjectionNorm));
	}
}
commandBandProjectionParams;


struct CommandCprimeParams : public Command
{
	CommandCprimeParams() : Command("Cprime-params", "jdftx/Output")
	{	
		format = "[<dk>=1E-4] [<degeneracyThreshold>=1E-6] [<vThreshold>=1E-4] [<realSpaceTruncated>=yes]";
		comments = "Control dC/dk calculation for L and Q output.\n"
			"Here, <dk> in a0^-1 controls the finite difference used for dC/dk,\n"
			"while <degeneracyThreshold> specifies the energy range within unitary\n"
			"rotations are accounted for in comparing wavefunctions between k.\n"
			"Within degenerate subspaces of energy, rotations are first revolved\n"
			"by the velocity operator and then obtained by best match between\n"
			"wavefunctions for sub-subspaces that have the same velocity within\n"
			"<vThreshold> (in Eh-a0 atomic units).\n"
			"If <realSpaceTruncated> is yes (default), then truncated directions\n"
			"are computed by a real space multiplication by r, instead of dC/dk.\n"
			"\n"
			"This is used for the calculation of orbital angular momenta L, output\n"
			"in the same binary format as the momenta and selected by 'dump End L'.\n"
			"\n"
			"Additionally, 'dump End Q' selects output of electric quadrupole matrix\n"
			"elements, defined as traceless symmetric tensor of (rj pk + pj rk).\n"
			"Specifically, the output is a nStates x 5 x nBands x nBands complex\n"
			"binary file, where the 5 components in order are xy, yz, zx, xxr, yyr\n"
			"(where xxr = xx - rr/3). Note that zzr = -(xxr + yyr) because the trace\n"
			"of the tensor is removed, and this redundant component is excluded.";
	}

	void process(ParamList& pl, Everything& e)
	{	e.dump.dumpCprime = std::make_shared<DumpCprime>();
		DumpCprime& dcp = *(e.dump.dumpCprime);
		pl.get(dcp.dk, 1E-4, "dk");
		pl.get(dcp.degeneracyThreshold, 1E-6, "degeneracyThreshold");
		pl.get(dcp.vThreshold, 1E-4, "vThreshold");
		pl.get(dcp.realSpaceTruncated, true, boolMap, "realSpaceTruncated");
	}

	void printStatus(Everything& e, int iRep)
	{	assert(e.dump.dumpCprime);
		const DumpCprime& dcp = *(e.dump.dumpCprime);
		logPrintf("%lg %lg %lg %s", dcp.dk, dcp.degeneracyThreshold, dcp.vThreshold, boolMap.getString(dcp.realSpaceTruncated));
	}
}
commandCprimeParams;
