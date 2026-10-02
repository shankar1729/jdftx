/* One-shot density-to-potential utility. Licensed under GPL v3 or later. */
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/DefectCoulomb.h>
#include <electronic/SpeciesInfo.h>
#include <commands/parser.h>
#include <core/Util.h>
#include <cmath>
#include <cstdlib>
#include <set>

namespace
{
	using Input = std::vector<std::pair<string,string>>;

	// Discard execution/restart/output controls, preserving the physical Hamiltonian.
	const std::set<string> discard = {
		"electronic-scf", "electronic-minimize", "ionic-minimize", "lattice-minimize", "ionic-dynamics",
		"initial-state", "wavefunction", "lcao-params", "elec-initial-fillings", "elec-initial-eigenvals",
		"fix-electron-density", "fix-electron-potential", "dump-only", "defect-coulomb", "defect-coulomb-center",
		"dump", "dump-name", "dump-interval", "density-of-states", "band-projection-params",
		"band-unfold", "Cprime-params", "polarizability", "electron-scattering", "slab-epsilon",
		"bulk-epsilon", "charged-defect", "charged-defect-correction", "bgw-params",
		"elec-ex-corr-compare", "elec-initial-state", "elec-density-of-states", "elec-eigenvalue-stats",
		"dump-Eresolved-density", "dump-Fermi-density"
	};

	void validateDensity(const Everything& e)
	{	for(const auto& channel: e.eVars.n)
		{	const double* data = channel->data();
			for(int i=0; i<e.gInfo.nr; i++)
				if(!std::isfinite(data[i])) die("Supplied density contains non-finite values.\n");
		}
	}

	void writeNscf(const Everything& e, const Input& physical, const string& prefix)
	{	if(!mpiWorld->isHead()) return;
		string filename = prefix + ".nscf.in";
		FILE* fp = fopen(filename.c_str(), "w");
		if(!fp) die("Cannot write NSCF input '%s'.\n", filename.c_str());
		fprintf(fp, "# Generated fixed-potential NSCF input.\n"
			"# Run from the original defect-input directory so relative species/files resolve.\n"
			"# Adjust the k points and elec-n-bands for your NSCF calculation.\n");
		for(const auto& entry: physical)
			if(entry.first!="symmetries" && entry.first!="symmetry-matrix" && entry.first!="fftbox" && entry.first!="ion-species")
				fprintf(fp, "%s %s\n", entry.first.c_str(), entry.second.c_str());
		// Pin the actual files used here so the NSCF executable need not have the
		// same build directory or pseudopotential-search environment as this utility.
		for(auto sp: e.iInfo.species)
		{	char* path = realpath(sp->potfilename.c_str(), nullptr);
			if(!path) die("Cannot resolve pseudopotential path '%s'.\n", sp->potfilename.c_str());
			fprintf(fp, "ion-species %s\n", path);
			free(path);
		}
		fprintf(fp, "symmetries none\nfftbox %d %d %d\n", e.gInfo.S[0], e.gInfo.S[1], e.gInfo.S[2]);
		fprintf(fp, "fix-electron-potential %s.$VAR\ndump-name %s.nscf.$VAR\ndump End State BandEigs\n",
			prefix.c_str(), prefix.c_str());
		if(fclose(fp)) die("Failed to close NSCF input '%s'.\n", filename.c_str());
	}
}

int main(int argc, char** argv)
{	InitParams ip("Construct a defect-only Coulomb NSCF potential without orbital solution; use scripts/defectDensityToVscloc.");
	initSystemCmdline(argc, argv, ip);
	{	Input physical;
		std::map<string,string> options;
		string inputCenter;
		const std::set<string> ldaGga = {
			"lda", "lda-PZ", "lda-PW", "lda-PW-prec", "lda-VWN", "lda-Teter",
			"gga", "gga-PBE", "gga-PBEsol", "gga-PW91"
		};
		for(auto entry: readInputFile(ip.inputFilename))
		{	if(entry.first.find("potential-build-")==0)
			{	trim(entry.second);
				if(options.count(entry.first)) die("Duplicate one-shot option '%s'.\n", entry.first.c_str());
				options[entry.first] = entry.second;
				continue;
			}
			if(entry.first=="defect-coulomb-center")
			{	if(inputCenter.length()) die("Duplicate defect-coulomb-center in source input.\n");
				inputCenter = entry.second;
				continue;
			}
			if(discard.count(entry.first))
			{	logPrintf("One-shot: ignoring %s execution/output control.\n", entry.first.c_str());
				continue;
			}
			if(entry.first=="add-U" || entry.first=="vibrations" || entry.first.find("perturb-")==0)
				die("One-shot density construction does not support %s.\n", entry.first.c_str());
			if(entry.first=="elec-ex-corr")
			{	istringstream args(entry.second); string xc; args >> xc;
				if(!ldaGga.count(xc)) die("One-shot construction requires an internal LDA/GGA functional; '%s' needs additional data or is unsupported.\n", xc.c_str());
			}
			physical.push_back(entry);
		}
		if(!options.count("potential-build-center") && inputCenter.length())
			options["potential-build-center"] = inputCenter;
		const std::set<string> validOptions = {"potential-build-mode", "potential-build-density", "potential-build-output", "potential-build-reference", "potential-build-tolerance", "potential-build-center"};
		for(const auto& option: options) if(!validOptions.count(option.first)) die("Unknown one-shot option '%s'.\n", option.first.c_str());
		string mode = options["potential-build-mode"], density = options["potential-build-density"], prefix = options["potential-build-output"];
		if((mode!="reference" && mode!="corrected") || !prefix.length() || density.find("$VAR")==string::npos)
			die("One-shot recipe requires mode reference/corrected, density pattern containing $VAR, and output prefix.\n");
		Input input = physical;
		input.emplace_back("fix-electron-density", density);
		input.emplace_back("dump-name", prefix + ".$VAR");
		input.emplace_back("dump", mode=="reference" ? "End DefectReference" : "End Vscloc DefectCoulomb");
		Everything e;
		parse(input, e, ip.printDefaults);
		e.eVars.skipWfnsInit = true;
		e.setup(); // Loads the supplied density; no orbitals are allocated or solved.
		validateDensity(e);
		if(e.exCorr.exxFactor() || e.exCorr.needsKEdensity() || e.exCorr.orbitalDep || e.eInfo.hasU)
			die("The density-only one-shot utility requires a local LDA/GGA Hamiltonian.\n");
		if(mode=="reference")
		{	DefectCoulomb::saveReference(e, prefix + ".defectReference");
			logPrintf("Wrote fixed clean reference from supplied density.\n");
		}
		else
		{	// Retain the original grid even if symmetry helped select it. Never symmetrize
			// the mixed potential using only the current ions: the reference can break that symmetry.
			e.symm.mode = SymmetriesNone;
			e.symm.setup(e); e.symm.setupMesh();
			auto correction = std::make_shared<DefectCoulomb>();
			correction->referenceFilename = options["potential-build-reference"];
			if(!correction->referenceFilename.length()) die("Corrected mode requires potential-build-reference.\n");
			if(options.count("potential-build-tolerance"))
			{	istringstream value(options["potential-build-tolerance"]);
				if(!(value >> correction->chargeTolerance) || !std::isfinite(correction->chargeTolerance) || correction->chargeTolerance<=0.)
					die("Invalid positive charge tolerance.\n");
			}
			if(options.count("potential-build-center"))
			{	istringstream value(options["potential-build-center"]);
				for(int k=0; k<3; k++)
					if(!(value >> correction->center[k]) || !std::isfinite(correction->center[k])) die("Invalid defect embedding center.\n");
				string extra; if(value >> extra) die("Defect embedding center requires exactly three coordinates.\n");
				if(e.iInfo.coordsType==CoordsCartesian) correction->center = inv(e.gInfo.R)*correction->center;
				correction->hasCenter = true;
			}
			e.eVars.defectCoulomb = correction;
			correction->setup(e, true);
			e.eVars.EdensityAndVscloc(e.ener);
			correction->report();
			e.dump(DumpFreq_End, 0); // Native Vscloc normalization and spin channel naming.
			writeNscf(e, physical, prefix);
			logPrintf("Constructed Vscloc once at supplied density; no SCF or band minimization was performed.\n");
		}
	}
	finalizeSystem();
	return 0;
}
