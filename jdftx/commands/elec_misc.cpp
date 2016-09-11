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

//! @file elec_misc.cpp Miscellaneous properties of the electronic system

struct CommandElecCutoff : public Command
{
	CommandElecCutoff() : Command("elec-cutoff", "jdftx/Electronic/Parameters")
	{
		format = "<Ecut> [<EcutRho>=0]";
		comments = "Electronic planewave cutoff in Hartree. Optionally specify charge density cutoff\n"
			"<EcutRho> in hartrees. If unspecified or zero, EcutRho is taken to be 4*Ecut.";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.cntrl.Ecut, 20., "Ecut");
		pl.get(e.cntrl.EcutRho, 0., "EcutRho");
		if(e.cntrl.EcutRho && e.cntrl.EcutRho < 4*e.cntrl.Ecut)
			throw string("<EcutRho> must be at least 4 <Ecut>");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg", e.cntrl.Ecut);
		if(e.cntrl.EcutRho)
			logPrintf(" %lg", e.cntrl.EcutRho);
	}
}
commandElecCutoff;

//-------------------------------------------------------------------------------------------------

struct CommandFftbox : public Command
{
	CommandFftbox() : Command("fftbox", "jdftx/Electronic/Parameters")
	{
		format = "<S0> <S1> <S2>";
		comments = "Manually override the real space grid dimensions used for scalar fields.\n"
		 "(The default values are calculated based on the EcutRho setting from elec-cutoff).";
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.gInfo.S[0], 0, "S0", true);
		pl.get(e.gInfo.S[1], 0, "S1", true);
		pl.get(e.gInfo.S[2], 0, "S2", true);
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%d %d %d", e.gInfo.S[0], e.gInfo.S[1], e.gInfo.S[2]);
	}
}
commandFftbox;

//-------------------------------------------------------------------------------------------------

struct CommandElecNbands : public Command
{
	CommandElecNbands() : Command("elec-n-bands", "jdftx/Electronic/Parameters")
	{
		format = "<n>";
		comments = "Manually specify the number of bands.\n\n"
			"(Default: set nBands assuming insulator, or in calculations with\n"
			"fermi-fillings, set equal to total number of atomic orbitals.)";
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.eInfo.nBands, 0, "n", true);
		if(e.eInfo.nBands<=0) throw string("<n> must be positive.\n");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%d", e.eInfo.nBands);
	}
}
commandElecNbands;

//-------------------------------------------------------------------------------------------------

struct CommandDavidsonBandRatio : public Command
{
	CommandDavidsonBandRatio() : Command("davidson-band-ratio", "jdftx/Electronic/Optimization")
	{
		format = "[<ratio>=1.1]";
		comments =
			"Ratio of number of bands in the Davidson working set to the\n"
			"number of actual bands in the calculation. Increasing this\n"
			"number should improve eigen-problem convergence at the\n"
			"expense of increased memory requirements.";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.cntrl.davidsonBandRatio, 1.1, "ratio");
		if(e.cntrl.davidsonBandRatio < 1.)
			throw string("<ratio> must be at least 1");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg", e.cntrl.davidsonBandRatio);
	}
}
commandDavidsonBandRatio;

//-------------------------------------------------------------------------------------------------

struct CommandLcaoParams : public Command
{
	CommandLcaoParams() : Command("lcao-params", "jdftx/Initialization")
	{
		format = "[<nIter>=-1] [<Ediff>=1e-6] [<smeaingWidth>=1e-3]";
		comments = "Control LCAO wavefunction initialization:\n"
			"+ <nIter>: maximum subspace iterations in LCAO (negative => auto-select)\n"
			"+ <Ediff>: energy-difference convergence threshold for subspace iteration\n"
			"+ <smearingWidth>: smearing width for the subspace iteration for constant fillings calculations.\n"
			"   If present, the smearing width from elec-smearing overrides this.\n";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.eVars.lcaoIter, -1, "nIter");
		pl.get(e.eVars.lcaoTol, 1e-6, "Ediff");
		pl.get(e.eInfo.smearingWidth, 1e-3, "smeaingWidth");
		if(e.eInfo.smearingWidth<=0) throw string("<kT> must be positive.\n");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%d %lg %lg", e.eVars.lcaoIter, e.eVars.lcaoTol, e.eInfo.smearingWidth);
	}
}
commandLcaoParams;

//-------------------------------------------------------------------------------------------------

EnumStringMap<SpinType> spinMap
(	SpinNone,   "no-spin",
	SpinZ,      "z-spin",
	SpinVector, "vector-spin",
	SpinOrbit,  "spin-orbit"
);

EnumStringMap<SpinType> spinDescMap
(	SpinNone,   "Unpolarized calculation",
	SpinZ,      "Spin-polarized calculation",
	SpinVector, "Non-collinear magnetism (supports spin-orbit coupling)",
	SpinOrbit,  "Non-collinear without magnetization, to allow for spin-orbit"
);

struct CommandSpintype : public Command
{
	CommandSpintype() : Command("spintype", "jdftx/Electronic/Parameters")
	{
		format = "<type>=" + spinMap.optionList();
		comments = "Select spin-polarization type:"
			+ addDescriptions(spinMap.optionList(), linkDescription(spinMap, spinDescMap));
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.eInfo.spinType, SpinNone, spinMap, "type");
	}

	void printStatus(Everything& e, int iRep)
	{	fputs(spinMap.getString(e.eInfo.spinType), globalLog);
	}
}
commandSpintype;

//-------------------------------------------------------------------------------------------------

//Base class for fix-electron-density and fix-electron-potential
struct CommandFixElectronHamiltonian : public Command
{	
    CommandFixElectronHamiltonian(string name) : Command("fix-electron-"+name, "jdftx/Electronic/Optimization")
	{
		format = "<filenamePattern>";
		comments = "Perform band structure calculations at fixed electron " + name + "\n"
			"(or spin " + name + ") read from the specified <filenamePattern>.\n"
			"This pattern must include $VAR which will be replaced by the appropriate\n"
			"variable names accounting for spin-polarization (same as used for dump).\n"
			"Meta-GGA calculations will also require the corresponding kinetic " + name + ".";
		
		require("spintype");
		forbid("elec-smearing");
		forbid("elec-ex-corr-compare");
		forbid("electronic-scf");
		forbid("vibrations");
		forbid("dump-only");
	}

	void processCommon(ParamList& pl, Everything& e, string& targetFilenamePattern)
	{	pl.get(targetFilenamePattern, string(), "filenamePattern", true);
		if(targetFilenamePattern.find("$VAR") == string::npos)
			throw string("<filenamePattern> must contain $VAR");
		e.cntrl.fixed_H = true;
	}

	void printStatusCommon(Everything& e, int iRep, const string& targetFilenamePattern)
	{	logPrintf("%s", targetFilenamePattern.c_str());
	}
};

struct CommandFixElectronDensity : public CommandFixElectronHamiltonian
{   CommandFixElectronDensity() : CommandFixElectronHamiltonian("density") { forbid("fix-electron-potential"); }
	void process(ParamList& pl, Everything& e) { processCommon(pl, e, e.eVars.nFilenamePattern); }
	void printStatus(Everything& e, int iRep) { printStatusCommon(e, iRep, e.eVars.nFilenamePattern); }
}
commandFixElectronDensity;

struct CommandFixElectronPotential : public CommandFixElectronHamiltonian
{   CommandFixElectronPotential() : CommandFixElectronHamiltonian("potential") { forbid("fix-electron-density"); }
	void process(ParamList& pl, Everything& e) { processCommon(pl, e, e.eVars.VFilenamePattern); }
	void printStatus(Everything& e, int iRep) { printStatusCommon(e, iRep, e.eVars.VFilenamePattern); }
}
commandFixElectronPotential;

//-------------------------------------------------------------------------------------------------

struct CommandConvergeEmptyStates : public Command
{
	CommandConvergeEmptyStates() : Command("converge-empty-states", "jdftx/Electronic/Optimization")
	{
		format = "yes|no";
		comments = "Whether to converge empty states after each electronic optimization (default no).\n"
			"Not required unless empty states are used in post-processing and need to be accurate.\n"
			"This is a shortcut to running a bandstructure calculation after a total energy\n"
			"or SCF calculation, and helps simplify workflow for DOS, polarizability, wannier\n"
			"and electron-phonon matrix element calculations.";
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.cntrl.convergeEmptyStates, false, boolMap, "shouldConverge", true);
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s", boolMap.getString(e.cntrl.convergeEmptyStates));
	}
}
commandConvergeEmptyStates;

//-------------------------------------------------------------------------------------------------

struct CommandWavefunctionDrag : public Command
{
	CommandWavefunctionDrag() : Command("wavefunction-drag", "jdftx/Ionic/Optimization")
	{
		format = "yes|no";
		comments =
			"Drag wavefunctions when ions are moved using atomic orbital projections (yes by default).";
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.cntrl.dragWavefunctions, true, boolMap, "shouldDrag", true);
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s", boolMap.getString(e.cntrl.dragWavefunctions));
	}
}
commandWavefunctionDrag;

//-------------------------------------------------------------------------------------------------

struct CommandCacheProjectors : public Command
{
	CommandCacheProjectors() : Command("cache-projectors", "jdftx/Miscellaneous")
	{
		format = "yes|no";
		comments =
			"Cache nonlocal-pseudopotential projectors (yes by default); turn off to save memory.";
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.cntrl.cacheProjectors, true, boolMap, "shouldCache", true);
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s", boolMap.getString(e.cntrl.cacheProjectors));
	}
}
commandCacheProjectors;

//-------------------------------------------------------------------------------------------------

struct CommandBasis : public Command
{
	CommandBasis() : Command("basis", "jdftx/Electronic/Parameters")
	{
		format = "<kdep>=" + kdepMap.optionList();
		comments = "Basis set at each k-point (default), or single basis set at gamma point";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.cntrl.basisKdep, BasisKpointDep, kdepMap, "kdep");
	}

	void printStatus(Everything& e, int iRep)
	{	fputs(kdepMap.getString(e.cntrl.basisKdep), globalLog);
	}
}
commandBasis;


//-------------------------------------------------------------------------------------------------

static EnumStringMap<ElecEigenAlgo> elecEigenMap(ElecEigenCG, "CG", ElecEigenDavidson, "Davidson");

struct CommandElecEigenAlgo : public Command
{
    CommandElecEigenAlgo() : Command("elec-eigen-algo", "jdftx/Electronic/Optimization")
	{
		format = "<algo>=" + elecEigenMap.optionList();
		comments = "Selects eigenvalue algorithm for band-structure calculations or inner loop of SCF.";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.cntrl.elecEigenAlgo, ElecEigenDavidson, elecEigenMap, "algo");
	}

	void printStatus(Everything& e, int iRep)
	{	fputs(elecEigenMap.getString(e.cntrl.elecEigenAlgo), globalLog);
	}
}
commandElecEigenAlgo;

//-------------------------------------------------------------------------------------------------

struct CommandRhoExternal : public Command
{
	CommandRhoExternal() : Command("rhoExternal", "jdftx/Coulomb interactions")
	{
		format = "<filename> [<includeSelfEnergy>=yes|no]";
		comments =
			"Include an external charge density [electrons/bohr^3] (real space binary)\n"
			"which interacts electrostatically with the electrons, nuclei and fluid.\n"
			"\n"
			"If <includeSelfEnergy>=yes (default no), then the Coulomb self-energy\n"
			"of rhoExternal is included in the output energy.";
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.eVars.rhoExternalFilename, string(), "filename", true);
		pl.get(e.eVars.rhoExternalSelfEnergy, false, boolMap, "includeSelfEnergy");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s %s", e.eVars.rhoExternalFilename.c_str(), boolMap.getString(e.eVars.rhoExternalSelfEnergy));
	}
}
commandRhoExternal;

//-------------------------------------------------------------------------------------------------

struct CommandVexternal : public Command
{
	CommandVexternal() : Command("Vexternal", "jdftx/Electronic/Parameters")
	{
		format = "<filename> | <filenameUp> <filenameDn>";
		comments =
			"Include an external potential (in hartrees) for the electrons\n"
			"(real space binary). Specify two files if V is spin-polarized.";
	}

	void process(ParamList& pl, Everything& e)
	{	e.eVars.VexternalFilename.resize(1);
		pl.get(e.eVars.VexternalFilename[0], string(), "filename", true);
		//Check if a second file has been specified:
		string filenameDn;
		pl.get(filenameDn, string(), "filenameDn");
		if(filenameDn.length())
			e.eVars.VexternalFilename.push_back(filenameDn);
	}

	void printStatus(Everything& e, int iRep)
	{	for(const string& filename: e.eVars.VexternalFilename)
			logPrintf("%s ", filename.c_str());
	}
}
commandVexternal;

//-------------------------------------------------------------------------------------------------

struct CommandBoxPotential : public Command
{
	CommandBoxPotential() : Command("box-potential", "jdftx/Electronic/Parameters")
	{
		format = "xmin xmax ymin ymax zmin zmax Vin Vout [<convolve_radius>=0.1]";
		comments =
			"Include an step-function shaped external potential (in hartrees) for the electrons";
	
		allowMultiple = true;
	}

	void process(ParamList& pl, Everything& e)
	{	ElecVars::BoxPotential bP;
		const char* dirNames[3] = { "x", "y", "z" };
		for(int k=0; k<3; k++)
		{	pl.get(bP.min[k], 0., dirNames[k]+string("min"), true);
			pl.get(bP.max[k], 0., dirNames[k]+string("max"), true);
			if(bP.max[k]<bP.min[k])
				throw(string("max must be smaller than min for each dimension"));
		}
		pl.get(bP.Vin, 0., "Vin", true);
		pl.get(bP.Vout, 0., "Vout", true);
		pl.get(bP.convolve_radius, 0.2, "convolve_radius", false);
		e.eVars.boxPot.push_back(bP);
	}

	void printStatus(Everything& e, int iRep)
	{	const ElecVars::BoxPotential& bP = e.eVars.boxPot[iRep];
		logPrintf("%.5g %.5g %.5g %.5g %.5g %.5g    %.5g %.5g  %.5g",
			 bP.min[0], bP.max[0], bP.min[1], bP.max[1], bP.min[2], bP.max[2],
			 bP.Vin, bP.Vout, bP.convolve_radius);
	}
}
commandBoxPot;

//-------------------------------------------------------------------------------------------------

struct CommandElectricField : public Command
{
	CommandElectricField() : Command("electric-field", "jdftx/Electronic/Parameters")
	{
		format = "<Ex> <Ey> <Ez>";
		comments = 
			"Apply an external electric field (in Cartesian coordinates, atomic\n"
			"units [Eh/e/a_0] and electron-is-positive sign convention).\n\n"
			"In truncated directions, the field will be applied as a ramp potential,\n"
			"and for periodic directions, it will be applied as a plane wave with\n"
			"the smallest commensurate wave vector and amplitude set by peak field.\n\n"
			"Coulomb truncation, if any, must be in embedded mode (command coulomb-truncation-embed).\n"
			"Symmetries will be automatically reduced to account for this field.";
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.coulombParams.Efield[0], 0., "Ex", true);
		pl.get(e.coulombParams.Efield[1], 0., "Ey", true);
		pl.get(e.coulombParams.Efield[2], 0., "Ez", true);
	}

	void printStatus(Everything& e, int iRep)
	{	for(int k=0; k<3; k++) logPrintf("%lg ", e.coulombParams.Efield[k]);
	}
}
commandElectricField;
