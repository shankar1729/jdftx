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

#include <commands/minimize.h>
#include <wannier/Wannier.h>
#include <core/LatticeUtils.h>

enum WannierMember
{	WM_localizationMeasure,
	WM_precond,
	WM_bStart,
	WM_outerWindow,
	WM_innerWindow,
	WM_frozenCenters,
	WM_saveWfns,
	WM_saveWfnsRealSpace,
	WM_saveMomenta,
	WM_loadRotations,
	WM_numericalOrbitals,
	WM_numericalOrbitalsOffset,
	WM_phononSup,
	WM_delim
};

EnumStringMap<WannierMember> wannierMemberMap
(	WM_localizationMeasure, "localizationMeasure",
	WM_precond, "precondition",
	WM_bStart, "bStart",
	WM_outerWindow, "outerWindow",
	WM_innerWindow, "innerWindow",
	WM_frozenCenters, "frozenCenters",
	WM_saveWfns, "saveWfns",
	WM_saveWfnsRealSpace, "saveWfnsRealSpace",
	WM_saveMomenta, "saveMomenta",
	WM_loadRotations, "loadRotations",
	WM_numericalOrbitals, "numericalOrbitals",
	WM_numericalOrbitalsOffset, "numericalOrbitalsOffset",
	WM_phononSup, "phononSupercell"
);

EnumStringMap<Wannier::LocalizationMeasure> localizationMeasureMap
(	Wannier::LM_FiniteDifference, "FiniteDifference",
	Wannier::LM_RealSpace, "RealSpace"
);

struct CommandWannier : public Command
{
	CommandWannier() : Command("wannier", "wannier")
	{
		format = "<key1> <args1...>  <key2> <args2...>  ...";
		comments =
			"Control wannier calculation and output. The possible <key>'s and their\n"
			"corresponding arguments are:\n"
			"\n+ localizationMeasure FiniteDifference | RealSpace\n\n"
			"   Controls how the localization of the Wannier functions is calculated.\n"
			"   The finite-difference reciprocal space measure of Marzari and Vanderbilt\n"
			"   (default) is inexpensive, but its error scales as Nkpoints^(-2./3).\n"
			"   The real space version is slower but its error scales as exp(-Nkpoints^(1/3)),\n"
			"   and is preferable for quantitative applications. Note that the real-space version\n"
			"   is not translationally invariant and wraps on a superlattice Wigner-Seitz cell\n"
			"   centered at the origin.\n"
			"\n+ precondition yes|no\n\n"
			"   Whether to use an inverse-Hemholtz preconditioner for the minimization.\n"
			"   Affects only the FiniteDifference localizationMeasure. (default: no)\n"
			"\n+ bStart <band>\n\n"
			"   For fixed band calculations, 0-based index of lowest band used.\n"
			"   The number of bands equals the number of wannier-centers specified.\n"
			"   Default: 0. Use in insulator calculations to ignore semi-core orbitals.\n"
			"\n+ outerWindow <eMin> <eMax>\n\n"
			"   Energy window within which bands contribute to wannier functions\n"
			"   bStart is ignored if outerWindow is specified.\n"
			"\n+ innerWindow <eMin> <eMax>\n\n"
			"   Inner energy window within which bands are used exactly.\n"
			"   Requires outerWindow, and innerWindow must be its subset.\n"
			"\n+ frozenCenters <nFrozen> <filename>\n\n"
			"   Include frozen Wannier centers imported as unitary rotations from <filename>\n"
			"   (.mlwfU output from another wannier run on the same jdftx state), and force\n"
			"   the current centers to be orthogonal to the imported ones. Note that all output\n"
			"   of this run will include the frozen as well as new centers for convenience.\n"
			"   \n"
			"   Requires innerWindow if outerWindow is present, and the frozen centers\n"
			"   must lie entirely within the inner window or fixed band range. (Effectively,\n" 
			"   the outer window used in generating the frozen centers must be a subset of\n"
			"   the present inner window / fixed band range.)\n"
			"\n+ saveWfns yes|no\n\n"
			"   Whether to write supercell wavefunctions in a spherical reciprocal-space basis.\n"
			"   Default: no.\n"
			"\n+ saveWfnsRealSpace yes|no\n\n"
			"   Whether to write supercell wavefunctions band-by-band in real space (can be enormous).\n"
			"   Default: no.\n"
			"\n+ saveMomenta yes|no\n\n"
			"   Whether to write momentum matrix elements in the same format as Hamiltonian.\n"
			"   The output is real and antisymmetric (drops the iota so as to half the output size).\n"
			"   Default: no.\n"
			"\n+ loadRotations yes|no\n\n"
			"   Whether to load rotations (.mlwU and .mlwfU2) from a previous Wannier run.\n"
			"   Default: no.\n"
			"\n+ numericalOrbitals <filename>\n\n"
			"   Load numerical orbitals from <filename> with basis described in <filename>.header\n"
			"   that can then be used as trial orbitals. The reciprocal space wavefunction output\n"
			"   of wannier is precisely in this format, so that supercell Wannier calculations can\n"
			"   be initialized from previously calculated bulk / unit cell Wannier functions.\n"
			"\n+ numericalOrbitalsOffset <x0> <x1> <x2>\n\n"
			"   Origin of the numerical orbitals in lattice coordinates of the input.\n"
			"   The default [ .5 .5 .5 ] (cell center) is appropriate for output from wannier.\n"
			"\n+ phononSupercell <N0> <N1> <N2>\n\n"
			"   If specified, wannier will read in output (phononHsub and phononBasis) of the\n"
			"   phonon code with this supercell and output Wannierized electron-phonon matrix\n"
			"   elements (mlwfHePh). This file will contain matrices for pairs of cells in the\n"
			"   order specified in the mlwfCellMapSqPh output file, with an outer loop over\n"
			"   the nuclear displacement modes in phononBasis.";
	}

	void process(ParamList& pl, Everything& e)
	{	Wannier& wannier = ((WannierEverything&)e).wannier;
		while(true)
		{	WannierMember key; pl.get(key, WM_delim, wannierMemberMap, "key");
			if(key==WM_delim) break;
			switch(key)
			{	case WM_localizationMeasure:
					pl.get(wannier.localizationMeasure, Wannier::LM_FiniteDifference,  localizationMeasureMap, "localizationMeasure", true);
					break;
				case WM_precond:
					pl.get(wannier.precond, false, boolMap, "precondition", true);
					break;
				case WM_bStart:
					pl.get(wannier.bStart, 0, "band", true);
					break;
				case WM_outerWindow:
					pl.get(wannier.eOuterMin, 0., "eMin", true);
					pl.get(wannier.eOuterMax, 0., "eMax", true);
					wannier.outerWindow = true;
					break;
				case WM_innerWindow:
					pl.get(wannier.eInnerMin, 0., "eMin", true);
					pl.get(wannier.eInnerMax, 0., "eMax", true);
					wannier.innerWindow = true;
					break;
				case WM_frozenCenters:
					pl.get(wannier.nFrozen, 0, "nFrozen", true);
					pl.get(wannier.frozenUfilename, string(), "filename", true);
					if(wannier.nFrozen <= 0) throw string("nFrozen must be > 0");
					break;
				case WM_saveWfns:
					pl.get(wannier.saveWfns, false, boolMap, "saveWfns", true);
					break;
				case WM_saveWfnsRealSpace:
					pl.get(wannier.saveWfnsRealSpace, false, boolMap, "saveWfnsRealSpace", true);
					break;
				case WM_saveMomenta:
					pl.get(wannier.saveMomenta, false, boolMap, "saveMomenta", true);
					break;
				case WM_loadRotations:
					pl.get(wannier.loadRotations, false, boolMap, "loadRotations", true);
					break;
				case WM_numericalOrbitals:
					pl.get(wannier.numericalOrbitalsFilename, string(), "filename", true);
					break;
				case WM_numericalOrbitalsOffset:
					pl.get(wannier.numericalOrbitalsOffset[0], 0., "x0", true);
					pl.get(wannier.numericalOrbitalsOffset[1], 0., "x1", true);
					pl.get(wannier.numericalOrbitalsOffset[2], 0., "x2", true);
					break;
				case WM_phononSup:
					pl.get(wannier.phononSup[0], 0, "N0", true);
					pl.get(wannier.phononSup[1], 0, "N1", true);
					pl.get(wannier.phononSup[2], 0, "N2", true);
					break;
				case WM_delim: //should never be encountered
					break;
			}
		}
	}

	void printStatus(Everything& e, int iRep)
	{	const Wannier& wannier = ((const WannierEverything&)e).wannier;
		logPrintf(" \\\n\tlocalizationMeasure %s", localizationMeasureMap.getString(wannier.localizationMeasure));
		logPrintf(" \\\n\tprecondition %s", boolMap.getString(wannier.precond));
		logPrintf(" \\\n\tsaveWfns %s", boolMap.getString(wannier.saveWfns));
		logPrintf(" \\\n\tsaveWfnsRealSpace %s", boolMap.getString(wannier.saveWfnsRealSpace));
		logPrintf(" \\\n\tsaveMomenta %s", boolMap.getString(wannier.saveMomenta));
		logPrintf(" \\\n\tloadRotations %s", boolMap.getString(wannier.loadRotations));
		if(wannier.outerWindow)
		{	logPrintf(" \\\n\touterWindow %lg %lg", wannier.eOuterMin, wannier.eOuterMax);
			if(wannier.innerWindow)
				logPrintf(" \\\n\tinnerWindow %lg %lg", wannier.eInnerMin, wannier.eInnerMax);
		}
		else
		{	logPrintf(" \\\n\tbStart %d", wannier.bStart);
		}
		if(wannier.nFrozen)
			logPrintf(" \\\n\tfrozenCenters %d %s", wannier.nFrozen, wannier.frozenUfilename.c_str());
		if(wannier.numericalOrbitalsFilename.length())
		{	logPrintf(" \\\n\tnumericalOrbitals %s", wannier.numericalOrbitalsFilename.c_str());
			const vector3<>& offs = wannier.numericalOrbitalsOffset;
			logPrintf(" \\\n\tnumericalOrbitalsOffset %lg %lg %lg", offs[0], offs[1], offs[2]);
		}
		if(wannier.phononSup.length_squared())
			logPrintf(" \\\n\tphononSupercell %d %d %d", wannier.phononSup[0], wannier.phononSup[1], wannier.phononSup[2]);
	}
}
commandWannier;


struct CommandWannierCenter : public Command
{
	CommandWannierCenter(string suffix=string()) : Command("wannier-center" + suffix, "wannier")
	{
		format = "<aorb1> [<aorb2> ...]";
		comments =
			"Specify trial orbital for a wannier function as a linear combination of\n"
			"atomic orbitals <aorb*>. The syntax for each <aorb> is:\n"
			"\n"
			"   <x0> <x1> <x2> [<a>=1|spName] [<orbDesc>=s] [<coeff>=1.0]\n"
			"\n"
			"which represents an atomic orbital centered at <x0>,<x1>,<x2>\n"
			"(coordinate system set by coords-type).\n"
			"\n"
			"If <a> is the name of one of the pseudopotentials, then its\n"
			"atomic orbitals will be used; otherwise a hydogenic orbital\n"
			"of decay length n*<a> bohrs, where n is the principal\n"
			"quantum number in <orbDesc> will be used.\n"
			"\n"
			"The orbital code <orbDesc> is as in command density-of-states.\n"
			"\n"
			"Specify <a>, <orbDesc> and <coeff> explicitly when using multiple\n"
			"orbitals; the defaults only apply to the single orbital case.\n"
			"\n"
			"Alternately, for using numerical trial orbitals that have been\n"
			"input using command wannier, use the syntax:\n"
			"\n"
			"   <x0> <x1> <x2> numerical <b> [<coeff>=1.0]\n"
			"\n"
			"where <b> is the 0-based index of the input orbital, and <coeff> may\n"
			"be used to linearly combine numerical orbitals with other numerical\n"
			"or atomic / hydrogenic orbitals as specified above.\n"
			"\n"
			"Specify this command once for each Wannier function.";
		
		allowMultiple = true;
		require("wannier-initial-state");
		require("wannier-dump-name");
		
		//Dependencies due to optional conversion from cartesian coords:
		require("latt-scale");
		require("coords-type");
		require("spintype");
	}

	void process(ParamList& pl, Everything& e)
	{	Wannier& wannier = ((WannierEverything&)e).wannier;
		Wannier::TrialOrbital t;
		while(true)
		{	Wannier::AtomicOrbital ao; string orbDesc;
			pl.get(ao.r[0], nan(""), "x0"); if(std::isnan(ao.r[0])) break;
			pl.get(ao.r[1], 0., "x1", true);
			pl.get(ao.r[2], 0., "x2", true);
			//Transform to lattice coordinates if necessary:
			if(e.iInfo.coordsType == CoordsCartesian)
				ao.r = inv(e.gInfo.R)*ao.r;
			//Determine trial orbital type:
			string aKey; pl.get(aKey, string(), "a");
			ao.numericalOrbIndex = -1;
			ao.sp = -1;
			ao.atom = -1;
			if(aKey == "numerical")
			{	pl.get(ao.numericalOrbIndex, -1, "b", true);
				if(ao.numericalOrbIndex<0) throw(string("Numerical orbital index must be non-negative"));
			}
			else
			{	for(int sp=0; sp<int(e.iInfo.species.size()); sp++)
					if(e.iInfo.species[sp]->name == aKey)
					{	ao.sp = sp;
						//Match to an atom if possible:
						const std::vector< vector3<> >& atpos = e.iInfo.species[sp]->atpos;
						for(size_t atom=0; atom<atpos.size(); atom++)
							if(circDistanceSquared(ao.r, atpos[atom]) < symmThresholdSq)
							{	ao.atom = atom;
								wannier.needAtomicOrbitals = true;
							}
						break;
					}
				if(ao.sp<0)
					ParamList(aKey).get(ao.a, 1., "a");
				pl.get(orbDesc, string(), "orbDesc");
				if(orbDesc.length())
				{	ao.orbitalDesc.parse(orbDesc);
					if(ao.orbitalDesc.m > ao.orbitalDesc.l)
						throw(string("Must specify a specific projection eg. px,py (not just p)"));
					if(ao.sp>=0 && ao.orbitalDesc.n + ao.orbitalDesc.l > 3)
						throw(string("Hydrogenic orbitals with n+l>4 not supported"));
				}
				else //default is nodeless s orbital
				{	ao.orbitalDesc.l = 0;
					ao.orbitalDesc.m = 0;
					ao.orbitalDesc.n = 0;
					ao.orbitalDesc.s = 0;
					ao.orbitalDesc.spinType = SpinNone;
				}
				if((e.eInfo.spinType==SpinOrbit || e.eInfo.spinType==SpinVector) && ao.orbitalDesc.spinType==SpinNone)
					throw string("Must specify an explicit spin projection for noncollinear modes");
				if((e.eInfo.spinType==SpinNone || e.eInfo.spinType==SpinZ) && ao.orbitalDesc.spinType!=SpinNone)
					throw string("Orbital projections must not specify spin for collinear modes");
				if(ao.orbitalDesc.spinType==SpinOrbit && ao.atom<0)
					throw string("Relativistic (l,j,mj) orbital projections must be centered on an atom");
			}
			pl.get(ao.coeff, 1., "coeff");
			t.push_back(ao);
		}
		if(!t.size()) throw(string("Trial orbital for each center must contain at least one atomic orbital"));
		wannier.trialOrbitals.push_back(t);
	}

	void printStatus(Everything& e, int iRep)
	{	const Wannier& wannier = ((const WannierEverything&)e).wannier;
		const Wannier::TrialOrbital& t = wannier.trialOrbitals[iRep];
		for(const Wannier::AtomicOrbital& ao: t)
		{	if(t.size()>1) logPrintf(" \\\n\t");
			vector3<> r = ao.r;
			if(e.iInfo.coordsType == CoordsCartesian)
				r = e.gInfo.R * r; //report cartesian positions
			logPrintf("%lg %lg %lg ", r[0], r[1], r[2]);
			if(ao.numericalOrbIndex >= 0)
				logPrintf("numerical %d", ao.numericalOrbIndex);
			else
			{	if(ao.sp>=0) logPrintf("%s", e.iInfo.species[ao.sp]->name.c_str());
				else logPrintf("%lg", ao.a);
				logPrintf(" %s", string(ao.orbitalDesc).c_str());
			}
			logPrintf("  %lg", ao.coeff);
		}
	}
}
commandWannierCenter;


struct CommandWannierCenterPinned : public CommandWannierCenter
{
	CommandWannierCenterPinned(string suffix=string()) : CommandWannierCenter("-pinned")
	{	comments =
			"Same as command wannier-center, except that the minimizer tries to\n"
			"center this orbital at the specified guess positions (weighted mean\n"
			"of the orbital centers, if usinga  linear combination). This can help\n"
			"improve the minimizer conditioning if the centers are known by symmetry.";
	}
	
	void process(ParamList& pl, Everything& e)
	{	Wannier& wannier = ((WannierEverything&)e).wannier;
		CommandWannierCenter::process(pl, e);
		Wannier::TrialOrbital& t = wannier.trialOrbitals.back();
		t.pinned = true;
		vector3<> r0 = t.front().r;
		vector3<> rSum; double wSum = 0.;
		for(const Wannier::AtomicOrbital& ao: t)
		{	//find position with minimum image convention from first position:
			vector3<> dr = ao.r - r0;
			for(int j=0; j<3; j++)
				dr[j] = floor(0.5 + dr[j]);
			vector3<> r = r0 + dr;
			//collect with weights:
			double weight = std::pow(ao.coeff, 2);
			rSum += weight * r;
			wSum += weight;
		}
		t.rCenter = (1./wSum) * rSum;
	}
}
commandWannierCenterPinned;


struct CommandWannierMinimize : public CommandMinimize
{	CommandWannierMinimize() : CommandMinimize("wannier", "wannier") {}
    MinimizeParams& target(Everything& e)
	{	Wannier& wannier = ((WannierEverything&)e).wannier;
		return wannier.minParams;
	}
    void process(ParamList& pl, Everything& e)
	{	Wannier& wannier = ((WannierEverything&)e).wannier;
		wannier.minParams.energyDiffThreshold = 1e-8; //override default value (0.) in MinimizeParams.h
		CommandMinimize::process(pl, e);
	}
}
commandWannierMinimize;


struct CommandWannierFilenames : public Command
{	virtual string& getTarget(Everything&)=0; //derived class determines where to save the file
	
	CommandWannierFilenames(string cmdSuffix, string comment) : Command("wannier-"+cmdSuffix, "wannier")
	{	format = "<format>";
		comments = 
			"Control the filename pattern for wannier " + comment + ", where <format> must contain\n"
			"a single instance of $VAR, which will be substituted by the name of each variable.";
	}

	void process(ParamList& pl, Everything& e)
	{	string& target = getTarget(e);
		pl.get(target, string("$STAMP.$VAR"), "format");
		if(target.find("$VAR")==string::npos)
			throw "<format> = " + target + " doesn't contain the pattern $VAR";
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s", getTarget(e).c_str());
	}
};

struct CommandWannierInitialState : public CommandWannierFilenames
{	CommandWannierInitialState() : CommandWannierFilenames("initial-state", "state initialization") {}
	string& getTarget(Everything& e) { return ((WannierEverything&)e).wannier.initFilename; }
}
commandWannierInitialState;

struct CommandWannierDumpName : public CommandWannierFilenames
{	CommandWannierDumpName() : CommandWannierFilenames("dump-name", "output") {}
	string& getTarget(Everything& e) { return ((WannierEverything&)e).wannier.dumpFilename; }
}
commandWannierDumpName;
