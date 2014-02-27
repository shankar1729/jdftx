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

//Wannier-specific commands affect this static object:
Wannier wannier;

enum WannierMember
{	WM_localizationMeasure,
	WM_bStart,
	WM_outerWindow,
	WM_innerWindow,
	WM_saveWfns,
	WM_saveWfnsRealSpace,
	WM_saveMomenta,
	WM_loadRotations,
	WM_numericalOrbitals,
	WM_numericalOrbitalsOffset,
	WM_delim
};

EnumStringMap<WannierMember> wannierMemberMap
(	WM_localizationMeasure, "localizationMeasure",
	WM_bStart, "bStart",
	WM_outerWindow, "outerWindow",
	WM_innerWindow, "innerWindow",
	WM_saveWfns, "saveWfns",
	WM_saveWfnsRealSpace, "saveWfnsRealSpace",
	WM_saveMomenta, "saveMomenta",
	WM_loadRotations, "loadRotations",
	WM_numericalOrbitals, "numericalOrbitals",
	WM_numericalOrbitalsOffset, "numericalOrbitalsOffset"
);

EnumStringMap<Wannier::LocalizationMeasure> localizationMeasureMap
(	Wannier::LM_FiniteDifference, "FiniteDifference",
	Wannier::LM_RealSpace, "RealSpace"
);

struct CommandWannier : public Command
{
	CommandWannier() : Command("wannier")
	{
		format = "<key1> <args1...>  <key2> <args2...>  ...";
		comments =
			"Control wannier calculation and output. The possible <key>'s and their\n"
			"corresponding arguments are:\n"
			"  localizationMeasure FiniteDifference | RealSpace\n"
			"    Controls how the localization of the Wannier functions is calculated.\n"
			"    The finite-difference reciprocal space measure of Marzari and Vanderbilt\n"
			"    (default) is inexpensive, but its error scales as Nkpoints^(-2./3).\n"
			"    The real space version is slower but its error scales as exp(-Nkpoints^(1/3)),\n"
			"    and is preferable for quantitative applications. Note that the real-space version\n"
			"    is not translationally invariant and wraps on a superlattice WIgner-Seitz cell\n"
			"    centered at the origin.\n"
			"  bStart <band>\n"
			"    For fixed band calculations, 0-based index of lowest band used.\n"
			"    The number of bands equals the number of wannier-centers specified.\n"
			"    Default: 0. Use in insulator calculations to ignore semi-core orbitals.\n"
			"  outerWindow <eMin> <eMax>\n"
			"    Energy window within which bands contribute to wannier functions\n"
			"    bStart is ignored if outerWindow is specified.\n"
			"  innerWindow <eMin> <eMax>\n"
			"    Inner energy window within which bands are used exactly.\n"
			"    Outer energy window must be specified to use this.\n"
			"  saveWfns yes|no\n"
			"    Whether to write supercell wavefunctions in a spherical reciprocal-space basis.\n"
			"    Default: no.\n"
			"  saveWfnsRealSpace yes|no\n"
			"    Whether to write supercell wavefunctions band-by-band in real space (can be enormous).\n"
			"    Default: no.\n"
			"  saveMomenta yes|no\n"
			"    Whether to write momentum matrix elements in the same format as Hamiltonian.\n"
			"    The output is real and antisymmetric (drops the iota so as to half the output size).\n"
			"    Default: no.\n"
			"  loadRotations yes|no\n"
			"    Whether to load rotations (.mlwU and .mlwfU2) from a previous Wannier run.\n"
			"    Default: no.\n"
			"  numericalOrbitals <filename>\n"
			"    Load numerical orbitals from <filename> with basis described in <filename>.header\n"
			"    that can then be used as trial orbitals. The reciprocal space wavefunction output\n"
			"    of wannier is precisely in this format, so that supercell Wannier calculations can\n"
			"    be initialized from previously calculated bulk / unit cell Wannier functions.\n"
			"  numericalOrbitalsOffset <x0> <x1> <x2>\n"
			"    Origin of the numerical orbitals in lattice coordinates of the input.\n"
			"    The default [ .5 .5 .5 ] (cell center) is appropriate for output from wannier.";
	}

	void process(ParamList& pl, Everything& e)
	{	while(true)
		{	WannierMember key; pl.get(key, WM_delim, wannierMemberMap, "key");
			if(key==WM_delim) break;
			switch(key)
			{	case WM_localizationMeasure:
					pl.get(wannier.localizationMeasure, Wannier::LM_FiniteDifference,  localizationMeasureMap, "localizationMeasure", true);
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
				case WM_delim: //should never be encountered
					break;
			}
		}
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf(" \\\n\tlocalizationMeasure %s", localizationMeasureMap.getString(wannier.localizationMeasure));
		logPrintf(" \\\n\tsaveWfns %s", boolMap.getString(wannier.saveWfns));
		logPrintf(" \\\n\tsaveWfnsRealSpace %s", boolMap.getString(wannier.saveWfnsRealSpace));
		logPrintf(" \\\n\tsaveMomenta %s", boolMap.getString(wannier.saveMomenta));
		logPrintf(" \\\n\tloadRotations %s", boolMap.getString(wannier.loadRotations));
		if(wannier.outerWindow) logPrintf(" \\\n\touterWindow %lg %lg", wannier.eOuterMin, wannier.eOuterMax);
		if(wannier.innerWindow) logPrintf(" \\\n\tinnerWindow %lg %lg", wannier.eInnerMin, wannier.eInnerMax);
		if(!(wannier.innerWindow || wannier.outerWindow))
			logPrintf(" \\\n\tbStart %d", wannier.bStart);
		if(wannier.numericalOrbitalsFilename.length())
		{	logPrintf(" \\\n\tnumericalOrbitals %s", wannier.numericalOrbitalsFilename.c_str());
			const vector3<>& offs = wannier.numericalOrbitalsOffset;
			logPrintf(" \\\n\tnumericalOrbitalsOffset %lg %lg %lg", offs[0], offs[1], offs[2]);
		}
	}
}
commandWannier;


struct CommandWannierCenter : public Command
{
	CommandWannierCenter() : Command("wannier-center")
	{
		format = "<aorb1> [<aorb2> ...]";
		comments =
			"Specify trial orbital for a wannier function as a linear combination of\n"
			"atomic orbitals <aorb*>. The syntax for each <horb> is:\n"
			"   <x0> <x1> <x2> [<a>=1|spName] [<orbDesc>=s] [<coeff>=1.0]\n"
			"which represents an atomic orbital centered at <x0>,<x1>,<x2>\n"
			"(coordinate system set by coords-type). If <a> is the name of one\n"
			"of the pseudopotentials, then its atomic orbitals will be used;\n"
			"otherwise a hydogenic orbital of decay length n*<a> bohrs, where\n"
			"n is the principal quantum number in <orbDesc> will be used.\n"
			"The orbital code <orbDesc> is as in command density-of-states.\n"
			"Specify a, orbDesc and coeff explicitly when using multiple\n"
			"orbitals; the defaults only apply to the single orbital case.\n"
			"   Alternately, for using numerical trial orbitals that have been\n"
			"input using command wannier, use the syntax:\n"
			"   <x0> <x1> <x2> numerical <b> [<coeff>=1.0]\n"
			"where <b> is the 0-based index of the input orbital, and coeff may\n"
			"be used to linearly combine numerical orbitals with other numerical\n"
			"or atomic / hydrogenic orbitals as specified above.\n"
			"   Specify this command once for each Wannier function.";
		allowMultiple = true;
		require("wannier-initial-state");
		require("wannier-dump-name");
		
		//Dependencies due to optional conversion from cartesian coords:
		require("latt-scale");
		require("coords-type");
	}

	void process(ParamList& pl, Everything& e)
	{	Wannier::TrialOrbital t;
		while(true)
		{	Wannier::AtomicOrbital ao; string orbDesc;
			pl.get(ao.r[0], nan(""), "x0"); if(isnan(ao.r[0])) break;
			pl.get(ao.r[1], 0., "x1", true);
			pl.get(ao.r[2], 0., "x2", true);
			string aKey; pl.get(aKey, string(), "a");
			ao.numericalOrbIndex = -1;
			ao.sp = -1;
			if(aKey == "numerical")
			{	pl.get(ao.numericalOrbIndex, -1, "b", true);
				if(ao.numericalOrbIndex<0) throw(string("Numerical orbital index must be non-negative"));
			}
			else
			{	for(int sp=0; sp<int(e.iInfo.species.size()); sp++)
					if(e.iInfo.species[sp]->name == aKey)
					{	ao.sp = sp;
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
				}
			}
			pl.get(ao.coeff, 1., "coeff");
			//Transform coordinates if necessary
			if(e.iInfo.coordsType == CoordsCartesian)
				ao.r = inv(e.gInfo.R)*ao.r;
			t.push_back(ao);
		}
		if(!t.size()) throw(string("Trial orbital for each center must contain at least one atomic orbital"));
		wannier.trialOrbitals.push_back(t);
	}

	void printStatus(Everything& e, int iRep)
	{	const Wannier::TrialOrbital& t = wannier.trialOrbitals[iRep];
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


struct CommandWannierMinimize : public CommandMinimize
{	CommandWannierMinimize() : CommandMinimize("wannier") {}
    MinimizeParams& target(Everything& e) { return wannier.minParams; }
    void process(ParamList& pl, Everything& e)
	{	wannier.minParams.energyDiffThreshold = 1e-8; //override default value (0.) in MinimizeParams.h
		CommandMinimize::process(pl, e);
	}
}
commandWannierMinimize;


struct CommandWannierFilenames : public Command
{	string& target;
	CommandWannierFilenames(string cmdSuffix, string comment, string& target) : Command("wannier-"+cmdSuffix), target(target)
	{	format = "<format>";
		comments = 
			"Control the filename pattern for wannier " + comment + ", where <format> must contain\n"
			"a single instance of $VAR, which will be substituted by the name of each variable.";
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(target, string("$STAMP.$VAR"), "format");
		if(target.find("$VAR")==string::npos)
			throw "<format> = " + target + " doesn't contain the pattern $VAR";
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s", target.c_str());
	}
};
CommandWannierFilenames commandWannierInitialState("initial-state", "state initialization", wannier.initFilename);
CommandWannierFilenames commandWannierDumpName("dump-name", "output", wannier.dumpFilename);
