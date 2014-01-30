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
{	WM_bStart,
	WM_outerWindow,
	WM_innerWindow,
	WM_saveWfns,
	WM_loadRotations,
	WM_delim
};

EnumStringMap<WannierMember> wannierMemberMap
(	WM_bStart, "bStart",
	WM_outerWindow, "outerWindow",
	WM_innerWindow, "innerWindow",
	WM_saveWfns, "saveWfns",
	WM_loadRotations, "loadRotations"
);

struct CommandWannier : public Command
{
	CommandWannier() : Command("wannier")
	{
		format = "<key1> <args1...>  <key2> <args2...>  ...";
		comments =
			"Control wannier calculation and output. The possible <key>'s and their\n"
			"corresponding arguments are:\n"
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
			"    Whether to write supercell wavefunctions (can be enormous).\n"
			"    Default: no.\n"
			"  loadRotations yes|no\n"
			"    Whether to load rotations (.mlwU and .mlwfU2) from a previous Wannier run.\n"
			"    Default: no.";
	}

	void process(ParamList& pl, Everything& e)
	{	while(true)
		{	WannierMember key; pl.get(key, WM_delim, wannierMemberMap, "key");
			if(key==WM_delim) break;
			switch(key)
			{	case WM_bStart:
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
				case WM_loadRotations:
					pl.get(wannier.loadRotations, false, boolMap, "loadRotations", true);
					break;
				case WM_delim: //should never be encountered
					break;
			}
		}
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("saveWfns %s loadRotations %s", boolMap.getString(wannier.saveWfns), boolMap.getString(wannier.loadRotations));
		if(wannier.outerWindow) logPrintf(" \\\n\touterWindow %lg %lg", wannier.eOuterMin, wannier.eOuterMax);
		if(wannier.innerWindow) logPrintf(" \\\n\tinnerWindow %lg %lg", wannier.eInnerMin, wannier.eInnerMax);
		if(!(wannier.innerWindow || wannier.outerWindow))
			logPrintf(" \\\n\tbStart %d", wannier.bStart);
	}
}
commandWannier;


struct CommandWannierCenter : public Command
{
	CommandWannierCenter() : Command("wannier-center")
	{
		format = "<horb1> [<horb2> ...]";
		comments =
			"Specify trial orbital for a wannier function as a linear combination of\n"
			"hydrogenic orbitals <horb*>. The syntax for each <horb> is:\n"
			"   <x0> <x1> <x2> [<a>=1] [<orbDesc>=s] [<coeff>=1.0]\n"
			"which represents a hydrogenic orbital centered at <x0>,<x1>,<x2>"
			"(coordinate system set by coords-type) with decay length <a> bohrs,\n"
			"and with orbital code <orbDesc> as in command density-of-states. Note that\n"
			"<a> sets the decay length of the nodeless orbital of specified angular\n"
			"momentum in <orbDesc>. Specify a, orbDesc and coeff explicitly when using\n"
			"multiple orbitals; the defaults only apply to the single orbital case.\n"
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
		{	Wannier::Hydrogenic h; string orbDesc;
			pl.get(h.r[0], nan(""), "x0"); if(isnan(h.r[0])) break;
			pl.get(h.r[1], 0., "x1", true);
			pl.get(h.r[2], 0., "x2", true);
			pl.get(h.a, 1., "a");
			pl.get(orbDesc, string(), "orbDesc");
			if(orbDesc.length())
			{	h.orbitalDesc.parse(orbDesc);
				if(h.orbitalDesc.m > h.orbitalDesc.l)
					throw(string("Must specify a specific projection eg. px,py (not just p)"));
				if(h.orbitalDesc.n + h.orbitalDesc.l > 3)
					throw(string("Hydrogenic orbitals with n+l>4 not supported"));
			}
			else //default is nodeless s orbital
			{	h.orbitalDesc.l = 0;
				h.orbitalDesc.m = 0;
				h.orbitalDesc.n = 0;
			}
			pl.get(h.coeff, 1., "coeff");
			//Transform coordinates if necessary
			if(e.iInfo.coordsType == CoordsCartesian)
				h.r = inv(e.gInfo.R)*h.r;
			t.push_back(h);
		}
		if(!t.size()) throw(string("Trial orbital for each center must contain at least one hydrogenic orbital"));
		wannier.trialOrbitals.push_back(t);
	}

	void printStatus(Everything& e, int iRep)
	{	const Wannier::TrialOrbital& t = wannier.trialOrbitals[iRep];
		for(const Wannier::Hydrogenic& h: t)
		{	if(t.size()>1) logPrintf(" \\\n\t");
			vector3<> r = h.r;
			if(e.iInfo.coordsType == CoordsCartesian)
				r = e.gInfo.R * r; //report cartesian positions
			logPrintf("%lg %lg %lg  %lg  %s  %lg", r[0], r[1], r[2], h.a, string(h.orbitalDesc).c_str(), h.coeff);
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
