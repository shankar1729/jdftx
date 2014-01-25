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
	WM_delim
};

EnumStringMap<WannierMember> wannierMemberMap
(	WM_bStart, "bStart",
	WM_outerWindow, "outerWindow",
	WM_innerWindow, "innerWindow",
	WM_saveWfns, "saveWfns"
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
				case WM_delim: //should never be encountered
					break;
			}
		}
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("saveWfns %s", boolMap.getString(wannier.saveWfns));
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
		format = "<x0> <x1> <x2> [<a>=1] [<orbDesc>=s]";
		comments =
			"Specify trial orbital for a wannier function as a hydrogenic orbital centered\n"
			"at <x0>,<x1>,<x2> (coordinate system set by coords-type) with decay length\n"
			"<a> bohrs, and with <orbDesc> as in command density-of-states. Note that\n"
			"<a> sets the decay length of the nodeless orbital of specified angular\n"
			"momentum in <orbDesc>. Specify the command once for each Wannier function.";
		allowMultiple = true;
		require("wannier-initial-state");
		require("wannier-dump-name");
		
		//Dependencies due to optional conversion from cartesian coords:
		require("latt-scale");
		require("coords-type");
	}

	void process(ParamList& pl, Everything& e)
	{	Wannier::Center center; string orbDesc;
		pl.get(center.r[0], 0., "x0", true);
		pl.get(center.r[1], 0., "x1", true);
		pl.get(center.r[2], 0., "x2", true);
		pl.get(center.a, 1., "a");
		pl.get(orbDesc, string(), "orbDesc");
		if(orbDesc.length())
		{	center.orbitalDesc.parse(orbDesc);
			if(center.orbitalDesc.m > center.orbitalDesc.l)
				throw(string("Must specify a specific projection eg. px,py (not just p)"));
			if(center.orbitalDesc.n + center.orbitalDesc.l > 3)
				throw(string("Hydrogenic orbitals with n+l>4 not supported"));
		}
		else //default is nodeless s orbital
		{	center.orbitalDesc.l = 0;
			center.orbitalDesc.m = 0;
			center.orbitalDesc.n = 0;
		}
		//Transform coordinates if necessary
		if(e.iInfo.coordsType == CoordsCartesian)
			center.r = inv(e.gInfo.R)*center.r;
		wannier.centers.push_back(center);
	}

	void printStatus(Everything& e, int iRep)
	{	const Wannier::Center& center = wannier.centers[iRep];
		vector3<> r = center.r;
		if(e.iInfo.coordsType == CoordsCartesian)
			r = e.gInfo.R * r; //report cartesian positions
		logPrintf("%lg %lg %lg  %lg  %s", r[0], r[1], r[2], center.a, string(center.orbitalDesc).c_str());
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
