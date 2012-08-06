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

struct CommandWannier : public Command
{
	CommandWannier() : Command("wannier")
	{
		format = "<center1> [<center2> ...] with each <center> = <band> <x0> <x1> <x2> <width>";
		comments =
			"Specify a group of bands to be combined into Maximally-Localized\n"
			"Wannier Functions (MLWF). There are as many centers as bands in\n"
			"each group. The center is selected by a guess gaussian centered\n"
			"at <x0>,<x1>,<x2> (coordinate system set by coords-type) with\n"
			"width <w> bohrs. The command may be specified multiple times,\n"
			"one for each group. Wannier functions will be saved in files\n"
			"<band>.mlwf (or <band>.mlwfUp/.mlwfDn with spin) at the end;\n"
			"the dump frequency cany be controlled using the dump command.";
		allowMultiple = true;
		require("wannier-supercell");
		
		//Dependencies due to optional conversion from cartesian coords:
		require("lattice");
		require("latt-scale");
		require("coords-type");
	}

	void process(ParamList& pl, Everything& e)
	{	std::vector<Wannier::Center> group;
		while(true)
		{	Wannier::Center center;
			pl.get(center.band, -1, "band");
			if(center.band==-1) break; //end of input
			pl.get(center.r[0], 0., "x0", true);
			pl.get(center.r[1], 0., "x1", true);
			pl.get(center.r[2], 0., "x2", true);
			pl.get(center.width, 0., "width", true);
			//Transform coordinates if necessary
			if(e.iInfo.coordsType == CoordsCartesian)
				center.r = inv(e.gInfo.R)*center.r;
			group.push_back(center);
		}
		if(!group.size())
			throw(string("Each wannier group must contain at least one center"));
		e.dump.wannier.group.push_back(group);
		e.dump.insert(std::make_pair(DumpFreq_End, DumpWannier)); //dump at end
	}

	void printStatus(Everything& e, int iRep)
	{	std::vector<Wannier::Center>& group = e.dump.wannier.group[iRep];
		for(size_t i=0; i<group.size(); i++)
		{	if(group.size()>1) logPrintf(" \\\n\t");
			vector3<> r = group[i].r;
			if(e.iInfo.coordsType == CoordsCartesian)
				r = e.gInfo.R * r; //report cartesian positions
			logPrintf("%d    %lg %lg %lg  %lg", group[i].band, r[0], r[1], r[2], group[i].width);
		}
	}
}
commandWannier;

EnumStringMap<bool> realMap
(	true, "Real",
	false, "Complex"
);

struct CommandWannierSupercell : public Command
{
	CommandWannierSupercell() : Command("wannier-supercell")
	{
		format = "<n0> <n1> <n2> [<scalar>="+realMap.optionList()+"]";
		comments = "Number of unit cells in wannier function output.\n"
			"The supercell spans from ceil(-<ni>/2) to ceil(<ni>/2)\n"
			"in lattice coordinates along each dimension i=0,1,2.\n"
			"<scalar>=Real outputs real double precision functions\n"
			"on the supercell after phase elimination (default)\n"
			"while Complex outputs the complex functions directly.";
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.dump.wannier.supercell[0], 0, "n0", true);
		pl.get(e.dump.wannier.supercell[1], 0, "n1", true);
		pl.get(e.dump.wannier.supercell[2], 0, "n2", true);
		pl.get(e.dump.wannier.convertReal, true, realMap, "scalar");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%d %d %d %s", e.dump.wannier.supercell[0],
			e.dump.wannier.supercell[1], e.dump.wannier.supercell[2],
			realMap.getString(e.dump.wannier.convertReal));
	}
}
commandWannierSupercell;
