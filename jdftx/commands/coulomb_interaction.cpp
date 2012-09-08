/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

EnumStringMap<CoulombTruncationParams::Type> truncationTypeMap
(	CoulombTruncationParams::Periodic,  "Periodic",
	CoulombTruncationParams::Slab,      "Slab",
	CoulombTruncationParams::Wire,      "Wire",
	CoulombTruncationParams::Isolated,  "Isolated",
	CoulombTruncationParams::Spherical, "Spherical"
);

EnumStringMap<int> truncationDirMap
(	0, "100",
	1, "010",
	2, "001"
);

struct CommandCoulombInteraction : public Command
{
	CommandCoulombInteraction() : Command("coulomb-interaction")
	{
		format = "<truncationType> [<args> ...]";
		comments =
			"Optionally truncate the coulomb interaction. The available truncation modes\n"
			"and the corresponding arguments are:\n"
			"   Periodic\n"
			"      Standard periodic (untruncated) coulomb interaction (Default)\n"
			"   Slab <dir>=" + truncationDirMap.optionList() + "\n"
			"      Truncate coulomb interaction along the specified lattice direction.\n"
			"      The other two lattice directions must be orthogonal to this one.\n"
			"      Useful for slab-like geometries.\n"
			"   Wire <dir>=" + truncationDirMap.optionList() + " <borderWidth> [<filename>]\n"
			"      Truncate coulomb interaction on the 2D Wigner-Seitz cell in the plane\n"
			"      perpendicular to <dir>. The other two lattice directions must be\n"
			"      orthogonal to this one. Useful for wire-like geometries.\n"
			"   Isolated <borderWidth> [<filename>]\n"
			"      Truncate coulomb interaction on the 3D Wigner-Seitz cell.\n"
			"   Spherical [<Rc>=0]\n"
			"      Truncate coulomb interaction on a sphere of radius <Rc> bohrs.\n"
			"      Rc=0 is understood to be the in-radius of the Wigner-Seitz cell.\n"
			"The Periodic, Slab and Spherical coulomb interaction kernels are computed\n"
			"analytically, whereas the Wire and Isolated use a numerically precomputed\n"
			"kernel using the O(NlogN) convolved Wigner-Seitz truncation scheme.\n"
			"This scheme smoothens the truncation boundary within a region of width\n"
			"<borderWidth>, which must be free of any charge density. Reducing this\n"
			"width increases the usable fraction of the unit cell, but increase the\n"
			"memory and time required to precompute the kernel, both of which scale as\n"
			"1/<borderWidth>^d with d=3 for Isolated and d=2 for Slab mode. These kernels\n"
			"can optionally be saved to <filename>. If such a file exists and matches\n"
			"the specified parameters, the kernel will be loaded from this file.\n"
			"    Note that for all the truncated modes, the charge density must be\n"
			"confined to a maximum separation of L/2-<borderWidth> in each truncated\n"
			"direction, where L is the length of the unit cell in that direction\n"
			"or 2 Rc for Spherical mode, and <borderWidth> is the truncation boundary\n"
			"width for the Slab and Isolated modes and zero for the rest. The center\n"
			"of the charge density is not important and may cross unit cell boundaries.";
		
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	CoulombTruncationParams& ctp = e.coulombTrunctaionParams;
		pl.get(ctp.type, CoulombTruncationParams::Periodic, truncationTypeMap, "truncationType");
		switch(ctp.type)
		{	case CoulombTruncationParams::Periodic:
				break;
			case CoulombTruncationParams::Slab:
				pl.get(ctp.iDir, 0, truncationDirMap, "dir", true);
				break;
			case CoulombTruncationParams::Wire:
			case CoulombTruncationParams::Isolated:
				if(ctp.type == CoulombTruncationParams::Wire)
					pl.get(ctp.iDir, 0, truncationDirMap, "dir", true);
				pl.get(ctp.borderWidth, 0., "borderWidth", true);
				pl.get(ctp.filename, string(), "filename");
				if(ctp.borderWidth <= 0.)
					throw string("Border width must be positive, recommend at least 3 bohrs.\n");
				break;
			case CoulombTruncationParams::Spherical:
				pl.get(ctp.Rc, 0., "Rc");
		}
	}

	void printStatus(Everything& e, int iRep)
	{	CoulombTruncationParams& ctp = e.coulombTrunctaionParams;
		logPrintf("%s", truncationTypeMap.getString(ctp.type));
		switch(ctp.type)
		{	case CoulombTruncationParams::Periodic:
				break;
			case CoulombTruncationParams::Slab:
				logPrintf(" %s", truncationDirMap.getString(ctp.iDir));
				break;
			case CoulombTruncationParams::Wire:
			case CoulombTruncationParams::Isolated:
				if(ctp.type == CoulombTruncationParams::Wire)
					logPrintf(" %s", truncationDirMap.getString(ctp.iDir));
				logPrintf(" %lg %s", ctp.borderWidth, ctp.filename.c_str());
				break;
			case CoulombTruncationParams::Spherical:
				if(ctp.Rc) logPrintf(" %lg", ctp.Rc);
		}
	}
}
commandCoulombInteraction;
