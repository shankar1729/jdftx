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

EnumStringMap<CoulombParams::ExchangeRegularization> exRegMethodMap
(	CoulombParams::None,                 "None",
	CoulombParams::AuxiliaryFunction,    "AuxiliaryFunction",
	CoulombParams::SphericalTruncated,   "SphericalTruncated",
	CoulombParams::WignerSeitzTruncated, "WignerSeitzTruncated"
);

EnumStringMap<CoulombParams::Geometry> truncationTypeMap
(	CoulombParams::Periodic,    "Periodic",
	CoulombParams::Slab,        "Slab",
	CoulombParams::Cylindrical, "Cylindrical",
	CoulombParams::Wire,        "Wire",
	CoulombParams::Isolated,    "Isolated",
	CoulombParams::Spherical,   "Spherical"
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
			"   Cylindrical <dir>=" + truncationDirMap.optionList() + " [<Rc>=0]\n"
			"      Truncate coulomb interaction on a cylinder of radius <Rc> bohrs\n"
			"      with axis along specified lattice direction. The other two lattice\n"
			"      directions must be orthogonal to this one. Rc=0 is understood to be\n"
			"      the in-radius of the 2D Wigner-Seitz cell perpendicular to <dir>.\n"
			"   Wire <dir>=" + truncationDirMap.optionList() + "\n"
			"      Truncate coulomb interaction on the 2D Wigner-Seitz cell in the plane\n"
			"      perpendicular to <dir>. The other two lattice directions must be\n"
			"      orthogonal to this one. Useful for wire-like geometries.\n"
			"   Isolated\n"
			"      Truncate coulomb interaction on the 3D Wigner-Seitz cell.\n"
			"   Spherical [<Rc>=0]\n"
			"      Truncate coulomb interaction on a sphere of radius <Rc> bohrs.\n"
			"      Rc=0 is understood to be the in-radius of the Wigner-Seitz cell.\n"
			"For all the truncated modes, the charge density must be confined to a\n"
			"maximum separation of L/2 in each truncated direction, where L is the\n"
			"length of the unit cell in that direction or 2 Rc for Spherical and\n"
			"Cylindrical modes. The center of the charge density is not important\n"
			"and may cross unit cell boundaries.";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	CoulombParams& cp = e.coulombParams;
		pl.get(cp.geometry, CoulombParams::Periodic, truncationTypeMap, "truncationType");
		if(cp.geometry==CoulombParams::Periodic) return; //no parameters
		//Get direction for the partially periodic modes:
		if(cp.geometry==CoulombParams::Slab
		|| cp.geometry==CoulombParams::Wire
		|| cp.geometry==CoulombParams::Cylindrical)
			pl.get(cp.iDir, 0, truncationDirMap, "dir", true);
		//Get optional radius for the cylindrical/spherical modes:
		if(cp.geometry==CoulombParams::Cylindrical
		|| cp.geometry==CoulombParams::Spherical)
			pl.get(cp.Rc, 0., "Rc");
	}

	void printStatus(Everything& e, int iRep)
	{	CoulombParams& cp = e.coulombParams;
		logPrintf("%s", truncationTypeMap.getString(cp.geometry));
		if(cp.geometry==CoulombParams::Periodic) return; //no parameters
		//Print direction for the partially periodic modes:
		if(cp.geometry==CoulombParams::Slab
		|| cp.geometry==CoulombParams::Wire
		|| cp.geometry==CoulombParams::Cylindrical)
			logPrintf(" %s", truncationDirMap.getString(cp.iDir));
		//Print optional radius for the cylindrical/spherical modes:
		if(cp.geometry==CoulombParams::Cylindrical
		|| cp.geometry==CoulombParams::Spherical)
			logPrintf(" %lg", cp.Rc);
	}
}
commandCoulombInteraction;


struct CommandExchangeRegularization : public Command
{
	CommandExchangeRegularization() : Command("exchange-regularization")
	{
		format = "<method>=" + exRegMethodMap.optionList();
		comments =
			"Regularization / singularity correction method for exact exchange.\n"
			"The allowed methods and defaults depend on the setting of <geometry>\n"
			"in command coulomb-interaction\n"
			"   None\n"
			"      No singularity correction. This is the default and only option\n"
			"      for non-periodic systems (<geometry> = Spherical / Isolated),\n"
			"      which have no G=0 singularity. This is allowed for 3D periodic\n"
			"      systems (<geometry> = Periodic), but is not recommended due to\n"
			"      extremely poor convergence with number of k-points.\n"
			"   AuxiliaryFunction\n"
			"      G=0 modification based on numerical integrals of an auxiliary\n"
			"      function, as described in P. Carrier et al, PRB 75, 205126 (2007).\n"
			"      This is allowed only for 3D periodic systems.\n"
			"   SphericalTruncated\n"
			"      Truncate exchange kernel on a sphere whose volume equals the k-point\n"
			"      sampled supercell, as in J. Spencer et al, PRB 77, 193110 (2008).\n"
			"      Allowed for any (partially) periodic <geometry>, but is recommended\n"
			"      only when the k-point sampled supercell is roughly isotropic.\n"
			"   WignerSeitzTruncated\n"
			"      Truncate exchange kernel on the Wigner-Seitz cell of the k-point\n"
			"      sampled supercell, as in R. Sundararaman et al (under preparation).\n"
			"      Default for any (partially) periodic <geometry>.";
		hasDefault = true;
		require("coulomb-interaction");
	};

	void process(ParamList& pl, Everything& e)
	{	CoulombParams& cp = e.coulombParams;
		//Select default regularization based on geometry:
		bool isIsolated = cp.geometry==CoulombParams::Isolated
				|| cp.geometry==CoulombParams::Spherical;
		cp.exchangeRegularization = isIsolated
			? CoulombParams::None
			: CoulombParams::WignerSeitzTruncated;
		pl.get(cp.exchangeRegularization, cp.exchangeRegularization, exRegMethodMap, "method");
		//Check compatibility of regularization with geometry:
		if(isIsolated && cp.exchangeRegularization!=CoulombParams::None)
			throw string("exchange-regularization <method> must be None for non-periodic"
				" coulomb-interaction <geometry> = Spherical or Isolated");
		if(cp.exchangeRegularization==CoulombParams::None
			&& !(isIsolated || cp.geometry==CoulombParams::Periodic))
			throw string("exchange-regularization <method> = None is supported only for"
				" non-periodic or 3D periodic values of coulomb-interaction <geometry>");
		if(cp.exchangeRegularization==CoulombParams::AuxiliaryFunction
			&& !(cp.geometry==CoulombParams::Periodic))
			throw string("exchange-regularization <method> = AuxiliaryFunction is supported"
				" only for coulomb-interaction <geometry> = Periodic");
	}
	
	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s", exRegMethodMap.getString(e.coulombParams.exchangeRegularization));
	}
}
commandCoulombParams;
