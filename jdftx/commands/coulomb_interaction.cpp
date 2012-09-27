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

EnumStringMap<ExchangeRegularization::Method> exRegMethodMap
(	ExchangeRegularization::None,                 "None",
	ExchangeRegularization::AuxiliaryFunction,    "AuxiliaryFunction",
	ExchangeRegularization::SphericalTruncated,   "SphericalTruncated",
	ExchangeRegularization::WignerSeitzTruncated, "WignerSeitzTruncated"
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
			"   Slab <dir>=" + truncationDirMap.optionList() + " [<ionMargin>=3]\n"
			"      Truncate coulomb interaction along the specified lattice direction.\n"
			"      The other two lattice directions must be orthogonal to this one.\n"
			"      Useful for slab-like geometries.\n"
			"   Cylindrical <dir>=" + truncationDirMap.optionList() + " [<Rc>=0] [<ionMargin>=3]\n"
			"      Truncate coulomb interaction on a cylinder of radius <Rc> bohrs\n"
			"      with axis along specified lattice direction. The other two lattice\n"
			"      directions must be orthogonal to this one. Rc=0 is understood to be\n"
			"      the in-radius of the 2D Wigner-Seitz cell perpendicular to <dir>.\n"
			"   Wire <dir>=" + truncationDirMap.optionList() + " <borderWidth> [<ionMargin>=3] [<filename>]\n"
			"      Truncate coulomb interaction on the 2D Wigner-Seitz cell in the plane\n"
			"      perpendicular to <dir>. The other two lattice directions must be\n"
			"      orthogonal to this one. Useful for wire-like geometries.\n"
			"   Isolated <borderWidth> [<ionMargin>=3] [<filename>]\n"
			"      Truncate coulomb interaction on the 3D Wigner-Seitz cell.\n"
			"   Spherical [<Rc>=0] [<ionMargin>=3]\n"
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
			"of the charge density is not important and may cross unit cell boundaries.\n"
			"    The <ionMargin> option helps ensure that the charge density satisfies\n"
			"the constraints mentioned above. This option checks that the ions (nuclei)\n"
			"satisfy the truncation constraints with an additional margin of <ionMargin>\n"
			"bohrs; the latter should therefore be set to a typical distance from nuclei\n"
			"where the electron density becomes negligible.";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	CoulombParams& cp = e.coulombParams;
		pl.get(cp.geometry, CoulombParams::Periodic, truncationTypeMap, "truncationType");
		switch(cp.geometry)
		{	case CoulombParams::Periodic:
				break;
			case CoulombParams::Slab:
				pl.get(cp.iDir, 0, truncationDirMap, "dir", true);
				pl.get(cp.ionMargin, 3., "ionMargin");
				break;
			case CoulombParams::Wire:
			case CoulombParams::Isolated:
				if(cp.geometry == CoulombParams::Wire)
					pl.get(cp.iDir, 0, truncationDirMap, "dir", true);
				pl.get(cp.borderWidth, 0., "borderWidth", true);
				pl.get(cp.ionMargin, 3., "ionMargin");
				pl.get(cp.filename, string(), "filename");
				if(cp.borderWidth <= 0.)
					throw string("Border width must be positive, recommend at least 3 bohrs.\n");
				break;
			case CoulombParams::Cylindrical:
			case CoulombParams::Spherical:
				if(cp.geometry == CoulombParams::Cylindrical)
					pl.get(cp.iDir, 0, truncationDirMap, "dir", true);
				pl.get(cp.Rc, 0., "Rc");
				pl.get(cp.ionMargin, 3., "ionMargin");
		}
		if(cp.geometry!=CoulombParams::Periodic && cp.ionMargin<=0.)
			throw string("Ion margin must be positive");
	}

	void printStatus(Everything& e, int iRep)
	{	CoulombParams& cp = e.coulombParams;
		logPrintf("%s", truncationTypeMap.getString(cp.geometry));
		switch(cp.geometry)
		{	case CoulombParams::Periodic:
				break;
			case CoulombParams::Slab:
				logPrintf(" %s %lg", truncationDirMap.getString(cp.iDir), cp.ionMargin);
				break;
			case CoulombParams::Wire:
			case CoulombParams::Isolated:
				if(cp.geometry == CoulombParams::Wire)
					logPrintf(" %s", truncationDirMap.getString(cp.iDir));
				logPrintf(" %lg  %lg %s", cp.borderWidth, cp.ionMargin, cp.filename.c_str());
				break;
			case CoulombParams::Cylindrical:
			case CoulombParams::Spherical:
				if(cp.geometry == CoulombParams::Cylindrical)
					logPrintf(" %s", truncationDirMap.getString(cp.iDir));
				logPrintf(" %lg %lg", cp.Rc, cp.ionMargin);
		}
	}
}
commandCoulombInteraction;


struct CommandExchangeRegularization : public Command
{
	CommandExchangeRegularization() : Command("exchange-regularization")
	{
		format = "<method>=" + exRegMethodMap.optionList() + " [<sigma>=0.33 <filename>]";
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
			"      Default for any (partially) periodic <geometry>.\n"
			"<sigma> is the gaussian smoothing width (in bohrs) for the Wigner-Seitz\n"
			"cell boundary and <filename>, if specified, will be used to cache the\n"
			"computed kernel for the WignerSeitzTruncated method.";
		
		hasDefault = true;
		require("coulomb-interaction");
	};

	void process(ParamList& pl, Everything& e)
	{	const CoulombParams& cp = e.coulombParams;
		ExchangeRegularization& exReg = e.coulombParams.exchangeRegularization;
		//Select default method based on geometry:
		bool isIsolated = cp.geometry==CoulombParams::Isolated
				|| cp.geometry==CoulombParams::Spherical;
		exReg.method = isIsolated
			? ExchangeRegularization::None
			: ExchangeRegularization::WignerSeitzTruncated;
		pl.get(exReg.method, exReg.method, exRegMethodMap, "method");
		//Check compatibility of method and geometry:
		if(isIsolated && exReg.method!=ExchangeRegularization::None)
			throw string("exchange-regularization <method> must be None for non-periodic"
				" coulomb-interaction <geometry> = Spherical or Isolated");
		if(exReg.method==ExchangeRegularization::None
			&& !(isIsolated || cp.geometry==CoulombParams::Periodic))
			throw string("exchange-regularization <method> = None is supported only for"
				" non-periodic or 3D periodic values of coulomb-interaction <geometry>");
		if(exReg.method==ExchangeRegularization::AuxiliaryFunction
			&& !(cp.geometry==CoulombParams::Periodic))
			throw string("exchange-regularization <method> = AuxiliaryFunction is supported"
				" only for coulomb-interaction <geometry> = Periodic");
		//Read optional parameters:
		if(exReg.method==ExchangeRegularization::WignerSeitzTruncated)
		{	pl.get(exReg.sigma, 0.33, "sigma");
			pl.get(exReg.filename, string(), "filename");
			if(exReg.sigma <= 0.)
				throw string("Wigner-seitz boundary gaussian width <sigma> must be positive.\n");
		}
	}
	
	void printStatus(Everything& e, int iRep)
	{	const ExchangeRegularization& exReg = e.coulombParams.exchangeRegularization;
		logPrintf("%s", exRegMethodMap.getString(exReg.method));
		if(exReg.method==ExchangeRegularization::WignerSeitzTruncated)
			logPrintf(" %lg %s", exReg.sigma, exReg.filename.c_str());
	}
}
commandExchangeRegularization;
