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
	CoulombParams::ProbeChargeEwald,     "ProbeChargeEwald",
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
	CommandCoulombInteraction() : Command("coulomb-interaction", "jdftx/Coulomb interactions")
	{
		format = "<truncationType> [<args> ...]";
		comments =
			"Optionally truncate the coulomb interaction. The available <truncationType>'s\n"
			"and the corresponding arguments are:\n"
			"\n+ Periodic\n\n"
			"    Standard periodic (untruncated) coulomb interaction (Default)\n"
			"\n+ Slab <dir>=" + truncationDirMap.optionList() + "\n\n"
			"    Truncate coulomb interaction along the specified lattice direction.\n"
			"    The other two lattice directions must be orthogonal to this one.\n"
			"    Useful for slab-like geometries.\n"
			"\n+ Cylindrical <dir>=" + truncationDirMap.optionList() + " [<Rc>=0]\n\n"
			"    Truncate coulomb interaction on a cylinder of radius <Rc> bohrs\n"
			"    with axis along specified lattice direction. The other two lattice\n"
			"    directions must be orthogonal to this one. Rc=0 is understood to be\n"
			"    the in-radius of the 2D Wigner-Seitz cell perpendicular to <dir>.\n"
			"\n+ Wire <dir>=" + truncationDirMap.optionList() + "\n\n"
			"    Truncate coulomb interaction on the 2D Wigner-Seitz cell in the plane\n"
			"    perpendicular to <dir>. The other two lattice directions must be\n"
			"    orthogonal to this one. Useful for wire-like geometries.\n"
			"\n+ Isolated\n\n"
			"    Truncate coulomb interaction on the 3D Wigner-Seitz cell.\n"
			"\n+ Spherical [<Rc>=0]\n\n"
			"    Truncate coulomb interaction on a sphere of radius <Rc> bohrs.\n"
			"    Rc=0 is understood to be the in-radius of the Wigner-Seitz cell.\n"
			"\n"
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


struct CommandCoulombTruncationEmbed : public Command
{
	CommandCoulombTruncationEmbed() : Command("coulomb-truncation-embed", "jdftx/Coulomb interactions")
	{
		format = "<c0> <c1> <c2>";
		comments =
			"Compute truncated Coulomb interaction in a double-sized box (doubled only\n"
			"along truncated directions). This relaxes the L/2 localization constraint\n"
			"otherwise required by truncated potentials (see command coulomb-interaction),\n"
			"but breaks translational invariance and requires the specification of a center.\n"
			"\n"
			"Coordinate system for center (<c0> <c1> <c2>) is as specified by coords-type.\n"
			"\n"
			"Default: not enabled; coulomb-interaction employs translationally invariant scheme";
		
		hasDefault = false;
		require("coulomb-interaction");
		//Dependencies due to coordinate system option:
		require("latt-scale");
		require("coords-type");
	}

	void process(ParamList& pl, Everything& e)
	{	e.coulombParams.embed = true;
		vector3<>& c = e.coulombParams.embedCenter;
		pl.get(c[0], 0., "c0", true);
		pl.get(c[1], 0., "c1", true);
		pl.get(c[2], 0., "c2", true);
		if(e.iInfo.coordsType==CoordsCartesian) c = inv(e.gInfo.R) * c; //Transform coordinates if necessary
		if(e.coulombParams.geometry==CoulombParams::Periodic)
			e.coulombParams.embed = false; //No embedding, just use to specify center
	}

	void printStatus(Everything& e, int iRep)
	{	vector3<> c = e.coulombParams.embedCenter;
		if(e.iInfo.coordsType==CoordsCartesian) c = e.gInfo.R * c; //Print in coordinate system chosen by user
		logPrintf("%lg %lg %lg", c[0], c[1], c[2]);
	}
}
commandCoulombTruncationEmbed;


struct CommandCoulombTruncationIonMargin : public Command
{
	CommandCoulombTruncationIonMargin() : Command("coulomb-truncation-ion-margin", "jdftx/Coulomb interactions")
	{
		format = "<margin>";
		comments =
			"Extra margin (in bohrs) around the ions, when checking localization constraints\n"
			"for truncated Coulomb potentials (see command coulomb-interaction). Set to a typical\n"
			"distance from nuclei where the electron density becomes negligible, so as to\n"
			"ensure the electron density satisfies those localization constraints.\n"
			"(Default: 5 bohrs, minimum allowed: 1 bohr)";
		hasDefault = false;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.coulombParams.ionMargin, 0., "margin", true);
		if(e.coulombParams.ionMargin < 1.) throw string("<margin> must be at least 1 bohr.");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg", e.coulombParams.ionMargin);
	}
}
commandCoulombTruncationIonMargin;


struct CommandExchangeRegularization : public Command
{
	CommandExchangeRegularization() : Command("exchange-regularization", "jdftx/Coulomb interactions")
	{
		format = "<method>=" + exRegMethodMap.optionList();
		comments =
			"Regularization / singularity correction method for exact exchange.\n"
			"The allowed methods and defaults depend on the setting of <geometry>\n"
			"in command coulomb-interaction\n"
			"\n+ None\n\n"
			"    No singularity correction: default and only option for non-periodic\n"
			"    systems with no G=0 singularity (<geometry> = Spherical / Isolated).\n"
			"    This is allowed for fully or partially periodic systems, but is not\n"
			"    recommended due to extremely poor convergence with number of k-points.\n"
			"\n+ AuxiliaryFunction\n\n"
			"    G=0 modification based on numerical integrals of an auxiliary\n"
			"    function, as described in Ref. \\cite AuxFunc-Carrier\n"
			"    Allowed for 3D/2D/1D periodic systems.\n"
			"\n+ ProbeChargeEwald\n\n"
			"    G=0 modification based on the Ewald sum of a single point charge\n"
			"    per k-point sampled supercell. Valid for 3D/2D/1D periodic systems.\n"
			"\n+ SphericalTruncated\n\n"
			"    Truncate exchange kernel on a sphere whose volume equals the k-point\n"
			"    sampled supercell, as in Ref. \\cite SphericalTruncation.\n"
			"    Allowed for any (partially) periodic <geometry>, but is recommended\n"
			"    only when the k-point sampled supercell is roughly isotropic.\n"
			"\n+ WignerSeitzTruncated\n\n"
			"    Truncate exchange kernel on the Wigner-Seitz cell of the k-point\n"
			"    sampled supercell, as in Ref. \\cite TruncatedEXX.\n"
			"    Default for any (partially) periodic <geometry>.";
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
	}
	
	void printStatus(Everything& e, int iRep)
	{	logPrintf("%s", exRegMethodMap.getString(e.coulombParams.exchangeRegularization));
	}
}
commandExchangeRegularization;


enum ExchangeParamsMember
{	EPM_blockSize,
	EPM_nOuterVxx,
	EPM_Delim
};
EnumStringMap<ExchangeParamsMember> epmMap
(	EPM_blockSize, "blockSize",
	EPM_nOuterVxx, "nOuterVxx"
);
EnumStringMap<ExchangeParamsMember> epmDescMap
(	EPM_blockSize, "Number of bands in blocks of FFTs used in exact-exchange calculation. Larger values are faster, but need more memory. (Default: 16)",
	EPM_nOuterVxx, "Maximum number of outer loop iterations to converge ACE exchange operator in SCF and band structure calculations. (Default: 20)"
);
struct CommandExchangeParams : public Command
{
	CommandExchangeParams() : Command("exchange-params", "jdftx/Coulomb interactions")
	{	
		format = "<key1> <value1> <key2> <value2> ...";
		comments = "Control exact exchange calculations. Possible keys and value types are:"
			+ addDescriptions(epmMap.optionList(), linkDescription(epmMap, epmDescMap))
			+ "\n\nAny number of these key-value pairs may be specified in any order.";
	}

	void process(ParamList& pl, Everything& e)
	{	while(true)
		{	ExchangeParamsMember key;
			pl.get(key, EPM_Delim, epmMap, "key");
			#define READ_AND_CHECK(param, target, op, val) \
				case EPM_##param: \
					pl.get(target, val, #param, true); \
					if(!(target op val)) throw string(#param " must be " #op " " #val); \
					break;
			switch(key)
			{	READ_AND_CHECK(blockSize, e.cntrl.exxBlockSize, >, 0)
				READ_AND_CHECK(nOuterVxx, e.cntrl.nOuterVxx, >, 0)
				case EPM_Delim: return; //end of input
			}
			#undef READ_AND_CHECK
		}
	}

	void printStatus(Everything& e, int iRep)
	{
		#define PRINT(param, target, format) logPrintf(" \\\n\t" #param " " format, target);
		PRINT(blockSize, e.cntrl.exxBlockSize, "%d")
		PRINT(nOuterVxx, e.cntrl.nOuterVxx, "%d")
		#undef PRINT
	}
}
commandExchangeParams;
