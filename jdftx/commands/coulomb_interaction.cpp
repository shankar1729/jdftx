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


/** \page CommandCoulombInteraction coulomb-interaction
Syntax
------
    coulomb-interaction <truncationType> [<args> ...]

Description
-----------
Optionally truncate the coulomb interaction. The available truncation modes
and the corresponding arguments are:

   + Periodic<br>
      Standard periodic (untruncated) coulomb interaction (Default)
   + Slab @<dir@>=100|010|001<br>
      Truncate coulomb interaction along the specified lattice direction.
      The other two lattice directions must be orthogonal to this one.
      Useful for slab-like geometries.
   + Cylindrical @<dir@>=100|010|001 [@<Rc@>=0]<br>
      Truncate coulomb interaction on a cylinder of radius @<Rc@> bohrs
      with axis along specified lattice direction. The other two lattice
      directions must be orthogonal to this one. Rc=0 is understood to be
      the in-radius of the 2D Wigner-Seitz cell perpendicular to @<dir@>.
   + Wire @<dir@>=100|010|001<br>
      Truncate coulomb interaction on the 2D Wigner-Seitz cell in the plane
      perpendicular to @<dir@>. The other two lattice directions must be
      orthogonal to this one. Useful for wire-like geometries.
   + Isolated<br>
      Truncate coulomb interaction on the 3D Wigner-Seitz cell.
   + Spherical [@<Rc@>=0]<br>
      Truncate coulomb interaction on a sphere of radius @<Rc@> bohrs.
      Rc=0 is understood to be the in-radius of the Wigner-Seitz cell.

For all the truncated modes, the charge density must be confined to a
maximum separation of L/2 in each truncated direction, where L is the
length of the unit cell in that direction or 2 Rc for Spherical and
Cylindrical modes. The center of the charge density is not important
and may cross unit cell boundaries.
*/
struct CommandCoulombInteraction : public Command
{
	CommandCoulombInteraction() : Command("coulomb-interaction", "Coulomb")
	{
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

/** \page CommandCoulombTruncationEmbed coulomb-truncation-embed
Syntax
------
    coulomb-truncation-embed <c0> <c1> <c2>

Description
-----------
Compute truncated Coulomb interaction in a double-sized box (doubled only
along truncated directions). This relaxes the L/2 localization constraint
otherwise required by truncated potentials (see command coulomb-interaction),
but breaks translational invariance and requires the specification of a center.
Coordinate system for center (@<c0@> @<c1@> @<c2@>) is as specified by coords-type.
(Default: not enabled i.e. employs translationally invariant scheme)

Requires
--------
\ref CommandCoulombInteraction
\ref CommandLattScale
\ref CommandCoordsType
*/
struct CommandCoulombTruncationEmbed : public Command
{
	CommandCoulombTruncationEmbed() : Command("coulomb-truncation-embed", "Coulomb")
	{
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
			throw string("coulomb-truncation-embed should only be specified for truncated geometries");
	}

	void printStatus(Everything& e, int iRep)
	{	vector3<> c = e.coulombParams.embedCenter;
		if(e.iInfo.coordsType==CoordsCartesian) c = e.gInfo.R * c; //Print in coordinate system chosen by user
		logPrintf("%lg %lg %lg", c[0], c[1], c[2]);
	}
}
commandCoulombTruncationEmbed;

/** \page CommandCoulombTruncationIonMargin coulomb-truncation-ion-margin
Syntax
------
    coulomb-truncation-ion-margin <margin>

Description
-----------
Extra margin (in bohrs) around the ions, when checking localization constraints
for truncated Coulomb potentials (see \ref CommandCoulombInteraction). Set to a typical
distance from nuclei where the electron density becomes negligible, so as to
ensure the electron density satisfies those localization constraints.
(Default: 5 bohrs, minimum allowed: 1 bohr)";
*/
struct CommandCoulombTruncationIonMargin : public Command
{
	CommandCoulombTruncationIonMargin() : Command("coulomb-truncation-ion-margin", "Coulomb")
	{
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


/** \page CommandExchangeRegularization exchange-regularization
Syntax
------
    exchange-regularization <method>

Description
-----------
Regularization / singularity correction method for exact exchange.
The allowed methods and defaults depend on the setting of @<geometry@>
in \ref CommandCoulombInteraction
 + None<br>
      No singularity correction: default and only option for non-periodic
      systems with no G=0 singularity (@<geometry@> = Spherical / Isolated).
      This is allowed for fully or partially periodic systems, but is not
      recommended due to extremely poor convergence with number of k-points.
 + AuxiliaryFunction<br>
      G=0 modification based on numerical integrals of an auxiliary
      function, as described in P. Carrier et al, PRB 75, 205126 (2007).
      Allowed for 3D/2D/1D periodic systems.
 + ProbeChargeEwald<br>
      G=0 modification based on the Ewald sum of a single point charge
      per k-point sampled supercell. Valid for 3D/2D/1D periodic systems.
 + SphericalTruncated<br>
      Truncate exchange kernel on a sphere whose volume equals the k-point
      sampled supercell, as in J. Spencer et al, PRB 77, 193110 (2008).
      Allowed for any (partially) periodic @<geometry@>, but is recommended
      only when the k-point sampled supercell is roughly isotropic.
 + WignerSeitzTruncated<br>
      Truncate exchange kernel on the Wigner-Seitz cell of the k-point
      sampled supercell, as in R. Sundararaman et al (under preparation).
      Default for any (partially) periodic @<geometry@>.
*/
struct CommandExchangeRegularization : public Command
{
	CommandExchangeRegularization() : Command("exchange-regularization", "Coulomb")
	{
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
commandCoulombParams;
