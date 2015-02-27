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

EnumStringMap<GridInfo::LatticeType> lattTypeMap
(	GridInfo::Triclinic, "Triclinic",
	GridInfo::Monoclinic, "Monoclinic",
	GridInfo::Orthorhombic, "Orthorhombic",
	GridInfo::Tetragonal, "Tetragonal",
	GridInfo::Rhombohedral, "Rhombohedral",
	GridInfo::Hexagonal, "Hexagonal",
	GridInfo::Cubic, "Cubic"
);

EnumStringMap<GridInfo::LatticeModification> lattModMap
(	GridInfo::BodyCentered, "Body-Centered",
	GridInfo::BaseCentered, "Base-Centered",
	GridInfo::FaceCentered, "Face-Centered"
);

struct CommandLattice : public Command
{
	CommandLattice() : Command("lattice")
	{
		format = " [<modification>] <lattice> <parameters...>\n"
			"\t| \\\n\t<R00> <R01> <R02> \\\n\t<R10> <R11> <R12> \\\n\t<R20> <R21> <R22>";
		comments = "Specify lattice by name and parameters, or manually lattice vectors (in columns).\n"
			"The options for the [<modification>] <lattice> <parameters...> scheme are:\n"
			"+ Triclinic <a> <b> <c> <alpha> <beta> <gamma>\n"
			"+ [Base-Centered] Monoclinic <a> <b> <c> <beta>\n"
			"+ [Base|Body|Face-Centered] Orthorhombic <a> <b> <c>\n"
			"+ [Body-Centered] Tetragonal <a> <c>\n"
			"+ Rhombohedral <a> <alpha>\n"
			"+ Hexagonal <a> <c>\n"
			"+ [Body|Face-Centered] Cubic <a>";
	}

	void process(ParamList& pl, Everything& e)
	{	GridInfo& g = e.gInfo;
		//Check for named lattice specification:
		try { pl.get(g.latticeModification, GridInfo::Simple, lattModMap, "modification", true); } catch(string) { pl.rewind(); } //no modification
		try { pl.get(g.latticeType, GridInfo::Manual, lattTypeMap, "type", true); } catch(string) { pl.rewind(); } //Not a named lattice
		//Read parameters:
		g.alpha = 90.;
		g.beta  = 90.;
		g.gamma = 90.;
		switch(g.latticeType)
		{
			case GridInfo::Manual:
			{	for(int j=0; j<3; j++) for(int k=0; k<3; k++)
				{	ostringstream oss; oss << "R" << j << k;
					try
					{	pl.get(g.R(j,k), 0.0, oss.str(), true);
					}
					catch(string)
					{	if(j==0 && k<2) throw string("First two parameters match neither <R00> <R01> ..., nor valid [<modification>] <lattice> ...");
						else throw; //propagate Rjk parse error
					}
				}
				break;
			}
			case GridInfo::Triclinic:
			{	pl.get(g.a, 0., "a", true);
				pl.get(g.b, 0., "b", true);
				pl.get(g.c, 0., "c", true);
				pl.get(g.alpha, 0., "alpha", true);
				pl.get(g.beta,  0., "beta", true);
				pl.get(g.gamma, 0., "gamma", true);
				break;
			}
			case GridInfo::Monoclinic:
			{	pl.get(g.a, 0., "a", true);
				pl.get(g.b, 0., "b", true);
				pl.get(g.c, 0., "c", true);
				pl.get(g.beta,  0., "beta", true);
				break;
			}
			case GridInfo::Orthorhombic:
			{	pl.get(g.a, 0., "a", true);
				pl.get(g.b, 0., "b", true);
				pl.get(g.c, 0., "c", true);
				break;
			}
			case GridInfo::Tetragonal:
			{	pl.get(g.a, 0., "a", true);
				g.b = g.a;
				pl.get(g.c, 0., "c", true);
				break;
			}
			case GridInfo::Rhombohedral:
			{	pl.get(g.a, 0., "a", true);
				g.b = g.a;
				g.c = g.a;
				pl.get(g.alpha, 0., "alpha", true);
				g.beta  = g.alpha;
				g.gamma = g.alpha;
				break;
			}
			case GridInfo::Hexagonal:
			{	pl.get(g.a, 0., "a", true);
				g.b = g.a;
				pl.get(g.c, 0., "c", true);
				g.gamma = 120.;
				break;
			}
			case GridInfo::Cubic:
			{	pl.get(g.a, 0., "a", true);
				g.b = g.a;
				g.c = g.a;
				break;
			}
		}
		//Check parameters and compute the lattice vectors:
		if(g.latticeType != GridInfo::Manual)
			g.setLatticeVectors();
	}

	void printStatus(Everything& e, int iRep)
	{	const GridInfo& g = e.gInfo;
		if(g.latticeType == GridInfo::Manual)
		{	matrix3<> Runscaled;
			for(int k=0; k<3; k++)
				Runscaled.set_col(k, (1./g.lattScale[k]) * g.R.column(k));
			for(int j=0; j<3; j++)
			{	logPrintf(" \\\n\t");
				for(int k=0; k<3; k++)
					logPrintf("%20.15lf ", Runscaled(j,k));
			}
		}
		else
		{	if(g.latticeModification != GridInfo::Simple)
				logPrintf("%s ", lattModMap.getString(g.latticeModification));
			logPrintf("%s ", lattTypeMap.getString(g.latticeType));
			switch(g.latticeType)
			{	case GridInfo::Manual: break; //never encountered
				case GridInfo::Triclinic: logPrintf("%lg %lg %lg %lg %lg %lg", g.a, g.b, g.c, g.alpha, g.beta, g.gamma); break;
				case GridInfo::Monoclinic: logPrintf("%lg %lg %lg %lg", g.a, g.b, g.c, g.beta); break;
				case GridInfo::Orthorhombic: logPrintf("%lg %lg %lg", g.a, g.b, g.c); break;
				case GridInfo::Tetragonal: logPrintf("%lg %lg", g.a, g.c); break;
				case GridInfo::Rhombohedral: logPrintf("%lg %lg", g.a, g.alpha); break;
				case GridInfo::Hexagonal: logPrintf("%lg %lg", g.a, g.c); break;
				case GridInfo::Cubic: logPrintf("%lg", g.a); break;
			}
		}
	}
}
commandLattice;


struct CommandLattScale : public Command
{
	CommandLattScale() : Command("latt-scale")
	{
		format = "<s0> <s1> <s2>";
		comments = "Scale lattice vector i by factor <si>";
		hasDefault = true;

		require("lattice");
	}

	void process(ParamList& pl, Everything& e)
	{	//check if any parameters have been specified:
		string key; pl.get(key, string(), ""); pl.rewind();
		if(!key.length()) { e.gInfo.lattScale = vector3<>(1.,1.,1.); return; }
		//Read parameters:
		for(int k=0; k<3; k++)
		{	ostringstream oss; oss << "s" << k;
			pl.get(e.gInfo.lattScale[k], 1., oss.str(), true);
			e.gInfo.R.set_col(k, e.gInfo.lattScale[k] * e.gInfo.R.column(k));
		}
	}

	void printStatus(Everything& e, int iRep)
	{	for(int k=0; k<3; k++) logPrintf("%lg ", e.gInfo.lattScale[k]);
	}
}
commandLattScale;


struct CommandLattMoveScale : public Command
{
	CommandLattMoveScale() : Command("latt-move-scale")
	{
		format = "<s0> <s1> <s2>";
		comments = "Preconditioning factor for each lattice vector (must be commensurate with symmetries)";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	vector3<>& s = e.cntrl.lattMoveScale;
		pl.get(s[0], 1., "s0");
		pl.get(s[1], 1., "s1");
		pl.get(s[2], 1., "s2");
	}

	void printStatus(Everything& e, int iRep)
	{	const vector3<>& s = e.cntrl.lattMoveScale;
		logPrintf("%lg %lg %lg", s[0], s[1], s[2]);
	}
}
commandLattMoveScale;

EnumStringMap<CoordsType> coordsMap(
	CoordsLattice, "lattice",
	CoordsCartesian, "cartesian" );

struct CommandCoordsType : public Command
{
	CommandCoordsType() : Command("coords-type")
	{
		format = "<coords>=" + coordsMap.optionList();
		comments = "Coordinate system used in specifying ion positions (default: lattice)";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.iInfo.coordsType, CoordsLattice, coordsMap, "coords");
	}

	void printStatus(Everything& e, int iRep)
	{	fputs(coordsMap.getString(e.iInfo.coordsType), globalLog);
	}
}
commandCoordsType;
