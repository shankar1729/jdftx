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

/*
EnumStringMap<S2quadType> quadDescMap
(	QuadEuler,       "Euler angles (<nBeta> x <nAlpha> x ceil(<nGamma>/Zn) nodes on S1/Zn)",
	QuadTetrahedron, "Tetrahedron rotation group (4 S2 nodes x ceil(3/n) on S1/Zn)", 
	QuadOctahedron,  "Octahedron rotation group (6 S2 nodes x ceil(4/n) on S1/Zn)",
	QuadIcosahedron, "Icosahedron rotation group (12 S2 nodes x ceil(5/n) on S1/Zn)",
	Quad7design_24,   "Spherical 7-design with 24 S2 nodes (x ceil(8/n) on S1/Zn)",
	Quad8design_36,   "Spherical 8-design with 36 S2 nodes (x ceil(9/n) on S1/Zn)",
	Quad9design_48,   "Spherical 9-design with 48 S2 nodes (x ceil(10/n) on S1/Zn)",
	Quad10design_60,  "Spherical 10-design with 60 S2 nodes (x ceil(11/n) on S1/Zn)",
	Quad11design_70,  "Spherical 11-design with 70 S2 nodes (x ceil(12/n) on S1/Zn)",
	Quad12design_84,  "Spherical 12-design with 84 S2 nodes (x ceil(13/n) on S1/Zn)",
	Quad13design_94,  "Spherical 13-design with 94 S2 nodes (x ceil(14/n) on S1/Zn)",
	Quad14design_108, "Spherical 14-design with 108 S2 nodes (x ceil(15/n) on S1/Zn)",
	Quad15design_120, "Spherical 15-design with 120 S2 nodes (x ceil(16/n) on S1/Zn)",
	Quad16design_144, "Spherical 16-design with 144 S2 nodes (x ceil(17/n) on S1/Zn)",
	Quad17design_156, "Spherical 17-design with 156 S2 nodes (x ceil(18/n) on S1/Zn)",
	Quad18design_180, "Spherical 18-design with 180 S2 nodes (x ceil(19/n) on S1/Zn)",
	Quad19design_204, "Spherical 19-design with 204 S2 nodes (x ceil(20/n) on S1/Zn)",
	Quad20design_216, "Spherical 20-design with 216 S2 nodes (x ceil(21/n) on S1/Zn)",
	Quad21design_240, "Spherical 21-design with 240 S2 nodes (x ceil(22/n) on S1/Zn)"
);


struct CommandFluidSO3quad : public Command
{
	CommandFluidSO3quad() : Command("fluid-SO3-quad")
	{
		format = "<quad> [ <nBeta> [<nAlpha> <nGamma>] ]";
		comments = "Specify SO(3) quadrature <quad> which may be one of:"
			+ addDescriptions(S2quadTypeMap.optionList(), linkDescription(S2quadTypeMap, quadDescMap))
			+ ".\nSample-count <nBeta> must be specified for the Euler quadrature\n"
			"while <nAlpha> and <nGamma> are optional (= 2<nBeta> if 0 or unspecified).";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	FluidSolverParams& fsp = e.eVars.fluidParams;
		pl.get(fsp.s2quadType, QuadOctahedron, S2quadTypeMap, "quad");
		if(fsp.s2quadType == QuadEuler)
		{	pl.get(fsp.quad_nBeta, 0u, "nBeta", true); //required
			pl.get(fsp.quad_nAlpha, 0u, "nAlpha"); //optional
			pl.get(fsp.quad_nGamma, 0u, "nGamma"); //optional
			if(!fsp.quad_nBeta) throw string("<nBeta> must be non-zero.");
		}
	}

	void printStatus(Everything& e, int iRep)
	{	const FluidSolverParams& fsp = e.eVars.fluidParams;
		logPrintf("%s", S2quadTypeMap.getString(fsp.s2quadType));
		if(fsp.s2quadType == QuadEuler)
		{	logPrintf(" %d", fsp.quad_nBeta);
			if(fsp.quad_nAlpha || fsp.quad_nGamma)
				logPrintf(" %d %d", fsp.quad_nAlpha, fsp.quad_nGamma);
		}
	}
}
commandFluidSO3quad;
*/
