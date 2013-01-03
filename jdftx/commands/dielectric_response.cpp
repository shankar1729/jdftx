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

#include <core/Units.h>
#include <commands/command.h>
#include <electronic/Everything.h>

struct CommandDielectricResponse : public Command
{
	CommandDielectricResponse() : Command("dielectric-response")
	{
		format = "[<epsilon-bulk>=80] [<Nbulk>=4.9383e-3] [<pMol>=0.944169] [<epsInf> = 1.77] [<linear>=no]";
		comments =
			"Options for the dielectric response of the fluid (defaults to water)\n"
			"\t<epsilon-bulk>: Bulk dielectric constant\n"
			"\t<Nbulk>: Bulk number density of molecules [in bohr^-3]\n"
			"\t<pMol>: Dipole moment of one molecule [in e/bohr^2]\n"
			"\t<epsInf>: High field (or frequency) dielectric constant of the fluid\n"
			"\t<linear>: Linearity of dielectric response = " + boolMap.optionList() + " (default: no)\n"
			"Note: <epsilon-bulk> only affects the Linear and Nonlinear 1.0 fluids\n"
			"      <linear> only affects the Nonlinear1.0 fluid (fluid 'Linear' is always linear!)\n"
			"      <Ndip> and <pMol> only affect the saturation behaviour of Nonlinear1.0.";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	FluidSolverParams& fsp = e.eVars.fluidParams;
		pl.get(fsp.epsilonBulk, 80.0, "epsilon-bulk");
		pl.get(fsp.Nbulk, 4.9383e-3, "Nbulk");
		pl.get(fsp.pMol, 0.92466, "pMol");
		pl.get(fsp.epsInf, 1.77, "epsInf");
		pl.get(fsp.linearDielectric, false, boolMap, "linear");
		
		if(fsp.epsInf < 1 or fsp.epsilonBulk < 1)
			throw string("Dielectric constant must be greater than 1.");
	}

	void printStatus(Everything& e, int iRep)
	{	const FluidSolverParams& fsp = e.eVars.fluidParams;
		logPrintf("%lg %lg %lg %lg %s", fsp.epsilonBulk, fsp.Nbulk, fsp.pMol, fsp.epsInf, boolMap.getString(fsp.linearDielectric));
	}
}
commandDielectricResponse;
