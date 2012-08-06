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
		format = "[<epsilon-bulk>=80] [<linear>=no]";
		comments =
			"Options for the dielectric response of the fluid\n"
			"\t<epsilon-bulk>: Bulk dielectric constant (default: 80)\n"
			"\t<linear>: Linearity of dielectric response = " + boolMap.optionList() + " (default: no)\n"
			"Note: <epsilon-bulk> only affects the Linear and Nonlinear 1.0 fluids\n"
			"      <linear> only affects the Nonlinear1.0 fluid (fluid 'Linear' is always linear!)";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	FluidSolverParams& fsp = e.eVars.fluidParams;
		pl.get(fsp.epsilonBulk, 80.0, "epsilon-bulk");
		pl.get(fsp.linearDielectric, false, boolMap, "linear");
	}

	void printStatus(Everything& e, int iRep)
	{	const FluidSolverParams& fsp = e.eVars.fluidParams;
		logPrintf("%lf %s", fsp.epsilonBulk, boolMap.getString(fsp.linearDielectric));
	}
}
commandDielectricResponse;
