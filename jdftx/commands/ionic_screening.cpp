/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver

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
#include <core/Units.h>


struct CommandIonicScreening : public Command
{
	CommandIonicScreening() : Command("ionic-screening")
	{
		format = "<concentration> <Zelectrolyte> [<linear>]";
		comments =
			"\t<concentration>: molar concentration of ions (default: 0.0, which turns off ionic screening)\n"
			"\t<Zelectrolyte>: magnitude of charge of the cations and anions (assumed equal)\n"
			"\t<linear>: linearity of screening = " + boolMap.optionList() + " (default: no)\n"
			" Note: <linear> only affects the Nonlinear1.0 fluid (fluid 'Linear' is always linear!)";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	FluidSolverParams& fsp = e.eVars.fluidParams;
		pl.get(fsp.ionicConcentration, 0.0, "concentration");
		pl.get(fsp.ionicZelectrolyte, 1, "Zelectrolyte");
		pl.get(fsp.linearScreening, false, boolMap, "linear");
		//convert to atomic units
		fsp.ionicConcentration *= mol/liter;
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg %d %s",
			e.eVars.fluidParams.ionicConcentration/(mol/liter), //report back in mol/liter
			e.eVars.fluidParams.ionicZelectrolyte,
			boolMap.getString(e.eVars.fluidParams.linearScreening));
	}
}
commandIonicScreening;
