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


struct CommandJdft1shape : public Command
{
	CommandJdft1shape() : Command("jdft1-shape")
	{
		format = "[<nc>=7e-4] [<sigma>=0.6]";
		comments = "Control the critical density <nc> and smoothing parameter <sigma> for the JDFT1 cavity";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.eVars.fluidParams.nc, 7e-4, "nc");
		pl.get(e.eVars.fluidParams.sigma, 0.6, "sigma");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg %lg", e.eVars.fluidParams.nc, e.eVars.fluidParams.sigma);
	}
}
commandJdft1shape;
