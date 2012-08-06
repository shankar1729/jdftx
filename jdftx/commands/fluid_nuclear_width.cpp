/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman, Kendra Letchworth Weaver

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

struct CommandFluidNuclearWidth : public Command
{
	CommandFluidNuclearWidth() : Command("fluid-nuclear-width")
	{
		format = "<width>";
		comments = "Manually specify width of gaussian representations of fluid nuclear sites (in bohr)";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.eVars.fluidParams.nuclearWidth, 0., "width");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg", e.eVars.fluidParams.nuclearWidth);
	}
}
commandFluidNuclearWidth;
