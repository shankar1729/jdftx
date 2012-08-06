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

struct CommandWavefunctionDrag : public Command
{
	CommandWavefunctionDrag() : Command("wavefunction-drag")
	{
		format = "[<radius>=0]";
		comments =
			"Typical extent of space (in bohrs) dragged by each atom as it moves.\n"
			"(used to guess wavefunctions after ionic steps). Set to 0 to disable.";
		hasDefault = true;
	}

	void process(ParamList& pl, Everything& e)
	{	pl.get(e.cntrl.dragRadius, 0.0, "radius");
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg", e.cntrl.dragRadius);
	}
}
commandWavefunctionDrag;
