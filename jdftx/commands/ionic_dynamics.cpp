/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Deniz Gunceler

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

#include <electronic/common.h>
#include <electronic/Everything.h>
#include <commands/command.h>
#include <commands/ParamList.h>
#include <core/VerletParams.h>
#include <core/Units.h>

struct CommandVerlet : public Command
{
	CommandVerlet() : Command("ionic-dynamics")
	{	format = "<time-step> <total-time> [<kT>=0.001] [<alpha>=0.0]";
		comments = "Applies Verlet molecular dynamics\n"
			   "Takes the time in femtoseconds and kT is in Hartree.\n"
			   "Give <alpha> if you want to the system to equilibrate with a heat bath at <kT>.\n"
			   "Also make sure to turn the symmetries off if the initial velocities don't satisfy the symmetry\n"
			   "conditions. Initial velocities will be assigned randomly if they are not given.";
		allowMultiple = false;
	}

	void process(ParamList& pl, Everything& e)
	{	VerletParams& vp = e.verletParams;
		double dt; pl.get(dt, 2.0, "time-step", true);
		double tMax; pl.get(tMax, 0.0, "max-time", true);
		vp.dt = dt*fs;
		vp.tMax = tMax*fs;
		pl.get(vp.kT, 0.001, "kT", false);
		pl.get(vp.alpha, 0.0, "alpha", false);
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg %lg %lg %lg",e.verletParams.dt/fs,e.verletParams.tMax/fs, e.verletParams.kT, e.verletParams.alpha);
	}
}
commandVerlet;
