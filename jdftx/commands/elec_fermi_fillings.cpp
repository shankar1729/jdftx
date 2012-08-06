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

struct CommandElecFermiFillings : public Command
{
	CommandElecFermiFillings() : Command("elec-fermi-fillings")
	{
		format = "<mixInterval> <kT> [<alpha>=0.5]";
		comments =
			"Fermi-Dirac fillings at a temperature <kT> (in Hartrees).\n"
			"  If <mixInterval> is zero, use the auxilliary hamiltonian method (recommended).\n"
			"  Else, mix fermi functions every <mixInterval> iterations with mixing fraction <alpha>.";
		
		require("elec-n-bands");
		forbid("fix-electron-density");
	}

	void process(ParamList& pl, Everything& e)
	{	ElecInfo& eInfo = e.eInfo;
		pl.get(eInfo.mixInterval, 0, "mixInterval", true);
		//Determine algorithm based on mixInterval:
		if(eInfo.mixInterval<0) throw string("<mixInterval> must be positive");
		eInfo.fillingsUpdate = eInfo.mixInterval ? ElecInfo::FermiFillingsMix : ElecInfo::FermiFillingsAux;
		pl.get(eInfo.kT, 0.0, "kT", true);
		if(eInfo.mixInterval) pl.get(eInfo.fillingMixFraction, 0.5, "alpha");
	}

	void printStatus(Everything& e, int iRep)
	{	const ElecInfo& eInfo = e.eInfo;
		logPrintf("%d %lg", eInfo.mixInterval, eInfo.kT);
		if(eInfo.mixInterval) logPrintf(" %lg", eInfo.fillingMixFraction);
	}
}
commandElecFermiFillings;
