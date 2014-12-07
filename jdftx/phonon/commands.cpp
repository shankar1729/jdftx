/*-------------------------------------------------------------------
Copyright 2014 Ravishankar Sundararaman

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

#include <phonon/Phonon.h>
#include <commands/command.h>
#include <core/Units.h>

Phonon phonon;


enum PhononMember
{	PM_sup,
	PM_dr,
 	PM_T,
	PM_delim
};

EnumStringMap<PhononMember> phononMemberMap
(	PM_sup, "supercell",
	PM_dr, "dr",
	PM_T, "T"
);

struct CommandPhonon : public Command
{
	CommandPhonon() : Command("phonon", "Phonon")
	{
		format = "<key1> <args1...>  <key2> <args2...>  ...";
		comments =
			"Control phonon calculation and output. The possible <key>'s and their\n"
			"corresponding arguments are:\n"
			"  supercell <N0> <N1> <N2>\n"
			"    Supercell for frozen phonon perturbation. Each entry must divide\n"
			"    the corresponding value in kpoint-folding. The k-point mesh must\n"
			"    be uniform and centered on Gamma.\n"
			"  dr <dr>\n"
			"    Amplitude (in bohrs) of frozen phonon perturbation (default 0.01).\n"
			"  T <T>\n"
			"    Temperature (in Kelvins) used for vibrational free energy estimation (default 298).\n";
		
		forbid("fix-electron-density");
		forbid("fix-electron-potential");
	}

	void process(ParamList& pl, Everything& e)
	{	while(true)
		{	PhononMember key; pl.get(key, PM_delim, phononMemberMap, "key");
			if(key==PM_delim) break;
			switch(key)
			{	case PM_sup:
					for(int j=0; j<3; j++)
					{	char paramName[8]; sprintf(paramName, "N%d", j);
						pl.get(phonon.sup[j], 0, paramName, true);
						if(phonon.sup[j]<=0)
							throw string("supercell dimensions must be positive");
					}
					break;
				case PM_dr:
					pl.get(phonon.dr, 0., "dr", true);
					break;
				case PM_T:
					pl.get(phonon.T, 0., "T", true);
					phonon.T *= Kelvin;
					break;
				case PM_delim: //should never be encountered
					break;
			}
		}
	}

	void printStatus(Everything& e, int iRep)
	{	logPrintf(" \\\n\tsupercell %d %d %d", phonon.sup[0], phonon.sup[1], phonon.sup[2]);
		logPrintf(" \\\n\tdr %lg", phonon.dr);
		logPrintf(" \\\n\tT %lg", phonon.T/Kelvin);
	}
}
commandPhonon;
