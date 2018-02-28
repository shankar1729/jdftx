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

enum PhononMember
{	PM_sup,
	PM_dr,
	PM_iPerturbation,
	PM_collectPerturbations,
	PM_saveHsub,
 	PM_T,
	PM_Fcut,
	PM_rSmooth,
	PM_delim
};

EnumStringMap<PhononMember> phononMemberMap
(	PM_sup, "supercell",
	PM_dr, "dr",
	PM_iPerturbation,"iPerturbation",
	PM_collectPerturbations, "collectPerturbations",
	PM_saveHsub, "saveHsub",
	PM_T, "T",
	PM_Fcut, "Fcut",
	PM_rSmooth, "rSmooth"
);

struct CommandPhonon : public Command
{
	CommandPhonon() : Command("phonon", "phonon")
	{
		format = "<key1> <args1...>  <key2> <args2...>  ...";
		comments =
			"Control phonon calculation and output. The possible <key>'s and their\n"
			"corresponding arguments are:\n"
			"\n+ supercell <N0> <N1> <N2>\n\n"
			"   Supercell for frozen phonon perturbation. Each entry must divide\n"
			"   the corresponding value in kpoint-folding. The k-point mesh must\n"
			"   be uniform and centered on Gamma.\n"
			"\n+ dr <dr>\n\n"
			"   Amplitude (in bohrs) of frozen phonon perturbation (default 0.1).\n"
			"\n+ iPerturbation <iPert>\n\n"
			"   Only run supercell calculation for perturbation number <iPert> (1-based).\n"
			"   Use the dry run to list out the number of perturbations, and the nStates for each\n"
			"   so that they can each be run individually with the optimum MPI configuration.\n"
			"   Use collectPerturbations below after all iPerturbation calculations complete.\n"
			"\n+ collectPerturbations\n\n"
			"   Collect results of previous individual supercell calculations.\n"
			"   Note that this requires all iPerturbation calculations (listed at\n"
			"   the end of the phonon dry run) to have already completed.\n"
			"\n+ saveHsub yes|no\n\n"
			"   Whether to compute / save phononHsub: the electron-phonon matrix elements.\n"
			"   Default: yes.\n"
			"\n+ T <T>\n\n"
			"   Temperature (in Kelvins) used for vibrational free energy estimation (default 298).\n"
			"\n+ Fcut <Fcut>\n\n"
			"   Fillings threshold to include in supercell calculation (default 1e-8).\n"
			"   The unit cell calculation may have extra bands for which matrix elements\n"
			"   are desired; this flag ensures that those extra bands do not affect the\n"
			"   performance or memory requirements of the supercell calculations.\n"
			"\n+ rSmooth <rSmooth>\n\n"
			"   Width in bohrs of the supercell boundary region over which matrix elements are smoothed.";
		
		forbid("fix-electron-density");
		forbid("fix-electron-potential");
	}

	void process(ParamList& pl, Everything& e)
	{	Phonon& phonon = ((PhononEverything&)e).phonon;
		while(true)
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
				case PM_iPerturbation:
					pl.get(phonon.iPerturbation, -1, "iPert", true);
					phonon.iPerturbation--; //internally store 0-based index
					if(phonon.iPerturbation<0)
						throw string("perturbation number must be positive");
					if(phonon.collectPerturbations)
						throw string("cannot use iPerturbation in the same calculation as collectPerturbations");
					break;
				case PM_collectPerturbations:
					phonon.collectPerturbations = true;
					if(phonon.iPerturbation>=0)
						throw string("cannot use iPerturbation in the same calculation as collectPerturbations");
					break;
				case PM_saveHsub:
					pl.get(phonon.saveHsub, true, boolMap, "saveHsub", true);
					break;
				case PM_T:
					pl.get(phonon.T, 0., "T", true);
					phonon.T *= Kelvin;
					break;
				case PM_Fcut:
					pl.get(phonon.Fcut, 0., "Fcut", true);
					if(phonon.Fcut < 0.) throw string("<Fcut> must be non-negative");
					break;
				case PM_rSmooth:
					pl.get(phonon.rSmooth, 1., "rSmooth", true);
					if(phonon.rSmooth <= 0.) throw string("<rSmooth> must be positive");
					break;
				case PM_delim: //should never be encountered
					break;
			}
		}
	}

	void printStatus(Everything& e, int iRep)
	{	const Phonon& phonon = ((const PhononEverything&)e).phonon;
		logPrintf(" \\\n\tsupercell %d %d %d", phonon.sup[0], phonon.sup[1], phonon.sup[2]);
		logPrintf(" \\\n\tdr %lg", phonon.dr);
		if(phonon.iPerturbation>=0) logPrintf(" \\\n\tiPerturbation %d", phonon.iPerturbation+1); //print 1-based index
		if(phonon.collectPerturbations) logPrintf(" \\\n\tcollectPerturbations");
		logPrintf(" \\\n\tsaveHsub %s", boolMap.getString(phonon.saveHsub));
		logPrintf(" \\\n\tT %lg", phonon.T/Kelvin);
		logPrintf(" \\\n\tFcut %lg", phonon.Fcut);
		logPrintf(" \\\n\trSmooth %lg", phonon.rSmooth);
	}
}
commandPhonon;
