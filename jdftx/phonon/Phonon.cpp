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
#include <core/Units.h>

Phonon::Phonon()
: dr(0.001), T(298*Kelvin)
{

}

void Phonon::setup()
{
	//Ensure phonon command specified:
	if(!sup.length())
		die("phonon supercell must be specified using the phonon command.\n");
	//Check symmetries
	if(e.symm.mode != SymmetriesNone)
		die("phonon does not support symmetries.\n");
	//Check kpoint and supercell compatibility:
	if(e.eInfo.qnums.size()>1 || e.eInfo.qnums[0].k.length_squared())
		die("phonon requires a Gamma-centered uniform kpoint mesh.\n");
	for(int j=0; j<3; j++)
	{	if(e.eInfo.kfold[j] % sup[j])
		{	die("kpoint folding %d is not a multiple of supercell count %d for lattice direction %d.\n",
				e.eInfo.kfold[j], sup[j], j);
		}
		eSup.eInfo.kfold[j] = e.eInfo.kfold[j] / sup[j];
	}
	//Grid:
	eSup.gInfo.S = Diag(sup) * e.gInfo.S; //for manual fftbox if any
	eSup.gInfo.R = e.gInfo.R * Diag(sup);
	int prodSup = sup[0] * sup[1] * sup[2];
	
	logPrintf("########### Initializing JDFTx for unit cell #############\n");
	SpeciesInfo::Constraint constraintFull;
	constraintFull.moveScale = 0;
	constraintFull.type = SpeciesInfo::Constraint::None;
	for(size_t sp=0; sp<e.iInfo.species.size(); sp++)
		e.iInfo.species[sp]->constraints.assign(e.iInfo.species[sp]->atpos.size(), constraintFull);
	e.setup();
	
	logPrintf("########### Initializing JDFTx for supercell #############\n");
	//Replicate atoms:
	for(size_t sp=0; sp<e.iInfo.species.size(); sp++)
	{	eSup.iInfo.species[sp]->atpos.clear();
		eSup.iInfo.species[sp]->initialMagneticMoments.clear();
		matrix3<> invSup = inv(Diag(vector3<>(sup)));
		vector3<int> iR;
		for(iR[0]=0; iR[0]<sup[0]; iR[0]++)
		for(iR[1]=0; iR[1]<sup[1]; iR[1]++)
		for(iR[2]=0; iR[2]<sup[2]; iR[2]++)
			for(vector3<> pos: e.iInfo.species[sp]->atpos)
				eSup.iInfo.species[sp]->atpos.push_back(invSup * (pos + iR));
		eSup.iInfo.species[sp]->constraints.assign(eSup.iInfo.species[sp]->atpos.size(), constraintFull);
	}
	//Remove initial state settings (incompatible with supercell):
	eSup.eVars.wfnsFilename.clear();
	eSup.eVars.HauxFilename.clear();
	eSup.eVars.eigsFilename.clear();
	eSup.eVars.fluidInitialStateFilename.clear();
	eSup.eInfo.initialFillingsFilename.clear();
	eSup.scfParams.historyFilename.clear();
	//ElecInfo:
	eSup.eInfo.nBands = e.eInfo.nBands * prodSup;
	//ElecVars:
	eSup.eVars.initLCAO = false;
	eSup.setup();
}
