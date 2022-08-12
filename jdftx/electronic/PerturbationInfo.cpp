/*
 * PerturbationInfo.cpp
 *
 *  Created on: Jul 25, 2022
 *      Author: brandon
 */

#include "PerturbationInfo.h"
#include <core/ScalarFieldIO.h>
#include <electronic/Everything.h>
#include <electronic/ElecVars.h>

PerturbationInfo::PerturbationInfo() {
}

PerturbationInfo::~PerturbationInfo() {
	// TODO Auto-generated destructor stub
}

void PerturbationInfo::setup(Everything *e, ElecVars *eVars) {
	if(dVexternalFilename.size())
	{	dVext.resize(dVexternalFilename.size());
		for(unsigned s=0; s<dVext.size(); s++)
		{	dVext[s] = ScalarFieldData::alloc(e->gInfo);
			logPrintf("Reading external potential from '%s'\n", dVexternalFilename[s].c_str());
			loadRawBinary(dVext[s], dVexternalFilename[s].c_str());
		}
		if(dVext.size()==1 && eVars->n.size()==2) //Replicate potential for second spin:
			dVext.push_back(dVext[0]->clone());
		if(e->iInfo.computeStress)
			die("\nStress calculation not supported with external potentials.\n\n");
	}
}
