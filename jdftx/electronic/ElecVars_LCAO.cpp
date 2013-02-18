/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman

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

#include <electronic/Everything.h>
#include <electronic/ElecVars.h>
#include <electronic/ColumnBundle.h>
#include <electronic/matrix.h>
#include <electronic/operators.h>

int ElecVars::LCAO()
{
	const ElecInfo &eInfo = e->eInfo;
	const IonInfo& iInfo = e->iInfo;
	
	//Count total atomic orbitals:
	int nAtomic = 0;
	for(auto sp: iInfo.species)
		nAtomic += sp->nAtomicOrbitals();
	if(nAtomic)
	{	logPrintf("linear combination of atomic orbitals\n");
		if(nAtomic < eInfo.nBands)
			logPrintf("Note: number of bands (%d) exceeds available atomic orbitals (%d)\n",
				eInfo.nBands, nAtomic);
	}
	else
	{	logPrintf("(No orbitals for LCAO)  ");
		return 0;
	}
	
	//Check exchange-correlation functional, and replace with PBE if not strictly (semi-)local
	ExCorr exCorrPBE; //defaults to gga-PBE
	const ExCorr* exCorr; //The functional to be used here: either e->exCorr or PBE
	if(e->exCorr.exxFactor() || e->exCorr.needsKEdensity()) //Hybrid or meta-GGA respectively
	{	logPrintf("Initializing semi-local functional for LCAO:\n");
		exCorrPBE.setup(*e);
		exCorr = &exCorrPBE;
	}
	else exCorr = &e->exCorr;
	
	//Get electron density obtained by adding those of the atoms:
	DataGptr nTilde;
	for(auto sp: e->iInfo.species)
		sp->accumulateAtomicDensity(nTilde);
	DataRptr V_LCAO; (*exCorr)(I(nTilde)+iInfo.nCore, &V_LCAO); //get exchange-correlation potential
	V_LCAO += I(iInfo.Vlocps + (*e->coulomb)(nTilde + rhoExternal), true); //add local-pseudopotential and mean-field electrostatic terms
	
	//Initialize one state at a time:
	for(int q=0; q<eInfo.nStates; q++)
	{	ColumnBundle psi(nAtomic, e->basis[q].nbasis, &e->basis[q], &eInfo.qnums[q], isGpuEnabled());
		//Fill with atomic orbitals:
		int iCol=0;
		for(auto sp: iInfo.species)
		{	sp->setAtomicOrbitals(psi, iCol);
			iCol += sp->nAtomicOrbitals();
		}
		psi = psi * invsqrt(psi^O(psi));
		//Compute the Hamiltonian:
		ColumnBundle Hpsi = -0.5*L(psi); //kinetic part
		iInfo.EnlAndGrad(eye(nAtomic), psi, Hpsi); //non-local pseudopotentials
		Hpsi += Idag_DiagV_I(psi, JdagOJ(V_LCAO + (Vexternal.size() ? Vexternal[eInfo.qnums[q].index()] : DataRptr())));
		//Diagonalize and set to eigenvectors:
		matrix evecs; diagMatrix eigs;
		(psi^Hpsi).diagonalize(evecs, eigs);
		psi = psi * evecs;
		//Pick the first nBands columns (which have the lowest eigenvalues:
		Y[q].setSub(0, psi);
	}
	return std::min(eInfo.nBands, nAtomic);
}
