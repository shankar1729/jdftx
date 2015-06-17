/*-------------------------------------------------------------------
Copyright 2012 Deniz Gunceler

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

#include <electronic/Dump_internal.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/operators.h>
#include <core/Operators.h>
#include <core/Util.h>

DumpSelfInteractionCorrection::DumpSelfInteractionCorrection(const Everything& everything)
{
	e = &everything;	
	DC.resize(3);
}

double DumpSelfInteractionCorrection::operator()(std::vector< diagMatrix >* correctedEigenvalues)
{
	// Loop over all quantum numbers (spin+kpoint) and bands; and corrects their eigenvalues
	double selfInteractionEnergy = 0;
	for(int q=e->eInfo.qStart; q<e->eInfo.qStop; q++)
	{	for(int iDir=0; iDir<3; iDir++)
			DC[iDir] = D(e->eVars.C[q], iDir);
		for(int n=0; n<e->eInfo.nBands; n++)
		{	double selfInteractionError = calcSelfInteractionError(q,n);
			if(correctedEigenvalues)
			{
				(*correctedEigenvalues)[q][n] = e->eVars.Hsub_eigs[q][n] - selfInteractionError;
			}
			selfInteractionEnergy += e->eVars.F[q][n]*e->eInfo.qnums[q].weight*selfInteractionError;
		}
	}
	
	// Frees DC
	for(int i=0; i<3; i++)
		DC[i].free();
	
	return selfInteractionEnergy;
	
}

double DumpSelfInteractionCorrection::calcSelfInteractionError(int q, int n)
{	
	// Get the real-space orbital density
	complexScalarField realSpaceOrbital = I(e->eVars.C[q].getColumn(n,0));
	ScalarField orbitalDensity; nullToZero(orbitalDensity, e->gInfo);
	callPref(eblas_accumNorm)(e->gInfo.nr, 1., realSpaceOrbital->dataPref(), orbitalDensity->dataPref());
	ScalarFieldTilde orbitalDensityTilde = J(orbitalDensity);
	
	// Calculate the Coulomb energy
	ScalarFieldTilde VorbitalTilde = (*e->coulomb)(orbitalDensityTilde);
	double coulombEnergy = 0.5*dot(orbitalDensityTilde, O(VorbitalTilde));
	
	// Sets up the spin dependent density for the orbital
	ScalarFieldArray orbitalSpinDensity(2);
	orbitalSpinDensity[0] = orbitalDensity;
	nullToZero(orbitalSpinDensity[1], e->gInfo);
	
	// Calculate the orbital KE if needed
	ScalarFieldArray KEdensity(2);
	if(e->exCorr.needsKEdensity())
	{	nullToZero(KEdensity[0], e->gInfo); nullToZero(KEdensity[1], e->gInfo);
		for(int iDir=0; iDir<3; iDir++)
		{	ScalarField tempKE; nullToZero(tempKE, e->gInfo);
			callPref(eblas_accumNorm)(e->gInfo.nr, 1., I(DC[iDir].getColumn(n,0))->dataPref(), tempKE->dataPref());
			KEdensity[0] += 0.5 * tempKE;
		}
	}
	
	double xcEnergy = e->exCorr(orbitalSpinDensity, 0, e->exCorr.needsKEdensity(), &KEdensity, 0);
	
	//logPrintf("\n%i\t%i\t%f\t%f\n", q, n, coulombEnergy, xcEnergy);
	
	return coulombEnergy+xcEnergy;
}

DumpSelfInteractionCorrection::~DumpSelfInteractionCorrection()
{
}

double DumpSelfInteractionCorrection::coulombExciton(int q1, int n1, int q2, int n2)
{	ScalarField density_1; nullToZero(density_1, e->gInfo);
	callPref(eblas_accumNorm)(e->gInfo.nr, 1., I(e->eVars.C[q1].getColumn(n1,0))->dataPref(), density_1->dataPref());
	ScalarField density_2; nullToZero(density_2, e->gInfo);
	callPref(eblas_accumNorm)(e->gInfo.nr, 1., I(e->eVars.C[q2].getColumn(n2,0))->dataPref(), density_2->dataPref());
	
	return 0.5*dot(J(density_2), O((*e->coulomb)(J(-density_1))));
}

void DumpSelfInteractionCorrection::dump(const char* filename)
{
	if(e->exCorr.exxFactor())
	{	logPrintf("WARNING: Perdew-Zunger self-interaction correction can't be used with exact exchange. (Skipping)\n");
		return;
	}
	if(e->eInfo.isNoncollinear())
	{	logPrintf("WARNING: Perdew-Zunger self-interaction correction can't be used with noncollinear spins. (Skipping)\n");
		return;
	}
	logPrintf("Dumping '%s'... ", filename);  logFlush();
	
	// Calculate correction
	std::vector< diagMatrix > correctedEigenvalues;
	correctedEigenvalues.resize(e->eInfo.qnums.size(), e->eInfo.nBands);
	(*this)(&correctedEigenvalues);
	
	e->eInfo.write(correctedEigenvalues, filename);
}
