/*-------------------------------------------------------------------#include <electronic/common.h>
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

#include <electronic/SelfInteractionCorrection.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/operators.h>
#include <core/Operators.h>
#include <core/Util.h>

SelfInteractionCorrection::SelfInteractionCorrection(const Everything& everything)
{
	e = &everything;	
	DC.resize(3);
}

double SelfInteractionCorrection::operator()(std::vector< diagMatrix >* correctedEigenvalues)
{
	// Loop over all quantum numbers (spin+kpoint) and bands; and corrects their eigenvalues
	double selfInteractionEnergy = 0;
	for(size_t q=0; q<e->eInfo.qnums.size(); q++)
	{	for(int iDir=0; iDir<3; iDir++)
			DC[iDir] = D(e->eVars.C[q], iDir);
		for(int n=0; n<e->eInfo.nBands; n++)
		{	double selfInteractionError = calcSelfInteractionError(q,n);
			if(correctedEigenvalues)
				correctedEigenvalues[q][n] = e->eVars.Hsub_eigs[q][n] - selfInteractionError;
			selfInteractionEnergy += e->eVars.F[q][n]*e->eInfo.qnums[q].weight*selfInteractionError;
		}
	}
	
	// Frees DC
	for(int i=0; i<3; i++)
		DC[i].free();
	
	return selfInteractionEnergy;
	
}

double SelfInteractionCorrection::calcSelfInteractionError(int q, int n)
{	
	// Get the real-space orbital density
	complexDataRptr realSpaceOrbital = I(e->eVars.C[q].getColumn(n));
	DataRptr orbitalDensity; nullToZero(orbitalDensity, e->gInfo);
	callPref(eblas_accumNorm)(e->gInfo.nr, 1., realSpaceOrbital->dataPref(), orbitalDensity->dataPref());
	DataGptr orbitalDensityTilde = J(orbitalDensity);
	
	// Calculate the Coulomb energy
	DataGptr VorbitalTilde = (*e->coulomb)(orbitalDensityTilde);
	double coulombEnergy = 0.5*dot(orbitalDensityTilde, O(VorbitalTilde));
	
	// Sets up the spin dependent density for the orbital
	DataRptrCollection orbitalSpinDensity(2);
	orbitalSpinDensity[0] = orbitalDensity;
	nullToZero(orbitalSpinDensity[1], e->gInfo);
	
	// Calculate the orbital KE if needed
	DataRptrCollection KEdensity(2);
	if(e->exCorr.needsKEdensity())
	{	nullToZero(KEdensity[0], e->gInfo); nullToZero(KEdensity[1], e->gInfo);
		for(int iDir=0; iDir<3; iDir++)
		{	DataRptr tempKE; nullToZero(tempKE, e->gInfo);
			callPref(eblas_accumNorm)(e->gInfo.nr, 1., I(DC[iDir].getColumn(n))->dataPref(), tempKE->dataPref());
			KEdensity[0] += 0.5 * tempKE;
		}
	}
	
	double xcEnergy = e->exCorr(orbitalSpinDensity, 0, e->exCorr.needsKEdensity(), &KEdensity, 0);
	
	//logPrintf("\n%i\t%i\t%f\t%f\n", q, n, coulombEnergy, xcEnergy);
	
	return coulombEnergy+xcEnergy;
}

SelfInteractionCorrection::~SelfInteractionCorrection()
{
}

double SelfInteractionCorrection::coulombExciton(int q1, int n1, int q2, int n2)
{	DataRptr density_1; nullToZero(density_1, e->gInfo);
	callPref(eblas_accumNorm)(e->gInfo.nr, 1., I(e->eVars.C[q1].getColumn(n1))->dataPref(), density_1->dataPref());
	DataRptr density_2; nullToZero(density_2, e->gInfo);
	callPref(eblas_accumNorm)(e->gInfo.nr, 1., I(e->eVars.C[q2].getColumn(n2))->dataPref(), density_2->dataPref());
	
	return 0.5*dot(J(density_2), O((*e->coulomb)(J(-density_1))));
}

void SelfInteractionCorrection::dump(const char* filename)
{
	if(e->exCorr.exxFactor())
	{	logPrintf("WARNING: Perdew-Zunger self-interaction correction can't be used with exact exchange!\n");
		return;
	}
	logPrintf("Dumping '%s'... ", filename);  logFlush();
	
	//Print to a buffer on each process:
	ostringstream oss;
	for(int q=e->eInfo.qStart; q<e->eInfo.qStop; q++)
	{	for(int iDir=0; iDir<3; iDir++)
			DC[iDir] = D(e->eVars.C[q], iDir);
		for(int n=0; n<e->eInfo.nBands; n++)
		{	double selfInteractionError = calcSelfInteractionError(q,n);
			oss << q
				<< '\t' << n
				<< '\t' << e->eVars.Hsub_eigs[q][n]
				<< '\t' << selfInteractionError
				<< '\t' << e->eVars.Hsub_eigs[q][n]-(e->eVars.F[q][n] ? 1. : 0.)*selfInteractionError
				<< '\t' << e->eVars.Hsub_eigs[q][n]-selfInteractionError
				<< '\n';
		}
	}
	if(!mpiUtil->isHead()) mpiUtil->send(oss.str(), 0, 0); //send to head for writing to file
	
	if(mpiUtil->isHead()) //Write to file
	{	FILE* fp = fopen(filename, "w");
		if(!fp) die("Error opening %s for writing.\n", filename);
		fprintf(fp, "WARNING: Self-Interaction Correction scheme is still experimental, take extreme care when using the numbers below!\n");
		fprintf(fp, "q\tn\t\"KS Eigenvalue\"\t\"Self-Interaction Error\"\t\"Corrected Eigenvalue\"\t\"Corrected Eigenvalue (no-fillings)\"\n");
		fputs(oss.str().c_str(), fp); //output own buffer
		for(int iSrc=1; iSrc<mpiUtil->nProcesses(); iSrc++)
		{	string buf;
			mpiUtil->recv(buf, iSrc, 0);
			fputs(buf.c_str(), fp); //output others' buffers in order
		}
		fclose(fp);
	}
	logPrintf("done\n"); logFlush();
}
