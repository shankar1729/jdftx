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
#include <core/Operators.h>
#include <core/Util.h>

DumpSelfInteractionCorrection::DumpSelfInteractionCorrection(const Everything& e) : e(e)
{
}

double DumpSelfInteractionCorrection::operator()(std::vector<diagMatrix>* correctedEigenvalues)
{
	// Loop over all quantum numbers (spin+kpoint) and bands; and corrects their eigenvalues
	double selfInteractionEnergy = 0;
	DC.resize(3);
	for(int q=0; q<e.eInfo.nStates; q++) //need to loop over all states since XC is MPI parallelized
	{	if(e.exCorr.needsKEdensity() && e.eInfo.isMine(q))
			for(int iDir=0; iDir<3; iDir++)
				DC[iDir] = D(e.eVars.C[q], iDir);
		if(e.eInfo.isMine(q) && correctedEigenvalues)
			(*correctedEigenvalues)[q].resize(e.eInfo.nBands);
		for(int n=0; n<e.eInfo.nBands; n++)
		{	double selfInteractionError = calcSelfInteractionError(q,n);
			if(e.eInfo.isMine(q))
			{	if(correctedEigenvalues)
					(*correctedEigenvalues)[q][n] = e.eVars.Hsub_eigs[q][n] - selfInteractionError;
				selfInteractionEnergy += e.eVars.F[q][n]*e.eInfo.qnums[q].weight*selfInteractionError;
			}
		}
	}
	DC.clear();
	mpiWorld->allReduce(selfInteractionEnergy, MPIUtil::ReduceSum);
	return selfInteractionEnergy;
	
}

double DumpSelfInteractionCorrection::calcSelfInteractionError(int q, int n)
{
	// Get the real-space orbital density
	ScalarFieldArray orbitalDensity; nullToZero(orbitalDensity, e.gInfo, 2);
	if(e.eInfo.isMine(q))
	{	QuantumNumber qnum; qnum.weight = 1.;
		diagMatrix Fqn = eye(1);
		ColumnBundle Cqn = e.eVars.C[q].getSub(n,n+1);
		std::vector<matrix> VdagCqn;
		for(const matrix& m: e.eVars.VdagC[q])
			VdagCqn.push_back(m(0, m.nRows(), n, n+1));
		orbitalDensity[0] = diagouterI(Fqn, Cqn, 1, &e.gInfo)[0];
		e.iInfo.augmentDensityInit();
		e.iInfo.augmentDensitySpherical(qnum, Fqn, VdagCqn); //pseudopotential contribution
		e.iInfo.augmentDensityGrid(orbitalDensity);
	}
	orbitalDensity[0]->bcastData(mpiWorld, e.eInfo.whose(q));
	ScalarFieldTilde orbitalDensityTilde = J(orbitalDensity[0]);
	
	// Calculate the Coulomb energy
	ScalarFieldTilde VorbitalTilde = (*e.coulomb)(orbitalDensityTilde);
	double coulombEnergy = 0.5*dot(orbitalDensityTilde, O(VorbitalTilde));
	
	// Calculate the orbital KE if needed
	ScalarFieldArray KEdensity(2);
	if(e.exCorr.needsKEdensity())
	{	nullToZero(KEdensity, e.gInfo);
		if(e.eInfo.isMine(q))
			for(int iDir=0; iDir<3; iDir++)
				KEdensity[0] += 0.5 * diagouterI(eye(1), DC[iDir].getSub(n,n+1), 1, &e.gInfo)[0];
		KEdensity[0]->bcastData(mpiWorld, e.eInfo.whose(q));
	}
	double xcEnergy = e.exCorr(orbitalDensity, 0, IncludeTXC(), &KEdensity, 0);
	return coulombEnergy + xcEnergy;
}

void DumpSelfInteractionCorrection::dump(const char* filename)
{
	if(e.exCorr.exxFactor())
	{	logPrintf("WARNING: Perdew-Zunger self-interaction correction can't be used with exact exchange. (Skipping)\n");
		return;
	}
	if(e.eInfo.isNoncollinear())
	{	logPrintf("WARNING: Perdew-Zunger self-interaction correction can't be used with noncollinear spins. (Skipping)\n");
		return;
	}
	logPrintf("Dumping '%s'... ", filename);  logFlush();
	
	// Calculate correction
	std::vector<diagMatrix> correctedEigenvalues(e.eInfo.qnums.size());
	double Esic = (*this)(&correctedEigenvalues);
	
	e.eInfo.write(correctedEigenvalues, filename);
	logPrintf("done\n");
	logPrintf("\tSelf-interaction energy: %.15lf\n", Esic);
}
