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

#include <electronic/ExCorr_OrbitalDep_GLLBsc.h>
#include <electronic/ColumnBundle.h>
#include <electronic/Everything.h>
#include <electronic/operators.h>

ExCorr_OrbitalDep_GLLBsc::ExCorr_OrbitalDep_GLLBsc(const Everything& e) : ExCorr::OrbitalDep(e)
{
}

DataRptrCollection ExCorr_OrbitalDep_GLLBsc::getPotential() const
{	int nSpins = e.eVars.n.size();
	if(!e.eVars.Hsub_eigs[e.eInfo.qStart].size()) return DataRptrCollection(nSpins); //no eigenvalues yet
	//Detect HOMO:
	std::vector<double> eHOMO(nSpins, -DBL_MAX);
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{	int s = e.eInfo.qnums[q].index();
		for(int b=0; b<e.eInfo.nBands; b++)
			if(e.eVars.F[q][b]>0.5)
				eHOMO[s] = std::max(eHOMO[s], e.eVars.Hsub_eigs[q][b]);
	}
	mpiUtil->allReduce(eHOMO.data(), eHOMO.size(), MPIUtil::ReduceMax);
	//Compute exchange response potential:
	const double Kx = 8*sqrt(2)/(3*M_PI*M_PI);
	DataRptrCollection Vx(nSpins);
	e.iInfo.augmentDensityInit();
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{	const QuantumNumber& qnum = e.eInfo.qnums[q];
		int s = qnum.index();
		diagMatrix Feff(e.eInfo.nBands);
		for(int b=0; b<e.eInfo.nBands; b++)
			Feff[b] = e.eVars.F[q][b] * Kx * sqrt(std::max(0., eHOMO[s]-e.eVars.Hsub_eigs[q][b]));
		Vx[s] += qnum.weight * diagouterI(Feff, e.eVars.C[q], &e.gInfo); //without the 1/n(r) denominator
		e.iInfo.augmentDensitySpherical(qnum, Feff, e.eVars.VdagC[q]); //ultrasoft contribution
	}
	e.iInfo.augmentDensityGrid(Vx);
	for(int s=0; s<nSpins; s++)
	{	nullToZero(Vx[s], e.gInfo);
		Vx[s]->allReduce(MPIUtil::ReduceSum);
		e.symm.symmetrize(Vx[s]);
		Vx[s] *= inv(e.eVars.n[s]); //denominator added here
	}
	return Vx;
}

void ExCorr_OrbitalDep_GLLBsc::dump() const
{
}