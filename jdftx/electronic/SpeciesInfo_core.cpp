/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

#include <electronic/SpeciesInfo.h>
#include <electronic/Everything.h>
#include <electronic/RadialSchrodinger.h>

RadialFunctionR getTau(const RadialFunctionR& n)
{	RadialFunctionR tau = n;
	for(size_t i=0; i<n.r.size(); i++)
	{	double iPlus = i+1<n.r.size() ? i+1 : i;
		double iMinus = i ? i-1 : i;
		double Dnbyn = log(n.f[iPlus]/n.f[iMinus]) / (n.r[i] * log(n.r[iPlus]/n.r[iMinus]));
		double tauTF = 0.3*n.f[i]*pow(3*M_PI*M_PI*n.f[i], 2./3); //Thomas-Fermi KE density
		double tauVW = 0.125 * Dnbyn * Dnbyn * n.f[i]; //von-Weisacker KE density
		tau.f[i] = tauTF + tauVW;
	}
	return tau;
}


void SpeciesInfo::setCore(RadialFunctionR& nCore)
{
	int nGridLoc = int(ceil(e->gInfo.GmaxGrid/dGloc))+5;
	
	//Truncate the radial grid where the integral of nCore is within double precision of its untruncated value:
	double nCoreTot = 0.;
	for(size_t i=0; i<nCore.r.size(); i++)
		nCoreTot += nCore.f[i] * (4*M_PI*nCore.r[i]*nCore.r[i]) * nCore.dr[i];
	double nCoreBack = 0.;
	for(size_t i=nCore.r.size()-1; ; i--)
	{	nCoreBack += nCore.f[i] * (4*M_PI*nCore.r[i]*nCore.r[i]) * nCore.dr[i];
		if(nCoreBack > 1e-16*nCoreTot)
		{	nCore.f.resize(i);
			nCore.r.resize(i);
			nCore.dr.resize(i);
			break;
		}
	}
	
	//Check if KE density is required:
	bool needTau = e->exCorr.needsKEdensity();
	for(auto ec: e->exCorrDiff)
		needTau |= ec->needsKEdensity();
	
	if(needTau)
	{	RadialFunctionR tauCore = getTau(nCore);
		logPrintf("  Transforming core KE density to a uniform radial grid of dG=%lg with %d points.\n",
			dGloc, nGridLoc);
		tauCore.transform(0, dGloc, nGridLoc, tauCoreRadial);
	}
	
	logPrintf("  Transforming core density to a uniform radial grid of dG=%lg with %d points.\n",
		dGloc, nGridLoc);
	nCore.transform(0, dGloc, nGridLoc, nCoreRadial);
}
