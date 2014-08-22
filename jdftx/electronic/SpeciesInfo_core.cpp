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

RadialFunctionR getTau(const RadialFunctionR& n, double rCut)
{	RadialFunctionR tau = n;
	size_t iMaxLast = 0; //location of last maxima
	for(size_t i=0; i<n.r.size(); i++)
	{	size_t iPlus = i+1<n.r.size() ? i+1 : i;
		size_t iMinus = i ? i-1 : i;
		double Dnbyn = log(n.f[iPlus]/n.f[iMinus]) / (n.r[i] * log(n.r[iPlus]/n.r[iMinus]));
		double tauTF = 0.3*n.f[i]*pow(3*M_PI*M_PI*n.f[i], 2./3); //Thomas-Fermi KE density
		double tauVW = 0.125 * Dnbyn * Dnbyn * n.f[i]; //von-Weisacker KE density
		tau.f[i] = tauTF + tauVW;
		if(i && (tau.f[i] > tau.f[i-1]) && (tau.f[i]>1e-3)) iMaxLast = i;
	}
	//Re-pseudize small r regions:
	if(!rCut) rCut = std::max(0.8, tau.r[iMaxLast] * 1.5);
	size_t iCut = 0; while(tau.r[iCut]<rCut) iCut++;
	rCut = tau.r[iCut]; //align to grid point
	double tauCut = tau.f[iCut];
	double tauPrimeCut = (tau.f[iCut+1] - tau.f[iCut-1]) / (tau.r[iCut] * log(tau.r[iCut+1]/tau.r[iCut-1]));
	double f2 = -tauPrimeCut/(2.*rCut); assert(f2 > 0);
	double f0 = tauCut + f2*rCut*rCut;
	for(size_t i=0; i<iCut; i++)
		tau.f[i] = f0 - f2*tau.r[i]*tau.r[i];
	return tau;
}


void SpeciesInfo::setCore(RadialFunctionR& nCore)
{
	if(e->exCorr.orbitalDep && e->exCorr.orbitalDep->ignore_nCore())
	{	logPrintf("  WARNING: Ignoring core density because that is not supported by the orbital-dependent functional.\n");
		return;
	}
	
	const double dG = e->gInfo.dGradial;
	int nGridLoc = int(ceil(e->gInfo.GmaxGrid/dG))+5;
	
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
	{	RadialFunctionR tauCore = getTau(nCore, tauCore_rCut);
		logPrintf("  Transforming core KE density to a uniform radial grid of dG=%lg with %d points.\n",
			dG, nGridLoc);
		tauCore.transform(0, dG, nGridLoc, tauCoreRadial);
		
		if(tauCorePlot)
		{	FILE* fp = fopen((name+".tauCoreRadial").c_str(), "w");
			fprintf(fp, "#  r      tauCore     nCore\n");
			for(size_t i=0; i<nCore.r.size(); i++)
				fprintf(fp, "%lg\t%le\t%le\n", nCore.r[i], tauCore.f[i], nCore.f[i]);
			fclose(fp);
		}
	}
	
	logPrintf("  Transforming core density to a uniform radial grid of dG=%lg with %d points.\n",
		dG, nGridLoc);
	nCore.transform(0, dG, nGridLoc, nCoreRadial);
}
