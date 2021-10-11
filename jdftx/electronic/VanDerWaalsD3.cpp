/*-------------------------------------------------------------------
Copyright 2021 Ravishankar Sundararaman.

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

#include <electronic/VanDerWaalsD3.h>
#include <electronic/VanDerWaalsD3_data.h>
#include <electronic/Everything.h>

namespace D3
{
	AtomParams getAtomParams(int Z)
	{	if((Z < 1) or (Z > D3::Zmax))
			die("\nAtomic number %d not supported in DFT-D3\n\n", Z);
		AtomParams ap;
		ap.Z = Z;
		ap.sqrtQ = sqrtQ[Z - 1]; //Note: per-atom arrays have 0-based index
		ap.k2Rcov = k2Rcov[Z - 1];
		for(int iCN=0; iCN<D3::C6::nCN; iCN++)
			if(Z == D3::C6::Z[iCN])
			{	ap.CN.push_back(D3::C6::CN[iCN]);
				ap.iCN.push_back(iCN);
			}
		return ap;
	}
	
	PairParams getPairParams(const AtomParams& ap1, const AtomParams& ap2)
	{	PairParams pp;
		//Radius:
		int iZmin = std::min(ap1.Z, ap2.Z) - 1; //index into arrays by lower of two atomic numbers
		int iZmax = std::max(ap1.Z, ap2.Z) - 1;
		int iZ12 = (iZmax*(iZmax+1))/2 + iZmin; //index into unique Z pairs
		pp.R0 = D3::R0_Angstrom[iZ12] * Angstrom;
		//C6 coefficients:
		pp.C6.init(ap1.CN.size(), ap2.CN.size());
		complex* C6data = pp.C6.data();
		for(int iCN1: ap1.iCN)
			for(int iCN2: ap2.iCN)
			{	int iCNmin = std::min(iCN1, iCN2);
				int iCNmax = std::max(iCN1, iCN2);
				int iCN12 = (iCNmax*(iCNmax+1))/2 + iCNmin; //index into unique iCN pairs
				*(C6data++) = D3::C6::coeff[iCN12];
			}
		return pp;
	}
}

//----- Implementation of class VanDerWaalsD3 -----

VanDerWaalsD3::VanDerWaalsD3(const Everything& e) : VanDerWaals(e)
{
	logPrintf("\nInitializing DFT-D3 calculator:\n");
	
	//Get parameters for exchange-correlation functional
	string xcName = e.exCorr.getName();
	D3::setXCscale(xcName, s6, sr6, s8, sr8);
	logPrintf("\tParameters set for %s functional\n", xcName.c_str());
	logPrintf("\ts6: %6.3lf  s_r6: %6.3lf\n", s6, sr6);
	logPrintf("\ts8: %6.3lf  s_r8: %6.3lf\n", s8, sr8);

	//Get per-atom parameters:
	logPrintf("\tPer-atom parameters loaded for:\n");
	for(size_t spIndex=0; spIndex<e.iInfo.species.size(); spIndex++)
	{	const auto& sp = e.iInfo.species[spIndex];
		assert(sp->atomicNumber);
		D3::AtomParams ap = D3::getAtomParams(sp->atomicNumber);
		atomParams.push_back(ap);
		logPrintf("\t%2s:  sqrtQ[a0]: %6.3f  Rcov[a0]: %6.3f  CN: [",
			sp->name.c_str(), ap.sqrtQ, ap.k2Rcov / D3::k2);
		for(double CNi: ap.CN) logPrintf(" %.2f", CNi);
		logPrintf(" ]\n");
	}

	//Get per-pair paramaters:
	pairParams.resize(atomParams.size());
	for(size_t iSp1=0; iSp1<atomParams.size(); iSp1++)
	{	pairParams[iSp1].resize(atomParams.size());
		for(size_t iSp2=0; iSp2<atomParams.size(); iSp2++)
			pairParams[iSp1][iSp2] = D3::getPairParams(atomParams[iSp1], atomParams[iSp2]);
	}

	Citations::add("DFT-D3 dispersion correction", "S. Grimme, J. Antony, S. Ehrlich and H. Krieg, J. Chem. Phys. 132, 154104 (2010)");
}


double VanDerWaalsD3::getScaleFactor(string exCorrName, double scaleOverride) const
{	return 0.;  //global scale factor not used in D3
}


double VanDerWaalsD3::energyAndGrad(std::vector<Atom>& atoms, const double scaleFac, matrix3<>* E_RRT) const
{	die("Not yet implemented.\n");
	return 0.;
}



