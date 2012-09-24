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

#include <electronic/VanDerWaals.h>
#include "Everything.h"

const static int maximumSupportedAtomicNumber = 54;

inline std::pair<double, double> calculateEnergyAndForceMagnitude(double Rij, double C6, double R0)
{
	//! This gets rid of the singularity in the Grimme vdw functional.
	//! All values less than the maximum of the Grimme vdw potential are set to the maximum
	if(Rij/R0 < 0.3000002494598603)
		return std::pair<double, double>(0.00114064201325433, 0.);
	
	double d = 20.;
	
	double exponential = exp(-d*(1.+Rij/R0));
	double fdamp = 1./(1. + exponential);
	double grad_fdamp = (d/R0) * pow(fdamp, 2) * exponential;
	
	double Rijn6 = 1./pow(Rij, 6);
	double grad_Rijn6 = -6.*Rijn6/Rij;
	
	double Energy = fdamp*C6*Rijn6;
	double ForceMag = C6*(fdamp*grad_Rijn6 + Rijn6*grad_fdamp);
	
	return std::pair<double, double>(Energy, ForceMag);
	
}


std::pair<double, double> VanDerWaals::getC6AndR0(int atomicNumber)
{
	if(atomicNumber > maximumSupportedAtomicNumber)
		die("\nAtomic numbers greater than %i are not yet supported!\n", maximumSupportedAtomicNumber);
	if(atomicNumber <= 0)
		die("\nYou can't have a negative or zero atomic number!\n")
		
	
	return std::pair<double, double>(vdwParams[atomicNumber].C6, vdwParams[atomicNumber].R0);

}

void VanDerWaals::setup(const Everything &everything){
  
	logPrintf("\nInitializing van der Waals corrections...\t");
	
	e = &everything;
	
	//! Checks for the pseudopotential format.  Only Fhi and USPP are allowed.
	for(unsigned sp=0; sp<e->iInfo.species.size(); sp++)
	{	
		SpeciesInfo& spInfo = *(e->iInfo.species[sp]);
		if(not (spInfo.getPSPFormat() != SpeciesInfo::Fhi or spInfo.getPSPFormat() != SpeciesInfo::Uspp))
			die("\nYou can use the Van der Waals corrections only for .fhi and .uspp pseudopotentials!\n");
	}
	
   //! Constructs the EXCorr -> scaling factor map
   scalingFactor["gga-PBE"] = 0.75;
   scalingFactor["hyb-gga-xc-b3lyp"] = 1.05;
   scalingFactor["mgga-TPSS"] = 1.;
   
	//! Checks whether the used EXCorr is supported
	if(scalingFactor.find(e->exCorr.getName()) == scalingFactor.end())
		die("\n%s Exchange-Correlation is not supported by Grimme Van der Waals corrections!\n", e->exCorr.getName().c_str());	
	
	//! Checks to see if the elec-xc-compare command has any unsupported EXCorr
	if(e->exCorrDiff.size())
	{
		string EXCorr;
		for(int counter = 0; counter < (int) e->exCorrDiff.size(); counter++){
			EXCorr = e->exCorrDiff[counter].get()->getName();
			if(scalingFactor.find(EXCorr) == scalingFactor.end())
				die("\nYou have chosen to compare the energy calculated with %s exchange-correlation\nVan der Waals Corrections are not supported for %s!\n", EXCorr.c_str(), EXCorr.c_str());
		}
	}
	
	//! Sets up the C6 and R0 parameters
   vdwParams.resize(maximumSupportedAtomicNumber+1);
//! vdwParams[0] = {C6, R0}
	vdwParams[1] = {0.14 , 1.001 };
	vdwParams[2] = {0.08 , 1.012 };
	vdwParams[3] = {1.61 , 0.825 };
	vdwParams[4] = {1.61 , 1.408 };
	vdwParams[5] = {3.13 , 1.485 };
	vdwParams[6] = {1.75 , 1.452 };
	vdwParams[7] = {1.23 , 1.397 };
	vdwParams[8] = {0.70 , 1.342 };
	vdwParams[9] = {0.75 , 1.287 };
	vdwParams[10] = {0.63 , 1.243 };
	vdwParams[11] = {5.71 , 1.144 };
	vdwParams[12] = {5.71 , 1.364 };
	vdwParams[13] = {10.79,1.639  };
	vdwParams[14] = {9.23 , 1.716 };
	vdwParams[15] = {7.84 , 1.705 };
	vdwParams[16] = {5.57 , 1.683 };
	vdwParams[17] = {5.07 , 1.639 };
	vdwParams[18] = {4.61 , 1.595 };
	vdwParams[19] = {10.80 , 1.485};
	vdwParams[20] = {10.80 , 1.474};
	vdwParams[21] = {10.80 , 1.562};
	vdwParams[22] = {10.80 , 1.562};
	vdwParams[23] = {10.80 , 1.562};
	vdwParams[24] = {10.80 , 1.562};
	vdwParams[25] = {10.80 , 1.562};
	vdwParams[26] = {10.80 , 1.562};
	vdwParams[27] = {10.80 , 1.562};
	vdwParams[28] = {10.80 , 1.562};
	vdwParams[29] = {10.80 , 1.562};
	vdwParams[30] = {10.80 , 1.562};
	vdwParams[31] = {16.99 , 1.650};
	vdwParams[32] = {17.10 , 1.727};
	vdwParams[33] = {16.37 , 1.760};
	vdwParams[34] = {12.64 , 1.771};
	vdwParams[35] = {12.47 , 1.749};
	vdwParams[36] = {12.01 , 1.727};
	vdwParams[37] = {24.67 , 1.628};
	vdwParams[38] = {24.67 , 1.606};
	vdwParams[39] = {24.67 , 1.639};
	vdwParams[40] = {24.67 , 1.639};
	vdwParams[41] = {24.67 , 1.639};
	vdwParams[42] = {24.67 , 1.639};
	vdwParams[43] = {24.67 , 1.639};
	vdwParams[44] = {24.67 , 1.639};
	vdwParams[45] = {24.67 , 1.639};
	vdwParams[46] = {24.67 , 1.639};
	vdwParams[47] = {24.67 , 1.639};
	vdwParams[48] = {24.67 , 1.639};
	vdwParams[49] = {37.32 , 1.672};
	vdwParams[50] = {38.71 , 1.804};
	vdwParams[51] = {38.44 , 1.881};
	vdwParams[52] = {31.74 , 1.892};
	vdwParams[53] = {31.50 , 1.892};
	vdwParams[54] = {29.99 , 1.881};
	
	logPrintf(" Done!\n");
}

double VanDerWaals::VDWEnergyAndGrad(std::vector<Atom>& atoms, string EXCorr)
{
	double Energy = 0.;  //! Total VDW Energy
	
	std::pair<double, double> temp, c1Params, c2Params;
	
	double Rij, C6, R0;
	vector3<> deltaR, unit; //! Unit vector from atom[c1] to atom[c2]
	
	int n[3]; //! How many units cells are summed in each direction
	n[0] = (int) 100./(e->gInfo.R * vector3<>(1.,0.,0.)).length(); // 100 Bohrs is approximately a picoHartree of energy
	n[1] = (int) 100./(e->gInfo.R * vector3<>(0.,1.,0.)).length();
	n[2] = (int) 100./(e->gInfo.R * vector3<>(0.,0.,1.)).length();
	
	//! Checks for Coulomb truncation
	switch(e->coulombParams.geometry)
	{
		case e->coulombParams.Slab:
			n[e->coulombParams.iDir] = 0;
		case e->coulombParams.Wire:
			n[(e->coulombParams.iDir+1)%3] = 0;
			n[(e->coulombParams.iDir+2)%3] = 0;
		default:
			break;
	}
	
	for(int c1 = 0; c1 < (int) atoms.size(); c1++)
	{
		c1Params = getC6AndR0(atoms[c1].atomicNumber);
		for(int c2 = 0; c2 < (int) atoms.size(); c2++)
		{
			if(c1 == c2) continue; // No self interaction		

			c2Params = getC6AndR0(atoms[c1].atomicNumber);
			C6 = sqrt(c1Params.first * c2Params.first);
			R0 = c1Params.second + c2Params.second;
			
			for(int t1 = -n[0]; t1<=n[0]; t1++)
			for(int t2 = -n[1]; t2<=n[1]; t2++)
			for(int t3 = -n[2]; t3<=n[2]; t3++)
			{
				deltaR = e->gInfo.R*(atoms[c2].pos + vector3<double>(t1, t2, t3) - atoms[c1].pos);
				Rij = deltaR.length();
				unit = deltaR/Rij;

				temp = calculateEnergyAndForceMagnitude(Rij, C6, R0);
				
				Energy -= temp.first;
				atoms[c1].force += temp.second*unit;
			}
		}
	}
	
	return 0.5*scalingFactor[EXCorr]*Energy;  // 0.5 is for double counting of interactions
}
