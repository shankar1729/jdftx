/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

This file is part of Fluid1D.

Fluid1D is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Fluid1D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Fluid1D.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#include <fluid1D/FluidMixture.h>
#include <fluid1D/IdealGasPomega.h>
#include <fluid1D/Fex_H2O_FittedCorrelations.h>
#include <fluid1D/Fex_H2O_ScalarEOS.h>
#include <fluid1D/Fex_H2O_BondedVoids.h>

#define ChosenFex ScalarEOS

#define FEX_NAME_quoter(FexSuffix) #FexSuffix
#define FEX_NAME(FexSuffix) FEX_NAME_quoter(FexSuffix)
string fexName = FEX_NAME(ChosenFex);

#define FEX_CLASS_paster(FexSuffix) Fex_H2O_##FexSuffix
#define FEX_CLASS(FexSuffix) FEX_CLASS_paster(FexSuffix)
#define FexClass FEX_CLASS(ChosenFex)


//Compute and return pair correlations for a specific temperature and bulk density
ScalarFieldCollection pairCorrelations(const GridInfo& gInfo, double T, double Nbulk)
{
	if(fexName=="FittedCorrelations" && fabs(T/Kelvin-298)>1)
	{	//provide NaN correlations for FittedCorrelations at invalid temperatures: (eases plot scripts)
		ScalarFieldCollection g; nullToZero(g, gInfo, 3);
		g *= NAN;
		return g;
	}
	
	logPrintf("\n\n------------- T=%.0lfK  Nbulk=%.2lebohr^-3 -------------\n", T/Kelvin, Nbulk);
	SO3quad quad(QuadTetrahedron, 2); //Not used and need it only to create an IdealGas object. Just picked the smallest quadrature
	TranslationOperatorLspline trans(gInfo); //Not used and need it only to create an IdealGas object.
	
	FluidMixture fluidMixture(gInfo, T);
	FexClass fex(fluidMixture); //Selected excess functionals using macros above
	IdealGasPomega idgas(&fex, 1.0, quad, trans); //Not used and needed only to register fex with fluidMixture
	
	//Compute and return pair correlations:
	std::vector<double> Nmol(1, Nbulk);
	return fluidMixture.getPairCorrelations(Nmol);
}

int main()
{
	//Read thermodynamic state points from a file (T in Kelvin, Nbulk in bohr^-3)
	std::vector< std::pair<double,double> > conditions;
	FILE* fp = fopen("epsr_conditions", "r");
	while(!feof(fp))
	{	double T, Nbulk;
		if(fscanf(fp, "%lf %lf", &T, &Nbulk) != 2) break;
		conditions.push_back(std::make_pair(T*Kelvin, Nbulk));
	}
	fclose(fp);
	
	GridInfo gInfo(GridInfo::Spherical, 1024, 0.0625); //Fine spherical grid to get good resolution in gXX
	ScalarFieldCollection gOO, gOH, gHH; //partial pair correlations for all conditions
	for(auto condition: conditions)
	{	double T = condition.first;
		double Nbulk = condition.second;
		ScalarFieldCollection gXX = pairCorrelations(gInfo, T, Nbulk);
		gOO.push_back(gXX[0]);
		gOH.push_back(gXX[1]);
		gHH.push_back(gXX[2]);
	}
	printToFile(gOO, (fexName+"/gOO").c_str());
	printToFile(gOH, (fexName+"/gOH").c_str());
	printToFile(gHH, (fexName+"/gHH").c_str());
}
