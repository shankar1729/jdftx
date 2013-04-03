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
#include <fluid1D/IdealGasMuEps.h>
#include <fluid1D/IdealGasPomega.h>
#include <fluid1D/Fex_H2O_FittedCorrelations.h>
#include <fluid1D/Fex_H2O_ScalarEOS.h>
#include <fluid1D/Fex_TM_ScalarEOS.h>
#include <fluid1D/Fex_H2O_BondedVoids.h>
#include <core1D/DataCollection.h>
#include <core1D/Operators.h>

int main(int argc, char** argv)
{	initSystem(argc, argv);

	//Setup simulation grid:
	GridInfo gInfo(GridInfo::Planar, 1, 0.125);
	
	FluidMixture fluidMixture(gInfo, 298*Kelvin);

	double Nguess = 5e-3; //default which works for water
	int Zn = 2; //default symmetry for water
	
	//----- Excess functional -----
	//Fex_H2O_FittedCorrelations fex(fluidMixture); string fexName = "FittedCorrelations";
	//Fex_H2O_ScalarEOS fex(fluidMixture); string fexName = "ScalarEOS";
	//Fex_H2O_BondedVoids fex(fluidMixture); string fexName = "BondedVoids";
	Fex_CHCl3_ScalarEOS fex(fluidMixture); string fexName = "CHCl3"; Nguess=1.1e-3; Zn = 3;
	//Fex_CCl4_ScalarEOS fex(fluidMixture); string fexName = "CCl4"; Nguess=0.9e-3; Zn = 3;
	
	//----- Setup quadrature for angular integration -----
	SO3quad quad(QuadEuler, Zn, 20, 1, 1); //use symmetries to construct extremely efficient specific quadrature
	
	//----- Translation operator -----
	TranslationOperatorLspline trans(gInfo);
	
	//----- Ideal gas -----
	IdealGasMuEps idgas(&fex, 1.0, quad, trans);
	//IdealGasPomega idgas(&fex, 1.0, quad, trans);

	double p = 1.01325*Bar;
	fluidMixture.setPressure(p, Nguess);

	//fluidMixture.verboseLog = true;
	
	MinimizeParams mp;
	mp.alphaTstart = 3e5;
	mp.nDim = gInfo.S * fluidMixture.get_nIndep();
	mp.nIterations=200;
	mp.nAlphaAdjustMax = 10;
	mp.dirUpdateScheme = MinimizeParams::FletcherReeves;
	
	FILE* fpEps = fopen((fexName + "/nonlineareps").c_str(), "w");
	double Dfield=1e-4;
	//double Dfield=6.4e-2;
	bool stateInitialized = false;
	for(; Dfield<6.5e-2; Dfield+=2e-3)
	{
		mp.knormThreshold = 1e-11;
		//mp.energyDiffThreshold = 1e-9 * gInfo.Volume() * pow(Dfield,2);
		fluidMixture.Eexternal = Dfield;

		if(!stateInitialized) //first iteration
		{	fluidMixture.initState(0.05); stateInitialized=true;
			mp.fdTest = true;
		}
		else mp.fdTest = false;
		
		fluidMixture.minimize(mp);
		fluidMixture.minimize(mp);

		ScalarFieldCollection N; double electricP;
		fluidMixture.getFreeEnergy(FluidMixture::Outputs(&N, &electricP));

		double nTyp = integral(N[0])/gInfo.Volume();
		double pTyp = electricP/gInfo.Volume();

		double epsilon = 1.0/(1.0 - 4*M_PI*pTyp/Dfield);
		double D_SI = Dfield/(eV/Angstrom); //Dfield in V/A
		printf("epsilon = %lf at D = %lf V/A\n", epsilon, D_SI);
		fprintf(fpEps, "%le\t%le\t%le\t%le\n", D_SI, epsilon,
			pTyp/(nTyp*fex.getMolecule()->get_dipole()), nTyp);
		fflush(fpEps);
		if(isnan(epsilon) || epsilon<0.0) break;
	}
	fclose(fpEps);
	return 0;
}
