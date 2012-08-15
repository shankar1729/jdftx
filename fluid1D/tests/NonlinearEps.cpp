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
#include <fluid1D/Fex_H2O_Lischner10.h>
#include <fluid1D/Fex_H2O_ScalarEOS.h>
#include <fluid1D/Fex_H2O_BondedVoids.h>
#include <core1D/DataCollection.h>
#include <core1D/Operators.h>

int main(int argc, char** argv)
{	initSystem(argc, argv);

	//Setup simulation grid:
	GridInfo gInfo(GridInfo::Planar, 1, 0.125);
	
	//----- Setup quadrature for angular integration -----
	const int Zn = 2; //Water molecule has Z2 symmetry about dipole axis
	SO3quad quad(QuadEuler, Zn, 20, 1, 1); //use symmetries to construct extremely efficient specific quadrature
	
	//----- Translation operator -----
	TranslationOperatorLspline trans(gInfo);
	
	FluidMixture fluidMixture(gInfo, 298*Kelvin);

	//----- Excess functional -----
	//Fex_H2O_Lischner10 fex(fluidMixture);
	Fex_H2O_ScalarEOS fex(fluidMixture);
	//Fex_H2O_BondedVoids fex(fluidMixture);

	//----- Ideal gas -----
	IdealGasMuEps idgas(&fex, 1.0, quad, trans);

	double p = 1.01325*Bar;
	fluidMixture.setPressure(p);

	MinimizeParams mp;
	mp.alphaTstart = 3e5;
	mp.nDim = gInfo.S * fluidMixture.get_nIndep();
	mp.nIterations=200;

	FILE* fpEps = fopen("nonlineareps", "w");
	double Dfield=1e-4;
	bool stateInitialized = false;
	for(; Dfield<5e-2; Dfield+=2e-3)
	{
		mp.energyDiffThreshold = 1e-9 * gInfo.Volume() * pow(Dfield,2);
		idgas.Eexternal = Dfield;

		if(!stateInitialized) //first iteration
		{	fluidMixture.initState(0.05); stateInitialized=true;
			mp.fdTest = true;
		}
		else mp.fdTest = false;
		
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
