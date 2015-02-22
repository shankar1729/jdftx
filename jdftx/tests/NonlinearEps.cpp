/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

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

#include <fluid/FluidMixture.h>
#include <fluid/IdealGasMuEps.h>
#include <fluid/Fex_H2O_FittedCorrelations.h>
#include <fluid/Fex_ScalarEOS.h>
#include <fluid/Fex_H2O_BondedVoids.h>
#include <core/ScalarFieldArray.h>
#include <core/Operators.h>

extern EnumStringMap<FluidComponent::Name> solventMap;

int main(int argc, char** argv)
{	initSystem(argc, argv);

	//Setup simulation grid:
	GridInfo gInfo;
	const double hGrid = 0.125;
	gInfo.S = vector3<int>(1, 1, 1);
	gInfo.R = Diag(hGrid * gInfo.S);
	gInfo.initialize();

	double T = 298*Kelvin;
	FluidComponent component(FluidComponent::CHCl3, T, FluidComponent::ScalarEOS);
	string solventName = solventMap.getString(component.name);
	component.s2quadType = QuadEuler;
	component.quad_nBeta = 11;
	component.quad_nAlpha = 2;
	component.quad_nGamma = 1;
	
	FluidMixture fluidMixture(gInfo, T);
	component.addToFluidMixture(&fluidMixture);
	double p = 1.01325*Bar;
	logPrintf("pV = %le\n", p*gInfo.detR);
	fluidMixture.initialize(p);

	MinimizeParams mp;
	mp.fpLog = globalLog;
	mp.nDim = gInfo.nr * fluidMixture.get_nIndep();
	mp.nAlphaAdjustMax = 10;
	mp.nIterations=200;

	FILE* fpEps = fopen((solventName+"/nonlineareps").c_str(), "w");
	double Dfield=1e-4;
	bool stateInitialized = false;
	for(; Dfield<6.5e-2; Dfield+=2e-3)
	{
		mp.energyDiffThreshold = 1e-9 * gInfo.detR * pow(Dfield,2);
		mp.knormThreshold = 1e-9 * gInfo.detR * Dfield;
		fluidMixture.Eexternal = vector3<>(0, 0, Dfield);

		if(!stateInitialized) //first iteration
		{	fluidMixture.initState(0.05); stateInitialized=true;
			mp.fdTest = true;
		}
		else mp.fdTest = false;
		
		fluidMixture.minimize(mp);

		ScalarFieldArray N; vector3<> electricP;
		fluidMixture.verboseLog = true;
		fluidMixture.getFreeEnergy(FluidMixture::Outputs(&N, &electricP));
		fluidMixture.verboseLog = false;
		
		double nTyp = integral(N[0])/gInfo.detR;
		double pTyp = electricP[2]/gInfo.detR;

		double epsilon = 1.0/(1.0 - 4*M_PI*pTyp/Dfield);
		double D_SI = Dfield/(eV/Angstrom); //Dfield in V/A
		logPrintf("epsilon = %lf at D = %lf V/A\n", epsilon, D_SI);
		fprintf(fpEps, "%le\t%le\t%le\t%le\n", D_SI, epsilon,
			pTyp/(nTyp*component.molecule.getDipole().length()), nTyp);
		fflush(fpEps);
		if(std::isnan(epsilon) || epsilon<0.0) break;
	}
	fclose(fpEps);
	
	finalizeSystem();
	return 0;
}
