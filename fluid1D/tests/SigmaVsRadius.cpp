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

#include <core/Minimize.h>
#include <fluid1D/FluidMixture.h>
#include <fluid1D/IdealGasPsiAlpha.h>
#include <fluid1D/IdealGasMuEps.h>
#include <fluid1D/IdealGasPomega.h>
#include <fluid1D/Fex_H2O_FittedCorrelations.h>
#include <fluid1D/Fex_H2O_ScalarEOS.h>
#include <fluid1D/Fex_H2O_BondedVoids.h>

int main(int argc, char** argv)
{	initSystem(argc, argv);

	//Setup simulation grid:
	GridInfo gInfo(GridInfo::Spherical, 512, 0.125);
	
	//----- Setup quadrature for angular integration -----
	const int Zn = 2; //Water molecule has Z2 symmetry about dipole axis
	SO3quad quad(QuadEuler, Zn, 20, 1); //Force nAlpha = 1

	//----- Translation operator -----
	TranslationOperatorLspline trans(gInfo);
	
	FluidMixture fluidMixture(gInfo, 298*Kelvin);

	//----- Excess functional -----
	Fex_H2O_FittedCorrelations fex(fluidMixture); string fexName = "FittedCorrelations";
	//Fex_H2O_ScalarEOS fex(fluidMixture); string fexName = "ScalarEOS";
	//Fex_H2O_BondedVoids fex(fluidMixture); string fexName = "BondedVoids";

	//----- Ideal gas -----
	IdealGasPomega idgas(&fex, 1.0, quad, trans);

	double p = 1.01325*Bar;
	fluidMixture.setPressure(p);
	
	//----- FDtest and CG parameters -----
	MinimizeParams mp;
	mp.alphaTstart = 3e4;
	mp.nDim = gInfo.S * fluidMixture.get_nIndep();
	mp.energyLabel = "Phi";
	mp.nIterations=100;
	mp.energyDiffThreshold=1e-11;
	
	FILE* fp = fopen((fexName + "/sigmavsradius").c_str(), "w");
	double Phi0 = 0.;
	for(int iRadius=0; iRadius<24; iRadius++)
	{
		//----- Initialize external potential -----
		double radius = iRadius ? 0.75*iRadius + 0.5*gInfo.rMax/gInfo.S : 0.;
		logPrintf("\n\nStarting solve #%d for radius %lf bohr:\n", iRadius+1, radius);
		nullToZero(idgas.V, gInfo);
		double* Vdata = idgas.V[0].data();
		for(int i=0; i<gInfo.S; i++)
			Vdata[i] = gInfo.r[i]<radius ? 1. : 0.;

		fluidMixture.initState(0.15);
		double Phi = fluidMixture.minimize(mp);
		if(!iRadius) { Phi0 = Phi; fprintf(fp, "%lf\t%le\n", 0., 0.); }
		else fprintf(fp, "%lf\t%le\n", radius, ((Phi-Phi0) - p * (4*M_PI*pow(radius,3)/3)) / (4*M_PI*pow(radius,2)));
	}
	fclose(fp);
}
