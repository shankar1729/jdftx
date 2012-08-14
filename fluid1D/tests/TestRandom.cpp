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

#include <fluid1D/FluidMixture.h>
#include <fluid1D/IdealGasMonoatomic.h>
#include <fluid1D/Fex_LJ.h>
#include <fluid1D/TranslationOperator.h>
#include <core/Units.h>
#include <cstdio>


int main(int argc, char** argv)
{	initSystem(argc, argv);

	//double T=298*Kelvin;
	//double T=130*Kelvin;
	double T=75*Kelvin;
	const double rMax = 64./pow(4*M_PI/3,1./3);
	const int S = 256;
	GridInfo gInfo(GridInfo::Spherical, S, rMax/S);
	
	FluidMixture fluidMixture(gInfo, T);
	//Argon component:
	Fex_LJ fexAr(fluidMixture, 119.8*Kelvin, 3.405*Angstrom, "Ar", +1*0); //params from SLOGwiki
	IdealGasMonoatomic idAr(&fexAr, 2.0);
	//Neon component:
	Fex_LJ fexNe(fluidMixture, 3.2135e-3*eV, 2.782*Angstrom, "Ne", -2*0); //params from SLOGwiki
	IdealGasMonoatomic idNe(&fexNe, 1.0);
	//Interaction between them:
	Fmix_LJ fmixLJ(fexAr, fexNe);
	//based on mole fractions, this should create a 2:1 Ar:Ne mixture

	fluidMixture.setPressure(1000*Bar);

	nullToZero(idAr.V[0], gInfo);
	double* Vdata = idAr.V[0].data();
	for(int i=0; i<gInfo.S; i++)
		Vdata[i] = gInfo.r[i]<6.0 ? 0.005 : 0.;
	idNe.V[0] = idAr.V[0];

	idNe.set_Nnorm(52);

	//fluidMixture.verboseLog = true;
	fluidMixture.initState(0.15);
	//fluidMixture.loadState("random.state");

	MinimizeParams mp;
	mp.alphaTstart = 3e1;
	mp.nDim = gInfo.S*fluidMixture.get_nIndep();
	mp.nIterations=100;
	mp.knormThreshold=5e-12;
	mp.fdTest = true;

	fluidMixture.minimize(mp);
	//fluidMixture.saveState("random.state");

	ScalarFieldCollection N;
	fluidMixture.getFreeEnergy(FluidMixture::Outputs(&N));

	printf("num(Ar) = %le\n", integral(N[0]));
	printf("num(Ne) = %le\n", integral(N[1]));

	N[0] *= 1.0/idAr.get_Nbulk();
	N[1] *= 1.0/idNe.get_Nbulk();
	printToFile(N, "random-ArNe");
}
