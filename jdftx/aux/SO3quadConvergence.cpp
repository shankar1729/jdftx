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
#include <fluid/IdealGasPomega.h>
#include <fluid/IdealGasPsiAlpha.h>
#include <fluid/Fex_ScalarEOS.h>
#include <core/ScalarFieldIO.h>

void initHardSphere(int i, vector3<> r, const vector3<>& r0, double radius, double height, double* phi)
{	phi[i] = ((r - r0).length() < radius ? height : 0.0);
}

int main(int argc, char** argv)
{	initSystem(argc, argv);
	
	//Parse command-line
	if(argc < 2) die("Usage: %s <quad> [<nBeta>] [INIT]\n", argv[0])
	S2quadType quadType = QuadOctahedron;
	if(!S2quadTypeMap.getEnum(argv[1], quadType))
		die("<quad> must be one of %s\n", S2quadTypeMap.optionList().c_str());
	int nBeta=0;
	if(quadType==QuadEuler)
	{	if(argc < 3) die("<nBeta> must be specified for Euler quadratures.\n")
		nBeta = atoi(argv[2]);
		if(nBeta <= 0) die("<nBeta> must be non-negative.")
	}
	int iArgInit = ((quadType==QuadEuler) ? 3 : 2);
	bool shouldInit = (iArgInit<argc && !strcasecmp(argv[iArgInit],"INIT"));
	
	//Setup grid:
	GridInfo gInfo;
	gInfo.S = vector3<int>(128, 128, 128);
	gInfo.R = Diag(0.25 * gInfo.S); //32 bohr^3 box
	gInfo.initialize();
	
	double T = 298*Kelvin;
	FluidComponent component(FluidComponent::H2O, T, FluidComponent::ScalarEOS);
	component.s2quadType = quadType;
	component.quad_nBeta = nBeta;
	component.representation = shouldInit ? FluidComponent::Pomega : FluidComponent::PsiAlpha;
	
	FluidMixture fluidMixture(gInfo, T);
	component.addToFluidMixture(&fluidMixture);
	double p = 1.01325*Bar;
	logPrintf("pV = %le\n", p*gInfo.detR);
	fluidMixture.initialize(p);

	//Initialize external potential (repel O from a 4A sphere)
	double radius = 4*Angstrom;
	nullToZero(component.idealGas->V, gInfo);
	applyFunc_r(gInfo, initHardSphere, gInfo.R*vector3<>(0.5,0.5,0.5), radius, 1.0, component.idealGas->V[0]->data());
	
	//----- Initialize state -----
	if(shouldInit) fluidMixture.initState(0.15);
	else fluidMixture.loadState("init.psiEff");

	//----- FDtest and CG -----
	MinimizeParams mp;
	mp.fpLog = globalLog;
	mp.nDim = gInfo.nr * fluidMixture.get_nIndep();
	mp.energyLabel = "Phi";
	mp.nIterations=500;
	mp.energyDiffThreshold=1e-16;
	//mp.fdTest = true;

	//------ Outputs ---------
	fluidMixture.minimize(mp);
	ScalarFieldArray N, psiEff;
	fluidMixture.getFreeEnergy(FluidMixture::Outputs(&N, 0, 0, shouldInit?&psiEff:0));
	if(shouldInit)
	{	saveDX(N[0], "init.NO");
		saveDX(N[1], "init.NH");
		saveToFile(psiEff, "init.psiEff");
	}
	
	finalizeSystem();
	return 0;
}


