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

#include <core/ScalarFieldIO.h>
#include <core/Minimize.h>
#include <fluid/FluidMixture.h>
#include <fluid/IdealGasPsiAlpha.h>
#include <fluid/IdealGasMuEps.h>
#include <fluid/IdealGasPomega.h>
#include <fluid/Fex_H2O_FittedCorrelations.h>
#include <fluid/Fex_ScalarEOS.h>
#include <fluid/Fex_H2O_BondedVoids.h>

void initHardSphere(int i, vector3<> r, const vector3<>& r0, double radius, double height, double* phi)
{	phi[i] = ((r - r0).length() < radius ? height : 0.0);
}

struct TestSphere
{	const GridInfo& gInfo;

	TestSphere(const GridInfo& gInfo):gInfo(gInfo) {}

	void test()
	{
		double T = 298*Kelvin;
		FluidComponent component(FluidComponent::H2O, T, FluidComponent::ScalarEOS);
		component.s2quadType = QuadOctahedron;
		component.representation = FluidComponent::MuEps;
		
		FluidMixture fluidMixture(gInfo, T);
		component.addToFluidMixture(&fluidMixture);
		double p = 1.01325*Bar;
		logPrintf("pV = %le\n", p*gInfo.detR);
		fluidMixture.initialize(p);

		//----- Initialize external potential -----
		double radius = 4.0*Angstrom;
		nullToZero(component.idealGas->V, gInfo);
		applyFunc_r(gInfo, initHardSphere, gInfo.R*vector3<>(0.5,0.5,0.5), radius, 1.0, component.idealGas->V[0]->data());

		//----- Initialize state -----
		bool loadState = false; //true;
		if(loadState) fluidMixture.loadState("TestSphere.psiEff");
		else fluidMixture.initState(0.15);

		//----- FDtest and CG -----
		MinimizeParams mp;
		mp.fpLog = globalLog;
		mp.nDim = gInfo.nr * fluidMixture.get_nIndep();
		mp.energyLabel = "Phi";
		mp.nIterations=100;
		mp.energyDiffThreshold=1e-11;
		//mp.fdTest = !loadState;

		logPrintf("Starting CG:\n");
		TIME("minimize", globalLog,
			fluidMixture.minimize(mp);
		);

		ScalarFieldArray N, psiEff;
		TIME("getFreeEnergy", globalLog,
			fluidMixture.getFreeEnergy(FluidMixture::Outputs(&N, 0, 0, &psiEff));
		);

		saveDX(N[0], "TestSphere.NO");
		saveDX(N[1], "TestSphere.NH");
		saveSphericalized(&N[0], N.size(), "TestSphere.N.spherical", 0.25);
		saveToFile(psiEff, "TestSphere.psiEff");
		
		//saveSphericalized(&N[0], N.size(), "Octahedron.Nspherical", 0.25);
	}
};


int main(int argc, char** argv)
{	initSystem(argc, argv);

	//Setup simulation grid:
	GridInfo gInfo;
	//gInfo.S = vector3<int>(32, 32, 32); double hGrid=1.0;
	gInfo.S = vector3<int>(64, 64, 64); double hGrid=0.5;
	//gInfo.S = vector3<int>(128, 128, 128); double hGrid=0.25;
	//gInfo.S = vector3<int>(256, 256, 256); double hGrid=0.125;
	//gInfo.S = vector3<int>(192, 192, 192); double hGrid=0.25;
	gInfo.R = Diag(hGrid * gInfo.S);
	gInfo.initialize();
	TestSphere(gInfo).test();
	
	finalizeSystem();
	return 0;
}
