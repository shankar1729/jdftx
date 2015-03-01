/*-------------------------------------------------------------------
Copyright 2011 Kendra Letchworth Weaver

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
#include <fluid/IdealGasMonoatomic.h>
#include <fluid/Fex_ScalarEOS.h>
#include <core/Units.h>
#include <core/ScalarFieldIO.h>

void initQuadraticWell(int i, vector3<> r, const vector3<>& r0, double width, double depth, double dr, double beta, double* phi)
{
	double sigma = width*sqrt(beta);
	double maxPsi = 0.000;
	double widthTrans = (width > 3.0*dr ? width/2.0 : 2.0*dr );
	double rTrans = sigma*sqrt(2.0*(maxPsi+depth))-2.0*width;

	double rmr0 = (r - r0).length(); //distance from center of well

	phi[i] = ((0.5*pow(rmr0/sigma,2)-depth)*(1.0+erf((rTrans-rmr0)/widthTrans))/2.0 + maxPsi*(1.0+erf((rmr0-rTrans)/widthTrans))/2.0);
	//phi[i] = (0.5*pow(rmr0/sigma,2)-depth);
}


struct TestGaussian
{	const GridInfo& gInfo;

	TestGaussian(const GridInfo& gInfo):gInfo(gInfo) {}

	void test()
	{
		double T = 298*Kelvin;
		FluidComponent component(FluidComponent::H2O, T, FluidComponent::ScalarEOS);
		component.s2quadType = Quad7design_24;
		
		FluidMixture fluidMixture(gInfo, T);
		component.addToFluidMixture(&fluidMixture);
		double p = 1.01325*Bar;
		logPrintf("pV = %le\n", p*gInfo.detR);
		fluidMixture.initialize(p);

		//Single configuration test:
		//#define SINGLE_CONFIG
		#ifdef SINGLE_CONFIG
		//Set external potential to zero
		idgas.V.clear();
		nullToZero(idgas.V, gInfo, fluidMixture.get_nDensities());
		//Set state to a kronecker delta in real space:
		fluidMixture.state.clear();
		nullToZero(fluidMixture.state, gInfo, fluidMixture.get_nDensities());
		fluidMixture.state[0] -= 100; //denisties ~ 1e-30
		fluidMixture.state[0]->data()[0] = 9;
		fluidMixture.verboseLog = true;
		//Evaluate:
		ScalarFieldArray N;
		double Phi = fluidMixture.getFreeEnergy(FluidMixture::Outputs(&N));
		logPrintf("Phi = %25.16lf, integral(N)=%25.16lf\n", Phi, integral(N[0]));
		saveRawBinary(N[0], "TestGaussian.N");
		return;
		#endif
		
		double widthStart = 0.2;
		int numWidths = 9;
		double deltaWidth = 0.1;

		double depthStart = 0.02;
		int numDepths = 1;
		double deltaDepth = 0.01;

		FILE* fplog = fopen("TestGaussian-Omegas-FixedTZ", "w");

		for(int i=0; i<numWidths; i++)
		{
			double width = widthStart + i*deltaWidth;
			for(int j=0; j<numDepths; j++)
			{
				double depth = depthStart + j*deltaDepth;

				//Initialize the OH potential:
				//formula assumes orthogonal lattice vectors for now dr^2=dx^2+dy^2+dz^2
				double dr = sqrt
					( (gInfo.R.column(0)/gInfo.S[0]).length_squared()
					+ (gInfo.R.column(1)/gInfo.S[1]).length_squared()
					+ (gInfo.R.column(2)/gInfo.S[2]).length_squared() );
				
				nullToZero(component.idealGas->V, gInfo, fluidMixture.get_nDensities());
				applyFunc_r(gInfo, initQuadraticWell, gInfo.R*vector3<>(.5,.5,.5),
					width,depth, dr, 1./T, component.idealGas->V[0]->data());

				//const char* stateFilename ="TestGaussian-EOSScalar/state.%c.bin";
				//Initial wavefunction
				fluidMixture.initState(0.0);
				//fluidMixture.loadState(stateFilename);

				MinimizeParams mp;
				mp.fpLog = globalLog;
				mp.dirUpdateScheme = MinimizeParams::HestenesStiefel;
				mp.nDim = gInfo.nr * fluidMixture.get_nIndep();
				mp.energyLabel = "Phi";
				mp.nIterations=200;
				mp.knormThreshold=1e-11;
				mp.fdTest = true;
				
				logPrintf("Starting CG:\n");
				TIME("minimize", globalLog, fluidMixture.minimize(mp); );

				ScalarFieldArray N;
				double Omega;
				TIME("getOmega calculation (with gradient)", globalLog,
					Omega = fluidMixture.getFreeEnergy(FluidMixture::Outputs(&N));
				);

				fprintf(fplog,"Width: %10.6le depth: %10.6le Omega: %10.6le\n",width, depth, Omega);

//				ScalarField saveR[2] = {n.O(), n.H()};
//				char filename[256];
//				sprintf(filename, "%s/testGaussian%s_nO", runFolder, runNameSuffix); saveDX(n.O(), filename);
//				sprintf(filename, "%s/testGaussian%s_nH", runFolder, runNameSuffix); saveDX(n.H(), filename);
//				sprintf(filename, "%s/testGaussian%s_n.spherical", runFolder, runNameSuffix); saveSphericalized(saveR, 2, filename, 0.25);
//				ScalarField saveR2[2] = {fex.phiOH.O(), fex.phiOH.H()};
//				sprintf(filename, "%s/testGaussian%s_phi.spherical", runFolder, runNameSuffix); saveSphericalized(saveR2, 2, filename, 0.25);
//				water.saveState(stateFilename);
			}
		}
		fclose(fplog);
	}
};


int main(int argc, char** argv)
{	initSystem(argc, argv);

	//Setup simulation grid:
	GridInfo gInfo;
	//gInfo.S = vector3<int>(256, 256, 256); double hGrid=0.0625;
	//gInfo.S = vector3<int>(128, 128, 128); double hGrid=0.125;
	gInfo.S = vector3<int>(128, 128, 128); double hGrid=0.25;
	gInfo.R = Diag(hGrid * gInfo.S);
	gInfo.initialize();

	TestGaussian(gInfo).test();

	finalizeSystem();
	return 0;
}
