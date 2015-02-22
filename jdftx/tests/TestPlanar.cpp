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
#include <fluid/IdealGasPsiAlpha.h>
#include <fluid/IdealGasMuEps.h>
#include <fluid/IdealGasPomega.h>
#include <fluid/Fex_ScalarEOS.h>

void setPhi(int i, vector3<> r, double* phiApplied, double* phiWall,
	double gridLength, double Dfield, double zWall)
{
	register double zEff;
	if(r[2]<0.5*gridLength) zEff = r[2] - 0.25*gridLength;
	else zEff = 0.75*gridLength - r[2];
	phiApplied[i] = -Dfield * zEff;
	phiWall[i] = ((fabs(zEff)>0.25*gridLength-zWall) ? 1.0 : 0.0);
}

int main(int argc, char** argv)
{	initSystem(argc, argv);

	//Setup simulation grid:
	GridInfo gInfo;
	gInfo.S = vector3<int>(1, 1, 1536);
	const double hGrid = 0.125;
	gInfo.R = Diag(vector3<>(0.5, 1., hGrid * gInfo.S[2])); //transverse dimensions so as to get energy per interface area
	gInfo.initialize();

	double T = 298*Kelvin;
	FluidComponent component(FluidComponent::CCl4, T, FluidComponent::ScalarEOS);
	component.s2quadType = QuadOctahedron;
	component.representation = FluidComponent::PsiAlpha;
	FluidMixture fluidMixture(gInfo, T);
	component.addToFluidMixture(&fluidMixture);
	
	std::vector<double> Nvap;
	double p = fluidMixture.getBoilingPressure(component.Nbulk, component.Pvap/T, &Nvap);
	fluidMixture.initialize(p);

	//----- Initialize geometry --------
	fluidMixture.state.clear();
	nullToZero(fluidMixture.state, gInfo, fluidMixture.get_nIndep());
	double L = gInfo.R(2,2);
	double zWall = 0.125 * L;
	double Rhs = pow(component.molecule.getVhs()*3./(4*M_PI), 1./3);
	double Nbulk = component.idealGas->get_Nbulk();
	double* stateData = fluidMixture.state[0]->data();
	for(int i=0; i<gInfo.S[2]; i++) //no potential, only initial state difference:
	{	double z = hGrid*i;
		z -= L * floor(z/L + 0.5);
		stateData[i] = log(Nvap[0]/Nbulk + 0.5*erfc((zWall-fabs(z))/Rhs)*(1.-Nvap[0]/Nbulk));
	}
	
	//----- FDtest and CG parameters -----
	MinimizeParams mp;
	mp.fpLog = globalLog;
	mp.nDim = gInfo.nr * fluidMixture.get_nIndep();
	mp.energyLabel = "sigma";
	mp.nIterations = 200;
	mp.energyDiffThreshold=1e-11;
	
	fluidMixture.minimize(mp);
	
	ScalarFieldArray N;
	double sigma = fluidMixture.getFreeEnergy(FluidMixture::Outputs(&N));
	
	double sigmaTarget = component.sigmaBulk;
	if(sigmaTarget) logPrintf("Error in sigma from target value = %le\n", sigma/sigmaTarget-1.);
	
	bool plotDensities = true;
	if(plotDensities)
	{	string fname = component.molecule.name + "/planar-psi-n";
		
		std::vector<const double*> Ndata(N.size());
		for(unsigned k=0; k<N.size(); k++) Ndata[k] = N[k]->data();
		FILE* fp = fopen(fname.c_str(), "w");
		fprintf(fp, "z");
		for(unsigned k=0; k<N.size(); k++)
			fprintf(fp, "\t%s", component.molecule.sites[k]->name.c_str());
		fprintf(fp, "\n");
		for(int i=0; i<gInfo.S[2]; i++)
		{	fprintf(fp, "%lg", i*hGrid);
			for(unsigned k=0; k<N.size(); k++)
				fprintf(fp, "\t%lg", Ndata[k][i]/(Nbulk*component.molecule.sites[k]->positions.size()));
			fprintf(fp, "\n");
		}
		fclose(fp);
		
		FILE* pp = popen("gnuplot -persist", "w");
		fprintf(pp, "set xlabel \"z [bohr]\"\n");
		fprintf(pp, "set ylabel \"N/Nbulk\"\n");
		fprintf(pp, "set key bottom right autotitle columnhead\n");
		fprintf(pp, "plot for [i=2:%lu] \"%s\" using 1:i w l\n", N.size()+1, fname.c_str());
		fclose(pp);
	}

	finalizeSystem();
	return 0;
}
