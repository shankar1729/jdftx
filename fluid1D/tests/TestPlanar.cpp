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
#include <fluid1D/Fex_H2O_Lischner10.h>
#include <fluid1D/Fex_H2O_ScalarEOS.h>
#include <fluid1D/Fex_H2O_BondedVoids.h>

#define TDEP
#define ChosenFex BondedVoids


#define FEX_NAME_quoter(FexSuffix) #FexSuffix
#define FEX_NAME(FexSuffix) FEX_NAME_quoter(FexSuffix)
string fexName = FEX_NAME(ChosenFex);

#define FEX_CLASS_paster(FexSuffix) Fex_H2O_##FexSuffix
#define FEX_CLASS(FexSuffix) FEX_CLASS_paster(FexSuffix)
#define FexClass FEX_CLASS(ChosenFex)


//Return planar liquid-vapor surface energy (and optionally plot density profiles)
double testPlanar(double T, double sigmaTarget=0., bool plotDensities=false)
{
	//Setup simulation grid:
	GridInfo gInfo(GridInfo::Planar, 768, 0.125);
	
	//----- Setup quadrature for angular integration -----
	const int Zn = 2; //Water molecule has Z2 symmetry about dipole axis
	SO3quad quad(QuadEuler, Zn, 20, 1); //Force nAlpha = 1

	//----- Translation operator -----
	TranslationOperatorLspline trans(gInfo);
	
	FluidMixture fluidMixture(gInfo, T);

	//----- Excess functional -----
	FexClass fex(fluidMixture);

	//----- Ideal gas -----
	IdealGasPsiAlpha idgas(&fex, 1.0, quad, trans);

	std::vector<double> Nvap;
	fluidMixture.setBoilingPressure(&Nvap);
	
	//----- Initialize geometry --------
	fluidMixture.state.clear();
	nullToZero(fluidMixture.state, gInfo, fluidMixture.get_nIndep());
	double rWall = 0.25 * gInfo.rMax;
	double psiVap = log(Nvap[0]/idgas.get_Nbulk());
	double* psiOdata = fluidMixture.state[0].data();
	for(int i=0; i<gInfo.S; i++) //no potential, only initial state difference:
		psiOdata[i] = gInfo.r[i]<rWall ? psiVap : 0.;
	
	//----- FDtest and CG parameters -----
	MinimizeParams mp;
	mp.alphaTstart = 3e4;
	mp.nDim = gInfo.S * fluidMixture.get_nIndep();
	mp.energyLabel = "sigma";
	mp.nIterations = 100;
	mp.energyDiffThreshold=1e-11;
	
	fluidMixture.minimize(mp);
		
	ScalarFieldCollection N;
	double sigma = fluidMixture.getFreeEnergy(FluidMixture::Outputs(&N));
	
	if(sigmaTarget) printf("Error in sigma from target value = %le\n", sigma/sigmaTarget-1.);
	
	if(plotDensities)
	{	string fname = fexName + "/planar-psi-n";
		N *= 1. / idgas.get_Nbulk();
		printToFile(N, fname.c_str());
		
		FILE* pp = popen("gnuplot -persist", "w");
		fprintf(pp, "set xlabel \"x [bohr]\"\n");
		fprintf(pp, "set ylabel \"n [bohr^-3]\"\n");
		fprintf(pp, "plot \"%s\" using 1:2 title \"NO / Nbulk\" w l, \"\" using 1:(0.5*$3) title \"NH / 2 Nbulk\" w l\n", fname.c_str());
		fprintf(pp, "set key bottom right\n");
		fclose(pp);
	}
	return sigma;
}

int main(int argc, char** argv)
{	initSystem(argc, argv);

	#ifdef TDEP
		FILE* fp = fopen((fexName+"/sigmavsT").c_str(), "w");
		for(double Tcel=0.0; Tcel<105.0; Tcel+=10.0)
		{	double T = (Tcel+273.16)*Kelvin;
			double sigmaSI = testPlanar(T) / (1e-3*Newton/meter);
			logPrintf("\n------------ T = %lf C, sigma = %le mN/m ---------------\n\n", Tcel, sigmaSI);
			fprintf(fp, "%lf %le\n", Tcel, sigmaSI);
		}
		fclose(fp);
	#else
		testPlanar(298*Kelvin, 71.98e-3 * Newton/meter, true);
	#endif
	return 0.;
}
