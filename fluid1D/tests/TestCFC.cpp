/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman

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
#include <fluid1D/IdealGasMonoatomic.h>
#include <fluid1D/TranslationOperator.h>
#include <fluid1D/Fex_TM_ScalarEOS.h>
#include <fluid1D/IdealGasPsiAlpha.h>
#include <fluid1D/IdealGasPomega.h>
#include <fluid1D/IdealGasMonoatomic.h>
#include <core/Units.h>
#include <cstdio>
#include <gsl/gsl_sf.h>
#include <list>

//#define CAVITATION_CDFT
//#define CAVITATION_MODEL
//#define PLANAR_TDEP
//#define PLANAR
#define CORRFUNC

//#define ChosenLiq CHCl3
#define ChosenLiq CCl4

const double sigmaCHCl3 = 26.67 * (1e-3*Newton/meter);
const double sigmaCCl4 = 26.17 * (1e-3*Newton/meter);

const double NbulkCHCl3 = 1.109e-3;
const double NbulkCCl4 = 0.9205e-3;

const double PvapCHCl3 = 26.0*KPascal;
const double PvapCCl4 = 15.128*KPascal;

const double RhsCHCl3 = 2.06*Angstrom;
const double RhsCCl4 = 2.17*Angstrom;


#define FEX_NAME_quoter(LiqName) #LiqName
#define FEX_NAME(LiqName) FEX_NAME_quoter(LiqName)
string fexName = FEX_NAME(ChosenLiq);

#define FEX_CLASS_paster(LiqName) Fex_##LiqName##_ScalarEOS
#define FEX_CLASS(LiqName) FEX_CLASS_paster(LiqName)
#define FexClass FEX_CLASS(ChosenLiq)

#define PASTE_TOKENS_helper(a, b) a##b
#define PASTE_TOKENS(a, b) PASTE_TOKENS_helper(a, b)


//Return planar liquid-vapor surface energy (and optionally plot density profiles)
double testPlanar(double T, double sigmaTarget=0., bool plotDensities=false)
{
	//Setup simulation grid:
	GridInfo gInfo(GridInfo::Planar, 768, 0.125);
	
	//----- Setup quadrature for angular integration -----
	const int Zn = 3; //Chloroform has Z3 symmetry about dipole axis
	SO3quad quad(QuadEuler, Zn, 18, 1); //Force nAlpha = 1

	//----- Translation operator -----
	TranslationOperatorLspline trans(gInfo);
	
	FluidMixture fluidMixture(gInfo, T);

	//----- Excess functional -----
	FexClass fex(fluidMixture);
	
	//----- Ideal gas -----
	IdealGasPsiAlpha idgas(&fex, 1.0, quad, trans);
	//IdealGasMonoatomic idgas(&fex, 1.0);
	
	std::vector<double> Nvap;
	fluidMixture.setBoilingPressure(&Nvap, 1e-8, PASTE_TOKENS(Nbulk,ChosenLiq), 1e-6);
	
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
	mp.fdTest = true;
	
	fluidMixture.minimize(mp);
	
	ScalarFieldCollection N;
	double sigma = fluidMixture.getFreeEnergy(FluidMixture::Outputs(&N));
	
	if(sigmaTarget) printf("Error in sigma from target value = %le\n", sigma/sigmaTarget-1.);
	
	if(plotDensities)
	{	string fname = fexName + "/planar-psi-n";
		N *= 1. / idgas.get_Nbulk();
		printToFile(N, fname.c_str());
		
		//double pMol = fex.getMolecule()->get_dipole();
		//logPrintf("Dipole moment = %lg Debye, bulk epsilon = %lg\n", pMol/0.393430307, 1 + 4*M_PI*idgas.get_Nbulk()*pMol*pMol/(3*T));
	
		FILE* pp = popen("gnuplot -persist", "w");
		fprintf(pp, "set xlabel \"x [bohr]\"\n");
		fprintf(pp, "set ylabel \"n [bohr^-3]\"\n");
		fprintf(pp, "plot \"%s\" using 1:2 title \"N / Nbulk\" w l\n", fname.c_str());
		fprintf(pp, "set key bottom right\n");
		fclose(pp);
	}
	return sigma;
}

void cavitation(GridInfo::CoordinateSystem coord, FILE* fp, bool model)
{
	const int nGrid = 512;
	const double hGrid = 0.125;
	const double T = 298*Kelvin;
	const double P = 101.3*KPascal;
	
	GridInfo gInfo(coord, nGrid, hGrid);
	SO3quad quad(QuadEuler, 3, 12, coord==GridInfo::Spherical ? 1 : 0);
	TranslationOperatorLspline trans(gInfo);
	FluidMixture fluidMixture(gInfo, T);
	FexClass fex(fluidMixture);
	IdealGasPomega idgas(&fex, 1.0, quad, trans);
	//IdealGasMonoatomic idgas(&fex, 1.0);
	fluidMixture.setPressure(P, PASTE_TOKENS(Nbulk, ChosenLiq));
	
	MinimizeParams mp;
	mp.alphaTstart = 3e4;
	mp.nDim = gInfo.S * fluidMixture.get_nIndep();
	mp.energyLabel = "Phi";
	mp.nIterations=500;
	mp.energyDiffThreshold=1e-11;
	mp.dirUpdateScheme = MinimizeParams::FletcherReeves;
	
	//Create list of sphere/cylinder radii equally spaced in curvature:
	std::list<double> Rlist;
	double Roffs = (coord==GridInfo::Spherical ? 0.5 : 0.25)*hGrid; //make sure the hard wall radii are midway between grid points for each sampling grid
	for(int i=0; i<40; i++) Rlist.push_back(round(36./((i+1)*hGrid))*hGrid + Roffs); //Equally spaced in curvature
	Rlist.sort(); Rlist.unique(); //remove duplicates
	std::vector<double> Rarr(Rlist.begin(), Rlist.end());

	//Prepare cavitation model:
	double Pvap = PASTE_TOKENS(Pvap, ChosenLiq);
	double Nbulk = PASTE_TOKENS(Nbulk, ChosenLiq);
	double sigma = PASTE_TOKENS(sigma, ChosenLiq);
	const double Gamma = log(Nbulk*T/Pvap)-1.;
	double d = 2.*fex.vdwRadius();
	const double Cp = 15. * (sigma/(d*Nbulk*T) - (1+Gamma)/6);
	printf("Cp = %lg  Gamma = %lg  d = %lg A\n", Gamma, Cp, d/Angstrom);

	SphericalKernel w(gInfo.S);
	for(int i=0; i<gInfo.S; i++)
		w[i] = gsl_sf_bessel_j0(gInfo.G[i]*d);

	for(int isDroplet=0; isDroplet<2; isDroplet++)
	{	for(double R: Rarr)
		{
			double Area = (coord==GridInfo::Spherical ? 4*M_PI*R*R : 2*M_PI*R);
			double Volume = (coord==GridInfo::Spherical ? 4*M_PI*R*R*R/3 : M_PI*R*R);
			
			//Setup shape function:
			ScalarField s(&gInfo);
			double* sData = s.data();
			for(int i=0; i<gInfo.S; i++)
				sData[i] = isDroplet
					? (gInfo.r[i]<R ? 1. : 0.)
					: (gInfo.r[i]<R ? 0. : 1.);
			
			if(model) //Compute model value:
			{	ScalarField sbar = I(w * J(s));
				double sigmaModel = integral(Diag(Nbulk*T*sbar)*(Diag(1-sbar)*(Gamma+Diag(sbar)*(1.-Gamma + Cp*(1-sbar))))) / Area;
				fprintf(fp, "%.8lf\t%.8le\n", (isDroplet?-1.:1.)/R, sigmaModel);
			}
			else //Compute classical DFT value:
			{	idgas.V[0] = 1. - s;
				idgas.set_Nnorm(isDroplet ? idgas.get_Nbulk()*Volume : 0.); //fixed N mode for droplet
				fluidMixture.initState(0.15);
				double sigmaCDFT = (fluidMixture.minimize(mp) + P * (isDroplet ? Volume-gInfo.Volume() : -Volume)) / Area;
				fprintf(fp, "%.8lf\t%.8le\n", (isDroplet?-1.:1.)/R, sigmaCDFT);
			}
		}
		fprintf(fp, "\n\n");
	}
}

//Compute and return pair correlations for a specific temperature and bulk density
ScalarFieldCollection pairCorrelations(const GridInfo& gInfo)
{
	SO3quad quad(QuadTetrahedron, 3); //Not used and need it only to create an IdealGas object. Just picked the smallest quadrature
	TranslationOperatorLspline trans(gInfo); //Not used and need it only to create an IdealGas object.
	
	double T = 298*Kelvin;
	FluidMixture fluidMixture(gInfo, T);
	FexClass fex(fluidMixture); //Selected excess functionals using macros above
	IdealGasPomega idgas(&fex, 1.0, quad, trans); //Not used and needed only to register fex with fluidMixture
	
	//Compute and return pair correlations:
	std::vector<double> Nmol(1, PASTE_TOKENS(Nbulk, ChosenLiq));
	return fluidMixture.getPairCorrelations(Nmol);
}

void printToFileCumulative(const ScalarField& g, const char* filename)
{	const double* gData = g.data();
	const std::vector<double>& r = g.gInfo->r;
	const std::vector<double>& w = g.gInfo->w;
	//Write file:
	FILE* fp = fopen(filename, "w");
	if(!fp) die("Could not open %s for writing.\n", filename)
	double gInt = 0.;
	for(int j=0; j<g.gInfo->S; j++)
	{	gInt += 0.5*gData[j]*w[j];
		fprintf(fp, "%le\t%le\t%le\n", r[j], gData[j], gInt);
		gInt += 0.5*gData[j]*w[j];
	}
	fclose(fp);
}

void saveCorrelationsCHCl3(const ScalarFieldCollection& gXX)
{	printToFileCumulative(gXX[0], (fexName+"/gCC").c_str());
	printToFileCumulative(gXX[1], (fexName+"/gCCl").c_str());
	printToFileCumulative(gXX[2], (fexName+"/gCH").c_str());
	printToFileCumulative(gXX[3], (fexName+"/gClCl").c_str());
	printToFileCumulative(gXX[4], (fexName+"/gClH").c_str());
}
void saveCorrelationsCCl4(const ScalarFieldCollection& gXX)
{	printToFileCumulative(gXX[0], (fexName+"/gCC").c_str());
	printToFileCumulative(gXX[1], (fexName+"/gCCl").c_str());
	printToFileCumulative(gXX[2], (fexName+"/gClCl").c_str());
}

int main(int argc, char** argv)
{	initSystem(argc, argv);

	#ifdef CAVITATION_CDFT
		FILE* fp = fopen((fexName+"/sigmaCDFT").c_str(), "w");
		cavitation(GridInfo::Spherical, fp, false);
		cavitation(GridInfo::Cylindrical, fp, false);
		fclose(fp);
	#endif
	#ifdef CAVITATION_MODEL
		FILE* fp = fopen((fexName+"/sigmaModel").c_str(), "w");
		cavitation(GridInfo::Spherical, fp, true);
		cavitation(GridInfo::Cylindrical, fp, true);
		fclose(fp);
	#endif

	#ifdef PLANAR_TDEP
		FILE* fp = fopen((fexName+"/sigmavsT").c_str(), "w");
		for(double Tcel=-10.0; Tcel<61.0; Tcel+=5.0)
		{	double T = (Tcel+273.16)*Kelvin;
			double sigmaSI = testPlanar(T) / (1e-3*Newton/meter);
			logPrintf("\n------------ T = %lf C, sigma = %le mN/m ---------------\n\n", Tcel, sigmaSI);
			fprintf(fp, "%lf %le\n", Tcel, sigmaSI);
		}
		fclose(fp);
	#endif
		
	#ifdef PLANAR
		testPlanar(298*Kelvin, PASTE_TOKENS(sigma, ChosenLiq), true);
	#endif
	
	#ifdef CORRFUNC
		GridInfo gInfo(GridInfo::Spherical, 1024, 0.0625); //Fine spherical grid to get good resolution in gXX
		PASTE_TOKENS(saveCorrelations, ChosenLiq)(pairCorrelations(gInfo));
	#endif
	
	finalizeSystem();
	return 0.;
}
