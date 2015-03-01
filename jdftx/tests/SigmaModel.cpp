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
#include <fluid/IdealGas.h>
#include <electronic/SphericalHarmonics.h>
#include <electronic/operators.h>

inline void initSphere(int i, vector3<> r, const vector3<>& rCenter, double radius, double valIn, double valOut, double* phi)
{	phi[i] = ((r - rCenter).length() < radius ? valIn : valOut);
}
inline void initCylinder(int i, vector3<> r, const vector3<>& rCenter, double radius, double valIn, double valOut, double* phi)
{	phi[i] = (hypot(r[1]-rCenter[1], r[2]-rCenter[2]) < radius ? valIn : valOut);
}

inline double wCavity_calc(double G, double d)
{	return bessel_jl(0, G*d); //corresponds to delta(d-r)
}

void printUsage()
{	logPrintf("Usage: SigmaModel <solventName> Sphere|Cylinder <radiusInBohrs> Droplet|Cavity Model|DFT\n");
}

extern EnumStringMap<FluidComponent::Name> solventMap;

int main(int argc, char** argv)
{	initSystem(argc, argv);
	
	//Parse commandline:
	if(argc != 6) { printUsage(); return 1; }
	//--- 1:
	FluidComponent::Name solventName;
	if(!solventMap.getEnum(argv[1], solventName)) { logPrintf("Unrecognized <solventName> '%s'\n", argv[1]); printUsage(); return 1; }
	//--- 2:
	bool cylindrical = false;
	if(string(argv[2])=="Cylinder") cylindrical = true;
	else if(string(argv[2])=="Sphere") cylindrical = false;
	else { logPrintf("Unrecognized geometry '%s'\n", argv[2]); printUsage(); return 1; }
	//--- 3:
	double radius = 0.;
	{	istringstream iss(argv[3]);
		iss >> radius;
		if(iss.fail()) { logPrintf("Could not parse <radiusInBohrs> from '%s'\n", argv[3]); printUsage(); return 1; }
		if(radius <= 0.) { logPrintf("<radiusInBohrs> must be > 0. (given %lg)\n", radius); printUsage(); return 1; }
	}
	//--- 4:
	bool droplet = false;
	if(string(argv[4])=="Droplet") droplet = true;
	else if(string(argv[4])=="Cavity") droplet = false;
	else { logPrintf("Unrecognized topology '%s', must be 'Droplet' or 'Cavity'\n", argv[4]); printUsage(); return 1; }
	//--- 5:
	bool model = false;
	if(string(argv[5])=="Model") model = true;
	else if(string(argv[5])=="DFT") model = false;
	else { logPrintf("Unrecognized mode '%s', must be 'Model' or 'DFT'\n", argv[5]); printUsage(); return 1; }
	
	//Create component:
	double T = 298*Kelvin;
	FluidComponent component(solventName, T, FluidComponent::ScalarEOS);
	component.s2quadType = QuadOctahedron;
	component.representation = FluidComponent::MuEps;
	
	//Setup simulation grid:
	double L = 2.*(radius + 3.*component.molecule.sites[0]->Rhs);
	GridInfo gInfo;
	if(cylindrical)
	{	gInfo.S = vector3<int>(1,512,512);
		gInfo.R = matrix3<>(1, L, L);
	}
	else
	{	gInfo.S = vector3<int>(128, 128, 128);
		gInfo.R = matrix3<>(L, L, L);
	}
	gInfo.initialize();
	
	//Initialize shape function:
	ScalarField s; nullToZero(s, gInfo);
	{	double valIn = droplet ? 1. : 0.;
		vector3<> rCenter = gInfo.R*vector3<>(0.5,0.5,0.5);
		if(cylindrical) applyFunc_r(gInfo, initCylinder, rCenter, radius, valIn, 1.-valIn, s->data());
		else  applyFunc_r(gInfo, initSphere, rCenter, radius, valIn, 1.-valIn, s->data());
	}
	double Area = (cylindrical ? 2*M_PI*radius : 4*M_PI*pow(radius,2));
	double Volume = (cylindrical ? M_PI*pow(radius,2) : (4*M_PI/3)*pow(radius,3));
	
	if(model) //Weighted-density model:
	{
		double NT = component.Nbulk * T;
		const double Gamma = log(NT/component.Pvap)-1.;
		double d = 2.*component.Rvdw;
		const double Cp = 15. * (component.sigmaBulk/(d*NT) - (1+Gamma)/6);
		RadialFunctionG w; w.init(0, gInfo.dGradial, gInfo.GmaxGrid, wCavity_calc, d);
		ScalarField sbar = I(w * J(s));
		double sigmaModel = integral(NT*sbar*(1-sbar)*(Gamma+sbar*(1.-Gamma + Cp*(1-sbar)))) / Area;
		logPrintf("\nSigmaModel_CurvatureAndSigma: %.8lf\t%.8le\n", (droplet?-1.:1.)/radius, sigmaModel);
		w.free();
	}
	else //Classical DFT:
	{
		FluidMixture fluidMixture(gInfo, T);
		component.addToFluidMixture(&fluidMixture);
		fluidMixture.initialize(0.);
		//--- Initialize external potential and state:
		component.idealGas->V[0] = 1.-s;
		nullToZero(component.idealGas->V, gInfo);
		fluidMixture.initState(0.15);
		if(droplet) component.Nnorm = component.idealGas->get_Nbulk() * Volume;
		//--- Minimize and print result:
		MinimizeParams mp;
		mp.fpLog = globalLog;
		mp.nDim = gInfo.nr * fluidMixture.get_nIndep();
		mp.energyLabel = "Phi";
		mp.nIterations = 100;
		mp.energyDiffThreshold = 1e-6 * (component.sigmaBulk * Area);
		double sigmaDFT = fluidMixture.minimize(mp) / Area;
		logPrintf("\nSigmaModel_CurvatureAndSigma: %.8lf\t%.8le\n", (droplet?-1.:1.)/radius, sigmaDFT);
	}
	
	finalizeSystem();
	return 0;
}
