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
#include <fluid/Fex_ScalarEOS.h>
#include <core/Units.h>
#include <core/ScalarFieldIO.h>

void initSphere(int i, vector3<> r, const vector3<>& r0, double radius, double hIn, double hOut, double* v)
{	v[i] = ((r - r0).length() < radius ? hIn : hOut);
}

//A wall with an attractive well just next to it (lets say a crude planarly averaged metal surface)
void initAttractiveWall(int i, vector3<> r, double xWall, double dxWall, double wallHeight, double wellDepth, double* v)
{	double& x = r[0];
	v[i] = x<xWall ? wallHeight : -wellDepth*exp(-pow((x-xWall)/dxWall,2));
}

//Add preconditioner to fluid mixture to boost low spatial frequency components with a cutoff length scale aTyp
class PreconditionedFluidMixture : public FluidMixture
{
	RealKernel K;
	static void setK(int i, double Gsq, double* K, double aTyp) { K[i] = 1.0/(1.0 + pow(aTyp,2)*Gsq); }
	
public:
	
	PreconditionedFluidMixture(const GridInfo& gInfo, double T, double aTyp)
	: FluidMixture(gInfo, T), K(gInfo)
	{	applyFuncGsq(gInfo, setK, K.data, aTyp);
	}

	ScalarFieldArray precondition(const ScalarFieldArray& grad)
	{	ScalarFieldArray Kgrad(grad.size());
		for(unsigned i=0; i<Kgrad.size(); i++) Kgrad[i] = I(K*J(grad[i]));
		return Kgrad;
	}
};


//################ Droplet on attractive wall ####################
int main(int argc, char** argv)
{	initSystem(argc, argv);

	GridInfo gInfo;
	gInfo.S = vector3<int>(64, 64, 64); double hGrid=1.0;
	//gInfo.S = vector3<int>(128, 128, 128); double hGrid=0.5;
	gInfo.R = Diag(gInfo.S * hGrid);
	gInfo.initialize();

	double T = 298*Kelvin;
	FluidComponent component(FluidComponent::H2O, T, FluidComponent::ScalarEOS);
	component.s2quadType = QuadOctahedron;
	component.Nnorm = 270;
	
	PreconditionedFluidMixture fluidMixture(gInfo, T, 1.0);
	component.addToFluidMixture(&fluidMixture);
	double p = 1.01325*Bar;
	logPrintf("pV = %le\n", p*gInfo.detR);
	fluidMixture.initialize(p);

	#define geomName "AttractiveWall-3bohr-3.0kT/drop_plane"

	bool loadState=false;
	const char stateFilename[] = "TestFixedN/" geomName "_state.bin";

	//Initialize potential: planar wall with attractive well:
	nullToZero(component.idealGas->V, gInfo);
	double xWall = 0.5*gInfo.R(0,0)-21.0;
	double dxWall = 2.0;
	applyFunc_r(gInfo, initAttractiveWall, xWall, dxWall, 100*T, 3.*T, component.idealGas->V[0]->data());

	if(loadState)
		fluidMixture.loadState(stateFilename);
	else
	{	//Initialize state biased towards center of cell
		const double muSet=-5.0;
		const double Rdroplet = 22.0; //guess droplet size:
		nullToZero(fluidMixture.state, gInfo, fluidMixture.get_nIndep());
		applyFunc_r(gInfo, initSphere, gInfo.R*vector3<>(0.5,0.5,0.5), Rdroplet, 0.0, muSet, fluidMixture.state[0]->data());

		RealKernel gauss(gInfo); initGaussianKernel(gauss, 2.0);
		fluidMixture.state[0] = I(gauss*J(fluidMixture.state[0]));

		ScalarField wallMask(ScalarFieldData::alloc(gInfo));
		applyFunc_r(gInfo, initAttractiveWall, xWall, 2*dxWall, 0.0, 1.0, wallMask->data());
		fluidMixture.state[0] *= (wallMask+1.0);
	}

	MinimizeParams mp;
	mp.fpLog = globalLog;
	mp.nDim = 2*gInfo.nr;
	mp.nIterations=10;
	mp.knormThreshold=1e-11;
	mp.dirUpdateScheme = MinimizeParams::FletcherReeves;
	//mp.dirUpdateScheme = MinimizeParams::SteepestDescent;
	//mp.updateTestStepSize = false;
	mp.fdTest = !loadState;

	int sysRet=system("mkdir -p TestFixedN/" geomName "_img");
	if(sysRet) { logPrintf("Error making image directory\n"); mpiUtil->exit(sysRet); }

	for(int loopCount=0; loopCount<100; loopCount++)
	{
		ScalarFieldArray N;
		TIME("getOmega calculation (with gradient)", globalLog,
			double omega = fluidMixture.getFreeEnergy(FluidMixture::Outputs(&N));
			if(std::isnan(omega)) break; //Don't save state after it has become nan
		);
		logPrintf("Ntot = %lf\n", gInfo.dV*sum(N[0]));

		logPrintf("Saving state:\n");
		fluidMixture.saveState(stateFilename);
		saveDX(N[0], "TestFixedN/" geomName "_nO");
		saveDX(N[1], "TestFixedN/" geomName "_nH");
		saveSphericalized(&N[0], 2, "TestFixedN/" geomName "_n.spherical", 0.25);
		//Invoke octave to create an image:
		FILE* pp = popen("octave -q", "w");
		fprintf(pp, "imwrite( waterSlice(\"TestFixedN/" geomName "_n%%s.bin\", [%d %d %d], 1, %d, 1e-2),", gInfo.S[0], gInfo.S[1], gInfo.S[2], gInfo.S[2]/2);
		fprintf(pp, " \"TestFixedN/" geomName "_img/img%04d.png\"); exit;\n", loopCount);
		fflush(pp); pclose(pp);

		logPrintf("Starting CG:\n");
		TIME("minimize", globalLog,
			fluidMixture.minimize(mp);
		);
	}
	
	finalizeSystem();
	return 0;
}


/*
//################ Spherical droplet ####################
int main(int argc, char** argv)
{	initSystem(argc, argv);

	#define geomName "sphere"

	double Nset;
	if(argc<2) die("Need to specify number of molecules as second argument")
	sscanf(argv[1], "%lg", &Nset);
	char logFilename[256]; sprintf(logFilename, "TestFixedN/" geomName "_N%s_cg.logPrintf", argv[1]);

	GridInfo gInfo;
	gInfo.S = vector3<int>(64, 64, 64); double hGrid=1.0;
	gInfo.R = Diag(hGrid * gInfo.S);
	gInfo.fpLog = fopen(logFilename, "w");
	gInfo.initialize();

	//========================  Setup quadrature for angular integration
	//Tetrahedron plato(0,4);
	Octahedron plato(1);
	//Icosahedron plato(0,6);
	//Tetrahedron plato(1);
	//Octahedron plato(1);
	SO3quad quad(&plato, gInfo.fpLog);

	//======================== Water molecule geometry
	RISMGeometry geom = rismGeometrySPCE();

	//======================== Excess functional
	//Fex_wdG_Cexp_STP fex(gInfo, quad, geom);
	//Fex_wdG_wpG_Cexp_STP fex(gInfo, quad, geom, 1.30*Angstrom);
	//Fex_wdG_Cexp_STP fex(gInfo, quad, geom, 1.30*Angstrom);
	//Fex_wdS_STP fex(gInfo, quad, geom);
	Fex_wdS fex(gInfo, quad, geom);
	//Fex_wdT_STP fex(gInfo, quad, geom);
	//Fex_wdT fex(gInfo, quad, geom);
	//Fex_wdT_wpT fex(gInfo, quad, geom);
	typedef decltype(fex) Fex;

	//======================== Representation
	RepresentationMuEps<Fex> water(gInfo, quad, geom, fex);
	typedef decltype(water.state) State;

	water.fixN(Nset);

	bool loadState=false;
	const char stateFilename[] = "TestFixedN/" geomName "_state.%c.bin";

	if(loadState)
		water.loadState(stateFilename);
	else
	{	//Initialize state biased towards center of cell, but no potential!
		const double muSet=-5.0;
		double Rdroplet = pow((Nset/fex.nl-gInfo.detR*exp(muSet))*3.0/(4*M_PI*(1-exp(muSet))), 1.0/3); //radius of droplet if no shell structure
		logPrintf("Ideal Rdroplet = %lf bohr\n", Rdroplet);
		nullToZero(water.state, gInfo);
		applyFuncGsq(gInfo, initSphere, 0.5*(gInfo.R.set_col(0]+gInfo.R.set_col(1]+gInfo.R.set_col(2]), Rdroplet, 0.0, muSet, water.state[0]->data());

		//FDtest:
		//State dpsi(gInfo); initRandomFlat(dpsi); dpsi *= (-1e-3/fex.T);
		//fdTest(water.state, dpsi, MinimizeParams(), water);
	}

	MinimizeParams mp;
	mp.fpLog = globalLog;
	mp.alphaMax = 1e300;
	mp.nDim = 2*gInfo.nr;
	mp.nIterations=1000;
	mp.knormThreshold=1e-11;
	mp.dirUpdateScheme = MinimizeParams::FletcherReeves;
	//mp.dirUpdateScheme = MinimizeParams::SteepestDescent;
	//mp.updateTestStepSize = false;

	LowfreqPrecon precon(gInfo, 1.0);

	logPrintf("Starting CG:\n");
	TIME("minimize", globalLog,
		ncgSolve(water.state, mp, water, precon);
	);

// 	ScalarFieldOH n;
// 	TIME("getOmega calculation (with gradient)", globalLog,
// 		water.getOmega(&n);
// 	);
// 	logPrintf("Ntot = %lf\n", gInfo.dV*sum(n.O()));
//
// 	logPrintf("Saving state:\n");
// 	water.saveState(stateFilename);
// 	ScalarField saveR[2] = {n.O(), n.H()};
// 	saveDX(n.O(), "TestFixedN/" geomName "_nO");
// 	saveDX(n.H(), "TestFixedN/" geomName "_nH");
// 	saveSphericalized(saveR, 2, "TestFixedN/" geomName "_n.spherical", 0.25);
	
	finalizeSystem();
	return 0;
}
*/
/*
//################ Dimer ####################
int main(int argc, char** argv)
{	initSystem(argc, argv);

	GridInfo gInfo;
	gInfo.S = vector3<int>(40,40,40); double hGrid=1.0;
	//gInfo.S = vector3<int>(128, 128, 128); double hGrid=0.5;
	gInfo.R = Diag(hGrid * gInfo.S);
	gInfo.initialize();

	//========================  Setup quadrature for angular integration
	//Tetrahedron plato(0,4);
	Octahedron plato(1);
	//Icosahedron plato(0,6);
	//Tetrahedron plato(1);
	//Octahedron plato(1);
	SO3quad quad(&plato);

	//======================== Water molecule geometry
	RISMGeometry geom = rismGeometrySPCE();

	//======================== Excess functional
	//Fex_wdG_Cexp_STP fex(gInfo, quad, geom);
	//Fex_wdG_wpG_Cexp_STP fex(gInfo, quad, geom, 1.30*Angstrom);
	//Fex_wdG_Cexp_STP fex(gInfo, quad, geom, 1.30*Angstrom);
	//Fex_wdS_STP fex(gInfo, quad, geom);
	Fex_wdS fex(gInfo, quad, geom);
	//Fex_wdT_STP fex(gInfo, quad, geom);
	//Fex_wdT fex(gInfo, quad, geom);
	//Fex_wdT_wpT fex(gInfo, quad, geom);
	typedef decltype(fex) Fex;

	//======================== Representation
	RepresentationMuEps<Fex> water(gInfo, quad, geom, fex);
	typedef decltype(water.state) State;

	#define geomName "dimer"

	bool loadState=false;
	const char stateFilename[] = "TestFixedN/" geomName "_state.%c.bin";

	water.fixN(2);

	if(loadState)
		water.loadState(stateFilename);
	else
	{	//Initialize state in two spheres separated by roughly the water-dimer equilibrium distance
		const double muOut=-5.0;
		const double muIn=2.0;
		const double Rdroplet = pow(1/(4*M_PI*fex.nl*exp(muIn)/3.0), 1.0/3);
		nullToZero(water.state, gInfo);
		applyFunc_r(gInfo, initSphere, 0.5*(gInfo.R.set_col(0]+gInfo.R.set_col(1]+gInfo.R.set_col(2]), Rdroplet, muIn-0.5*muOut, 0.5*muOut, water.state[0]->data());
		ScalarFieldTilde trans1(ScalarFieldTildeData::alloc(gInfo)); initTranslation(trans1, vector3<>(-2.532,0.017,-0.039));
		ScalarFieldTilde trans2(ScalarFieldTildeData::alloc(gInfo)); initTranslation(trans2, vector3<>(+2.567,0.087,0.013)); //slightly assymetric w.r.t grid to help symmetry breaking, if any
		water.state[0] = I((trans1+trans2)*J(water.state[0]));

		//FDtest:
		//State dpsi(gInfo); initRandomFlat(dpsi); dpsi *= (-1e-3/fex.T);
		//fdTest(water.state, dpsi, MinimizeParams(), water);
	}

	MinimizeParams mp;
	mp.fpLog = globalLog;
	mp.alphaMax = 1e300;
	mp.nDim = 2*gInfo.nr;
	mp.nIterations=20;
	mp.knormThreshold=1e-11;
	mp.dirUpdateScheme = MinimizeParams::FletcherReeves;
	//mp.dirUpdateScheme = MinimizeParams::SteepestDescent;
	//mp.updateTestStepSize = false;

	LowfreqPrecon precon(gInfo, 1.0);

	//int sysRet=system("mkdir -p TestFixedN/" geomName "_img");
	//if(sysRet) { logPrintf("Error making image directory\n"); mpiUtil->exit(sysRet); }

	for(int loopCount=0; loopCount<20; loopCount++)
	{
		ScalarFieldOH n;
		TIME("getOmega calculation (with gradient)", globalLog,
			double omega = water.getOmega(&n);
			if(std::isnan(omega)) break; //Don't save state after it has become nan
		);
		logPrintf("Ntot = %lf\n", gInfo.dV*sum(n.O()));

		logPrintf("Saving state:\n");
		water.saveState(stateFilename);
		ScalarField saveR[2] = {n.O(), n.H()};
		saveDX(n.O(), "TestFixedN/" geomName "_nO");
		saveDX(n.H(), "TestFixedN/" geomName "_nH");
		saveSphericalized(saveR, 2, "TestFixedN/" geomName "_n.spherical", 0.25);

		//Invoke octave to create an image:
// 		FILE* pp = popen("octave", "w");
// 		fprintf(pp, "imwrite( waterSlice(\"TestFixedN/" geomName "_n%%s.bin\", [%d %d %d], 1, %d, 1e-2), \"TestFixedN/" geomName "_img/yz%04d.png\"); ", gInfo.S[0], gInfo.S[1], gInfo.S[2], gInfo.S[2]/2, loopCount);
// 		fprintf(pp, "imwrite( waterSlice(\"TestFixedN/" geomName "_n%%s.bin\", [%d %d %d], 2, %d, 1e-2), \"TestFixedN/" geomName "_img/zx%04d.png\"); ", gInfo.S[0], gInfo.S[1], gInfo.S[2], gInfo.S[1]/2, loopCount);
// 		fprintf(pp, "imwrite( waterSlice(\"TestFixedN/" geomName "_n%%s.bin\", [%d %d %d], 3, %d, 1e-2), \"TestFixedN/" geomName "_img/xy%04d.png\"); ", gInfo.S[0], gInfo.S[1], gInfo.S[2], gInfo.S[0]/2, loopCount);
// 		fprintf(pp, "exit;\n");
// 		fflush(pp); pclose(pp);

		logPrintf("Starting CG:\n");
		TIME("minimize", globalLog,
			ncgSolve(water.state, mp, water, precon);
		);
	}
	
	finalizeSystem();
	return 0;
}
*/
