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
#include <fluid/IdealGasMonoatomic.h>
#include <fluid/Fex_LJ.h>
#include <fluid/TranslationOperator.h>
#include <fluid/MixedFMT_internal.h>
#include <core/Units.h>
#include <core/DataIO.h>
#include <cstdio>


void initHardSphere(int i, vector3<> r, const vector3<>& r0, double radius, double height, double* phi)
{	phi[i] = ((r - r0).length() < radius ? height : 0.0);
	//phi[i] = height/(1+pow((r - r0).length()/radius,6));
}


std::vector<double*> pArr(double* x, int N)
{	std::vector<double*> v(N);
	for(int i=0; i<N; i++) v[i]=x+i;
	return v;
}
std::vector<const double*> const_pArr(const double* x, int N)
{	std::vector<const double*> v(N);
	for(int i=0; i<N; i++) v[i]=x+i;
	return v;
}

void fdtestFMT()
{	double n0=5e-3, n1=1e-2, n2=0.1, n3=0.6, n1v[]={2e-3,-3e-4,5e-3}, n2v[]={1e-2,-5e2,2e-2}, n2m[]={3e-2,2e-2,-2e-2,3e-2,5e-2};
	double grad_n0=0, grad_n1=0, grad_n2=0, grad_n3=0, grad_n1v[]={0,0,0}, grad_n2v[]={0,0,0,}, grad_n2m[]={0,0,0,0,0};
	phiFMT_calc(0, &n0, &n1, &n2, &n3, const_pArr(n1v,3), const_pArr(n2v,3), const_pArr(n2m,5),
		&grad_n0, &grad_n1, &grad_n2, &grad_n3, pArr(grad_n1v,3), pArr(grad_n2v,3), pArr(grad_n2m,5));
	#define FDTEST(dir) \
		printf("FDtesting along " #dir "\n"); \
		for(double h=1e-9; h<1e-1; h*=10) \
		{	double grad_n0P=0, grad_n1P=0, grad_n2P=0, grad_n3P=0, grad_n1vP[]={0,0,0}, grad_n2vP[]={0,0,0,}, grad_n2mP[]={0,0,0,0,0}; \
			dir += h; \
			double phiP = phiFMT_calc(0, &n0, &n1, &n2, &n3, const_pArr(n1v,3), const_pArr(n2v,3), const_pArr(n2m,5), \
					&grad_n0P, &grad_n1P, &grad_n2P, &grad_n3P, pArr(grad_n1vP,3), pArr(grad_n2vP,3), pArr(grad_n2mP,5)); \
			dir -= h; \
			double grad_n0M=0, grad_n1M=0, grad_n2M=0, grad_n3M=0, grad_n1vM[]={0,0,0}, grad_n2vM[]={0,0,0,}, grad_n2mM[]={0,0,0,0,0}; \
			dir -= h; \
			double phiM = phiFMT_calc(0, &n0, &n1, &n2, &n3, const_pArr(n1v,3), const_pArr(n2v,3), const_pArr(n2m,5), \
					&grad_n0M, &grad_n1M, &grad_n2M, &grad_n3M, pArr(grad_n1vM,3), pArr(grad_n2vM,3), pArr(grad_n2mM,5)); \
			dir += h; \
			printf("\th=%le: Ratio=%.15lf\n", h, (phiP-phiM)/(2*h*grad_##dir)); \
		}
	FDTEST(n0)
	FDTEST(n1)
	FDTEST(n2)
	FDTEST(n3)
	FDTEST(n1v[0])
	FDTEST(n1v[1])
	FDTEST(n1v[2])
	FDTEST(n2v[0])
	FDTEST(n2v[1])
	FDTEST(n2v[2])
	FDTEST(n2m[0])
	FDTEST(n2m[1])
	FDTEST(n2m[2])
	FDTEST(n2m[3])
	FDTEST(n2m[4])
	#undef FDTEST
}

void checkSymmetry(DataRptr a, DataRptr b, const TranslationOperator& T, const vector3<> t)
{	DataRptr c=0; T.taxpy(t, 1.0, a, c); double n1 = dot(c, b);
	DataRptr d=0; T.taxpy(-t, 1.0, b, d); double n2 = dot(d, a);
	printf("Symmetry error = %le\n", fabs(n1-n2)/fabs(n1));
}

int main(int argc, char** argv)
{	initSystem(argc, argv);

	//fdtestFMT(); return 0;

	/*
	SiteProperties Oprop(gInfo, 1.4*Angstrom,0.0*Angstrom, -0.8476, 0);
	SiteProperties Hprop(gInfo, 0.0*Angstrom,0.0*Angstrom, +0.4238, 0);
	SiteProperties Vprop(gInfo, 0.3*Angstrom,0.0*Angstrom, 0, 0);
	double rOH=1*Angstrom, rOV=Oprop.sphereRadius+Vprop.sphereRadius;

	Molecule voidWater("H2O",
		&Oprop, vector3<>(0,0,0),
		&Hprop, vector3<>(-1,-1,1)*(rOH/sqrt(3)), vector3<>(1,1,1)*(rOH/sqrt(3)),
		&Vprop, vector3<>(1,1,-1)*(rOV/sqrt(3)), vector3<>(-1,-1,-1)*(rOV/sqrt(3)),
		&Vprop, vector3<>(1,-1,1)*(rOV/sqrt(3)), vector3<>(-1,1,1)*(rOV/sqrt(3)) );

	printf("nSites = %d\n", voidWater.nSites);
	printf("nIndices = %d\n", voidWater.nIndices);
	for(int i=0; i<voidWater.nSites; i++)
	{	const Site& s = voidWater.site[i];
		printf("Site %d: index=%d, radius=%.1lfA, |pos|=%lfA\n",
			i, s.index, s.prop->sphereRadius/Angstrom, length(s.pos)/Angstrom);
	}
	*/
	//double T=298*Kelvin;
	//double T=130*Kelvin;
	double T=75*Kelvin;
	GridInfo gInfo;
	//gInfo.S = vector3<int>(64, 64, 64); double hGrid = 0.5;
	gInfo.S = vector3<int>(128, 128, 128); double hGrid = 0.5;
	gInfo.R = Diag(hGrid * gInfo.S);
	gInfo.initialize();
	vector3<> rCenter = gInfo.R * vector3<>(0.5,0.5,0.5);

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
	applyFunc_r(gInfo, initHardSphere, rCenter, 6.0, 0.005, idAr.V[0]->data());
	idNe.V[0] = idAr.V[0];

	DataRptr rhoExternal(DataR::alloc(gInfo));
	applyFunc_r(gInfo, initHardSphere, rCenter, 6.0, 1.016177/(4*M_PI*pow(6.0,3)/3), rhoExternal->data());
	//fluidMixture.rhoExternal = J(rhoExternal);

	idNe.set_Nnorm(52);

	//fluidMixture.verboseLog = true;
	fluidMixture.initState(0.15);
	//fluidMixture.loadState("random.state");

	MinimizeParams mp;
	mp.alphaTstart = 3e1;
	mp.nDim = gInfo.nr*fluidMixture.get_nIndep();
	mp.nIterations=100;
	mp.knormThreshold=5e-12;
	mp.dirUpdateScheme=MinimizeParams::HestenesStiefel;
	mp.fdTest = true;

	fluidMixture.minimize(mp);
	//fluidMixture.saveState("random.state");

	DataRptrCollection N;
	fluidMixture.getFreeEnergy(FluidMixture::Outputs(&N));

	printf("Qext    = %lf\n", integral(rhoExternal));
	printf("num(Ar) = %lf\n", integral(N[0]));
	printf("num(Ne) = %lf\n", integral(N[1]));

	N[0] *= 1.0/idAr.get_Nbulk();
	N[1] *= 1.0/idNe.get_Nbulk();
	saveSphericalized(&N[0], N.size(), "random-ArNe", 0.25);

	//-------- Test translation operators
	//vector3<> offset = 1.5*gInfo.h[0]-0.5*gInfo.h[1]+0.5*gInfo.h[2];
	vector3<> offset(0.567, 0.243, -0.7578);
	//vector3<> offset(1.5, 2.0, 8.0);
	TranslationOperatorFourier opF(gInfo);
	TranslationOperatorSpline opC(gInfo, TranslationOperatorSpline::Constant);
	TranslationOperatorSpline opL(gInfo, TranslationOperatorSpline::Linear);
	DataRptr test0(DataR::alloc(gInfo,isGpuEnabled())); initRandom(test0);
	DataRptr test1(DataR::alloc(gInfo,isGpuEnabled())); initRandom(test1);
	RealKernel gauss(gInfo); initGaussianKernel(gauss, 0.0);
	test0 = I(gauss*J(test0));
	test1 = I(gauss*J(test1));
	checkSymmetry(test0, test1, opF, offset);
	checkSymmetry(test0, test1, opC, offset);
	checkSymmetry(test0, test1, opL, offset);

	printStats(N[0], "N0");
	printStats(N[1], "N1");
	DataRptrCollection NoffsF(N.size()), NoffsC(N.size()), NoffsL(N.size());
	TIME("TranslationFourier", stdout,
		opF.taxpy(offset, 1.0, N[0], NoffsF[0]); ) printStats(NoffsF[0], "NF0");
		opF.taxpy(offset, 1.0, N[1], NoffsF[1]);   printStats(NoffsF[1], "NF1");
	TIME("TranslationConstSpline", stdout,
		opC.taxpy(offset, 1.0, N[0], NoffsC[0]); ) printStats(NoffsC[0], "NC0");
		opC.taxpy(offset, 1.0, N[1], NoffsC[1]);   printStats(NoffsC[1], "NC1");
	TIME("TranslationLinearSpline", stdout,
		opL.taxpy(offset, 1.0, N[0], NoffsL[0]); ) printStats(NoffsL[0], "NL0");
		opL.taxpy(offset, 1.0, N[1], NoffsL[1]);   printStats(NoffsL[1], "NL1");
	vector3<> offsetCenter = rCenter + offset;
	saveSphericalized(&NoffsF[0], NoffsF.size(), "random-ArNe-TF", 0.25, &offsetCenter);
	saveSphericalized(&NoffsC[0], NoffsC.size(), "random-ArNe-TC", 0.25, &offsetCenter);
	saveSphericalized(&NoffsL[0], NoffsL.size(), "random-ArNe-TL", 0.25, &offsetCenter);
}
