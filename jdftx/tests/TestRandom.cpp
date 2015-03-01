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
#include <core/ScalarFieldIO.h>
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
	double Phi_n0=0, Phi_n1=0, Phi_n2=0, Phi_n3=0, Phi_n1v[]={0,0,0}, Phi_n2v[]={0,0,0,}, Phi_n2m[]={0,0,0,0,0};
	phiFMT_calc(0, &n0, &n1, &n2, &n3, const_pArr(n1v,3), const_pArr(n2v,3), const_pArr(n2m,5),
		&Phi_n0, &Phi_n1, &Phi_n2, &Phi_n3, pArr(Phi_n1v,3), pArr(Phi_n2v,3), pArr(Phi_n2m,5));
	#define FDTEST(dir) \
		logPrintf("FDtesting along " #dir "\n"); \
		for(double h=1e-9; h<1e-1; h*=10) \
		{	double Phi_n0P=0, Phi_n1P=0, Phi_n2P=0, Phi_n3P=0, Phi_n1vP[]={0,0,0}, Phi_n2vP[]={0,0,0,}, Phi_n2mP[]={0,0,0,0,0}; \
			dir += h; \
			double phiP = phiFMT_calc(0, &n0, &n1, &n2, &n3, const_pArr(n1v,3), const_pArr(n2v,3), const_pArr(n2m,5), \
					&Phi_n0P, &Phi_n1P, &Phi_n2P, &Phi_n3P, pArr(Phi_n1vP,3), pArr(Phi_n2vP,3), pArr(Phi_n2mP,5)); \
			dir -= h; \
			double Phi_n0M=0, Phi_n1M=0, Phi_n2M=0, Phi_n3M=0, Phi_n1vM[]={0,0,0}, Phi_n2vM[]={0,0,0,}, Phi_n2mM[]={0,0,0,0,0}; \
			dir -= h; \
			double phiM = phiFMT_calc(0, &n0, &n1, &n2, &n3, const_pArr(n1v,3), const_pArr(n2v,3), const_pArr(n2m,5), \
					&Phi_n0M, &Phi_n1M, &Phi_n2M, &Phi_n3M, pArr(Phi_n1vM,3), pArr(Phi_n2vM,3), pArr(Phi_n2mM,5)); \
			dir += h; \
			logPrintf("\th=%le: Ratio=%.15lf\n", h, (phiP-phiM)/(2*h*Phi_##dir)); \
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

void checkSymmetry(ScalarField a, ScalarField b, const TranslationOperator& T, const vector3<> t)
{	ScalarField c=0; T.taxpy(t, 1.0, a, c); double n1 = dot(c, b);
	ScalarField d=0; T.taxpy(-t, 1.0, b, d); double n2 = dot(d, a);
	logPrintf("Symmetry error = %le\n", fabs(n1-n2)/fabs(n1));
}

void testTranslationOperators(const ScalarFieldArray& N)
{	const GridInfo& gInfo = N[0]->gInfo;
	vector3<> rCenter = gInfo.R * vector3<>(0.5,0.5,0.5);

	//vector3<> offset = 1.5*gInfo.h[0]-0.5*gInfo.h[1]+0.5*gInfo.h[2];
	vector3<> offset(0.567, 0.243, -0.7578);
	//vector3<> offset(1.5, 2.0, 8.0);
	TranslationOperatorFourier opF(gInfo);
	TranslationOperatorSpline opC(gInfo, TranslationOperatorSpline::Constant);
	TranslationOperatorSpline opL(gInfo, TranslationOperatorSpline::Linear);
	ScalarField test0(ScalarFieldData::alloc(gInfo,isGpuEnabled())); initRandom(test0);
	ScalarField test1(ScalarFieldData::alloc(gInfo,isGpuEnabled())); initRandom(test1);
	RealKernel gauss(gInfo); initGaussianKernel(gauss, 0.0);
	test0 = I(gauss*J(test0));
	test1 = I(gauss*J(test1));
	checkSymmetry(test0, test1, opF, offset);
	checkSymmetry(test0, test1, opC, offset);
	checkSymmetry(test0, test1, opL, offset);

	printStats(N[0], "N0");
	printStats(N[1], "N1");
	ScalarFieldArray NoffsF(N.size()), NoffsC(N.size()), NoffsL(N.size());
	TIME("TranslationFourier", globalLog,
		opF.taxpy(offset, 1.0, N[0], NoffsF[0]); ) printStats(NoffsF[0], "NF0");
		opF.taxpy(offset, 1.0, N[1], NoffsF[1]);   printStats(NoffsF[1], "NF1");
	TIME("TranslationConstSpline", globalLog,
		opC.taxpy(offset, 1.0, N[0], NoffsC[0]); ) printStats(NoffsC[0], "NC0");
		opC.taxpy(offset, 1.0, N[1], NoffsC[1]);   printStats(NoffsC[1], "NC1");
	TIME("TranslationLinearSpline", globalLog,
		opL.taxpy(offset, 1.0, N[0], NoffsL[0]); ) printStats(NoffsL[0], "NL0");
		opL.taxpy(offset, 1.0, N[1], NoffsL[1]);   printStats(NoffsL[1], "NL1");
	vector3<> offsetCenter = rCenter + offset;
	saveSphericalized(&NoffsF[0], NoffsF.size(), "random-ArNe-TF", 0.25, &offsetCenter);
	saveSphericalized(&NoffsC[0], NoffsC.size(), "random-ArNe-TC", 0.25, &offsetCenter);
	saveSphericalized(&NoffsL[0], NoffsL.size(), "random-ArNe-TL", 0.25, &offsetCenter);

}

int main(int argc, char** argv)
{	initSystem(argc, argv);

	//fdtestFMT(); return 0;

	//double T=298*Kelvin;
	//double T=130*Kelvin;
	double T=75*Kelvin;
	GridInfo gInfo;
	//gInfo.S = vector3<int>(64, 64, 64); double hGrid = 0.5;
	gInfo.S = vector3<int>(128, 128, 128); double hGrid = 0.5;
	gInfo.R = Diag(hGrid * gInfo.S);
	gInfo.initialize();
	vector3<> rCenter = gInfo.R * vector3<>(0.5,0.5,0.5);

	double RhsAr = 1.702*Angstrom, epsAr = 119.8*Kelvin;
	double RhsNe = 1.391*Angstrom, epsNe = 3.2135e-3*eV;
	std::shared_ptr<FluidComponent> componentAr,componentNe;
	componentAr = std::make_shared<FluidComponent>(FluidComponent::CustomAnion, T, FluidComponent::MeanFieldLJ);	
	componentAr->molecule.setModelMonoatomic("Ar", +1., RhsAr);
	componentAr->epsLJ = epsAr;
	componentAr->Nbulk = 2e-3;
	componentNe = std::make_shared<FluidComponent>(FluidComponent::CustomCation, T, FluidComponent::MeanFieldLJ);
	componentNe->molecule.setModelMonoatomic("Ne", -2., RhsNe);
	componentNe->epsLJ = epsNe;
	componentNe->Nbulk = 1e-3;
	
	FluidMixture fluidMixture(gInfo, T);
	componentAr->addToFluidMixture(&fluidMixture);
	componentNe->addToFluidMixture(&fluidMixture);
	Fmix_LJ fmixLJ(&fluidMixture, componentAr, componentNe, sqrt(epsAr*epsNe), RhsAr+RhsNe);

	fluidMixture.initialize(1000*Bar);

	const double Radius = 3.;
	nullToZero(componentAr->idealGas->V[0], gInfo);
	applyFunc_r(gInfo, initHardSphere, rCenter, Radius, 1, componentAr->idealGas->V[0]->data());
	componentNe->idealGas->V[0] = componentAr->idealGas->V[0];

	ScalarField rhoExternal(ScalarFieldData::alloc(gInfo));
	applyFunc_r(gInfo, initHardSphere, rCenter, Radius, 1.016177/(4*M_PI*pow(Radius,3)/3), rhoExternal->data());
	fluidMixture.rhoExternal = J(rhoExternal);

	//componentNe.Nnorm = 52;

	//fluidMixture.verboseLog = true;
	fluidMixture.initState(0.15*0);
	//fluidMixture.loadState("random.state");

	MinimizeParams mp;
	mp.fpLog = globalLog;
	mp.nDim = gInfo.nr*fluidMixture.get_nIndep();
	mp.nIterations=200;
	mp.knormThreshold=5e-12;
	mp.dirUpdateScheme=MinimizeParams::HestenesStiefel;
	mp.fdTest = true;

	fluidMixture.minimize(mp);
	fluidMixture.saveState("random.state");

	ScalarFieldArray N; ScalarFieldTilde Phi_rhoExternalTilde;
	fluidMixture.getFreeEnergy(FluidMixture::Outputs(&N, 0, &Phi_rhoExternalTilde));

	logPrintf("Qext    = %lf\n", integral(rhoExternal));
	logPrintf("num(Ar) = %lf\n", integral(N[0]));
	logPrintf("num(Ne) = %lf\n", integral(N[1]));

	N[0] *= 1.0/componentAr->idealGas->get_Nbulk();
	N[1] *= 1.0/componentNe->idealGas->get_Nbulk();
	saveSphericalized(&N[0], N.size(), "random-ArNe", 0.25);
	//ScalarField dTot = I(Phi_rhoExternalTilde - 4*M_PI*Linv(O(J(rhoExternal))));
	//saveSphericalized(&dTot, 1, "random-dTot", 0.25);

	//testTranslationOperators(N);
	
	finalizeSystem();
	return 0;
}
