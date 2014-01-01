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

#include <core/BlasExtra.h>
#include <core/Util.h>
#include <core/Operators.h>
#include <core/GridInfo.h>
#include <core/DataIO.h>
#include <core/Coulomb.h>
#include <core/LatticeUtils.h>
#include <core/Random.h>
#include <fluid/SO3quad.h>
#include <electronic/SphericalHarmonics.h>
#include <electronic/ExCorr_internal_GGA.h>
#include <electronic/operators.h>
#include <gsl/gsl_sf.h>
#include <stdlib.h>


class OperatorTest
{	const GridInfo& gInfo;

public:
	OperatorTest(const GridInfo& gInfo):gInfo(gInfo) {}
	void test()
	{
		{	puts("Testing memory usage");
			DataGptr gtemp(DataG::alloc(gInfo)); initTranslation(gtemp, vector3<>(1.0,2.0,3.0));
			gtemp *= (1.0/nrm2(gtemp));
			puts("Before:");
			printf("\tnorm(gtemp) = %lf\n", nrm2(gtemp));
			DataGptr gtemp2 = (5.0 * (gtemp * 3.0)) * 7.0;
			puts("After:");
			printf("\tnorm(gtemp) = %lf\n", nrm2(gtemp));
			printf("\tnorm(gtemp2) = %lf\n", nrm2(gtemp2));
			O(O(O(O(gtemp * O(O(O(O(O(O(gtemp)))))*gtemp)))));
		}

		{	puts("\nTest 1: Norm, dot product");
			DataRptr r1(DataR::alloc(gInfo));
			initRandom(r1);
			printf("\tNorm of random vector = %le (expectation value: %le)\n", nrm2(r1), sqrt(gInfo.nr));
			printf("\tSelf dot product = %le (expectation value: %le)\n", dot(r1, r1), double(gInfo.nr));

			puts("\nTest 2: Linear combines");
			DataRptr r2 = 3.0*r1 - r1;
			DataRptr r3 = r1 + r2;
			r2 -= 2.0*((r3 - 0.5*r1*2.0) - r1);
			printf("\tLinear combination norm = %le (should be 0 within roundoff)\n",  nrm2(r2));

			puts("\nTest 3: Transform tests");
			DataGptr g1 = J(r1);
			printf("\tTransform invertibility norm(A - I(J(A))) = %le (should be 0 within roundoff)\n", nrm2(I(g1)-r1));
			printf("\tRepeated to check c2r input integrity: error = %le \n", nrm2(I(g1)-r1));
			printf("\tParseval check: sqrt(N) norm(J(A))- norm(A) = %le\n", sqrt(gInfo.nr) * nrm2(g1) - nrm2(r1));
			printf("\tParseval check: sqrt(N J(A).J(A)) - norm(A) = %le\n", sqrt(gInfo.nr * dot(g1, g1)) - nrm2(r1));

		}

		{	puts("\nTest 4: Poisson solver (from DFT mini course)");
			const double sigma1 = 0.75, sigma2 = 0.5;
			DataRptr gauss1, gauss2;
			{	DataGptr translateCenter(DataG::alloc(gInfo)); RealKernel gaussian(gInfo);
				initTranslation(translateCenter, 0.5*(gInfo.R.column(0)+gInfo.R.column(1)+gInfo.R.column(2)));
				gauss1 = I(gaussConvolve(translateCenter*(1./gInfo.detR), sigma1));
				gauss2 = I(gaussConvolve(translateCenter*(1./gInfo.detR), sigma2));
			}
			DataRptr n = gauss2 - gauss1;
			CoulombParams cp; cp.geometry = CoulombParams::Periodic;
			DataRptr phi = I((*cp.createCoulomb(gInfo))(J(n)));
			//DataRptr phi = I(Linv(-4*M_PI * O(J(n))));
			printf("\tNormalization check on g1: %20.16lf\n", sum(gauss1)*gInfo.detR/gInfo.nr);
			printf("\tNormalization check on g2: %20.16lf\n", sum(gauss2)*gInfo.detR/gInfo.nr);
			printf("\tTotal charge check: %20.16lf\n", sum(n)*gInfo.detR/gInfo.nr);
			double Unum = 0.5*dot(J(phi),O(J(n)));
			double Uanal=((1/sigma1+1/sigma2)/2-sqrt(2)/sqrt(pow(sigma1,2)+pow(sigma2,2)))/sqrt(M_PI);
			printf("\tNumeric, analytic Coulomb energy: %20.16lf,%20.16lf\n", Unum, Uanal);
			saveDX(n, "poissontest_n");
			saveDX(phi, "poissontest_phi");
			DataRptr saveR[2] = {n, phi};
			saveSphericalized(saveR, 2, "poissontest.spherical", 0.25);
		}
	}

	void timeParallel()
	{	puts("\nTiming templated parallelization:");
		DataRptr in(DataR::alloc(gInfo)), out;
		initZero(in); in+=1.5; in*=2;

		TIME("inv", stdout,
			for(int i=0; i<100; i++) out = inv(in);
		)
		printf("Output relative error = %le\n", nrm2(out)/sqrt(out->nElem*pow(1.0/3,2))-1);

		TIME("exp", stdout,
			for(int i=0; i<100; i++) out = exp(in);
		)
		printf("Output relative error = %le\n", nrm2(out)/sqrt(out->nElem*pow(exp(3.0),2))-1);
	}
};

void sync()
{
	#ifdef GPU_ENABLED
	cudaThreadSynchronize();
	#endif
}

void timeEblas3(const GridInfo& gInfo)
{
	//A not-so-untypical wavefunction size:
	int colLength = 65145;
	int nCols = 213;

	//A couple of column bundles:
	Data cb1(gInfo, nCols*colLength, 2, isGpuEnabled());
	Data cb2(gInfo, nCols*colLength, 2, isGpuEnabled());
	//An nBandsxnBands matrix:
	Data mat(gInfo, nCols*nCols, 2, isGpuEnabled());

	sync();

	TIME("Row-major operator^", stdout,
		 for(int i=0; i<10; i++)
			 callPref(eblas_zgemm)(CblasNoTrans, CblasConjTrans, nCols, nCols, colLength,
				1.0, (complex*)cb1.dataPref(), nCols, (complex*)cb2.dataPref(), nCols,
				0.0, (complex*)mat.dataPref(), nCols);
		sync();
	)

	TIME("Column-major operator^", stdout,
		 for(int i=0; i<10; i++)
			 callPref(eblas_zgemm)(CblasConjTrans, CblasNoTrans, nCols, nCols, colLength,
				1.0, (complex*)cb1.dataPref(), colLength, (complex*)cb2.dataPref(), colLength,
				0.0, (complex*)mat.dataPref(), nCols);
		sync();
	)

	TIME("Row-major cb*mat", stdout,
		 for(int i=0; i<10; i++)
			callPref(eblas_zgemm)(CblasNoTrans, CblasNoTrans, nCols, colLength, nCols,
			1.0, (complex*)mat.dataPref(), nCols, (complex*)cb1.dataPref(), nCols,
			0.0, (complex*)cb2.dataPref(), nCols);
		 sync();
	)

	TIME("Column-major cb*mat", stdout,
		 for(int i=0; i<10; i++)
			callPref(eblas_zgemm)(CblasNoTrans, CblasNoTrans, colLength, nCols, nCols,
			1.0, (complex*)cb1.dataPref(), colLength, (complex*)mat.dataPref(), nCols,
			0.0, (complex*)cb2.dataPref(), colLength);
		 sync();
	)

	mpiUtil->exit(0);
}


void testHarmonics()
{	//Spherical bessel functions:
	for(int l=0; l<=6; l++)
	{	double relErrSum=0., relErrMax=0.;
		double absErrSum=0., absErrMax=0.;
		 int errCount=0;
		for(double x=1e-8; x<1e3; x*=1.1)
		{	double jl = bessel_jl(l,x);
			double jlRef = gsl_sf_bessel_jl(l,x);
			double absErr = fabs(jl-jlRef);
			absErrSum += absErr;
			if(absErr>absErrMax) absErrMax=absErr;
			double relErr = fabs(jl-jlRef)/std::max(1e-15,fabs(jlRef));
			relErrSum += relErr;
			if(relErr>relErrMax) relErrMax=relErr;
			errCount++;
		}
		printf("j%d relError: mean=%le max=%le  absError: mean=%le max=%le\n",
			l, relErrSum/errCount, relErrMax, absErrSum/errCount, absErrMax);
	}
}

void testYlmProd()
{	vector3<> qHat(Random::normal(),Random::normal(),Random::normal()); qHat *= 1./qHat.length(); //random test unit vector
	const int lMax = 3;
	for(int l1=0; l1<=lMax; l1++) for(int m1=-l1; m1<=l1; m1++)
	{	double Ylm1 = Ylm(l1, m1, qHat);
		for(int l2=0; l2<=lMax; l2++) for(int m2=-l2; m2<=l2; m2++)
		{	double Ylm2 = Ylm(l2, m2, qHat);
			std::vector<YlmProdTerm> expansion = expandYlmProd(l1,m1, l2,m2);
			double YlmProd = 0.;
			for(const auto& term: expansion) YlmProd += term.coeff * Ylm(term.l, term.m, qHat);
			logPrintf("(%d,%+d)*(%d,%+d): %le\n", l1,m1, l2,m2, YlmProd-Ylm1*Ylm2);
		}
	}
}

void fdtest2D(double F(double,double,double&,double&), const char* name)
{	double A = 1.23856;
	double B = 3.104262;
	double F_A, F_B; F(A, B, F_A, F_B);
	double dumpA, dumpB;
	printf("FD testing %s:\n", name); 
	for(double d=1e-1; d>1e-12; d*=0.1)
	{	double FnumA = (0.5/d)*(F(A+d,B,dumpA,dumpB) - F(A-d,B,dumpA,dumpB));
		double FnumB = (0.5/d)*(F(A,B+d,dumpA,dumpB) - F(A,B-d,dumpA,dumpB));
		printf("\td: %le  RelErrA: %le  RelErrB: %le\n", d, FnumA/F_A-1, FnumB/F_B-1);
	}
}
#define FDtest2D(x) fdtest2D(x, #x)

void fdtestGGAs()
{
// 	FDtest2D(integralErfcGaussian<1>);
// 	FDtest2D(integralErfcGaussian<2>);
// 	FDtest2D(integralErfcGaussian<3>);
// 	FDtest2D(integralErfcGaussian<5>);
	FDtest2D(GGA_eval<GGA_X_wPBE_SR>);
}

void timePointGroupOps(const GridInfo& gInfo)
{
	std::vector<matrix3<int>> sym = getSymmetries(gInfo.R);
	matrix3<int> m = sym[3]; //some random symmetry matrix
	m.print(globalLog, " %2d ");
	matrix3<int> mInv = det(m) * adjugate(m); //since |det(m)| = 0
	
	DataRptr x(DataR::alloc(gInfo)); initRandom(x);
	logPrintf("Rel. error in implementations = %le\n",
		nrm2(pointGroupGather(x, m) - pointGroupScatter(x, mInv)) / nrm2(x));
	
	int nRepetitions = int(1e8/gInfo.nr);
	logPrintf("Timing %d repetitions.\n", nRepetitions);
	DataRptr y;
	
	TIME("Gather", globalLog,
		for(int iRep=0; iRep<nRepetitions; iRep++)
			 y = pointGroupGather(x, m);
	)
	TIME("Scatter", globalLog,
		for(int iRep=0; iRep<nRepetitions; iRep++)
			 y = pointGroupScatter(x, mInv);
	)
}

void print(const DataGptr& x)
{	complex* xData = x->data();
	const GridInfo& g = x->gInfo;
	vector3<int> iv;
	for(iv[0]=0; iv[0]<g.S[0]; iv[0]++)
	{	for(iv[1]=0; iv[1]<g.S[1]; iv[1]++)
		{	for(iv[2]=0; iv[2]<1+g.S[2]/2; iv[2]++)
				logPrintf("\t%2lg", xData[g.halfGindex(g.wrapGcoords(iv))].real());
			logPrintf("\n");
		}
		logPrintf("\n");
	}
}
void testChangeGrid()
{	GridInfo g1, g2;
	g1.R = g2.R = matrix3<>(1,1,1);
	g1.S = vector3<int>(3,3,3);
	g2.S = vector3<int>(4,4,4);
	g1.initialize();
	g2.initialize();
	DataGptr x1, x2, x3;
	nullToZero(x1, g1);
	complex* x1data = x1->data();
	vector3<int> iv; int i=0;
	for(iv[0]=0; iv[0]<g1.S[0]; iv[0]++)
	for(iv[1]=0; iv[1]<g1.S[1]; iv[1]++)
	for(iv[2]=0; iv[2]<1+g1.S[2]/2; iv[2]++)
		x1data[g1.halfGindex(g1.wrapGcoords(iv))] = i++;
	
	x2 = J(changeGrid(I(x1,true), g2));
	x3 = J(changeGrid(I(x2,true), g1));
	logPrintf("\n--------- x1 --------\n"); print(x1);
	logPrintf("\n--------- x2 --------\n"); print(x2);
	logPrintf("\n--------- x3 --------\n"); print(x3);
}

int main(int argc, char** argv)
{	initSystem(argc, argv);
	testChangeGrid(); return 0;
	//testYlmProd(); return 0;
	//testHarmonics(); return 0;
	//fdtestGGAs(); return 0;

// 	const int Zn = 2;
// 	SO3quad quad(QuadEuler, Zn, 11); //quad.print();
// 	SO3quad quad2(Quad10design_60, Zn); // quad.print();
// 	return 0;
	
	GridInfo gInfo;
	gInfo.S = vector3<int>(128, 128, 128);
	gInfo.R.set_col(0, vector3<>(0.0, 12.0, 12.0));
	gInfo.R.set_col(1, vector3<>(12.0, 0.0, 12.0));
	gInfo.R.set_col(2, vector3<>(12.0, 12.0, 0.0));
	gInfo.initialize();
	//testCavitation(gInfo); return 0;
	timePointGroupOps(gInfo); return 0;
	//timeEblas3(gInfo);
	
	OperatorTest op(gInfo);
	op.test();
	op.timeParallel();
	
	finalizeSystem();
	return 0;
}
