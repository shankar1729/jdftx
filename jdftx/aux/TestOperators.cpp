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

#include <electronic/SpeciesInfo.h>
#include <electronic/ColumnBundle.h>
#include <electronic/ExCorr_internal_GGA.h>
#include <core/Util.h>
#include <core/Operators.h>
#include <core/ScalarFieldIO.h>
#include <core/Coulomb.h>
#include <core/LatticeUtils.h>
#include <core/Random.h>
#include <core/SphericalHarmonics.h>
#include <core/Blip.h>
#include <fluid/SO3quad.h>
#include <gsl/gsl_sf.h>
#include <stdlib.h>

class OperatorTest
{	const GridInfo& gInfo;

public:
	OperatorTest(const GridInfo& gInfo):gInfo(gInfo) {}
	void test()
	{
		{	logPrintf("Testing memory usage\n");
			ScalarFieldTilde gtemp(ScalarFieldTildeData::alloc(gInfo)); initTranslation(gtemp, vector3<>(1.0,2.0,3.0));
			gtemp *= (1.0/nrm2(gtemp));
			logPrintf("Before:\n");
			logPrintf("\tnorm(gtemp) = %lf\n", nrm2(gtemp));
			ScalarFieldTilde gtemp2 = (5.0 * (gtemp * 3.0)) * 7.0;
			logPrintf("After:\n");
			logPrintf("\tnorm(gtemp) = %lf\n", nrm2(gtemp));
			logPrintf("\tnorm(gtemp2) = %lf\n", nrm2(gtemp2));
			O(O(O(O(gtemp * O(O(O(O(O(O(gtemp)))))*gtemp)))));
		}

		{	logPrintf("\nTest 1: Norm, dot product\n");
			ScalarField r1(ScalarFieldData::alloc(gInfo));
			initRandom(r1);
			logPrintf("\tNorm of random vector = %le (expectation value: %le)\n", nrm2(r1), sqrt(gInfo.nr));
			logPrintf("\tSelf dot product = %le (expectation value: %le)\n", dot(r1, r1), double(gInfo.nr));

			logPrintf("\nTest 2: Linear combines\n");
			ScalarField r2 = 3.0*r1 - r1;
			ScalarField r3 = r1 + r2;
			r2 -= 2.0*((r3 - 0.5*r1*2.0) - r1);
			logPrintf("\tLinear combination norm = %le (should be 0 within roundoff)\n",  nrm2(r2));

			logPrintf("\nTest 3: Transform tests\n");
			ScalarFieldTilde g1 = J(r1);
			logPrintf("\tTransform invertibility norm(A - I(J(A))) = %le (should be 0 within roundoff)\n", nrm2(I(g1)-r1));
			logPrintf("\tRepeated to check c2r input integrity: error = %le \n", nrm2(I(g1)-r1));
			logPrintf("\tParseval check: sqrt(N) norm(J(A))- norm(A) = %le\n", sqrt(gInfo.nr) * nrm2(g1) - nrm2(r1));
			logPrintf("\tParseval check: sqrt(N J(A).J(A)) - norm(A) = %le\n", sqrt(gInfo.nr * dot(g1, g1)) - nrm2(r1));

		}

		{	logPrintf("\nTest 4: Poisson solver (from DFT mini course)\n");
			const double sigma1 = 0.75, sigma2 = 0.5;
			ScalarField gauss1, gauss2;
			{	ScalarFieldTilde translateCenter(ScalarFieldTildeData::alloc(gInfo)); RealKernel gaussian(gInfo);
				initTranslation(translateCenter, 0.5*(gInfo.R.column(0)+gInfo.R.column(1)+gInfo.R.column(2)));
				gauss1 = I(gaussConvolve(translateCenter*(1./gInfo.detR), sigma1));
				gauss2 = I(gaussConvolve(translateCenter*(1./gInfo.detR), sigma2));
			}
			ScalarField n = gauss2 - gauss1;
			CoulombParams cp; cp.geometry = CoulombParams::Periodic;
			ScalarField phi = I((*cp.createCoulomb(gInfo))(J(n)));
			//ScalarField phi = I(Linv(-4*M_PI * O(J(n))));
			logPrintf("\tNormalization check on g1: %20.16lf\n", sum(gauss1)*gInfo.detR/gInfo.nr);
			logPrintf("\tNormalization check on g2: %20.16lf\n", sum(gauss2)*gInfo.detR/gInfo.nr);
			logPrintf("\tTotal charge check: %20.16lf\n", sum(n)*gInfo.detR/gInfo.nr);
			double Unum = 0.5*dot(J(phi),O(J(n)));
			double Uanal=((1/sigma1+1/sigma2)/2-sqrt(2)/sqrt(pow(sigma1,2)+pow(sigma2,2)))/sqrt(M_PI);
			logPrintf("\tNumeric, analytic Coulomb energy: %20.16lf,%20.16lf\n", Unum, Uanal);
			saveDX(n, "poissontest_n");
			saveDX(phi, "poissontest_phi");
			ScalarField saveR[2] = {n, phi};
			saveSphericalized(saveR, 2, "poissontest.spherical", 0.25);
		}
		
		{	logPrintf("\nTest 5: Asymmetric Fourier transform test\n");
			ScalarField r1(ScalarFieldData::alloc(gInfo)); initRandom(r1);
			ScalarField r2(ScalarFieldData::alloc(gInfo)); initRandom(r2);
			logPrintf("\tAnti-hermiticity of D: %le\n", fabs(dot(r1,gradient(r2)[0])/dot(r2,gradient(r1)[0]) + 1.));
		}
	}

	void timeParallel()
	{	logPrintf("\nTiming templated parallelization:\n");
		ScalarField in(ScalarFieldData::alloc(gInfo)), out;
		initZero(in); in+=1.5; in*=2;

		TIME("inv", globalLog,
			for(int i=0; i<100; i++) out = inv(in);
		)
		logPrintf("Output relative error = %le\n", nrm2(out)/sqrt(out->nElem*pow(1.0/3,2))-1);

		TIME("exp", globalLog,
			for(int i=0; i<100; i++) out = exp(in);
		)
		logPrintf("Output relative error = %le\n", nrm2(out)/sqrt(out->nElem*pow(exp(3.0),2))-1);
	}
};

void sync()
{
	#ifdef GPU_ENABLED
	cudaDeviceSynchronize();
	#endif
}

void timeEblas3(const GridInfo& gInfo)
{
	//A not-so-untypical wavefunction size:
	int colLength = 65145;
	int nCols = 213;

	//A couple of column bundles:
	ColumnBundle cb1(nCols, colLength, 0, 0, isGpuEnabled());
	ColumnBundle cb2(nCols, colLength, 0, 0, isGpuEnabled());
	//An nBandsxnBands matrix:
	matrix mat(nCols, nCols, isGpuEnabled());

	sync();

	TIME("Row-major operator^", globalLog,
		 for(int i=0; i<10; i++)
			 callPref(eblas_zgemm)(CblasNoTrans, CblasConjTrans, nCols, nCols, colLength,
				1.0, (complex*)cb1.dataPref(), nCols, (complex*)cb2.dataPref(), nCols,
				0.0, (complex*)mat.dataPref(), nCols);
		sync();
	)

	TIME("Column-major operator^", globalLog,
		 for(int i=0; i<10; i++)
			 callPref(eblas_zgemm)(CblasConjTrans, CblasNoTrans, nCols, nCols, colLength,
				1.0, (complex*)cb1.dataPref(), colLength, (complex*)cb2.dataPref(), colLength,
				0.0, (complex*)mat.dataPref(), nCols);
		sync();
	)

	TIME("Row-major cb*mat", globalLog,
		 for(int i=0; i<10; i++)
			callPref(eblas_zgemm)(CblasNoTrans, CblasNoTrans, nCols, colLength, nCols,
			1.0, (complex*)mat.dataPref(), nCols, (complex*)cb1.dataPref(), nCols,
			0.0, (complex*)cb2.dataPref(), nCols);
		 sync();
	)

	TIME("Column-major cb*mat", globalLog,
		 for(int i=0; i<10; i++)
			callPref(eblas_zgemm)(CblasNoTrans, CblasNoTrans, colLength, nCols, nCols,
			1.0, (complex*)cb1.dataPref(), colLength, (complex*)mat.dataPref(), nCols,
			0.0, (complex*)cb2.dataPref(), colLength);
		 sync();
	)

	mpiWorld->exit(0);
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
		logPrintf("j%d relError: mean=%le max=%le  absError: mean=%le max=%le\n",
			l, relErrSum/errCount, relErrMax, absErrSum/errCount, absErrMax);
		
		//Test the Y transformation matrix inside SpeciesInfo::getYlmToSpinAngleMatrix (need to modify that code to return the intermediate Y)
// 		matrix Y = SpeciesInfo::getYlmToSpinAngleMatrix(l, 2*l+1);
// 		int nTest = 20;
// 		matrix Yreal(nTest, 2*l+1), Ycomp(nTest, 2*l+1);
// 		for(int iTest=0; iTest<nTest; iTest++)
// 		{	double theta = Random::uniform(0, M_PI); double ct, st; sincos(theta, &st, &ct);
// 			double phi = Random::uniform(0, 2*M_PI); double cp, sp; sincos(phi, &sp, &cp);
// 			vector3<> qhat(st*cp, st*sp, ct);
// 			for(int m=-l; m<=l; m++)
// 			{	Yreal.data()[Yreal.index(iTest,l+m)] = Ylm(l, m, qhat);
// 				Ycomp.data()[Ycomp.index(iTest,l+m)] = gsl_sf_legendre_sphPlm(l, abs(m), ct) * cis(m*phi) * (m<0 ? pow(-1,m) : 1);
// 			}
// 		}
// 		logPrintf("Y%d relError: %le\n", l, nrm2(Ycomp - Yreal*Y)/nrm2(Ycomp));
// 		logPrintf("\n"); Ycomp.print(globalLog, "%+9.5f%+9.5fi\n");
// 		logPrintf("\n"); (Yreal*Y).print(globalLog, "%+9.5f%+9.5fi\n");
		
		//Test SpeciesInfo::getYlmToSpinAngleMatrix
		for(int j2=2*l-1; j2<=2*l+1; j2+=2) if(j2>0)
		{	logPrintf("\n\n------------ l=%d  j=%lg ----------\n\n", l, 0.5*j2);
			matrix flj = SpeciesInfo::getYlmOverlapMatrix(l, j2);
			flj.print(globalLog, "%+9.5f%+9.5fi ");
		}
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
			double err = YlmProd-Ylm1*Ylm2;
			logPrintf("(%d,%+d)*(%d,%+d): %le%s\n", l1,m1, l2,m2, err, (fabs(err)>1e-16 ? " ERR" : ""));
		}
	}
}

/*
double YlmPrime(int l, int m, int iDir, const vector3<>& qHat)
{	const double alpha = sqrt((2*l+1.)/(2*l-1));
	if(m == 0)
	{	switch(iDir)
		{	case 0: return alpha * (l>1 ? -sqrt(0.5*l*(l-1))*Ylm(l-1,+1,qHat) : 0.);
			case 1: return alpha * (l>1 ? -sqrt(0.5*l*(l-1))*Ylm(l-1,-1,qHat) : 0.);
			case 2: return alpha * (l>0 ? l*Ylm(l-1,0,qHat) : 0.);
		}
	}
	else if(m > 0)
	{	switch(iDir)
		{	case 0: return alpha*0.5 *
				( (l>m+1 ? -sqrt((l-m)*(l-m-1))*Ylm(l-1,m+1,qHat) : 0.)
				+ (l>=m  ?  sqrt((l+m)*(l+m-1))*Ylm(l-1,m-1,qHat) * (m==1 ? sqrt(2.) : 1.) : 0.) );
			case 1: return alpha*0.5 *
				( (l>m+1 ? -sqrt((l-m)*(l-m-1))*Ylm(l-1,-(m+1),qHat) : 0.)
				+ (l>=m  ? -sqrt((l+m)*(l+m-1))*Ylm(l-1,-(m-1),qHat) * (m==1 ? 0. : 1.) : 0.) );
			case 2: return alpha * (l>0 ? sqrt(l*l-m*m)*Ylm(l-1,m,qHat) : 0.);
		}
	}
	else //m < 0
	{	switch(iDir)
		{	case 0: return alpha*0.5 *
				( (-l<=m  ?  sqrt((l-m)*(l-m-1))*Ylm(l-1,m+1,qHat) * (m==-1 ? 0. : 1.) : 0.)
				+ (-l<m-1 ? -sqrt((l+m)*(l+m-1))*Ylm(l-1,m-1,qHat) : 0.) );
			case 1: return alpha*0.5 *
				( (-l<=m  ?  sqrt((l-m)*(l-m-1))*Ylm(l-1,-(m+1),qHat) * (m==-1 ? sqrt(2.) : 1.) : 0.)
				+ (-l<m-1 ?  sqrt((l+m)*(l+m-1))*Ylm(l-1,-(m-1),qHat) : 0.) );
			case 2: return alpha * (l>0 ? sqrt(l*l-m*m)*Ylm(l-1,m,qHat) : 0.);
		}
	}
}

double YlmPrime(int l, int m, int iDir, const vector3<>& qHat)
{	const double alpha = sqrt((2*l+1.)/(2*l-1));
	if(iDir == 2)
		return (l>0 ? (alpha*sqrt(l*l-m*m)) * Ylm(l-1,m,qHat) : 0.);
	else
	{	int dSign = (iDir==0 ? +1 : -1); //+1 for x, -1 for y
		if(m == 0)
			return (l>1 ? (-alpha * sqrt(0.5*l*(l-1))) * Ylm(l-1,dSign,qHat) : 0.);
		else
		{	double result = 0.;
			int mSign = (m<0 ? -1 : +1);
			for(int s=-1; s<=+1; s+=2)
			{	int mOut = m - s;
				if(abs(mOut) <= l-1)
					result += ((0.5*alpha) * (iDir==0 ? mSign*s : -mSign)
						* sqrt((l+m*s)*(l+m*s-1)*(mOut ? 1 : 1+dSign*mSign)))
						* Ylm(l-1, dSign*mOut, qHat);
			}
			return result;
		}
	}
}
*/

void testYlmDeriv()
{	vector3<> qHat(Random::normal(),Random::normal(),Random::normal()); qHat *= 1./qHat.length(); //random test unit vector
	const int lMax = 6;
	const double dq = 1e-6;
	for(int l=0; l<=lMax; l++)
		for(int m=-l; m<=l; m++)
		{	vector3<> Yprime = YlmPrime(l, m, qHat); //analyticla derivative
			for(int iDir=0; iDir<3; iDir++)
			{	//Numerical derivative:
				vector3<> qHatP = qHat; qHatP[iDir] += dq;
				vector3<> qHatM = qHat; qHatM[iDir] -= dq;
				double YlmPrimeNum = (0.5/dq) * (Ylm(l,m,qHatP) - Ylm(l,m,qHatM));
				double err = Yprime[iDir] - YlmPrimeNum;
				logPrintf("Y(%d,%+d)_%c: %le%s\n", l,m, "xyz"[iDir], err, (fabs(err)>1e-10 ? " ERR" : ""));
			}
		}
}

void testVnlPrime()
{
	//Forward declaration from SpeciesInfo_internal.h:
	void Vnl(int nbasis, int atomStride, int nAtoms, int l, int m, const vector3<> k, const vector3<int>* iGarr,
		const matrix3<> G, const vector3<>* pos, const RadialFunctionG& VnlRadial, complex* Vnl,
		const vector3<>* derivDir=0, const int stressDir=-1);
	
	vector3<> k(0.25, 0.1, 0.3);
	vector3<int> iG(1, -1, 0);
	matrix3<> R(3.,1.,0., -2.,4.,1., 0.,2.,5.);
	matrix3<> G = 2*M_PI*inv(R);
	vector3<> pos(-0.2, 0.3, 0.1);
	vector3<> derivDir(0.7, 0.2, -0.5);
	nProcsAvailable = 1;
	RadialFunctionG Vradial;
	logPrintf("|kpG|: %le\n", ((k+iG)*G).length());
	for(int l=0; l<=3; l++)
	{	Vradial.init(l, 0.02, 20., RadialFunctionG::cusplessExpTilde, -8., 0.1);
		for(int m=-l; m<=l; m++)
		{	//Numerical derivative:
			complex Vp, Vm; double alpha = 1e-4; vector3<> dkLat = alpha*derivDir*inv(G);
			Vnl(1, 1, 1, l, m, k+dkLat, &iG, G, &pos, Vradial, &Vp);
			Vnl(1, 1, 1, l, m, k-dkLat, &iG, G, &pos, Vradial, &Vm);
			complex VprimeNum = (0.5/alpha)*(Vp - Vm);
			//Analytic derivative:
			complex Vprime;
			Vnl(1, 1, 1, l, m, k, &iG, G, &pos, Vradial, &Vprime, &derivDir);
			double err = (Vprime - VprimeNum).abs() / Vprime.abs();
			logPrintf("Vprime(%d,%+d): %le%s\n", l,m, err, (fabs(err)>1e-8 ? " ERR" : ""));
		}
		Vradial.free();
	}
}

void fdtest2D(double F(double,double,double&,double&), const char* name)
{	double A = 1.23856;
	double B = 3.104262;
	double F_A, F_B; F(A, B, F_A, F_B);
	double dumpA, dumpB;
	logPrintf("FD testing %s:\n", name); 
	for(double d=1e-1; d>1e-12; d*=0.1)
	{	double FnumA = (0.5/d)*(F(A+d,B,dumpA,dumpB) - F(A-d,B,dumpA,dumpB));
		double FnumB = (0.5/d)*(F(A,B+d,dumpA,dumpB) - F(A,B-d,dumpA,dumpB));
		logPrintf("\td: %le  RelErrA: %le  RelErrB: %le\n", d, FnumA/F_A-1, FnumB/F_B-1);
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

void print(const ScalarFieldTilde& x)
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
	ScalarFieldTilde x1, x2, x3;
	nullToZero(x1, g1);
	complex* x1data = x1->data();
	vector3<int> iv; int i=0;
	for(iv[0]=0; iv[0]<g1.S[0]; iv[0]++)
	for(iv[1]=0; iv[1]<g1.S[1]; iv[1]++)
	for(iv[2]=0; iv[2]<1+g1.S[2]/2; iv[2]++)
		x1data[g1.halfGindex(g1.wrapGcoords(iv))] = i++;
	
	x2 = J(changeGrid(I(x1), g2));
	x3 = J(changeGrid(I(x2), g1));
	logPrintf("\n--------- x1 --------\n"); print(x1);
	logPrintf("\n--------- x2 --------\n"); print(x2);
	logPrintf("\n--------- x3 --------\n"); print(x3);
}

void testHugeFileIO()
{	matrix M(15000,15000);
	logPrintf("Testing huge file I/O with %lg GB.\n", pow(0.5,30)*(M.nData()*sizeof(complex)));
	complex* Mdata = M.data();
	for(size_t i=0; i<M.nData(); i++) Mdata[i] = i;
	const char* fname = "testHugeFileIO.dat";
	//Write:
	logPrintf("Writing ... "); logFlush();
	MPIUtil::File fp; mpiWorld->fopenWrite(fp, fname);
	mpiWorld->fwriteData(M, fp);
	mpiWorld->fclose(fp);
	logPrintf("done.\n"); logFlush();
	//Read back:
	logPrintf("Reading ... "); logFlush();
	M.zero();
	mpiWorld->fopenRead(fp, fname, M.nData()*sizeof(complex));
	mpiWorld->freadData(M, fp);
	mpiWorld->fclose(fp);
	logPrintf("done.\n"); logFlush();
	//Check results:
	double rmsErr = 0.;
	for(size_t i=0; i<M.nData(); i++) rmsErr += (Mdata[i]-i).norm();
	rmsErr = sqrt(rmsErr/M.nData());
	logPrintf("rmsErr = %le\n", rmsErr);
}

void testResample()
{	GridInfo gInfo1;
	gInfo1.S = vector3<int>(1,1,2) * 32;
	gInfo1.R.set_col(0, vector3<>(0,5,5));
	gInfo1.R.set_col(1, vector3<>(5,0,5));
	gInfo1.R.set_col(2, vector3<>(9,9,0));
	gInfo1.initialize();
	GridInfo gInfo2;
	gInfo2.R = matrix3<>(1,1,1) * 20;
	gInfo2.S = vector3<int>(1,1,1) * 128;
	gInfo2.initialize();
	
	//Create a gaussian on grid1:
	ScalarFieldTilde delta(ScalarFieldTildeData::alloc(gInfo1));
	initTranslation(delta, vector3<>(0,0,0));
	ScalarFieldTilde gauss1 = gaussConvolve(delta*(1./gInfo1.detR), 0.5);
	saveRawBinary(I(gauss1), "testResample.n1");
	
	//Convert to other grid:
	BlipResampler resample(gInfo1, gInfo2);
	ScalarField gauss2 = resample(gauss1);
	saveRawBinary(gauss2, "testResample.n2");
}

void testMatrixLinalg()
{	//Make a hermitian pos-def matrix:
	int N = 128;
	matrix A(N, N);
	randomize(A);
	A = A * dagger(A);
	//Make a random rhs matrix:
	int Nrhs = 64;
	matrix b(N, Nrhs);
	randomize(b);
	//Test invApply:
	matrix x = invApply(A, b);
	logPrintf("Relative error in invApply = %le\n", nrm2(A * x - b)/nrm2(b));
	std::vector<bool> upperArr = {{false, true}};
	for(bool upper: upperArr)
	{	const char* upperStr = upper ? "upper" : "lower";
		//Test Cholesky factorization:
		matrix L = cholesky(A, upper);
		matrix LL = upper ? (dagger(L) * L) : (L * dagger(L)); //appropriate product
		logPrintf("Relative error in cholesky(%s) = %le\n", upperStr, nrm2(LL - A)/nrm2(A));
		//Test triangular matrix inversion:
		matrix invL = invTriangular(L, upper);
		logPrintf("Relative error in invTriangular(%s) = %le\n", upperStr, nrm2(L * invL - eye(N))/sqrt(N));
	}
	//Test orthoMatrix:
	matrix U = orthoMatrix(A);
	logPrintf("Relative error in orthoMatrix = %le\n", nrm2(dagger(U) * A * U - eye(N)));
}

int main(int argc, char** argv)
{	initSystem(argc, argv);
	//testHarmonics(); return 0;
	//testYlmProd(); return 0;
	//testYlmDeriv(); return 0;
	//testVnlPrime(); return 0;
	//fdtestGGAs(); return 0;
	//testChangeGrid(); return 0;
	//testHugeFileIO(); return 0;
	//testResample(); return 0;
	testMatrixLinalg(); return 0;
	
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
	//timeEblas3(gInfo);
	
	OperatorTest op(gInfo);
	op.test();
	op.timeParallel();
	
	finalizeSystem();
	return 0;
}
