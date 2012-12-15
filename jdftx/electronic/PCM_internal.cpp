/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman, Deniz Gunceler

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

#include <electronic/PCM_internal.h>
#include <core/Operators.h>
#include <core/DataMultiplet.h>

//----------------------- The JDFT `shape function' and gradient ------------------

void pcmShapeFunc(int N, const double* nCavity, double* shape, const double nc, const double sigma)
{	threadedLoop(pcmShapeFunc_calc, N, nCavity, shape, nc, sigma);
}
void pcmShapeFunc_grad(int N, const double* nCavity, const double* grad_shape, double* grad_nCavity, const double nc, const double sigma)
{	threadedLoop(pcmShapeFunc_grad_calc, N, nCavity, grad_shape, grad_nCavity, nc, sigma);
}
#ifdef GPU_ENABLED
void pcmShapeFunc_gpu(int N, const double* nCavity, double* shape, const double nc, const double sigma);
void pcmShapeFunc_grad_gpu(int N, const double* nCavity, const double* grad_shape, double* grad_nCavity, const double nc, const double sigma);
#endif
void pcmShapeFunc(const DataRptr& nCavity, DataRptr& shape, const double nc, const double sigma)
{	nullToZero(shape, nCavity->gInfo);
	callPref(pcmShapeFunc)(nCavity->gInfo.nr, nCavity->dataPref(), shape->dataPref(), nc, sigma);
}
void pcmShapeFunc_grad(const DataRptr& nCavity, const DataRptr& grad_shape, DataRptr& grad_nCavity, const double nc, const double sigma)
{	nullToZero(grad_nCavity, nCavity->gInfo);
	callPref(pcmShapeFunc_grad)(nCavity->gInfo.nr, nCavity->dataPref(), grad_shape->dataPref(), grad_nCavity->dataPref(), nc, sigma);
}

double cavitationEnergyAndGrad(const DataRptr& shape, DataRptr& Acavity_shape, double cavityTension, double cavityPressure)
{
	DataRptrVec shape_x = gradient(shape);
	DataRptr surfaceDensity = sqrt(shape_x[0]*shape_x[0] + shape_x[1]*shape_x[1] + shape_x[2]*shape_x[2]);
	double surfaceArea = integral(surfaceDensity);
	double volume = integral(1.-shape);

	DataRptr invSurfaceDensity = inv(surfaceDensity);
	Acavity_shape = (-cavityTension)*divergence(shape_x*invSurfaceDensity); // Surface term
	Acavity_shape += (-cavityPressure); // Volume term
	
	return surfaceArea*cavityTension + volume*cavityPressure;
}

//------------- Helper classes for NonlinearPCM  -------------
namespace NonlinearPCMeval
{
	Screening::Screening(bool linear, double T, double Nion, double Zion)
	: linear(linear), NT2(Nion*T*2.), ZbyT(Zion/T), m2NZ(-2.*Nion*Zion)
	{
	}
	
	void ScreeningFreeEnergy_sub(size_t iStart, size_t iStop, const double* phi, const double* s, double* rho, double* A, double* A_phi, double* A_s, const Screening& eval)
	{	for(size_t i=iStart; i<iStop; i++) eval.freeEnergy_calc(i, phi, s, rho, A, A_phi, A_s);
	}
	void Screening::freeEnergy(size_t N, const double* phi, const double* s, double* rho, double* A, double* A_phi, double* A_s) const
	{	threadLaunch(ScreeningFreeEnergy_sub, N, phi, s, rho, A, A_phi, A_s, *this);
	}
	
	void ScreeningConvertDerivative_sub(size_t iStart, size_t iStop, const double* phi, const double* s, const double* A_rho, double* A_phi, double* A_s, const Screening& eval)
	{	for(size_t i=iStart; i<iStop; i++) eval.convertDerivative_calc(i, phi, s, A_rho, A_phi, A_s);
	}
	void Screening::convertDerivative(size_t N, const double* phi, const double* s, const double* A_rho, double* A_phi, double* A_s) const
	{	threadLaunch(ScreeningConvertDerivative_sub, N, phi, s, A_rho, A_phi, A_s, *this);
	}
	
	
	Dielectric::Dielectric(bool linear, double T, double Nmol, double pMol, double epsBulk, double epsInf)
	: linear(linear), Np(Nmol * pMol), pByT(pMol / T), NT(Nmol * T),
		Nchi((epsInf-1.)/(4*M_PI)), alpha(3 - 4*M_PI*Np*pByT/(epsBulk-epsInf))
	{	//Check parameter validity:
		if(alpha < 0.)
			die("\nCurrent Nonlinear PCM parameters imply negative correlation for dipole rotations.\n"
				"\tHINT: Reduce solvent molecule dipole or epsInf to fix this");
		
		logPrintf("alpha=%lf  Np=%le  p/T=%le  NT=%le  Nchi=%lf   4pi(Npp/T(3-alpha)+Nchi)=%lf\n",
			alpha, Np, pByT, NT, Nchi, 4*M_PI*(Np*pByT/(3.-alpha)+Nchi));
		
		double x = 2.9869354;
		double F, F_x, Chi, Chi_x;
		compute(x, F, F_x, Chi, Chi_x);
		printf("FD testing Dielectric::Compute:\n"); 
		for(double d=1e-1; d>1e-12; d*=0.1)
		{	//+
			double Fplus, Fplus_x, ChiPlus, ChiPlus_x;
			compute(x+d, Fplus, Fplus_x, ChiPlus, ChiPlus_x);
			//-
			double Fminus, Fminus_x, ChiMinus, ChiMinus_x;
			compute(x-d, Fminus, Fminus_x, ChiMinus, ChiMinus_x);
			//report
			double F_xNum = (0.5/d)*(Fplus-Fminus);
			double Chi_xNum = (0.5/d)*(ChiPlus-ChiMinus);
			printf("\td: %le  RelErr(F): %le  RelErr(Chi): %le\n", d, F_xNum/F_x-1, Chi_xNum/Chi_x-1);
		}
	}
	
	void Dielectric::initLookupTables()
	{
	}
	void Dielectric::freeLookupTables()
	{
	}
	
	void DielectricFreeEnergy_sub(size_t iStart, size_t iStop, vector3<const double*> gradPhi, const double* s, vector3<double*> p, double* A, vector3<double*> A_gradPhi, double* A_s, const Dielectric& eval)
	{	for(size_t i=iStart; i<iStop; i++) eval.freeEnergy_calc(i, gradPhi, s, p, A, A_gradPhi, A_s);
	}
	void Dielectric::freeEnergy(size_t N,
		vector3<const double*> gradPhi, const double* s, vector3<double*> p, double* A, vector3<double*> A_gradPhi, double* A_s) const
	{	threadLaunch(DielectricFreeEnergy_sub, N, gradPhi, s, p, A, A_gradPhi, A_s, *this);
	}
	
	void DielectricConvertDerivative_sub(size_t iStart, size_t iStop, vector3<const double*> gradPhi, const double* s, vector3<const double*> A_p, vector3<double*> A_gradPhi, double* A_s, const Dielectric& eval)
	{	for(size_t i=iStart; i<iStop; i++) eval.convertDerivative_calc(i, gradPhi, s, A_p, A_gradPhi, A_s);
	}
	void Dielectric::convertDerivative(size_t N, vector3<const double*> gradPhi, const double* s, vector3<const double*> A_p, vector3<double*> A_gradPhi, double* A_s) const
	{	threadLaunch(DielectricConvertDerivative_sub, N, gradPhi, s, A_p, A_gradPhi, A_s, *this);
	}
}
