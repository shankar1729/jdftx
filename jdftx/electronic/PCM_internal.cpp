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
#include <core/Units.h>

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
	Screening::Screening(bool linear, double T, double Nion, double Zion, double Rplus, double Rminus, double epsBulk)
	: linear(linear), NT(Nion*T), NZ(Nion*Zion),
	x0plus(Nion * (4.*M_PI/3)*pow(Rplus,3)),
	x0minus(Nion * (4.*M_PI/3)*pow(Rminus,3)),
	x0(x0plus + x0minus)
	{
		if(x0 >= 1.) die("Bulk ionic concentration exceeds hard sphere limit = %lg mol/liter.\n", (Nion/x0) / (mol/liter));
		
		double screenLength = sqrt(T*epsBulk/(8*M_PI*Nion*Zion*Zion));
		if(linear) logPrintf("   Linear ions with screening length = %lg bohrs.\n", screenLength);
		else logPrintf("   Nonlinear ions with screening length = %lg bohrs and Z = %lg at T = %lg K.\n", screenLength, Zion, T/Kelvin);
	}
	
	void ScreeningFreeEnergy_sub(size_t iStart, size_t iStop, double mu0, const double* muPlus, const double* muMinus, const double* s, double* rho, double* A, double* A_muPlus, double* A_muMinus, double* A_s, const Screening& eval)
	{	for(size_t i=iStart; i<iStop; i++) eval.freeEnergy_calc(i, mu0, muPlus, muMinus, s, rho, A, A_muPlus, A_muMinus, A_s);
	}
	void Screening::freeEnergy(size_t N, double mu0, const double* muPlus, const double* muMinus, const double* s, double* rho, double* A, double* A_muPlus, double* A_muMinus, double* A_s) const
	{	threadLaunch(ScreeningFreeEnergy_sub, N, mu0, muPlus, muMinus, s, rho, A, A_muPlus, A_muMinus, A_s, *this);
	}
	
	void ScreeningConvertDerivative_sub(size_t iStart, size_t iStop, double mu0, const double* muPlus, const double* muMinus, const double* s, const double* A_rho, double* A_muPlus, double* A_muMinus, double* A_s, const Screening& eval)
	{	for(size_t i=iStart; i<iStop; i++) eval.convertDerivative_calc(i, mu0, muPlus, muMinus, s, A_rho, A_muPlus, A_muMinus, A_s);
	}
	void Screening::convertDerivative(size_t N, double mu0, const double* muPlus, const double* muMinus, const double* s, const double* A_rho, double* A_muPlus, double* A_muMinus, double* A_s) const
	{	threadLaunch(ScreeningConvertDerivative_sub, N, mu0, muPlus, muMinus, s, A_rho, A_muPlus, A_muMinus, A_s, *this);
	}
	
	
	Dielectric::Dielectric(bool linear, double T, double Nmol, double pMol, double epsBulk, double epsInf)
	: linear(linear), Np(Nmol * pMol), NT(Nmol * T),
		alpha(3 - 4*M_PI*Np*pMol/(T*(epsBulk-epsInf))),
		X((epsInf-1.)*T/(4*M_PI*Np*pMol))
	{	//Check parameter validity:
		if(alpha < 0.)
			die("\nCurrent Nonlinear PCM parameters imply negative correlation for dipole rotations.\n"
				"\tHINT: Reduce solvent molecule dipole or epsInf to fix this.\n");
		if(linear) logPrintf("   Linear dielectric with epsBulk = %lg.\n", epsBulk);
		else logPrintf("   Nonlinear dielectric with epsBulk = %lg and epsInf = %lg at T = %lg K.\n", epsBulk, epsInf, T/Kelvin);
	}
	
	void DielectricFreeEnergy_sub(size_t iStart, size_t iStop, vector3<const double*> eps, const double* s, vector3<double*> p, double* A, vector3<double*> A_eps, double* A_s, const Dielectric& eval)
	{	for(size_t i=iStart; i<iStop; i++) eval.freeEnergy_calc(i, eps, s, p, A, A_eps, A_s);
	}
	void Dielectric::freeEnergy(size_t N,
		vector3<const double*> eps, const double* s, vector3<double*> p, double* A, vector3<double*> A_eps, double* A_s) const
	{	threadLaunch(DielectricFreeEnergy_sub, N, eps, s, p, A, A_eps, A_s, *this);
	}
	
	void DielectricConvertDerivative_sub(size_t iStart, size_t iStop, vector3<const double*> eps, const double* s, vector3<const double*> A_p, vector3<double*> A_eps, double* A_s, const Dielectric& eval)
	{	for(size_t i=iStart; i<iStop; i++) eval.convertDerivative_calc(i, eps, s, A_p, A_eps, A_s);
	}
	void Dielectric::convertDerivative(size_t N, vector3<const double*> eps, const double* s, vector3<const double*> A_p, vector3<double*> A_eps, double* A_s) const
	{	threadLaunch(DielectricConvertDerivative_sub, N, eps, s, A_p, A_eps, A_s, *this);
	}
}
