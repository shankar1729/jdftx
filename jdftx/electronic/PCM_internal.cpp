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
#include <electronic/Everything.h>
#include <electronic/operators.h>
#include <core/Operators.h>
#include <core/DataMultiplet.h>
#include <core/Units.h>
#include <core/Util.h>


namespace ShapeFunction
{
	void compute(int N, const double* n, double* shape, const double nc, const double sigma)
	{	threadedLoop(compute_calc, N, n, shape, nc, sigma);
	}
	void propagateGradient(int N, const double* n, const double* grad_shape, double* grad_n, const double nc, const double sigma)
	{	threadedLoop(propagateGradient_calc, N, n, grad_shape, grad_n, nc, sigma);
	}
	#ifdef GPU_ENABLED
	void compute_gpu(int N, const double* n, double* shape, const double nc, const double sigma);
	void propagateGradient_gpu(int N, const double* n, const double* grad_shape, double* grad_n, const double nc, const double sigma);
	#endif
	void compute(const DataRptr& n, DataRptr& shape, const double nc, const double sigma)
	{	nullToZero(shape, n->gInfo);
		callPref(compute)(n->gInfo.nr, n->dataPref(), shape->dataPref(), nc, sigma);
	}
	void propagateGradient(const DataRptr& n, const DataRptr& grad_shape, DataRptr& grad_n, const double nc, const double sigma)
	{	nullToZero(grad_n, n->gInfo);
		callPref(propagateGradient)(n->gInfo.nr, n->dataPref(), grad_shape->dataPref(), grad_n->dataPref(), nc, sigma);
	}
	
	void expandDensityHelper(int N, double alpha, const double* nBar, const double* DnBarSq, double* nEx, double* nEx_nBar, double* nEx_DnBarSq)
	{	threadedLoop(expandDensity_calc, N, alpha, nBar, DnBarSq, nEx, nEx_nBar, nEx_DnBarSq);
	}
	#ifdef GPU_ENABLED
	void expandDensityHelper_gpu(int N, double alpha, const double* nBar, const double* DnBarSq, double* nEx, double* nEx_nBar, double* nEx_DnBarSq);
	#endif
	void expandDensity(const RadialFunctionG& w, double R, const DataRptr& n, DataRptr& nEx, const DataRptr* A_nEx, DataRptr* A_n)
	{	//Compute weighted densities:
		DataGptr nBarTilde = w * J(n);
		DataRptr nBar = I(nBarTilde);
		DataRptrVec DnBar = I(gradient(R * nBarTilde));
		DataRptr DnBarSq = lengthSquared(DnBar);
		//Compute the elementwise function and optionally its derivatives:
		nullToZero(nEx, n->gInfo);
		DataRptr nEx_nBar, nEx_DnBarSq;
		if(A_n)
		{	assert(A_nEx);
			nullToZero(nEx_nBar, n->gInfo);
			nullToZero(nEx_DnBarSq, n->gInfo);
		}
		callPref(expandDensityHelper)(n->gInfo.nr, R*R*R, nBar->dataPref(), DnBarSq->dataPref(), nEx->dataPref(),
			(A_n ? nEx_nBar->dataPref() : 0), (A_n ? nEx_DnBarSq->dataPref() : 0));
		//Propagate gradients if necessary:
		if(A_n)
		{	DataGptr A_nBarTilde = Idag((*A_nEx) * nEx_nBar); //contribution through nBar
			DataRptr A_DnBarSq = (*A_nEx) * nEx_DnBarSq;
			A_nBarTilde -= (2.*R) * divergence(Idag(A_DnBarSq * DnBar)); //contribution through DnBar
			(*A_n) += Jdag(w * A_nBarTilde);
		}
	}
}


//------------- Helper classes for NonlinearPCM  -------------
namespace NonlinearPCMeval
{
	Screening::Screening(bool linear, double T, double Nion, double Zion, double VhsPlus, double VhsMinus, double epsBulk)
	: linear(linear), NT(Nion*T), NZ(Nion*Zion),
	x0plus(Nion * VhsPlus),
	x0minus(Nion * VhsMinus),
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
		if(!pMol) die("\nNonlinearPCM shluld only be used for polar solvents with non-zero dipole moment.\n");
		if(alpha < 0.)
			die("\nCurrent Nonlinear PCM parameters imply negative correlation for dipole rotations.\n"
				"\tHINT: Reduce solvent molecule dipole or epsInf to fix this.\n");
		if(linear) logPrintf("   Linear dielectric with epsBulk = %lg.\n", epsBulk);
		else logPrintf("   Nonlinear dielectric with epsBulk = %lg and epsInf = %lg with density Nmol = %lg of dipoles pMol = %lg at T = %lg K.\n", epsBulk, epsInf, Nmol, pMol, T/Kelvin);
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
