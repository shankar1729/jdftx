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
#include <core/Operators.h>
#include <core/DataMultiplet.h>
#include <core/Units.h>
#include <core/Util.h>

PCM::PCM(const Everything& e, const FluidSolverParams& fsp): FluidSolver(e), params(fsp)
{	k2factor = (8*M_PI/params.T) * params.ionicConcentration * pow(params.ionicZelectrolyte,2);
}

void PCM::dumpDebug(FILE* fp) const
{	DataRptrVec shape_x = gradient(shape);
	DataRptr surfaceDensity = sqrt(shape_x[0]*shape_x[0] + shape_x[1]*shape_x[1] + shape_x[2]*shape_x[2]);
	
	fprintf(fp, "Cavity volume = %f\n", integral(1.-shape));
	fprintf(fp, "Cavity surface Area = %f\n", integral(surfaceDensity));

	fprintf(fp, "\nComponents of Adiel:\n");
	Adiel.print(fp, true, "   %13s = %25.16lf\n");	
}

//----------------------- The JDFT `shape function' and gradient ------------------

namespace ShapeFunction
{
	void compute(int N, const double* nCavity, double* shape, const double nc, const double sigma)
	{	threadedLoop(compute_calc, N, nCavity, shape, nc, sigma);
	}
	void propagateGradient(int N, const double* nCavity, const double* grad_shape, double* grad_nCavity, const double nc, const double sigma)
	{	threadedLoop(propagateGradient_calc, N, nCavity, grad_shape, grad_nCavity, nc, sigma);
	}
	#ifdef GPU_ENABLED
	void compute_gpu(int N, const double* nCavity, double* shape, const double nc, const double sigma);
	void propagateGradient_gpu(int N, const double* nCavity, const double* grad_shape, double* grad_nCavity, const double nc, const double sigma);
	#endif
	void compute(const DataRptr& nCavity, DataRptr& shape, const FluidSolverParams& fsp)
	{	nullToZero(shape, nCavity->gInfo);
		callPref(compute)(nCavity->gInfo.nr, nCavity->dataPref(), shape->dataPref(), fsp.nc, fsp.sigma);
	}
	void propagateGradient(const DataRptr& nCavity, const DataRptr& grad_shape, DataRptr& grad_nCavity, const FluidSolverParams& fsp)
	{	nullToZero(grad_nCavity, nCavity->gInfo);
		callPref(propagateGradient)(nCavity->gInfo.nr, nCavity->dataPref(), grad_shape->dataPref(), grad_nCavity->dataPref(), fsp.nc, fsp.sigma);
	}
}

void citePCM(const FluidSolverParams& fsp)
{
	switch(fsp.pcmVariant)
	{	case PCM_SGA13:
			Citations::add("Linear/nonlinear dielectric/ionic fluid model with weighted-density cavitation and dispersion",
				"R. Sundararaman, D. Gunceler, and T.A. Arias, (under preparation)");
			break;
		case PCM_GLSSA13:
			Citations::add("Linear/nonlinear dielectric/ionic fluid model with effective cavity tension",
				"D. Gunceler, K. Letchworth-Weaver, R. Sundararaman, K.A. Schwarz and T.A. Arias, arXiv:1301.6189");
			break;
		case PCM_LA12:
		case PCM_PRA05:
			if(fsp.ionicConcentration)
				Citations::add("Linear dielectric fluid model with ionic screening",
					"K. Letchworth-Weaver and T.A. Arias, Phys. Rev. B 86, 075140 (2012)");
			else
				Citations::add("Linear dielectric fluid model",
					"S.A. Petrosyan SA, A.A. Rigos and T.A. Arias, J Phys Chem B. 109, 15436 (2005)");
			break;
	}
}


namespace Cavitation
{
	void energyAndGrad(EnergyComponents& E, const DataRptr& shape, DataRptr& E_shape, const FluidSolverParams& fsp)
	{	switch(fsp.pcmVariant)
		{	case PCM_SGA13:
				energyAndGradWDA(E, shape, E_shape, fsp);
				break;
			case PCM_GLSSA13:
				energyAndGradEffectiveTension(E, shape, E_shape, fsp);
				break;
			case PCM_LA12:
			case PCM_PRA05:
			default:
				break; //no contribution
		}
	}
	
	void energyAndGradEffectiveTension(EnergyComponents& E, const DataRptr& shape, DataRptr& E_shape, const FluidSolverParams& fsp)
	{	DataRptrVec shape_x = gradient(shape);
		DataRptr surfaceDensity = sqrt(shape_x[0]*shape_x[0] + shape_x[1]*shape_x[1] + shape_x[2]*shape_x[2]);
		double surfaceArea = integral(surfaceDensity);
		DataRptr invSurfaceDensity = inv(surfaceDensity);
		E_shape += (-fsp.cavityTension)*divergence(shape_x*invSurfaceDensity); // Surface term
		E["CavityTension"] = surfaceArea*fsp.cavityTension;
	}
	
	void energyAndGradWDA(EnergyComponents& E, const DataRptr& shape, DataRptr& E_shape, const FluidSolverParams& fsp)
	{	die("Not yet implemented.\n");
	}
	
	void print(const FluidSolverParams& fsp)
	{	switch(fsp.pcmVariant)
		{	case PCM_SGA13:
				logPrintf("   Weighted density cavitation model constrained by Nbulk: %lg bohr^-3, Pvap: %lg kPa and sigmaBulk: %lg Eh/bohr^2 at T: %lg K.\n",
					fsp.Nbulk, fsp.Pvap/KPascal, fsp.sigmaBulk, fsp.T/Kelvin);
				logPrintf("   Weighted density dispersion model using vdW pair potentials.\n");
				break;
			case PCM_GLSSA13:
				logPrintf("   Effective cavity tension: %lg Eh/bohr^2 to account for cavitation and dispersion.\n", fsp.cavityTension);
				break;
			case PCM_LA12:
			case PCM_PRA05:
			default:
				logPrintf("   No cavitation model.\n");
				break;
		}
	}
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
