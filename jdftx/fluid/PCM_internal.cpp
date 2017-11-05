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

#include <fluid/PCM_internal.h>
#include <electronic/Everything.h>
#include <core/Operators.h>
#include <core/VectorField.h>
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
	void compute(const ScalarField& n, ScalarField& shape, const double nc, const double sigma)
	{	nullToZero(shape, n->gInfo);
		callPref(compute)(n->gInfo.nr, n->dataPref(), shape->dataPref(), nc, sigma);
	}
	void propagateGradient(const ScalarField& n, const ScalarField& grad_shape, ScalarField& grad_n, const double nc, const double sigma)
	{	nullToZero(grad_n, n->gInfo);
		callPref(propagateGradient)(n->gInfo.nr, n->dataPref(), grad_shape->dataPref(), grad_n->dataPref(), nc, sigma);
	}
}

namespace ShapeFunctionCANDLE
{
	void compute_or_grad(int N, bool grad,
		const double* n, vector3<const double*> Dn, vector3<const double*> Dphi, double* shape,
		const double* A_shape, double* A_n, vector3<double*> A_Dn, vector3<double*> A_Dphi, double* A_pCavity,
		const double nc, const double invSigmaSqrt2, const double pCavity)
	{	threadedLoop(compute_or_grad_calc, N, grad, n, Dn, Dphi, shape, A_shape, A_n, A_Dn, A_Dphi, A_pCavity, nc, invSigmaSqrt2, pCavity);
	}
	#ifdef GPU_ENABLED
	void compute_or_grad_gpu(int N, bool grad,
		const double* n, vector3<const double*> Dn, vector3<const double*> Dphi, double* shape,
		const double* A_shape, double* A_n, vector3<double*> A_Dn, vector3<double*> A_Dphi, double* A_pCavity,
		const double nc, const double invSigmaSqrt2, const double pCavity);
	#endif
	void compute(const ScalarField& n, const ScalarFieldTilde& phi,
		ScalarField& shape, double nc, double sigma, double pCavity)
	{	VectorField Dn = gradient(n);
		VectorField Dphi = I(gradient(phi));
		nullToZero(shape, n->gInfo);
		callPref(compute_or_grad)(n->gInfo.nr, false,
			n->dataPref(), Dn.const_dataPref(), Dphi.const_dataPref(), shape->dataPref(),
			0, 0, vector3<double*>(), vector3<double*>(), 0,
			nc, sqrt(0.5)/sigma, pCavity);
	}
	void propagateGradient(const ScalarField& n, const ScalarFieldTilde& phi, const ScalarField& E_shape,
		ScalarField& E_n, ScalarFieldTilde& E_phi, double& E_pCavity, double nc, double sigma, double pCavity)
	{	VectorField Dn = gradient(n);
		VectorField Dphi = I(gradient(phi));
		nullToZero(E_n, n->gInfo);
		VectorField E_Dn; nullToZero(E_Dn, n->gInfo);
		VectorField E_Dphi; nullToZero(E_Dphi, n->gInfo);
		ScalarField E_pCavityArr; nullToZero(E_pCavityArr, n->gInfo);
		callPref(compute_or_grad)(n->gInfo.nr, true,
			n->dataPref(), Dn.const_dataPref(), Dphi.const_dataPref(), 0,
			E_shape->dataPref(), E_n->dataPref(), E_Dn.dataPref(), E_Dphi.dataPref(), E_pCavityArr->dataPref(),
			nc, sqrt(0.5)/sigma, pCavity);
		Dn=0; Dphi=0; //free memory
		E_n -= divergence(E_Dn);
		E_phi -= divergence(J(E_Dphi));
		E_pCavity += integral(E_pCavityArr);
	}
}

namespace ShapeFunctionSGA13
{
	void expandDensityHelper(int N, double alpha, const double* nBar, const double* DnBarSq, double* nEx, double* nEx_nBar, double* nEx_DnBarSq)
	{	threadedLoop(expandDensity_calc, N, alpha, nBar, DnBarSq, nEx, nEx_nBar, nEx_DnBarSq);
	}
	#ifdef GPU_ENABLED
	void expandDensityHelper_gpu(int N, double alpha, const double* nBar, const double* DnBarSq, double* nEx, double* nEx_nBar, double* nEx_DnBarSq);
	#endif
	void expandDensity(const RadialFunctionG& w, double R, const ScalarField& n, ScalarField& nEx, const ScalarField* A_nEx, ScalarField* A_n)
	{	//Compute weighted densities:
		ScalarFieldTilde nBarTilde = w * J(n);
		ScalarField nBar = I(nBarTilde);
		VectorField DnBar = I(gradient(R * nBarTilde));
		ScalarField DnBarSq = lengthSquared(DnBar);
		//Compute the elementwise function and optionally its derivatives:
		nullToZero(nEx, n->gInfo);
		ScalarField nEx_nBar, nEx_DnBarSq;
		if(A_n)
		{	assert(A_nEx);
			nullToZero(nEx_nBar, n->gInfo);
			nullToZero(nEx_DnBarSq, n->gInfo);
		}
		callPref(expandDensityHelper)(n->gInfo.nr, R*R*R, nBar->dataPref(), DnBarSq->dataPref(), nEx->dataPref(),
			(A_n ? nEx_nBar->dataPref() : 0), (A_n ? nEx_DnBarSq->dataPref() : 0));
		//Propagate gradients if necessary:
		if(A_n)
		{	ScalarFieldTilde A_nBarTilde = Idag((*A_nEx) * nEx_nBar); //contribution through nBar
			ScalarField A_DnBarSq = (*A_nEx) * nEx_DnBarSq;
			A_nBarTilde -= (2.*R) * divergence(Idag(A_DnBarSq * DnBar)); //contribution through DnBar
			(*A_n) += Jdag(w * A_nBarTilde);
		}
	}
}


namespace ShapeFunctionSoftSphere
{
	void compute_thread(size_t iStart, size_t iStop, const vector3<int>& S, const matrix3<>& RTR,
		int nAtoms, const vector3<>* x, int nReps, const vector3<int>* reps, const double* radius, double* shape, double sigmaInv)
	{	vector3<> Sinv(1./S[0], 1./S[1], 1./S[2]);
		THREAD_rLoop( compute_calc(i, iv, Sinv, RTR, nAtoms, x, nReps, reps, radius, shape, sigmaInv); )
	}
	void compute(const vector3<int>& S, const matrix3<>& RTR,
		int nAtoms, const vector3<>* x, int nReps, const vector3<int>* reps, const double* radius, double* shape, double sigmaInv)
	{	threadLaunch(compute_thread, S[0]*S[1]*S[2], S, RTR, nAtoms, x, nReps, reps, radius, shape, sigmaInv);
	}
	#ifdef GPU_ENABLED
	void compute_gpu(const vector3<int>& S, const matrix3<>& RTR,
		int nAtoms, const vector3<>* x, int nReps, const vector3<int>* reps, const double* radius, double* shape, double sigmaInv);
	#endif
	void compute(const std::vector<vector3<>>& x, const std::vector<vector3<int>>& reps, const std::vector<double>& radius, ScalarField& shape, double sigma)
	{	const GridInfo& gInfo = shape->gInfo;
		ManagedArray<vector3<>> xManaged(x);
		ManagedArray<vector3<int>> repsManaged(reps);
		ManagedArray<double> radiusManaged(radius);
		callPref(compute)(gInfo.S, gInfo.RTR, x.size(), xManaged.dataPref(), reps.size(), repsManaged.dataPref(),
			radiusManaged.dataPref(), shape->dataPref(), 1./sigma);
	}

	void propagateGradient_thread(size_t iStart, size_t iStop, const vector3<int>& S, const matrix3<>& RTR,
		const vector3<>& x, int nReps, const vector3<int>* reps, double radius,
		const double* shape, const double* E_shape, vector3<double*> E_x, double* E_radius, double sigmaInv)
	{	vector3<> Sinv(1./S[0], 1./S[1], 1./S[2]);
		THREAD_rLoop( propagateGradient_calc(i, iv, Sinv, RTR, x, nReps, reps, radius, shape, E_shape, E_x, E_radius, sigmaInv); )
	}
	void propagateGradient(const vector3<int>& S, const matrix3<>& RTR,
		const vector3<>& x, int nReps, const vector3<int>* reps, double radius,
		const double* shape, const double* E_shape, vector3<double*> E_x, double* E_radius, double sigmaInv)
	{	threadLaunch(propagateGradient_thread, S[0]*S[1]*S[2], S, RTR, x, nReps, reps, radius, shape, E_shape, E_x, E_radius, sigmaInv);
	}
	#ifdef GPU_ENABLED
	void propagateGradient_gpu(const vector3<int>& S, const matrix3<>& RTR,
		const vector3<>& x, int nReps, const vector3<int>* reps, double radius,
		const double* shape, const double* E_shape, vector3<double*> E_x, double* E_radius, double sigmaInv);
	#endif
	void propagateGradient(const std::vector<vector3<>>& x, const std::vector<vector3<int>>& reps, const std::vector<double>& radius,
		const ScalarField& shape, const ScalarField& E_shape, std::vector<vector3<>>& E_x, std::vector<double>& E_radius, double sigma)
	{	const GridInfo& gInfo = shape->gInfo;
		ManagedArray<vector3<int>> repsManaged(reps);
		VectorField E_xField; nullToZero(E_xField, gInfo);
		ScalarField E_radiusField; nullToZero(E_radiusField, gInfo);
		E_x.resize(x.size());
		E_radius.resize(x.size());
		for(int iAtom=0; iAtom<int(x.size()); iAtom++)
		{	callPref(propagateGradient)(gInfo.S, gInfo.RTR, x[iAtom], reps.size(), repsManaged.dataPref(), radius[iAtom],
				shape->dataPref(), E_shape->dataPref(), E_xField.dataPref(), E_radiusField->dataPref(), 1./sigma);
			E_x[iAtom] += gInfo.dV * sumComponents(E_xField);
			E_radius[iAtom] += integral(E_radiusField);
		}
	}
}


namespace ShapeFunctionSCCS
{
	void compute(int N, const double* n, double* shape, const double rhoMin, const double rhoMax, const double epsBulk)
	{	threadedLoop(compute_calc, N, n, shape, rhoMin, rhoMax, epsBulk);
	}
	void propagateGradient(int N, const double* n, const double* grad_shape, double* grad_n, const double rhoMin, const double rhoMax, const double epsBulk)
	{	threadedLoop(propagateGradient_calc, N, n, grad_shape, grad_n, rhoMin, rhoMax, epsBulk);
	}
	#ifdef GPU_ENABLED
	void compute_gpu(int N, const double* n, double* shape, const double rhoMin, const double rhoMax, const double epsBulk);
	void propagateGradient_gpu(int N, const double* n, const double* grad_shape, double* grad_n, const double rhoMin, const double rhoMax, const double epsBulk);
	#endif
	void compute(const ScalarField& n, ScalarField& shape, const double rhoMin, const double rhoMax, const double epsBulk)
	{	nullToZero(shape, n->gInfo);
		callPref(compute)(n->gInfo.nr, n->dataPref(), shape->dataPref(), rhoMin, rhoMax, epsBulk);
	}
	void propagateGradient(const ScalarField& n, const ScalarField& grad_shape, ScalarField& grad_n, const double rhoMin, const double rhoMax, const double epsBulk)
	{	nullToZero(grad_n, n->gInfo);
		callPref(propagateGradient)(n->gInfo.nr, n->dataPref(), grad_shape->dataPref(), grad_n->dataPref(), rhoMin, rhoMax, epsBulk);
	}
}

//------------- Helper classes for NonlinearPCM  -------------
namespace NonlinearPCMeval
{
	Screening::Screening(bool linear, double T, double Nion, double Zion, double VhsPlus, double VhsMinus, double epsBulk)
	: linear(linear), NT(Nion*T), ZbyT(Zion/T), NZ(Nion*Zion),
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
	
	void ScreeningPhiToState_sub(size_t iStart, size_t iStop, const double* phi, const double* s, const RadialFunctionG& xLookup, bool setState, double* muPlus, double* muMinus, double* kappaSq, const Screening& eval)
	{	for(size_t i=iStart; i<iStop; i++) eval.phiToState_calc(i, phi, s, xLookup, setState, muPlus, muMinus, kappaSq);
	}
	void Screening::phiToState(size_t N, const double* phi, const double* s, const RadialFunctionG& xLookup, bool setState, double* muPlus, double* muMinus, double* kappaSq) const
	{	threadLaunch(ScreeningPhiToState_sub, N, phi, s, xLookup, setState, muPlus, muMinus, kappaSq, *this);
	}
	
	
	Dielectric::Dielectric(bool linear, double T, double Nmol, double pMol, double epsBulk, double epsInf)
	: linear(linear), Np(Nmol * pMol), pByT(pMol/T), NT(Nmol * T),
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
	
	void DielectricPhiToState_sub(size_t iStart, size_t iStop, vector3<const double*> Dphi, const double* s, const RadialFunctionG& gLookup, bool setState, vector3<double*> eps, double* epsilon, const Dielectric& eval)
	{	for(size_t i=iStart; i<iStop; i++) eval.phiToState_calc(i, Dphi, s, gLookup, setState, eps, epsilon);
	}
	void Dielectric::phiToState(size_t N, vector3<const double*> Dphi, const double* s, const RadialFunctionG& gLookup, bool setState, vector3<double*> eps, double* epsilon) const
	{	threadLaunch(DielectricPhiToState_sub, N, Dphi, s, gLookup, setState, eps, epsilon, *this);
	}

}
