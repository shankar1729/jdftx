/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

This file is part of Fluid1D.

Fluid1D is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Fluid1D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Fluid1D.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#include <fluid1D/IdealGasMonoatomic.h>
#include <core/BlasExtra.h>

IdealGasMonoatomic::IdealGasMonoatomic(Fex* fex, double xBulk): IdealGas(1,fex,xBulk)
{	assert(molecule->nSites==1); //IdealGasMonoatomic must be used only with single site molecules.
	assert(molecule->site[0].pos==vector3<>(0,0,0)); //The single site for IdealGasMonoatomic must be the origin
}

void IdealGasMonoatomic::initState(const ScalarField* Vex, ScalarField* psi, double scale, double Elo, double Ehi) const
{	ScalarField Veff; nullToZero(Veff, gInfo);
	Veff += V[0];
	Veff += Vex[0];
	//Statistics and min/max capping
	double Emin, Emax, Emean = integral(Veff)/gInfo.Volume();
	eblas_capMinMax(gInfo.S, Veff.data(), Emin, Emax, Elo, Ehi);
	logPrintf("\tIdealGasMonoatomic[%s] single molecule energy: min = %le, max = %le, mean = %le\n",
		   molecule->name.c_str(), Emin, Emax, Emean);
	//Set state:
	psi[0] = (-scale/T)*Veff;
}

void IdealGasMonoatomic::getDensities(const ScalarField* psi, ScalarField* N, double& P) const
{	N[0] = Nbulk * exp(psi[0]);
}

double IdealGasMonoatomic::compute(const ScalarField* psi, const ScalarField* N, ScalarField* grad_N,
	const double& P, double& grad_P, const double Nscale, double& grad_Nscale) const
{	ScalarField PhiNI_N = DiagJdagOJ1(T*psi[0] - mu);
	if(V[0]) PhiNI_N += JdagOJ(V[0]);
	grad_N[0] += PhiNI_N;
	return dot(N[0], PhiNI_N) - T*integral(N[0]);
}

void IdealGasMonoatomic::convertGradients(const ScalarField* psi, const ScalarField* N,
	const ScalarField* grad_N, double grad_P, ScalarField* grad_psi, const double Nscale) const
{	grad_psi[0] = Diag(N[0])*grad_N[0];
}

