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

#include <fluid/IdealGasMonoatomic.h>
#include <core/BlasExtra.h>

IdealGasMonoatomic::IdealGasMonoatomic(const FluidMixture* fluidMixture, const FluidComponent* comp): IdealGas(1,fluidMixture,comp)
{	assert(molecule.isMonoatomic()); //IdealGasMonoatomic must be used only with single site molecules.
}

void IdealGasMonoatomic::initState(const DataRptr* Vex, DataRptr* psi, double scale, double Elo, double Ehi) const
{	DataRptr Veff; nullToZero(Veff, gInfo);
	Veff += V[0];
	Veff += Vex[0];
	//Statistics and min/max capping
	double Emin, Emax, Emean = sum(Veff)/gInfo.nr;
	callPref(eblas_capMinMax)(gInfo.nr, Veff->dataPref(), Emin, Emax, Elo, Ehi);
	logPrintf("\tIdealGasMonoatomic[%s] single molecule energy: min = %le, max = %le, mean = %le\n",
		   molecule.name.c_str(), Emin, Emax, Emean);
	//Set state:
	psi[0] = (-scale/T)*Veff;
}

void IdealGasMonoatomic::getDensities(const DataRptr* psi, DataRptr* N, vector3<>& P0) const
{	N[0] = Nbulk * exp(psi[0]);
	P0 = vector3<>();
}

double IdealGasMonoatomic::compute(const DataRptr* psi, const DataRptr* N, DataRptr* Phi_N, const double Nscale, double& Phi_Nscale) const
{	DataRptr PhiNI_N = T*psi[0] + V[0] - (mu + T);
	Phi_N[0] += PhiNI_N;
	Phi_N[0] += T;
	return gInfo.dV*dot(N[0], PhiNI_N);
}

void IdealGasMonoatomic::convertGradients(const DataRptr* psi, const DataRptr* N, const DataRptr* Phi_N, const vector3<>& Phi_P0, DataRptr* Phi_psi, const double Nscale) const
{	Phi_psi[0] = N[0]*Phi_N[0];
}


