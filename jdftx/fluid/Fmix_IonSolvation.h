/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver

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

#ifndef JDFTX_FLUID_FMIX_IONSOLVATION_H
#define JDFTX_FLUID_FMIX_IONSOLVATION_H

#include <fluid/FluidMixture.h>
#include <fluid/Fex_HardSphereIon.h>

class Fmix_IonSolvation : public Fmix
{
public:
    Fmix_IonSolvation(FluidMixture& fluidMixture, const Fex& fexH2O, const Fex& fexIon, double Esolv, double Rsolv)
    : Fmix(fluidMixture), fexH2O(fexH2O), fexIon(fexIon), Ksolv(gInfo)
    {	initGaussianKernel(Ksolv, Rsolv);
		Kmul = -Esolv*(4*M_PI*pow(Rsolv,3))/3;
	}
	Fmix_IonSolvation(FluidMixture& fluidMixture, const Fex& fexH2O, const HardSphereIon* Ion)
    : Fmix(fluidMixture), fexH2O(fexH2O), fexIon(*(Ion->fex)), Ksolv(gInfo)
    {	initGaussianKernel(Ksolv, Ion->Rsolv);
		Kmul = -Ion->Esolv*(4*M_PI*pow(Ion->Rsolv,3))/3;
	}
    string getName() const { return fexIon.getMolecule()->name + "+H2O"; }
    double compute(const DataGptrCollection& Ntilde, DataGptrCollection& grad_Ntilde) const
	{	unsigned iIon = fluidMixture.get_offsetDensity(&fexIon);
		unsigned iH2O = fluidMixture.get_offsetDensity(&fexH2O);
		DataGptr V_Ion = gInfo.nr * Kmul*(Ksolv * Ntilde[iIon]);
		DataGptr V_H2O = gInfo.nr * Kmul*(Ksolv * Ntilde[iH2O]);
		grad_Ntilde[iIon] += V_H2O;
		grad_Ntilde[iH2O] += V_Ion;
		return gInfo.dV*dot(V_Ion,Ntilde[iH2O]);
	}
	double computeUniform(const std::vector<double>& N, std::vector<double>& grad_N) const
	{	unsigned iIon = fluidMixture.get_offsetDensity(&fexIon);
		unsigned iH2O = fluidMixture.get_offsetDensity(&fexH2O);
		grad_N[iIon] += Kmul*Ksolv.data[0]*N[iH2O];
		grad_N[iH2O] += Kmul*Ksolv.data[0]*N[iIon];
		return N[iIon]*Kmul*Ksolv.data[0]*N[iH2O];
	}
private:
	const Fex &fexH2O, &fexIon;
	RealKernel Ksolv; double Kmul; //shape function and prefactor
};

#endif //  JDFTX_FLUID_FMIX_IONSOLVATION_H
