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

#ifndef JDFTX_FLUID_FEX_HARDSPHEREION_H
#define JDFTX_FLUID_FEX_HARDSPHEREION_H

#include <fluid/FluidMixture.h>

struct HardSphereIon
{
	string name; //!< name of the hard sphere ion  
	double Z; //!< site charge on the ion in electrons. Total electronic charge is Znuc + Z 
	double rHS; //!< hard sphere radius 
	double Concentration; //!< concentration in Mol/L
	double Esolv; //!< mixing parameter with H2O: depth of interaction potential in hartree
	double Rsolv; //!< mixing parameter with H2O; width of gaussian kernel interaction potential in bohr
	bool MixFunc; //!< True if hard sphere mixes with water
	double Znuc; //!< Nuclear charge of ion in coupling functional electron density model
	double CouplingWidth; //!< Exponential width of electron density model in convolution coupling
	string CouplingFilename; //!< Filename to specify electron density model for convolution coupling
	ConvolutionCouplingSiteModel convCouplingModel; //!<True if hard sphere couples to electrons through convolution coupling
	Fex* fex;
	IdealGas* idgas;
	Fmix* fmix;
};

class Fex_HardSphereIon : public Fex
{
public:
    Fex_HardSphereIon(FluidMixture& fluidMixture, double rHS, double Z, string name)
    : Fex(fluidMixture), Qkernel(gInfo), prop(gInfo, rHS,0, Z,&Qkernel), molecule(name, &prop, vector3<>(0,0,0))
    {
		initGaussianKernel(Qkernel, rHS);
	}
	
	Fex_HardSphereIon(FluidMixture& fluidMixture, const HardSphereIon* Ion)
    : Fex(fluidMixture), Qkernel(gInfo), prop(gInfo, Ion->rHS, 0, Ion->Z, &Qkernel), molecule(Ion->name, &prop, vector3<>(0,0,0))
    {
		initGaussianKernel(Qkernel, Ion->rHS);
	}

    const Molecule* getMolecule() const { return &molecule; }
    double get_aDiel() const {return 1.0; }
    double compute(const DataGptr* Ntilde, DataGptr* grad_Ntilde) const { return 0.0; }
    double computeUniform(const double* N, double* grad_N) const { return 0.0; }

private:
	RealKernel Qkernel;
	SiteProperties prop;
	Molecule molecule;
};

#endif //JDFTX_FLUID_FEX_HARDSPHEREION_H