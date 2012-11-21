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

#ifndef JDFTX_FLUID_FEX_H2O_SCALAREOS_H
#define JDFTX_FLUID_FEX_H2O_SCALAREOS_H
#include <fluid/Fex.h>

class Fex_H2O_ScalarEOS : public Fex
{
public:
	//! Create water with the ScalarEOS functional (can choose soft or hard sphere version)
	Fex_H2O_ScalarEOS(FluidMixture& fluidMixture);

	const Molecule* getMolecule() const { return &molecule; }
	double get_aDiel() const;
	double compute(const DataGptr* Ntilde, DataGptr* grad_Ntilde) const;
	double computeUniform(const double* N, double* grad_N) const;

protected: //Kendra: NOTE: changed these to protected to allow Fex_H2O_Custom to access these variables (but not modify them)
		   //If preferable, I can write get functions, but that seemed a little unnecessary.
	RealKernel fex_LJatt, siteChargeKernel;
	struct ScalarEOS_eval* eval;
	
private:
	SiteProperties propO;
	SiteProperties propH;
	Molecule molecule;
};

class Fex_H2O_Custom : public Fex_H2O_ScalarEOS
{
public:
	//! Create customizable water with the ScalarEOS functional (can choose soft or hard sphere version)
	Fex_H2O_Custom(FluidMixture& fluidMixture, std::vector<H2OSite>& H2OSites);
	
	const Molecule* getMolecule() const { return &customMolecule; }
	double get_aDiel() const;
	
	~Fex_H2O_Custom()
	{
		for (uint i=0; i < prop.size(); i++)
		{
			delete prop[i];
		}
	}

	
private:
	std::vector<SiteProperties*> prop;
	std::vector<std::vector<vector3<>>> pos;
	Molecule customMolecule;
	double aDielFactor;
};


#endif // JDFTX_FLUID_FEX_H2O_SCALAREOS_H
