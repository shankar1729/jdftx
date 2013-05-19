/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman

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

#ifndef FLUID1D_FLUID1D_FEX_TM_SCALAREOS_H
#define FLUID1D_FLUID1D_FEX_TM_SCALAREOS_H

#include <fluid1D/Fex.h>

//! Base class of functionals for relatively nonpolar fluids based on the Tao-Mason equation of state
class Fex_TM_ScalarEOS : public Fex
{
public:
	//!Fluid EOS determined by critical point (Tc,Pc) and acentricity (omega)
	//!Hard sphere radius is used to separate out the repulsive part to be treated by FMT
	//!Note: assumes that only a single site of the molecule has a hard sphere radius which is equal to the value the supplied value
	Fex_TM_ScalarEOS(FluidMixture& fluidMixture, double Tc, double Pc, double omega, double sphereRadius);
	virtual ~Fex_TM_ScalarEOS();
	
	double compute(const ScalarFieldTilde* Ntilde, ScalarFieldTilde* grad_Ntilde) const;
	double computeUniform(const double* N, double* grad_N) const;
	void directCorrelations(const double* N, ScalarFieldTildeCollection& C) const;
	
	double vdwRadius() const; //!< get vdW radius corresponding to equation of state
private:
	SphericalKernel fex_LJatt;
protected:
	SphericalKernel siteChargeKernel;
	struct TaoMasonEOS_eval* eval;
};

class Fex_CHCl3_ScalarEOS : public Fex_TM_ScalarEOS
{
public:
	Fex_CHCl3_ScalarEOS(FluidMixture& fluidMixture);
	const Molecule* getMolecule() const { return &molecule; }
	double get_aDiel() const { return 1.; }
private:
	SiteProperties propC;
	SiteProperties propH;
	SiteProperties propCl;
	Molecule molecule;
};

class Fex_CCl4_ScalarEOS : public Fex_TM_ScalarEOS
{
public:
	Fex_CCl4_ScalarEOS(FluidMixture& fluidMixture);
	const Molecule* getMolecule() const { return &molecule; }
	double get_aDiel() const { return 1.; }
private:
	SiteProperties propC;
	SiteProperties propCl;
	Molecule molecule;
};

#endif //FLUID1D_FLUID1D_FEX_TM_SCALAREOS_H
