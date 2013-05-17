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

#include <fluid/IdealGasPsiAlpha.h>
#include <fluid/FluidComponent.h>
#include <fluid/Euler.h>


IdealGasPsiAlpha::IdealGasPsiAlpha(const FluidMixture* fluidMixture, const FluidComponent* comp, const SO3quad& quad, const TranslationOperator& trans)
: IdealGasPomega(fluidMixture, comp, quad, trans, comp->molecule.sites.size())
{
}

string IdealGasPsiAlpha::representationName() const
{    return "PsiAlpha";
}

void IdealGasPsiAlpha::initState_o(int o, const matrix3<>& rot, double scale, const DataRptr& Eo, DataRptr* indep) const
{	//InitState loop unused, overriden below:
}

void IdealGasPsiAlpha::initState(const DataRptr* Vex, DataRptr* psi, double scale, double Elo, double Ehi) const
{	IdealGasPomega::initState(Vex, psi, scale, Elo, Ehi);
	//Initialize the state (simply a constant factor times the potential):
	for(unsigned i=0; i<molecule.sites.size(); i++)
	{	DataRptr Veff_i; nullToZero(Veff_i, gInfo);
		Veff_i += V[i];
		Veff_i += Vex[i];
		psi[i] = (-scale/T)*Veff_i;
	}
}

void IdealGasPsiAlpha::getDensities_o(int o, const matrix3<>& rot, const DataRptr* psi, DataRptr& logPomega_o) const
{	for(unsigned i=0; i<molecule.sites.size(); i++)
		for(vector3<> pos: molecule.sites[i]->positions)
			trans.taxpy(-rot*pos, 1., psi[i], logPomega_o);
}

void IdealGasPsiAlpha::convertGradients_o(int o, const matrix3<>& rot, const DataRptr& Phi_logPomega_o, DataRptr* Phi_psi) const
{	for(unsigned i=0; i<molecule.sites.size(); i++)
		for(vector3<> pos: molecule.sites[i]->positions)
			trans.taxpy(rot*pos, 1., Phi_logPomega_o, Phi_psi[i]);
}
