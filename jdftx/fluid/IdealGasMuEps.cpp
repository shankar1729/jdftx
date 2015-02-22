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

#include <fluid/IdealGasMuEps.h>
#include <fluid/Euler.h>

IdealGasMuEps::IdealGasMuEps(const FluidMixture* fluidMixture, const FluidComponent* comp,  const SO3quad& quad, const TranslationOperator& trans)
: IdealGasPomega(fluidMixture, comp, quad, trans, 4)
{
}

string IdealGasMuEps::representationName() const
{    return "MuEps";
}

void IdealGasMuEps::initState_o(int o, const matrix3<>& rot, double scale, const ScalarField& Eo, ScalarField* mueps) const
{	vector3<> pVec = rot * pMol;
	mueps[0] += (-quad.weight(o)*scale/T) * Eo;
	for(int k=0; k<3; k++)
		mueps[k+1] += (-pVec[k]*quad.weight(o)*scale/T) * Eo;
}

void IdealGasMuEps::getDensities_o(int o, const matrix3<>& rot, const ScalarField* mueps, ScalarField& logPomega_o) const
{	vector3<> pVec = rot * pMol;
	logPomega_o += mueps[0];
	for(int k=0; k<3; k++)
		logPomega_o += pVec[k]*mueps[k+1];
}

void IdealGasMuEps::convertGradients_o(int o, const matrix3<>& rot, const ScalarField& Phi_logPomega_o, ScalarField* Phi_mueps) const
{	vector3<> pVec = rot * pMol;
	Phi_mueps[0] += Phi_logPomega_o;
	for(int k=0; k<3; k++)
		Phi_mueps[k+1] += pVec[k] * Phi_logPomega_o;
}
