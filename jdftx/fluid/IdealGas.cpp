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

#include <fluid/IdealGas.h>
#include <fluid/FluidMixture.h>

IdealGas::IdealGas(int nIndep, const FluidMixture* fluidMixture, const FluidComponent* comp)
: nIndep(nIndep), molecule(comp->molecule), gInfo(fluidMixture->gInfo), T(fluidMixture->T),
V(molecule.sites.size()), Nbulk(0), mu(0)
{
}

double IdealGas::get_Nbulk()
{	return Nbulk;
}

void IdealGas::overrideBulk(double Nbulk, double mu)
{	this->Nbulk = Nbulk;
	this->mu = mu;
}
