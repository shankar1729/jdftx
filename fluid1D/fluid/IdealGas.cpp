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

#include <fluid1D/IdealGas.h>
#include <fluid1D/FluidMixture.h>

IdealGas::IdealGas(int nIndep, Fex* fex, double xBulk)
: nIndep(nIndep), molecule(fex->getMolecule()),
fluidMixture(fex->fluidMixture), gInfo(fluidMixture.gInfo), T(fluidMixture.T), V(molecule->nIndices),
xBulk(xBulk), Nnorm(0), Nbulk(0), mu(0)
{	fluidMixture.addComponent(this, fex);
}

double IdealGas::get_Nbulk()
{	return Nbulk;
}

void IdealGas::set_Nnorm(double Nnorm)
{	if(Nnorm>0) this->Nnorm=Nnorm;
}

void IdealGas::overrideBulk(double Nbulk, double mu)
{	this->Nbulk = Nbulk;
	this->mu = mu;
}
