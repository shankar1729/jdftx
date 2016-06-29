/*-------------------------------------------------------------------
Copyright 2013 Deniz Gunceler

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

#include <electronic/BandMinimizer.h>
#include <electronic/Everything.h>

BandMinimizer::BandMinimizer(Everything& e, int qActive, bool precond):
qActive(qActive), e(e), eVars(e.eVars),
precond(precond)
{	assert(e.cntrl.fixed_H); // Check whether the electron Hamiltonian is fixed
	e.elecMinParams.energyLabel = relevantFreeEnergyName(e);
}

void BandMinimizer::step(const ColumnBundle& dir, double alpha)
{	assert(dir.nCols() == e.eVars.C[qActive].nCols());
	axpy(alpha, dir, e.eVars.C[qActive]);
}

double BandMinimizer::compute(ColumnBundle* grad)
{	return e.eVars.bandEnergyAndGrad(qActive, e.ener, grad, &Kgrad);
}

ColumnBundle BandMinimizer::precondition(const ColumnBundle& grad)
{	return precond ? Kgrad : grad;
}

void BandMinimizer::constrain(ColumnBundle& dir)
{	
}
