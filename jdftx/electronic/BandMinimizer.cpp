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
#include <electronic/ColumnBundle.h>

BandMinimizer::BandMinimizer(Everything& e, int q): e(e), eVars(e.eVars), eInfo(e.eInfo), q(q)
{	assert(e.cntrl.fixed_H); // Check whether the electron Hamiltonian is fixed
	e.elecMinParams.energyLabel = relevantFreeEnergyName(e);
}

void BandMinimizer::step(const ColumnBundle& dir, double alpha)
{	assert(dir.nCols() == eVars.C[q].nCols());
	axpy(alpha, dir, eVars.C[q]);
	eVars.orthonormalize(q);
}

double BandMinimizer::compute(ColumnBundle* grad, ColumnBundle* Kgrad)
{	if(grad) grad->free();
	diagMatrix Fq = eye(eInfo.nBands);
	const QuantumNumber& qnum = eInfo.qnums[q];
	ColumnBundle Hq;
	double KEq = eVars.applyHamiltonian(q, Fq, Hq, e.ener, true);
	if(grad)
	{	double KErollover = 2.*KEq/(qnum.weight*eInfo.nBands);
		Hq -=  O(eVars.C[q])*eVars.Hsub[q]; //orthonormality contribution
		Hq *= qnum.weight;
		std::swap(*grad, Hq);
		if(Kgrad)
		{	*Kgrad = *grad;
			precond_inv_kinetic(*Kgrad, KErollover); //apply precondition in place
		}
	}
	return qnum.weight * trace(eVars.Hsub[q]).real();
}

void BandMinimizer::constrain(ColumnBundle& dir)
{	dir -= eVars.C[q] * (eVars.C[q]^O(dir));
}
