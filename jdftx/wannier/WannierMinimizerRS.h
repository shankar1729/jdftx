/*-------------------------------------------------------------------
Copyright 2014 Ravishankar Sundararaman

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

#ifndef JDFTX_WANNIER_WANNIERMINIMIZERRS_H
#define JDFTX_WANNIER_WANNIERMINIMIZERRS_H

#include <wannier/WannierMinimizer.h>

class WannierMinimizerRS : public WannierMinimizer
{
public:
	WannierMinimizerRS(const Everything& e, const Wannier& wannier);

	void initialize(int iSpin);
	double getOmega(bool grad);
	double getOmegaI(bool grad);
	
private:
	VectorField r; ScalarField rSq; //r and r^2 wrapped on the Wigner-Seitz cell
	int iSpin; //!< spin channel currently being minimized
	int nStart, nStop; //!< center MPI division
	double getOmega(bool grad, bool invariant); //compute Omega or OmegaI depending on invariant
};

#endif //  JDFTX_WANNIER_WANNIERMINIMIZERRS_H