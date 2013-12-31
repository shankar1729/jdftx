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

#ifndef JDFTX_ELECTRONIC_RESIDUALMINIMIZERPARAMS_H
#define JDFTX_ELECTRONIC_RESIDUALMINIMIZERPARAMS_H

#include <core/Util.h>
#include <vector>

enum MixedVariable
{
	density, //! Mix electron density (n) and kinetic energy density (tau)
	potential //! Mix the local electronic potential (Vscloc) and the kinetic energy potential (Vtau)
};

enum VectorExtrapolation
{
	plain,  //! No vector extrapolation, just half-mixes the new density
	DIIS   //! Direct Inversion in the Iterative Subspace (DIIS), mixes all cached densities to minimize the residual (reminiscient of Krylov subspace methods)
};

struct EigenShift
{	int q; //! Quantum number
	int n; //! Band index
	double shift; //! Energy shift
	
	bool fromHOMO; //! Whether n-indexing is done from HOMO (or not)
	
	EigenShift(int q, int n, double shift, bool fromHOMO): q(q), n(n), shift(shift), fromHOMO(fromHOMO){}
};

struct ResidualMinimizerParams
{
	int nIterations; //! maximum iterations (single point calculation if 0)
	double energyDiffThreshold; //! convergence threshold for energy difference between successive iterations
	MixedVariable mixedVariable; //! Whether we are mixing the density or the potential
	VectorExtrapolation vectorExtrapolation; //! Vector extrapolation method used to construct the new density
	int history; //! Number of past residuals and vectors are kept cached and used in DIIS
	bool verbose; //! Whether the inner eigensolver will print process
	double damping;  //! Fraction of the old variable that will be mixed with the new one

	std::vector<EigenShift> eigenShifts; //! A list of all eigenshifts, used for non-ground-state calculations
	
	ResidualMinimizerParams()
	{
		nIterations = 10;
		energyDiffThreshold = 1e-6;
		mixedVariable = potential;
		vectorExtrapolation = DIIS;
		verbose = false;
		damping = 0.5;
		history = 15;
	}
	
};

#endif // JDFTX_ELECTRONIC_RESIDUALMINIMIZERPARAMS_H