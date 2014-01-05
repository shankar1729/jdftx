/*-------------------------------------------------------------------
Copyright 2013 Deniz Gunceler, Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_SCFPARAMS_H
#define JDFTX_ELECTRONIC_SCFPARAMS_H

#include <core/Util.h>
#include <vector>

struct SCFparams
{
	int nIterations; //!< maximum iterations (single point calculation if 0)
	double energyDiffThreshold; //!< convergence threshold for energy difference between successive iterations
	double eigDiffThreshold; //!< convergence threshold on the RMS change of eigenvalues
	double residualThreshold; //!< convergence threshold on the residual
	
	enum MixedVariable
	{	MV_Density, //! Mix electron density (n) and kinetic energy density (tau)
		MV_Potential //! Mix the local electronic potential (Vscloc) and the kinetic energy potential (Vtau)
	}
	mixedVariable; //! Whether we are mixing the density or the potential
	
	enum VectorExtrapolation
	{	VE_Plain,  //!< No vector extrapolation, just mixes the new density with the damping fraction
		VE_DIIS   //!< Direct Inversion in the Iterative Subspace (DIIS), mixes all cached densities to minimize the residual (reminiscient of Krylov subspace methods)
	}
	vectorExtrapolation; //!< Vector extrapolation method used to construct the new density
	
	int history; //!< Number of past residuals and vectors are kept cached and used in DIIS
	bool verbose; //!< Whether the inner eigensolver will print progress
	double mixFraction;  //!< Maximum fraction of the new variable that will be mixed with the old one
	double qKerker; //!< Wavevector controlling Kerker preconditioning (if negative, auto-set to Gmin)
	double qMetric; //!< Wavevector controlling the DIIS metric (if negative, auto-set to Gmin)
	
	struct EigenShift
	{	int q; //!< Quantum number
		int n; //!< Band index
		double shift; //!< Energy shift
		bool fromHOMO; //!< Whether n-indexing is done from HOMO (or not)
	};
	std::vector<EigenShift> eigenShifts; //! A list of all eigenshifts, used for non-ground-state calculations
	
	SCFparams()
	{	nIterations = 20;
		energyDiffThreshold = 1e-8;
		eigDiffThreshold = 1e-8;
		residualThreshold = 1e-7;
		mixedVariable = MV_Potential;
		vectorExtrapolation = VE_DIIS;
		verbose = false;
		mixFraction = 0.5;
		qKerker = -1.; //i.e. auto-set
		qMetric = -1.; //i.e. auto-set
		history = 15;
	}
};

#endif // JDFTX_ELECTRONIC_SCFPARAMS_H
