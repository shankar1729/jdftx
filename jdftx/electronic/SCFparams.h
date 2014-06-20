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
	int nEigSteps; //!< number of steps of the eigenvalue solver per iteration (use elecMinParams.nIterations if 0)
	double energyDiffThreshold; //!< convergence threshold for energy difference between successive iterations
	double eigDiffThreshold; //!< convergence threshold on the RMS change of eigenvalues
	double residualThreshold; //!< convergence threshold on the residual

	string historyFilename; //!< Read SCF history in order to resume a previous run
	
	enum MixedVariable
	{	MV_Density, //! Mix electron density (n) and kinetic energy density (tau)
		MV_Potential //! Mix the local electronic potential (Vscloc) and the kinetic energy potential (Vtau)
	}
	mixedVariable; //! Whether we are mixing the density or the potential
	
	int history; //!< Number of past residuals and vectors are kept cached and used in DIIS
	bool verbose; //!< Whether the inner eigensolver will print progress
	double mixFraction;  //!< Mixing fraction for total density / potential
	double mixFractionMag;  //!< Mixing fraction for magnetization density / potential
	double qKerker; //!< Wavevector controlling Kerker preconditioning (if negative, auto-set to Gmin)
	double qMetric; //!< Wavevector controlling the DIIS metric (if negative, auto-set to Gmin)
	double sp_constraint; //! Whether to ensure that the exchange is exact in the single particle limit
	bool MOMenabled; //! Whether to use the maximum-overlap method or not
	
	struct EigenShift
	{	int q; //!< Quantum number
		int n; //!< Band index
		double shift; //!< Energy shift
		bool fromHOMO; //!< Whether n-indexing is done from HOMO (or not)
	};
	std::vector<EigenShift> eigenShifts; //! A list of all eigenshifts, used for non-ground-state calculations
	
	SCFparams()
	{	nIterations = 50;
		nEigSteps = 2; //for Davidson; the default for CG is 40 (and set by the command)
		energyDiffThreshold = 1e-8;
		eigDiffThreshold = 1e-8;
		residualThreshold = 1e-7;
		mixedVariable = MV_Potential;
		verbose = false;
		mixFraction = 0.5;
		mixFractionMag = 1.5;
		qKerker = 0.8;
		qMetric = 0.8;
		history = 10;
		sp_constraint = 0.;
		MOMenabled = false;
	}
};

#endif // JDFTX_ELECTRONIC_SCFPARAMS_H
