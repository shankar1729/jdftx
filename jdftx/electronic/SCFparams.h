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
#include <core/PulayParams.h>

//! @addtogroup ElecSystem
//! @{

//! Parameters controlling SCF iteration
struct SCFparams : public PulayParams
{
	int nEigSteps; //!< number of steps of the eigenvalue solver per iteration (use elecMinParams.nIterations if 0)
	double eigDiffThreshold; //!< convergence threshold on the RMS change of eigenvalues

	string historyFilename; //!< Read SCF history in order to resume a previous run
	
	enum MixedVariable
	{	MV_Density, //!< Mix electron density (n) and kinetic energy density (tau)
		MV_Potential //!< Mix the local electronic potential (Vscloc) and the kinetic energy potential (Vtau)
	}
	mixedVariable; //!< Whether we are mixing the density or the potential
	
	double qKerker; //!< Wavevector controlling Kerker preconditioning
	double qKappa; //!< wavevector controlling long-range damping (if negative, auto-set to zero or fluid Debye wave-vector as appropriate)
	
	bool verbose; //!< Whether the inner eigensolver will print progress
	double mixFractionMag;  //!< Mixing fraction for magnetization density / potential
	
	SCFparams()
	{	nEigSteps = 2; //for Davidson; the default for CG is 40 (and set by the command)
		eigDiffThreshold = 1e-8;
		mixedVariable = MV_Density;
		qKerker = 0.8;
		qKappa = -1.;
		verbose = false;
		mixFractionMag = 1.5;
	}
};

//! @}
#endif // JDFTX_ELECTRONIC_SCFPARAMS_H
