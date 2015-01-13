/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman
Copyright 1996-2003 Sohrab Ismail-Beigi

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

#ifndef JDFTX_ELECTRONIC_CONTROL_H
#define JDFTX_ELECTRONIC_CONTROL_H

#include <electronic/common.h>
#include <core/vector3.h>

//! K-point dependence of basis
enum BasisKdep { BasisKpointDep, BasisKpointIndep } ; 
static EnumStringMap<BasisKdep> kdepMap(BasisKpointDep, "kpoint-dependent", BasisKpointIndep, "single" );

//! Electronic eigenvalue method
enum ElecEigenAlgo { ElecEigenCG, ElecEigenDavidson };

class Control
{
public:
	bool fixed_H; //!< fixed Hamiltonian (band structure) mode for electronic sector
	bool fixOccupied; //!< whether to hold occupied orbitals fixed in band structure calculations
	double occupiedThreshold; //!< fillings threshold for occupied states
	bool cacheProjectors; //!< whether to cache nonlocal projectors
	double davidsonBandRatio; //!< ratio of number of Davidson working bands to actual bands in system (>= 1)
	
	ElecEigenAlgo elecEigenAlgo; //!< Eigenvalue algorithm
	BasisKdep basisKdep; //!< k-dependence of basis
	double Ecut, EcutRho; //!< energy cutoff for electrons and charge density grid (EcutRho=0 => EcutRho = 4 Ecut)
	
	bool dragWavefunctions; //!< whether to drag wavefunctions using atomic orbital projections on ionic steps
	vector3<> lattMoveScale; //!< preconditioning factor for each lattice vector during lattice minimization
	
	int fluidGummel_nIterations; //!< max iterations of the fluid<->electron self-consistency loop
	double fluidGummel_Atol; //!< stopping free-energy tolerance for the fluid<->electron self-consistency loop

	double overlapConditionThreshold; //!< Threshold for overlap condition number at which wavefunctions are re-orthogonalized
	int overlapCheckInterval; //!< Number of electronic steps between overlap condition checks
	
	bool shouldPrintEigsFillings; //!< whether eigenvalues and fillings should be printed at each iteration
	bool shouldPrintEcomponents; //!< whether energy components should be printed at each iteration
	bool shouldPrintMuSearch; //!< whether mu bisection progress should be printed
	bool shouldPrintKpointsBasis; //!< whether individual kpoint and basis details should be printed at the beginning
	
	bool invertKS; //!< Kohn-Sham inversion (sequence of band structure solves to find optimum potential)
	bool invertKS_nonlocal; //!< whether to retain non-local portions of pseudopotential for Kohn-Sham inversion
	double invertKS_sigma; //!< bandwidth cutoff for the external potential
	string invertKS_chiGuessFilename; //!< filename pattern of variables (wfns/fillings/eigenvals) generating guess chi
	
	bool scf; //! whether SCF iteration or total energy minimizer will be called
	bool convergeEmptyStates; //! whether to converge empty states after every electronic minimization
	bool dumpOnly; //! run a single-electronic-point energy evaluation and process the end dump
	
	Control()
	:	fixed_H(false),
		fixOccupied(false), occupiedThreshold(0), cacheProjectors(true), davidsonBandRatio(1.1),
		elecEigenAlgo(ElecEigenDavidson), basisKdep(BasisKpointDep), Ecut(0), EcutRho(0), dragWavefunctions(true),
		fluidGummel_nIterations(10), fluidGummel_Atol(1e-5),
		overlapConditionThreshold(1.5), overlapCheckInterval(20),
		shouldPrintEigsFillings(false), shouldPrintEcomponents(false), shouldPrintMuSearch(false), shouldPrintKpointsBasis(false),
		invertKS(false), invertKS_nonlocal(true), invertKS_sigma(0.), scf(false), convergeEmptyStates(false), dumpOnly(false)
	{
	}
};
#endif // JDFTX_ELECTRONIC_CONTROL_H
