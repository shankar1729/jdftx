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


class Control
{
public:
	bool fixed_n; //!< fixed density (band structure) mode for electronic sector
	bool fixOccupied; //!< whether to hold occupied orbitals fixed in band structure calculations
	double occupiedThreshold; //!< fillings threshold for occupied states
	
	BasisKdep basisKdep; //!< k-dependence of basis
	double Ecut; //!< energy cutoff

	double dragRadius; //!< typical region of space around each atom "dragged" on each ionic step (0 to disable)
	
	int fluidGummel_nIterations; //!< max iterations of the fluid<->electron self-consistency loop
	double fluidGummel_Atol; //!< stopping free-energy tolerance for the fluid<->electron self-consistency loop

	bool shouldPrintEigsFillings; //!< whether eigenvalues and fillings should be printed at each iteration
	bool shouldPrintEcomponents; //!< whether energy components should be printed at each iteration
	bool shouldPrintMuSearch; //!< whether mu bisection progress should be printed

	bool invertKS; //!< Kohn-Sham inversion (sequence of band structure solves to find optimum potential)
	bool invertKS_nonlocal; //!< whether to retain non-local portions of pseudopotential for Kohn-Sham inversion
	double invertKS_sigma; //!< bandwidth cutoff for the external potential
	string invertKS_chiGuessFilename; //!< filename pattern of variables (wfns/fillings/eigenvals) generating guess chi

	Control()
	:	fixed_n(false),
		fixOccupied(false), occupiedThreshold(0),
		basisKdep(BasisKpointDep), Ecut(0), dragRadius(0),
		fluidGummel_nIterations(10), fluidGummel_Atol(1e-5),
		shouldPrintEigsFillings(false), shouldPrintEcomponents(false), shouldPrintMuSearch(false),
		invertKS(false), invertKS_nonlocal(true), invertKS_sigma(0.)
	{
	}
};
#endif // JDFTX_ELECTRONIC_CONTROL_H
