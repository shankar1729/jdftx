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

#ifndef JDFTX_ELECTRONIC_EVERYTHING_H
#define JDFTX_ELECTRONIC_EVERYTHING_H

#include <core/GridInfo.h>
#include <core/MinimizeParams.h>
#include <core/Coulomb.h>
#include <electronic/Control.h>
#include <electronic/Basis.h>
#include <electronic/IonInfo.h>
#include <electronic/Symmetries.h>
#include <electronic/ElecInfo.h>
#include <electronic/ElecVars.h>
#include <electronic/Energies.h>
#include <electronic/ExCorr.h>
#include <electronic/Dump.h>
#include <electronic/SCFparams.h>
#include <electronic/IonDynamicsParams.h>
#include <memory>

//! @addtogroup ElectronicDFT
//! @{

//! Container class that contains, well, everything.
class Everything
{
public:
	Control cntrl;  //!< control variables
	Dump dump;      //!< output options
	GridInfo gInfo; //!< main grid descriptor
	std::shared_ptr<GridInfo> gInfoWfns; //!< tighter grid sufficient for wavefunction operations
	std::vector<Basis> basis; //!< wavefunction basis for all k points
	IonInfo iInfo;   //!< ionic system
	Symmetries symm; //!< symmetries
	Symmetries symmUnperturbed; //!< symmetries of unperturbed system in vibration calculations (symm is set to mode=SymmNone in these calculations)
	ExCorr exCorr; //!< Exchange and correlation functional
	std::vector<std::shared_ptr<ExCorr> > exCorrDiff; //!< Other exchange and correlation functionals for comparison
	std::shared_ptr<class ExactExchange> exx; //!< Exact exchange
	ElecInfo eInfo; //!< Auxiliary electronic information
	ElecVars eVars; //!< Electroic variables
	Energies ener;  //!< Energy components
	
	MinimizeParams elecMinParams; //!< electronic minimization parameters
	MinimizeParams ionicMinParams; //!< ionic minimization parameters
	MinimizeParams fluidMinParams; //!< fluid minimization parameters
	MinimizeParams latticeMinParams; //!< lattice minimization parameters
	MinimizeParams inverseKSminParams; //!< Inverse Kohn-sham minimization parameters
	IonDynamicsParams ionDynamicsParams; //!< Molecular dynamics parameters
	SCFparams scfParams; //!< Self-consistent field mixing parameters
	
	CoulombParams coulombParams; //!< Coulomb truncation parameters
	std::shared_ptr<Coulomb> coulomb; //!< Coulomb interaction (optionally truncated)

	std::shared_ptr<VanDerWaals> vanDerWaals; //! Pair potential for vdw correction
	std::shared_ptr<class Vibrations> vibrations; //! Vibrational mode calculator

	//! Call the setup/initialize routines of all the above in the necessray order
	void setup();
	void updateSupercell(bool force=false); //!< (re-)initialize coulombParams.supercell if necessary (or if forced)
};

//! @}
#endif // JDFTX_ELECTRONIC_EVERYTHING_H
