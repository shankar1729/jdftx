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
#include <electronic/common.h>
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
#include <memory>

class Everything
{
public:
	Control cntrl;
	Dump dump;
	GridInfo gInfo;
	std::shared_ptr<GridInfo> gInfoWfns; //tighter grid sufficient for wavefunction operations
	std::vector<Basis> basis;
	IonInfo iInfo;
	Symmetries symm;
	Symmetries symmUnperturbed; //symmetries of unperturbed system in vibration calculations (symm is set to mode=SymmNone in these calculations)
	ExCorr exCorr; //!< Exchange and correlation functional
	std::vector<std::shared_ptr<ExCorr> > exCorrDiff; //!< Other exchange and correlation functionals for comparison
	std::shared_ptr<ExactExchange> exx; //!< Exact exchange
	ElecInfo eInfo;
	ElecVars eVars;
	Energies ener;
	
	MinimizeParams elecMinParams; //!< electronic minimization parameters
	MinimizeParams ionicMinParams; //!< ionic minimization parameters
	MinimizeParams fluidMinParams; //!< fluid minimization parameters
	MinimizeParams latticeMinParams; //!< lattice minimization parameters
	MinimizeParams inverseKSminParams; //!< Inverse Kohn-sham minimization parameters
	SCFparams scfParams; //!< Self-consistent field mixing parameters
	
	CoulombParams coulombParams; //!< Coulomb truncation parameters
	std::shared_ptr<Coulomb> coulomb; //!< Coulomb interaction (optionally truncated)

	std::shared_ptr<VanDerWaals> vanDerWaals; //! Pair potential for vdw correction
	std::shared_ptr<Vibrations> vibrations; //! Vibrational mode calculator

	//! Call the setup/initialize routines of all the above in the necessray order
	void setup();
	void updateSupercell(bool force=false); //!< (re-)initialize coulombParams.supercell if necessary (or if forced)
};

#endif // JDFTX_ELECTRONIC_EVERYTHING_H
