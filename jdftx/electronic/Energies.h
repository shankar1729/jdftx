/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_ENERGIES_H
#define JDFTX_ELECTRONIC_ENERGIES_H

#include <electronic/common.h>

struct Energies
{
	double KE; //!<!< Kinetic energy
	double Enl; //!< Non-local pseudopotential energy
	double Eloc; //!< Local pseudopotential energy (includes electron-nuclear electrostatic piece)
	double EH; //!< Hartree energy (mean field electron-electron electrostatic interaction)
	double Eewald; //!< Ewald sum  (nuclear-nuclear electrostatic interaction)
	double A_diel; //!< Free energy of fluid + coupling
	double Eexternal; //!< Energy due to coupling with the external potential and charge
	double EXX; //!< exact exchange energy
	double Exc; //!< Exchange-correlation energy
	
	double Exc_core; //!< Exchange-correlation energy subtraction for partial cores
	double Epulay; //!< Pulay correction energy
	
	double Etot; //!< Total: sum of all above
	
	double TS; //!< Fillings entropy
	double F; //!< Helmholtz energy (Etot-TS)
	
	double muN; //!< Fixed electron chemical potential Lagrange multiplier
	double G; //!< Grand free energy (F-muN)
	
	double Eband; //!< band structure energy (tr Hsub)
	
	Energies(); //!< initialize all to 0.
	void updateTotals(); //!< Recompute total energies (Etot, F and G)
	void print(FILE* fp=globalLog) const; //!< print energies to logPrintf / stream
};

double relevantFreeEnergy(const Everything&); //!< Return the variational free-energy that should be considered for a given system
const char* relevantFreeEnergyName(const Everything&); //!< return the name of the relevant free energy

void print_Hsub_eigs(const Everything&); //!< Print eigenvalues to logPrintf


#endif // JDFTX_ELECTRONIC_ENERGIES_H
