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
#include <core/EnergyComponents.h>

struct Energies
{
	EnergyComponents E; //!< All components of the internal energy (excluding TS and muN)
	
	double TS; //!< Fillings entropy
	double F() const { return double(E)-TS; } //!< Helmholtz energy (Etot-TS)
	
	double muN; //!< Fixed electron chemical potential Lagrange multiplier
	double G() const { return F()-muN; } //!< Grand free energy (F-muN)
	
	double Eband; //!< band structure energy (tr Hsub)
	
	Energies(); //!< initialize all to 0.
	void print(FILE* fp=globalLog) const; //!< print energies to logPrintf / stream
};

double relevantFreeEnergy(const Everything&); //!< Return the variational free-energy that should be considered for a given system
const char* relevantFreeEnergyName(const Everything&); //!< return the name of the relevant free energy

void print_Hsub_eigs(const Everything&); //!< Print eigenvalues to logPrintf


#endif // JDFTX_ELECTRONIC_ENERGIES_H
