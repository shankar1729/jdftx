/*-------------------------------------------------------------------
Copyright 2015 Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_ELECTRONSCATTERING_H
#define JDFTX_ELECTRONIC_ELECTRONSCATTERING_H

#include <electronic/common.h>

struct ElectronScattering
{
	double eta; //!< frequency resolution and half the imaginary part ascribed to probe frequency (set to eInfo.kT if unspecified)
	double Ecut; //!< energy cut-off for dielectric matrices (set to cntrl.Ecut if unspecified)
	double fCut; //!< threshold for considering states fully occupied / unoccupied (default: 1e-6)
	double omegaMax; //!< maximum energy transfer to account for and hence maximum frequency in dielectric grid (if zero, autodetermine from available eigenvalues)
	
	ElectronScattering();
	void dump(const Everything& e); //!< compute and dump Im(Sigma_ee) for each eigenstate
};

#endif //JDFTX_ELECTRONIC_ELECTRONSCATTERING_H
