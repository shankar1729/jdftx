/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_VIBRATIONS_H
#define JDFTX_ELECTRONIC_VIBRATIONS_H

#include <electronic/common.h>
#include <core/vector3.h>
#include <core/VectorField.h>

class Vibrations
{
public:
	double dr; //!< Perturbation amplitude
	bool centralDiff; //!< whether to use central difference derivatives
	bool useConstraints; //!< whether to use ion constraints to restrict vibrational modes
	bool translationSym; //!< whether to use translation symmetry to optimize force calculations
	bool rotationSym; //!< whether to project out rotational modes: valid only for molecules
	double omegaMin; //!< frequency cutoff for free energy calculation and detailed mode print out
	double T; //!< ionic temperature used for entropy and free energy estimation
	double omegaResolution; //!< frequency resolution used for identifying and reporting degeneracies
	
	Vibrations();
	void setup(Everything* e);
	void calculate();
	
private:
	Everything* e;
	vector3<> getSplit() const; //get optimum latttice coordinates for splitting periodicity in a molecular geometry
	IonicGradient getCMcoords() const; //get cartesian coordinates of all atoms relative to molecule center of mass
	VectorField Ptest; //vector field that measures dipole moment in lattice coordinates
	vector3<> getPel() const; //get electronic dipole moment at current state in cartesian coordinates
};

#endif //JDFTX_ELECTRONIC_VIBRATIONS_H
