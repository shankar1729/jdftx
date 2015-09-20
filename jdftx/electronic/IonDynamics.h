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

#ifndef JDFTX_ELECTRONIC_IONDYNAMICS_H
#define JDFTX_ELECTRONIC_IONDYNAMICS_H

#include <electronic/common.h>
#include <electronic/IonicMinimizer.h>
#include <core/vector3.h>
#include <core/matrix3.h>
#include <vector>

struct IonDynamics
{	Everything& e;
	double initialPotentialEnergy;
	double kineticEnergy, potentialEnergy;
	double pressure, totalMomentumNorm;
	double totalMass; int numberOfAtoms;
	vector3<> totalMomentum;
	
	IonicMinimizer imin; //Just to be able to call IonicMinimizer::step(). Doesn't minimize anything.

	IonDynamics(Everything& e) : e(e), totalMass(0.0), numberOfAtoms(0), imin(IonicMinimizer(e)){};
	//Virtual functions from Minimizable:
	void step(const IonicGradient&,const double&);
	void velocitiesInit();
	double computeAcceleration(IonicGradient&); //!< Write acceleration into `accel` in cartesian coordinates and return relevant energy.
	void computeMomentum(); //!< Calculate totalMomentum and totalMomentumNorm.
	void computePressure(); //!< Updates pressure.
	void computeKineticEnergy();//!< Updates `kineticEnergy`
	bool report(double t);
	
	void run(); //!< Run the simulation
	
};

#endif //JDFTX_ELECTRONIC_IONDYNAMICS_H
