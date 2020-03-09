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

#ifndef JDFTX_ELECTRONIC_IONICDYNAMICS_H
#define JDFTX_ELECTRONIC_IONICDYNAMICS_H

#include <electronic/IonicMinimizer.h>
#include <core/matrix3.h>

//! @addtogroup IonicSystem
//! @{

//! Ionic dynamics (AIMD)
class IonicDynamics
{
public:
	IonicDynamics(Everything& e); //!< Initialize dynamics
	void run(); //!< Run the simulation
private:
	Everything& e;
	int nAtomsTot; //!< total number of atoms in system
	double Mtot; //!< total mass of system
	double kineticEnergy; //!< current kinetic energy
	double potentialEnergy; //!< current potential energy
	double pressure; //!< current pressure
	
	IonicMinimizer imin; //Just to be able to call IonicMinimizer::step(). Doesn't minimize anything.

	//similar to the virtual functions of Minimizable:
	void step(const IonicGradient&, const double&);   //!< Given the acceleration, take a time step. Scale the velocities if heat bath exists
	double computeAcceleration(IonicGradient& accel); //!< Write acceleration into 'accel' in cartesian coordinates and return relevant energy.
	bool report(double t);
	
	//Utility functions
	void velocitiesInit();       //!< Randomly initialize the velocities assuming potential energy is zero (E_tot = KE = 3NkT at t=0)
	void computePressure();      //!< Updates pressure.
	void computeKineticEnergy(); //!< Updates `kineticEnergy`
};

//! @}
#endif //JDFTX_ELECTRONIC_IONICDYNAMICS_H
