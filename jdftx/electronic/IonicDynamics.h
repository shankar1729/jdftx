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
	int nAtomsTot; //!< total number of atoms
	int nDOF; //!< number of degrees of freedom
	double Mtot; //!< total mass of system
	IonicMinimizer imin; //!< Just to be able to call IonicMinimizer::step(). Doesn't minimize anything.
	bool nAccumNeeded; //!< Whether accumulated electron density is needed
	
	//Current thermodynamic properties:
	double KE; //!< current kinetic energy
	double PE; //!< current potential energy
	double T; //!< current temperature
	double p; //!< current pressure
	matrix3<> stress; //!< current stress tensor (includes kinetic stress whereas IonInfo::stress does not)
	
	//Utility functions
	void initializeVelocities(); //!< Initialize Maxwell-Boltzmann distribution of velocities
	IonicGradient getVelocities(); //!< Get Cartesian velocities from SpeciesInfo lattice-coordinate versions
	void setVelocities(const IonicGradient&); //!< Set SpeciesInfo velocities from Cartesian-coordinate version
	void computePressure(); //!< Update pressure and stress
	void computeKE(); //!< Update kinetic energy and temperature
	void computePE(IonicGradient& accel); //!< Update potential energy and set acceleration
	bool report(int iter, double t); //!< Report properties at current step
};

//! @}
#endif //JDFTX_ELECTRONIC_IONICDYNAMICS_H
