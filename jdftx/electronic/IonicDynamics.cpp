/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth-Weaver

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

#include <electronic/IonicMinimizer.h>
#include <electronic/IonInfo.h>
#include <electronic/Symmetries.h>
#include <electronic/Everything.h>
#include <electronic/ElecMinimizer.h>
#include <electronic/Dump.h>
#include <electronic/LatticeMinimizer.h>
#include <electronic/IonicDynamics.h>
#include <core/Random.h>
#include <core/BlasExtra.h>

IonicDynamics::IonicDynamics(Everything& e)
: e(e), imin(IonicMinimizer(e))
{
	logPrintf("---------- Ionic Dynamics -----------\n");
	
	//Initialize mass, atom count and check velocities:
	bool vInitNeeded = false;
	Mtot = 0.;
	nAtomsTot = 0;
	for(const auto& sp: e.iInfo.species)
	{	Mtot += (sp->mass * amu) * sp->atpos.size();
		nAtomsTot += sp->atpos.size();
		for(const vector3<>& vel: sp->velocities)
			if(isnan(vel.length_squared()))
		      vInitNeeded = true; //If any velocity missing, randomize them all
	}
	
	//Initialize velocities if necessary:
	if(vInitNeeded)
		initializeVelocities();
	
	computeKineticEnergy();
}

void IonicDynamics::initializeVelocities()
{
	//Initialize random velocities with |v|^2 ~ 1/M
	vector3<> pTot; //total momentum
	for(auto& sp: e.iInfo.species)
	{	double invsqrtM = 1./sqrt(sp->mass);
		for(vector3<>& vel: sp->velocities)
		{	for(int iDir=0; iDir<3; iDir++)
				vel[iDir] = Random::normal() * invsqrtM;
			pTot += sp->mass * vel;
		}
	}
	
	//Remove center of mass momentum:
	for(auto& sp: e.iInfo.species)
		for(vector3<>& vel: sp->velocities)
			vel -= pTot / (sp->mass * nAtomsTot);
	
	//Rescale to current temperature:
	computeKineticEnergy();
	double keRatio = (1.5*e.ionicDynParams.T0*nAtomsTot)/kineticEnergy;
	double velocityScaleFactor = sqrt(keRatio);
	for(auto& sp: e.iInfo.species)
		for(vector3<>& vel: sp->velocities)
			vel *= velocityScaleFactor;
}


double IonicDynamics::computeAcceleration(IonicGradient& accel)
{	if(not e.iInfo.checkPositions())
	{	e.iInfo.printPositions(globalLog);
		die("\nMD step caused pseudopotential core overlaps.\n");
	}
	
	//Initialize ion-dependent quantities at this position:
	e.iInfo.update(e.ener);

	//Minimize the system:
	elecFluidMinimize(e);
	
	//Calculate forces
	e.iInfo.ionicEnergyAndGrad(); //compute forces in lattice coordinates
	accel = e.gInfo.invRT * e.iInfo.forces; //forces in cartesian coordinates (not accel. yet)
	
	//Compute net acceleration for subtraction:
	vector3<> netAccel;
	for(unsigned sp=0; sp<accel.size(); sp++)
		for(const vector3<>& a: accel[sp])
			netAccel += a;
	netAccel *= (1./Mtot); // This is the net acceleration of the unit cell.

	for(unsigned sp=0; sp<accel.size(); sp++)
	{	double invM = 1./(e.iInfo.species[sp]->mass * amu);
		for(vector3<>& a: accel[sp])
			a = invM*a - netAccel; //convert to acceleration with mean removed
	}
	return relevantFreeEnergy(e);
}

void IonicDynamics::computeKineticEnergy()
{	kineticEnergy = 0.;
	for(const auto& sp: e.iInfo.species)
		for(vector3<>& vel: sp->velocities)
			kineticEnergy += (0.5 * sp->mass*amu) * e.gInfo.RTR.metric_length_squared(vel);
}

void IonicDynamics::computePressure()
{	pressure = e.iInfo.computeStress
		? (-1./3)*trace(e.iInfo.stress)
		: NAN; //pressure not available
}

bool IonicDynamics::report(int iter, double t)
{	logPrintf("\nVerletMD: Step: %3d  tMD[fs]: %9.2lf  Ekin: %8.4lf  Epot: %8.4lf  Etot: %8.4lf  P[Bar]: %8.4le  t[s]: %9.2lf\n",
		iter, t/fs, kineticEnergy, potentialEnergy, kineticEnergy + potentialEnergy, pressure/Bar, clock_sec());
	return imin.report(iter);
}

void IonicDynamics::step(const IonicGradient& accel, const double& dt)
{	const IonicDynamicsParams& idp = e.ionicDynParams;
	IonicGradient dpos;
	dpos.init(e.iInfo);
	//Rescale the velocities to track the temperature
	//Assumes that kinetic energy gives approximately the input temperature
	double averageKineticEnergy = 1.5 * nAtomsTot * idp.T0;
	double alpha = idp.dt / idp.tDampT;
	double scaleFactor = 1.0 + 2.0 * alpha*(averageKineticEnergy-kineticEnergy)/ kineticEnergy;
	// Prevent scaling from being too aggressive
	if (scaleFactor < 0.5) scaleFactor = 0.5;
	if (scaleFactor > 1.5) scaleFactor = 1.5;
	
	for(unsigned sp=0; sp < e.iInfo.species.size(); sp++) //initialize dpos with last step.
	{	SpeciesInfo& spInfo = *(e.iInfo.species[sp]);
		for(unsigned atom=0; atom<spInfo.atpos.size(); atom++)
		{	spInfo.velocities[atom] *= scaleFactor;
			dpos[sp][atom] = spInfo.velocities[atom]*dt;
		}
	}
	//given the previous dpos (v_{n-1}*dt) calculate new dpos (v_n*dt) 
	dpos = dpos + e.gInfo.invR*accel*(dt*dt); //accel is in cartesian, dpos in lattice
	
	//call the IonicMinimizer::step() that moves the wavefunctions and the nuclei
	imin.step(e.gInfo.R * dpos, 1.0); //takes its argument in cartesian coordinates
	//Update the velocities
	for(unsigned sp=0; sp < e.iInfo.species.size(); sp++)
	{	SpeciesInfo& spInfo = *(e.iInfo.species[sp]);
		for(unsigned atom=0; atom<spInfo.atpos.size(); atom++)
			spInfo.velocities[atom] = dpos[sp][atom]/dt;
		mpiWorld->bcastData(spInfo.atpos);
		spInfo.sync_atpos();
	}
}

void IonicDynamics::run()
{	const IonicDynamicsParams& idp = e.ionicDynParams;
	IonicGradient accel; //in cartesian coordinates
	
	accel.init(e.iInfo);
	nullToZero(e.eVars.nAccum, e.gInfo);
	
	for(int iter=0; iter<idp.nSteps; iter++)
	{	double t = iter*idp.dt;
		potentialEnergy = computeAcceleration(accel);
		computeKineticEnergy();
		computePressure();
			
		report(iter, t);
		step(accel, e.ionicDynParams.dt);
		
		//Accumulate the averaged electronic density over the trajectory
		for(unsigned s=0; s<e.eVars.nAccum.size(); s++)
		{	double fracNew = idp.dt/(t+idp.dt);
			e.eVars.nAccum[s] += fracNew*(e.eVars.n[s] - e.eVars.nAccum[s]);
		}
	}
}
