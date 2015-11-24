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

#include <electronic/IonicMinimizer.h>
#include <electronic/IonInfo.h>
#include <electronic/Symmetries.h>
#include <electronic/Everything.h>
#include <electronic/ElecMinimizer.h>
#include <electronic/Dump.h>
#include <electronic/LatticeMinimizer.h>
#include <electronic/operators.h>
#include <electronic/IonDynamics.h>
#include <core/Random.h>
#include <core/BlasExtra.h>

void IonDynamics::velocitiesInit()
{	bool velocitiesGiven=true;
	IonInfo& iInfo = e.iInfo;
	for(unsigned sp=0; sp < iInfo.species.size(); sp++) 
	{	SpeciesInfo& spInfo = *(iInfo.species[sp]);
		totalMass += spInfo.mass*amu*spInfo.atpos.size();
		numberOfAtoms += spInfo.atpos.size();
		if (!spInfo.velocities.size())
		      velocitiesGiven=false;	//If any of the velocities was missing then randomize them all.
	}
	if (velocitiesGiven)
	{	computeMomentum(); 
		computeKineticEnergy();
		return;
	}
	double dt = e.verletParams.dt, kT = e.verletParams.kT;
	double v,theta,phi;
	for(unsigned sp=0; sp<e.iInfo.species.size(); sp++) // Initialize random velocities
	{	SpeciesInfo& spInfo = *(iInfo.species[sp]);
		for(unsigned atom=0; atom<spInfo.atpos.size(); atom++)
		{	v = Random::uniform(0.0,0.1);
			theta=Random::uniform(0,M_PI), phi = Random::uniform(0,2*M_PI);
			vector3<> vel(v*dt*sin(theta)*sin(phi),v*dt*sin(theta)*cos(phi),v*dt*cos(theta));
			spInfo.velocities.push_back(vel);
		}
	}
	computeMomentum();
	logPrintf("----------Ion Dynamics-----------\ndensity = %lg (in atomic units)\n",totalMass/e.gInfo.detR);
	vector3<> averageDrift = totalMomentum / totalMass;
	//Subtract average momentum from the individual momentums
	for(unsigned sp=0; sp<e.iInfo.species.size(); sp++)
	{	SpeciesInfo& spInfo = *(iInfo.species[sp]);
		for(unsigned atom=0; atom<spInfo.velocities.size(); atom++)
			spInfo.velocities[atom] -= averageDrift;
	}
	
	//Now our lattice does not have an overall momentum
	//We can scale the speeds to be give us the right temperature.
	computeKineticEnergy();
	double energyRatio = (3.0*kT)/(kineticEnergy/numberOfAtoms);
	double velocityScaleFactor = sqrt(energyRatio);
	for(unsigned sp=0; sp<e.iInfo.species.size(); sp++) // Scale the velocities
	{	SpeciesInfo& spInfo = *(iInfo.species[sp]);
		for(unsigned atom=0; atom<spInfo.velocities.size(); atom++)
			spInfo.velocities[atom] *= velocityScaleFactor;
	}
	computeKineticEnergy();
	//check if the momentum sums up to zero and we have the correct energy
	computeMomentum();
	if( e.gInfo.RTR.metric_length_squared(totalMomentum)>1.0e-8)
	{	logPrintf("\nPROBLEM! Lattice is drifting");
		die("\nBye!");
		return;
	}
	if( std::abs(kineticEnergy-3.0*kT*numberOfAtoms)>1.0e-8)
	{	logPrintf("\nPROBLEM! Energy is not initialized correctly.");
		die("\nBye!");
		return;
	}
}

double IonDynamics::computeAcceleration(IonicGradient& accel)
{	if(not e.iInfo.checkPositions())
	{	e.iInfo.printPositions(globalLog);
		die("\nMD step caused pseudopotential core overlaps.\n");
	}
	
	//Initialize ion-dependent quantities at this position:
	e.iInfo.update(e.ener);

	//Minimize the electronic system:
	elecFluidMinimize(e);
	
	//Calculate forces
	e.iInfo.ionicEnergyAndGrad(e.iInfo.forces); //compute forces in lattice coordinates
	accel = e.gInfo.invRT * e.iInfo.forces; //forces in cartesian coordinates (not accel. yet)
	vector3<> netAccel;
	for(unsigned sp=0; sp<accel.size(); sp++)
	{	for(unsigned atom=0; atom<accel[sp].size(); atom++)
			netAccel += accel[sp][atom];
	}
	netAccel *= (1.0/totalMass); // This is the net acceleration of the unit cell.

	for(unsigned sp=0; sp<accel.size(); sp++)
	{	SpeciesInfo& spInfo = *(e.iInfo.species[sp]);
		for(unsigned atom=0; atom<accel[sp].size(); atom++)
		{	accel[sp][atom] *= 1.0/(spInfo.mass*amu);  //Divide force by mass to get the acceleration
			accel[sp][atom] -= netAccel; // Subtract the net acceleration to cancel the drift.
		}
	}
	return relevantFreeEnergy(e);
}

void IonDynamics::computeMomentum()
{	vector3<> p(0.0,0.0,0.0);
	double mass;
	for(unsigned sp=0; sp<e.iInfo.species.size(); sp++)
	{	mass = (e.iInfo.species[sp])->mass*amu;
		SpeciesInfo& spInfo = *(e.iInfo.species[sp]);
		for(unsigned atom=0; atom<spInfo.velocities.size(); atom++)
			p = p + spInfo.velocities[atom]*mass;
	}
	totalMomentumNorm = sqrt( e.gInfo.RTR.metric_length_squared(p) ); //in cartesian
	totalMomentum = p; //in lattice
}

void IonDynamics::computeKineticEnergy()
{	kineticEnergy=0.0;
	double mass;
	for(unsigned sp=0; sp<e.iInfo.species.size(); sp++)
	{	mass = (e.iInfo.species[sp])->mass*amu;
		SpeciesInfo& spInfo = *(e.iInfo.species[sp]);
		for(unsigned atom=0; atom<spInfo.velocities.size(); atom++)
			kineticEnergy += mass*e.gInfo.RTR.metric_length_squared(spInfo.velocities[atom])/2.0; // KE += m * v^2 / 2
	}
}

void IonDynamics::computePressure() // <!  A very similar code can be found in LatticeMinimize.cpp
{	double h = 1.0e-5; //Magic number from LatticeMinimize
	matrix3<> Rorig = e.gInfo.R;
	double V_0 = e.gInfo.detR;
	matrix3<> direction(1.0,1.0,1.0); // identity
	e.gInfo.R = Rorig + Rorig*(-2*h*direction);
	LatticeMinimizer::updateLatticeDependent(e);
	const double En2h = relevantFreeEnergy(e);
	
	e.gInfo.R = Rorig + Rorig*(-h*direction);
	LatticeMinimizer::updateLatticeDependent(e);
	const double Enh = relevantFreeEnergy(e);
	
	e.gInfo.R = Rorig + Rorig*(h*direction);
	LatticeMinimizer::updateLatticeDependent(e);
	const double Eph = relevantFreeEnergy(e);
	
	e.gInfo.R = Rorig + Rorig*(2*h*direction);
	LatticeMinimizer::updateLatticeDependent(e);
	const double Ep2h = relevantFreeEnergy(e);
	
	double centralDifference = (1./(12.*h))*(En2h - 8.*Enh + 8.*Eph - Ep2h); 
	/*
	We want dE/dV but with this centralDifference we have
	dE/dx as E=E(V(x)) and V(x) = V_0*(1+x)^3
		
		dE/dV = (dE/dx)  / (dV/dx)
		
		dV/dx = 3*V_0*(1+x)^2
		
	we are computing the derivative at x=0 so:
		dE/dV = (dE/dx) / (3*V_0)
	and
		P = -dE/dV
	*/
	pressure = -centralDifference/(3*V_0);
	e.gInfo.R = Rorig; // Reset R
	LatticeMinimizer::updateLatticeDependent(e);
}

bool IonDynamics::report(double t)
{	logPrintf("\nVerletMD t = %f fs (dt = %f in atomic units)",t/fs,e.verletParams.dt);
	logPrintf("\nE_Kin = %lg \t E_pot = %lg \t E_tot = %lg \t pressure(in Bar) = %lg",
		  kineticEnergy,potentialEnergy - initialPotentialEnergy, 
		  kineticEnergy + potentialEnergy - initialPotentialEnergy, pressure/Bar);
	logPrintf("\n"); e.iInfo.printPositions(globalLog);
	logPrintf("\n"); e.iInfo.forces.print(e, globalLog);
	logPrintf("# Energy components:\n"); e.ener.print(); logPrintf("\n");
	return false;
}

void IonDynamics::step(const IonicGradient& accel,const double& dt)
{	IonicGradient dpos;
	dpos.init(e.iInfo);
	
	//Rescale the velocities to track the temperature
	//Assumes that kinetic energy gives approximately the input temperature
	double averageKineticEnergy = 1.5 * numberOfAtoms * e.verletParams.kT;
	double scaleFactor = 1.0 + 2.0 * e.verletParams.alpha*(averageKineticEnergy-kineticEnergy)/
								kineticEnergy;
	if (scaleFactor < 0.5)
		scaleFactor = 0.5;
	else if(scaleFactor > 1.5)
		scaleFactor = 1.5;
	
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
		
		mpiUtil->bcast((double*)spInfo.atpos.data(), 3*spInfo.atpos.size());
		spInfo.sync_atpos();
	}
}

void IonDynamics::run()
{	IonicGradient accel; //in lattice coordinates x_n-x_{n-1}
	//Initialize old positions according to the temperature
	velocitiesInit();
	accel.init(e.iInfo);
	initialPotentialEnergy = (double)NAN; // ground state potential
	
	for(double t=0.0; t<e.verletParams.tMax; t+=e.verletParams.dt)
	{	potentialEnergy = computeAcceleration(accel);
		computeMomentum();computeKineticEnergy();computePressure();
			
		if (std::isnan(initialPotentialEnergy ))
			initialPotentialEnergy = potentialEnergy;
		report(t);
		step(accel, e.verletParams.dt);	
	}
}
