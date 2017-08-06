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
#include <electronic/IonDynamics.h>
#include <core/Random.h>
#include <core/BlasExtra.h>

// Sigmoid-like (Fermi function upside down) confining gravitational acceleration, 
//   zero around origin constant far away
// Positive since it is the magnitude
// ps[0] -> g_far    ps[1] -> s     ps[2] -> r_0
// g(r) = g_far / (1 + e^(-u) where u=(r-r_0)/s
inline double g_smoothLinear(double r, const std::vector<double>& ps)
{	return ps[0] / (1.0 + exp(-(r - ps[2]) / ps[1]));}

// integral of sigmoid function (not potential energy, just potential)
inline double V_smoothLinear(double r, const std::vector<double>& ps)
{	return ps[0] * ps[1] * log(1.0 + exp((r-ps[2])/ps[1]));}

void IonDynamics::velocitiesInit()
{	bool velocitiesGiven=true;
	IonInfo& iInfo = e.iInfo;
	IonDynamicsParams& idp = e.ionDynamicsParams;
	for(unsigned sp=0; sp < iInfo.species.size(); sp++) 
	{	SpeciesInfo& spInfo = *(iInfo.species[sp]);
		totalMass += spInfo.mass*amu*spInfo.atpos.size();
		numberOfAtoms += spInfo.atpos.size();
		if (!spInfo.velocities.size())
		      velocitiesGiven=false;	//If any of the velocities was missing then randomize them all.
	}
	if (velocitiesGiven)
	{	computeMomentum(); 
		switch(idp.driftType)
		{  
		  case DriftVelocity: removeNetDriftVelocity(); break;
		  case DriftNone: break;
		  case DriftMomentum: 
	          default: removeNetAvgMomentum(); break;
		}		
		computeKineticEnergy();
		return;
	}
	double dt = idp.dt, kT = idp.kT;
	double v,theta,phi;
	for(unsigned sp=0; sp<iInfo.species.size(); sp++) // Initialize random velocities
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
	 //OK to always remove the net momentum of randomly initialized velocities using the default or specified scheme
	switch(idp.driftType)
		{  
		  case DriftVelocity: removeNetDriftVelocity(); break;
		  case DriftNone: 
		  case DriftMomentum: 
	          default: removeNetAvgMomentum(); break;
		}
	
	//Now our lattice should not have an overall momentum
	//We can scale the speeds to give us the right temperature.
	computeKineticEnergy();
	double energyRatio = (3.0*kT)/(kineticEnergy/numberOfAtoms);
	double velocityScaleFactor = sqrt(energyRatio);
	for(unsigned sp=0; sp<iInfo.species.size(); sp++) // Scale the velocities
	{	SpeciesInfo& spInfo = *(iInfo.species[sp]);
		for(unsigned atom=0; atom<spInfo.velocities.size(); atom++)
			spInfo.velocities[atom] *= velocityScaleFactor;
	}
	computeKineticEnergy();
	//check if the momentum sums up to zero and we have the correct energy
	computeMomentum();
	if(idp.driftType!=DriftNone) assert(e.gInfo.RTR.metric_length_squared(totalMomentum) < 1.0e-8);
	assert( std::abs(kineticEnergy-3.0*kT*numberOfAtoms) < 1.0e-8);
}

double IonDynamics::computeAcceleration(IonicGradient& accel)
{	if(not e.iInfo.checkPositions())
	{	e.iInfo.printPositions(globalLog);
		die("\nMD step caused pseudopotential core overlaps.\n");
	}
	
	//Initialize ion-dependent quantities at this position:
	e.iInfo.update(e.ener);

	//Minimize the system:
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
	//accel is the acceleration in cartesian coordinates now.
	
	//add virtual confining gravitation-like acceleration, calculation done in cartesian coordinates
	double virtualPotentialEnergy = 0.0;
	const auto& confineType = e.ionDynamicsParams.confineType;
	if (confineType != ConfineNone)
	{	const auto& ps = e.ionDynamicsParams.confineParameters;
		auto spuriousPressure = 0.0;
		bool isPolynomial = (confineType == ConfineLinear    || 
				     confineType == ConfineQuadratic || 
				     confineType == ConfineCubic);
		auto V0 = 0.0;	if (confineType == ConfineSmoothLinear) V0 = V_smoothLinear(0.0, ps);
		for(unsigned sp=0; sp<accel.size(); sp++)
		{	SpeciesInfo& spInfo = *(e.iInfo.species[sp]);
			double mass = spInfo.mass*amu;
			for(unsigned atom=0; atom<accel[sp].size(); atom++)
			{	auto r = e.gInfo.R * spInfo.atpos[atom];
				auto rNorm = r.length();
				auto rHat = r/rNorm;
				auto force = 0.0;
				assert(fabs(rHat.length_squared() - 1.0) < 1e-7);
				if (isPolynomial)
				{	auto V = 0.0;
					auto g = 0.0, gtmp = 0.0; // Ftmp is the helper that holds the terms while taking the derivative inside for loop if polynomial
					for (auto p: ps)
					{	V += p; V *= rNorm;
						g += gtmp; g *= rNorm;  g += p; gtmp *= rNorm; gtmp += p; //derivative of V
					}
					virtualPotentialEnergy += mass*V;
					force = mass*g;
					accel[sp][atom] += (-g*rHat);
				}
				else if (confineType == ConfineSmoothLinear)
				{	virtualPotentialEnergy += mass*(V_smoothLinear(rNorm, ps) - V0);
					auto a = g_smoothLinear(rNorm, ps);
					force = mass * a;
					accel[sp][atom] += (-a*rHat);
				}
				spuriousPressure += force/(4*M_PI*rNorm*rNorm);
				logPrintf("SpuriousInwardForce %s %19.15lf at radius %19.15lf of mass %19.15lf\n", 
					  spInfo.name.c_str(), force, rNorm, mass);
			}
		}
		
		logPrintf("ConfiningPotentialEnergy = %lg   SpuriousPressure = %lg\n", virtualPotentialEnergy, spuriousPressure);
	}
	
	return relevantFreeEnergy(e) + virtualPotentialEnergy;
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
{	int iter = lrint(t/e.ionDynamicsParams.dt); //round to int to get iteration number
	//Dump:
	e.dump(DumpFreq_Dynamics, iter);
	
	logPrintf("\nVerletMD t = %f fs (dt = %f in atomic units) Iter: %d",t/fs,e.ionDynamicsParams.dt, iter);
	logPrintf("\nE_Kin = %lg \t E_pot = %lg \t E_tot = %lg \t pressure(in Bar) = %lg Momentum = %lg",
		  kineticEnergy, potentialEnergy - initialPotentialEnergy, 
		  kineticEnergy + potentialEnergy - initialPotentialEnergy, pressure/Bar, totalMomentumNorm);
	logPrintf("\n"); e.iInfo.printPositions(globalLog);
	logPrintf("\n"); e.iInfo.forces.print(e, globalLog);
	logPrintf("# Energy components:\n"); e.ener.print(); logPrintf("\n");
	return false;
}

void IonDynamics::step(const IonicGradient& accel, const double& dt)
{	IonicGradient dpos;
	dpos.init(e.iInfo);
	//Rescale the velocities to track the temperature
	//Assumes that kinetic energy gives approximately the input temperature
	double averageKineticEnergy = 1.5 * numberOfAtoms * e.ionDynamicsParams.kT;
	double scaleFactor = 1.0 + 2.0 * e.ionDynamicsParams.alpha*(averageKineticEnergy-kineticEnergy)/
								kineticEnergy;
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
		
		mpiUtil->bcast((double*)spInfo.atpos.data(), 3*spInfo.atpos.size());
		spInfo.sync_atpos();
	}
}

void IonDynamics::run()
{	IonicGradient accel; //in cartesian coordinates
	//Initialize old positions according to the temperature
	velocitiesInit();
	
	if (e.ionDynamicsParams.confineType != ConfineNone)
	{	centerOfMassToOrigin();
	}
	
	accel.init(e.iInfo);
	initialPotentialEnergy = (double)NAN; // ground state potential
	nullToZero(e.eVars.nAccumulated,e.gInfo);
	
	for(double t=0.0; t<e.ionDynamicsParams.tMax; t+=e.ionDynamicsParams.dt)
	{	potentialEnergy = computeAcceleration(accel);
		computeMomentum();computeKineticEnergy();computePressure();
			
		if (e.ionDynamicsParams.confineType == ConfineNone) assert(totalMomentumNorm<1e-7);
		
		if (std::isnan(initialPotentialEnergy ))
			initialPotentialEnergy = potentialEnergy;
		report(t);
		step(accel, e.ionDynamicsParams.dt);
		//Accumulate the averaged electronic density over the trajectory
		for(unsigned s=0; s<e.eVars.nAccumulated.size(); s++)
		{
		  e.eVars.nAccumulated[s]=(e.eVars.nAccumulated[s]*t+e.eVars.n[s]*e.ionDynamicsParams.dt)*(1.0/(t+e.ionDynamicsParams.dt));
		}
	}
}

void IonDynamics::removeNetDriftVelocity()  
{	vector3<> averageDrift = totalMomentum / totalMass;
	//Subtract average drift velocity of center of mass from the individual velocities
	for(auto& spInfo : e.iInfo.species)
	{	for(auto& v : spInfo->velocities)
			v -= averageDrift;
	}
}

void IonDynamics::removeNetAvgMomentum()
{	vector3<> averageMomentum = totalMomentum / numberOfAtoms;
	//Subtract average momentum from the individual momentums
	for(auto& spInfo : e.iInfo.species)
	{	for(unsigned atom=0; atom<spInfo->velocities.size(); atom++)
			spInfo->velocities[atom] -= averageMomentum / (spInfo->mass*amu);
	}

}


void IonDynamics::centerOfMassToOrigin()
{	//Calculate the weighted average of positions
	vector3<> centerOfMass;
	for(const auto& spInfo : e.iInfo.species)
	{	auto w = (spInfo->mass) * amu / totalMass;
		for(const auto& pos : spInfo->atpos)
			centerOfMass += w * pos;
	}
	
	//Move everything
	for(auto& spInfo : e.iInfo.species)
	{	for(auto& pos : spInfo->atpos)
			pos -= centerOfMass;
		mpiUtil->bcast((double*)spInfo->atpos.data(), 3*spInfo->atpos.size());
		spInfo->sync_atpos();
	}
	
	
	logPrintf("Atoms are shifted to make the center of mass at origin\n");
}