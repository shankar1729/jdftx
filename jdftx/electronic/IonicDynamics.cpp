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
: e(e),
	statT(e.ionicDynParams.statMethod!=IonicDynamicsParams::StatNone),
	statP(statT and (not (std::isnan)(e.ionicDynParams.P0))),
	statStress(statT and (not (std::isnan)(trace(e.ionicDynParams.stress0)))),
	lmin(LatticeMinimizer(e, true, statP, statStress)), nAccumNeeded(false)
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
			if(std::isnan(vel.length_squared()))
		      vInitNeeded = true; //If any velocity missing, randomize them all
	}
	nDOF = e.ionicMinParams.nDim;
	if(nDOF == 3*nAtomsTot) nDOF -= 3; //discount overall translation
	if(not nDOF)
		die("No degrees of freedom for IonicDynamics.\n\n");
	
	//Thermostat initialization:
	const IonicDynamicsParams& idp = e.ionicDynParams;
	const int nDimThermo = (idp.statMethod==IonicDynamicsParams::NoseHoover)
		? (idp.chainLengthT + ((statP or statStress) ? idp.chainLengthP : 0))
		: 0;
	e.iInfo.thermostat.resize(nDimThermo, 0.);
	if(statStress) stressTarget = idp.stress0;
	if(statP) stressTarget = -idp.P0 * matrix3<>(1,1,1);
	assert(not (statStress and statP));
	
	//Initialize velocities if necessary:
	if(vInitNeeded)
		initializeVelocities();
	
	//Whether nAccum is needed:
	for(auto dumpPair: e.dump)
		if(dumpPair.second == DumpElecDensityAccum)
			nAccumNeeded = true;
}

void IonicDynamics::initializeVelocities()
{	//Initialize random velocities with |v|^2 ~ 1/M
	LatticeGradient vel; vel.init(e.iInfo);
	randomize(vel); //unit normal distribution
	for(size_t sp=0; sp<e.iInfo.species.size(); sp++)
	{	double invsqrtM = 1./sqrt(e.iInfo.species[sp]->mass);
		for(vector3<>& v: vel.ionic[sp]) v *= invsqrtM;
	}
	vel.thermostat = e.iInfo.thermostat;
	vel.barostat = e.iInfo.barostat;
	
	//Apply constraints:
	lmin.constrain(vel);
	setVelocities(vel);
	
	//Rescale to current temperature:
	computeKE();
	double keRatio = (0.5 * nDOF * e.ionicDynParams.T0)/KE;
	double velocityScaleFactor = sqrt(keRatio);
	for(auto& sp: e.iInfo.species)
		for(vector3<>& vel: sp->velocities)
			vel *= velocityScaleFactor;
}

//Get Cartesian velocities from SpeciesInfo lattice-coordinate versions
LatticeGradient IonicDynamics::getVelocities()
{	LatticeGradient v; v.init(e.iInfo);
	for(size_t sp=0; sp<e.iInfo.species.size(); sp++)
	{	const SpeciesInfo& spInfo = *(e.iInfo.species[sp]);
		for(size_t at=0; at<spInfo.atpos.size(); at++)
			v.ionic[sp][at] = e.gInfo.R * spInfo.velocities[at];
	}
	v.thermostat = e.iInfo.thermostat;
	v.barostat = e.iInfo.barostat;
	//Berendsen barostat (sets the strain rate directly):
	const IonicDynamicsParams& idp = e.ionicDynParams;
	if(idp.statMethod == IonicDynamicsParams::Berendsen and (statP or statStress))
		v.lattice = (1./(idp.tDampP*idp.B0)) * (stressTarget - stress);
	return v;
}

//Set SpeciesInfo velocities from Cartesian-coordinate version
void IonicDynamics::setVelocities(const LatticeGradient& v)
{	for(size_t sp=0; sp<e.iInfo.species.size(); sp++)
	{	SpeciesInfo& spInfo = *(e.iInfo.species[sp]);
		for(size_t at=0; at<spInfo.atpos.size(); at++)
			spInfo.velocities[at] = e.gInfo.invR * v.ionic[sp][at];
	}
	e.iInfo.thermostat = v.thermostat;
	e.iInfo.barostat = v.barostat;	
}


void IonicDynamics::computePressure()
{	if(e.iInfo.computeStress)
	{	stress = e.iInfo.stress;
		//Add kinetic stress
		for(size_t sp=0; sp<e.iInfo.species.size(); sp++)
		{	const SpeciesInfo& spInfo = *(e.iInfo.species[sp]);
			double prefac = (spInfo.mass*amu) / e.gInfo.detR;
			for(const vector3<>& vel: spInfo.velocities)
			{	vector3<> vCart = e.gInfo.R * vel;
				stress -= prefac * outer(vCart,vCart);
			}
		}
		//Compute pressure:
		p = (-1./3)*trace(stress);
	}
	else //stress and pressure not available
	{	stress = matrix3<>(NAN, NAN, NAN);
		p = NAN;
	}
}

void IonicDynamics::computeKE()
{	KE = 0.;
	for(const auto& sp: e.iInfo.species)
		for(vector3<>& vel: sp->velocities)
			KE += (0.5 * sp->mass*amu) * e.gInfo.RTR.metric_length_squared(vel);
	T = 2*KE/nDOF;
}

LatticeGradient IonicDynamics::computePE()
{	LatticeGradient gradUnused, accel; accel.init(e.iInfo);
	PE = lmin.compute(&gradUnused, &accel); //In dynamicsMode, Lattice/IonicMinimizer::compute replaces Kgrad with acceleration
	if(std::isnan(PE))
		die("\nIonicDynamics: step caused pseudopotential core overlap (try core-overlap-check none).\n\n");
	return accel;
}

LatticeGradient IonicDynamics::thermostat()
{	const IonicDynamicsParams& idp = e.ionicDynParams;
	//Update KE and pressure first:
	computeKE();
	computePressure();
	//Compute velocity-dependent forces (non-zero if thermostatting):
	LatticeGradient accelV; accelV.init(e.iInfo);
	if(statT)
	{	double relErrKE = KE/(0.5*nDOF*idp.T0) - 1.;
		double omegaDamp = 1./idp.tDampT;
		//Compute damping term depending on thermostat:
		double minusGammaDamp = 0.;
		switch(idp.statMethod)
		{	case IonicDynamicsParams::Berendsen:
			{	minusGammaDamp = -0.5*omegaDamp*relErrKE;
				break;
			}
			case IonicDynamicsParams::NoseHoover:
			{	minusGammaDamp = -e.iInfo.thermostat[0];
				double omegaDampSq = omegaDamp*omegaDamp;
				for(int j=0; j<idp.chainLengthT; j++)
				{	//Coupling to system / previous thermostat DOF:
					accelV.thermostat[j] = (j==0)
						? omegaDampSq * relErrKE
						: (j==1 ? nDOF : 1) * std::pow(e.iInfo.thermostat[j-1],2) - omegaDampSq;
					//Coupling to next thermostat DOF:
					if(j+1 < idp.chainLengthT)
						accelV.thermostat[j] -= e.iInfo.thermostat[j] * e.iInfo.thermostat[j+1];
				}
				break;
			}
			case IonicDynamicsParams::StatNone: break; //Never reached (just to suppress compiler warning)
		}
		//Add damping term to atom velocities
		for(size_t sp=0; sp<e.iInfo.species.size(); sp++)
		{	const SpeciesInfo& spInfo = *(e.iInfo.species[sp]);
			for(size_t at=0; at<spInfo.atpos.size(); at++)
			{	vector3<> vCart = e.gInfo.R * spInfo.velocities[at]; //current Cartesian velocity
				accelV.ionic[sp][at] = minusGammaDamp * vCart;
			}
		}
		lmin.constrain(accelV);
	}
	return accelV;
}


bool IonicDynamics::report(int iter, double t)
{	logPrintf("\nIonicDynamics: Step: %3d  PE: %10.6lf  KE: %10.6lf  T[K]: %8.3lf  P[Bar]: %8.4le  tMD[fs]: %9.2lf  t[s]: %9.2lf\n",
		iter, PE, KE, T/Kelvin, p/Bar, t/fs, clock_sec());
	if(e.iInfo.computeStress)
	{	logPrintf("\n# Stress tensor including kinetic terms in Cartesian coordinates [Eh/a0^3]:\n");
		stress.print(globalLog, "%12lg ", true, 1e-14);
	}
	return lmin.report(iter);
}

void IonicDynamics::run()
{	const IonicDynamicsParams& idp = e.ionicDynParams;
	
	//Initial energies and forces
	if(nAccumNeeded) nullToZero(e.eVars.nAccum, e.gInfo);
	LatticeGradient accel = computePE(), accelV = thermostat(); //in Cartesian coordinates
	
	for(int iter=0; iter<=idp.nSteps; iter++)
	{	double t = iter*idp.dt;
		report(iter, t);
		if(iter==idp.nSteps) break;
		
		//Velocity Verlet step:
		//--- velocity update: first half step
		LatticeGradient vel = getVelocities();
		axpy(0.5*idp.dt, accel+accelV, vel);
		//--- position and position-dependent acceleration update:
		lmin.step(vel, idp.dt);
		accel = computePE();
		//--- velocity update: second half step estimator
		axpy(0.5*idp.dt, accel+accelV, vel); //note second-order error here due to first-order error in accelV
		setVelocities(vel);
		//--- velocity update: second half step corrector
		LatticeGradient accelVnew = thermostat();
		axpy(0.5*idp.dt, accelVnew-accelV, vel); //corrects second-order error introduced in estimator step
		setVelocities(vel);
		accelV = thermostat();
		
		//Accumulate the averaged electronic density over the trajectory
		if(nAccumNeeded and (not e.iInfo.ljOverride))
			for(unsigned s=0; s<e.eVars.nAccum.size(); s++)
			{	double fracNew = idp.dt/(t+idp.dt);
				e.eVars.nAccum[s] += fracNew*(e.eVars.n[s] - e.eVars.nAccum[s]);
			}
	}
}
