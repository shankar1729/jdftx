/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Deniz Gunceler, Kendra Letchworth-Weaver

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

#include <commands/command.h>
#include <commands/ParamList.h>
#include <electronic/Everything.h>
#include <electronic/IonDynamicsParams.h>
#include <core/Units.h>

struct CommandIonicDynamics : public Command
{
	CommandIonicDynamics() : Command("ionic-dynamics", "jdftx/Ionic/Dynamics")
	{	format = "<time-step> <total-time> [<kT>=0.001] [<alpha>=0.0]";
		comments = "Applies Verlet molecular dynamics\n"
			   "Takes the time in femtoseconds and kT is in Hartree.\n"
			   "Give <alpha> if you want to the system to equilibrate with a heat bath at <kT>.\n"
			   "Also make sure to turn the symmetries off if the initial velocities don't satisfy the symmetry\n"
			   "conditions. Initial velocities will be assigned randomly if they are not given.";
		allowMultiple = false;
	}

	void process(ParamList& pl, Everything& e)
	{	IonDynamicsParams& idp = e.ionDynamicsParams;
		double dt; pl.get(dt, 2.0, "time-step", true);
		double tMax; pl.get(tMax, 0.0, "max-time", true);
		idp.dt = dt*fs;
		idp.tMax = tMax*fs;
		pl.get(idp.kT, 0.001, "kT", false);
		pl.get(idp.alpha, 0.0, "alpha", false);
	}

	void printStatus(Everything& e, int iRep)
	{	IonDynamicsParams& idp = e.ionDynamicsParams;
		logPrintf("%lg %lg %lg %lg", idp.dt/fs, idp.tMax/fs, idp.kT, idp.alpha);
	}
}
commandIonicDynamics;

EnumStringMap<ConfiningPotentialType> confiningPotentialTypeMap
(	ConfineNone, "None",
	ConfineLinear, "Linear",
	ConfineQuadratic, "Quadratic",
	ConfineCubic, "Cubic",
	ConfineSmoothLinear, "SmoothLinear"
);

struct CommandConfineDynamics : public Command
{
	CommandConfineDynamics() : Command("confine-dynamics", "jdftx/Ionic/Dynamics")
	{	format = "<type> <parameter0> [<parameter1> ...]";
		comments = "Apply an attractive gravitation-like potential of <type>=None|Linear|Quadratic|Cubic|SmoothLinear\n"
			   "to every atom. Use this command to keep small clusters of molecules together during MD.\n"
			   "\n"
		           "Force is proportional to the mass and pulls towards the origin.\n"
			   "For polynomial potentials the parameters start from the coefficient\n"
			   "of the largest degree term. The constant term does not exist.\n"
			   "\n"
			   "<type>=SmoothLinear is the integral of a smoothened step function (radially) for the force.\n"
			   "Its three parameters are\n"
			   "    <parameter0> acceleration far from the origin (in atomic units)\n"
			   "                 tiny values are typical ~ 1e-8\n"
			   "    <parameter1> width of the smoothing in Bohr\n"
			   "    <parameter2> the activation radius, at which the force turns on\n"
			   "                 roughly the radius of the cluster of molecules.\n";
		allowMultiple = false;
		hasDefault = false;
	}

	void process(ParamList& pl, Everything& e)
	{	IonDynamicsParams& idp = e.ionDynamicsParams;
		pl.get(idp.confineType, ConfineNone, confiningPotentialTypeMap, "type", true);
		int numOfParams;
		switch (idp.confineType)
		{	case ConfineLinear: numOfParams = 1; break;
			case ConfineQuadratic: numOfParams = 2; break;
			case ConfineCubic: numOfParams = 3; break;
			case ConfineSmoothLinear: numOfParams = 3; break;
			default: numOfParams = 0; break;
		}
		double p;
		for(int i = 0; i < numOfParams; ++i)
		{	ostringstream oss; oss << "parameter" << i;
			pl.get(p, 0.0, oss.str(), true);
			idp.confineParameters.push_back(p);
		}
	}

	void printStatus(Everything& e, int iRep)
	{	IonDynamicsParams& idp = e.ionDynamicsParams;
		logPrintf("%s", confiningPotentialTypeMap.getString(idp.confineType));
		for (auto p: idp.confineParameters)
		{	logPrintf(" %lg", p);
		}
	}
}
commandConfineDynamic;

EnumStringMap<DriftRemovalType> driftRemovalTypeMap
(	DriftNone, "None",
	DriftVelocity, "DriftVelocity",
	DriftMomentum, "AvgMomentum"
);

struct CommandNetDriftRemoval : public Command
{	
        CommandNetDriftRemoval() : Command("net-drift-removal", "jdftx/Ionic/Dynamics")
	{	format = "<scheme>=AvgMomentum";
		comments = "Deal with the net momentum per unit cell when initializing ionic dynamics.\n"
		           "Select <scheme>=AvgMomentum|DriftVelocity|None \n"
			   "DriftVelocity removes the net drift velocity of the center of mass from each atom in the unit cell, \n"
			   "AvgMomentum removes the average net momentum per atom from each atom in the unit cell, \n"
			   "None does not remove the momentum from the unit cell. \n"
		           "CAUTION: None may cause mean positions of your atoms to drift over time! \n";
		allowMultiple = false;
		hasDefault = false;
	}

	void process(ParamList& pl, Everything& e)
	{	IonDynamicsParams& idp = e.ionDynamicsParams;
	        pl.get(idp.driftType, DriftMomentum, driftRemovalTypeMap, "scheme", true);
	}

	void printStatus(Everything& e, int iRep)
	{	
		IonDynamicsParams& idp = e.ionDynamicsParams;
		logPrintf("%s", driftRemovalTypeMap.getString(idp.driftType));
	}
}
CommandNetDriftRemoval;
