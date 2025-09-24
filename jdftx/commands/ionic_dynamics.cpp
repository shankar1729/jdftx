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
#include <electronic/IonicDynamicsParams.h>
#include <electronic/IonicGaussianPotential.h>
#include <core/Units.h>


EnumStringMap<IonicDynamicsParams::StatMethod> statMethodMap
(	IonicDynamicsParams::StatNone, "None",
	IonicDynamicsParams::Berendsen, "Berendsen",
	IonicDynamicsParams::NoseHoover, "NoseHoover"
);

//An enum entry for each configurable option of IonicDynamicsParams
enum IonicDynamicsParamsMember
{	IDPM_dt,
	IDPM_nSteps,
	IDPM_statMethod,
	IDPM_T0,
	IDPM_P0,
	IDPM_stress0,
	IDPM_tDampT,
	IDPM_tDampP,
	IDPM_chainLengthT,
	IDPM_chainLengthP,
	IDPM_B0,
	IDPM_Delim //!< delimiter to detect end of input
};

EnumStringMap<IonicDynamicsParamsMember> idpmMap
(	IDPM_dt, "dt",
	IDPM_nSteps, "nSteps",
	IDPM_statMethod, "statMethod",
	IDPM_T0, "T0",
	IDPM_P0, "P0",
	IDPM_stress0, "stress0",
	IDPM_tDampT, "tDampT",
	IDPM_tDampP, "tDampP",
	IDPM_chainLengthT, "chainLengthT",
	IDPM_chainLengthP, "chainLengthP",
	IDPM_B0, "B0"
);

EnumStringMap<IonicDynamicsParamsMember> idpmDescMap
(	IDPM_dt, "time step [fs]",
	IDPM_nSteps, "number of molecular dynamics steps",
	IDPM_statMethod, statMethodMap.optionList() + " (method for thermostat and/or barostat)",
	IDPM_T0, "thermostat temperature or temperature for initial velocities [Kelvin]",
	IDPM_P0, "barostat pressure [bar] (default NAN => no hydrostatic barostat)",
	IDPM_stress0, "barostat stress components xx yy zz yz zx xy [bar] (default NANs => no anisotropic barostat)",
	IDPM_tDampT, "thermostat damping time [fs]",
	IDPM_tDampP, "barostat damping time [fs]",
	IDPM_chainLengthT, "Nose-Hoover chain length for thermostat",
	IDPM_chainLengthP, "Nose-Hoover chain length for barostat",
	IDPM_B0, "Characteristic bulk modulus [bar] for Berendsen barostat (damping ~ B0 * tDampP)"
);

struct CommandIonicDynamics : public Command
{
	CommandIonicDynamics() : Command("ionic-dynamics", "jdftx/Ionic/Optimization")
	{	format = "<key1> <value1> <key2> <value2> ...";
		comments = "Born-Oppenheimer molecular dynamics, controlled by keys:"
		+ addDescriptions(idpmMap.optionList(), linkDescription(idpmMap, idpmDescMap))
		+ "\n\nAny number of these key-value pairs may be specified in any order.\n\n"
			"Note that nSteps must be non-zero to activate dynamics.\n"
			"Default mode is NVE; specify statMethod to add a thermostat or barostat";
	}

	void process(ParamList& pl, Everything& e)
	{	IonicDynamicsParams& idp = e.ionicDynParams;
		const double nanVal = NAN;
		while(true)
		{	IonicDynamicsParamsMember key;
			pl.get(key, IDPM_Delim, idpmMap, "key");
			switch(key)
			{	case IDPM_dt: pl.get(idp.dt, 1., "dt", true); idp.dt *= fs; break;
				case IDPM_nSteps: pl.get(idp.nSteps, 0, "nSteps", true); break;
				case IDPM_statMethod: pl.get(idp.statMethod, IonicDynamicsParams::StatNone, statMethodMap, "statMethod", true); break;
				case IDPM_T0: pl.get(idp.T0, nanVal, "T0", true); idp.T0 *= Kelvin; break;
				case IDPM_P0: pl.get(idp.P0, nanVal, "P0", true); idp.P0 *= Bar; break;
				case IDPM_stress0:
				{	pl.get(idp.stress0(0,0), nanVal, "stress0_xx", true);
					pl.get(idp.stress0(1,1), nanVal, "stress0_yy", true);
					pl.get(idp.stress0(2,2), nanVal, "stress0_zz", true);
					pl.get(idp.stress0(1,2), nanVal, "stress0_yz", true);
					pl.get(idp.stress0(2,0), nanVal, "stress0_zx", true);
					pl.get(idp.stress0(0,1), nanVal, "stress0_xy", true);
					idp.stress0(2,1) = idp.stress0(1,2);
					idp.stress0(0,2) = idp.stress0(2,0);
					idp.stress0(1,0) = idp.stress0(0,1);
					idp.stress0 *= Bar;
					break;
				}
				case IDPM_tDampT: pl.get(idp.tDampT, 50., "tDampT", true); idp.tDampT *= fs; break;
				case IDPM_tDampP: pl.get(idp.tDampP, 100., "tDampP", true); idp.tDampP *= fs; break;
				case IDPM_chainLengthT: pl.get(idp.chainLengthT, 3, "chainLengthT", true); break;
				case IDPM_chainLengthP: pl.get(idp.chainLengthP, 3, "chainLengthP", true); break;
				case IDPM_B0: pl.get(idp.B0, nanVal, "B0", true); idp.B0 *= Bar; break;
				case IDPM_Delim: 
					if((not std::isnan(idp.P0)) and (not std::isnan(trace(idp.stress0))))
						throw(string("Cannot specify both P0 (hydrostatic) and stress0 (anisotropic) barostats"));
					return; //end of input
			}
		}
	}

	void printStatus(Everything& e, int iRep)
	{	const IonicDynamicsParams& idp = e.ionicDynParams;
		logPrintf(" \\\n\tdt         %lg", idp.dt/fs);
		logPrintf(" \\\n\tnSteps     %d", idp.nSteps);
		logPrintf(" \\\n\tstatMethod %s", statMethodMap.getString(idp.statMethod));
		logPrintf(" \\\n\tT0         %lg", idp.T0/Kelvin);
		logPrintf(" \\\n\tP0         %lg", idp.P0/Bar);
		logPrintf(" \\\n\tstress0 %lg %lg %lg  %lg %lg %lg", 
			idp.stress0(0,0)/Bar, idp.stress0(1,1)/Bar, idp.stress0(2,2)/Bar,
			idp.stress0(1,2)/Bar, idp.stress0(2,0)/Bar, idp.stress0(0,1)/Bar);
		logPrintf(" \\\n\ttDampT       %lg", idp.tDampT/fs);
		logPrintf(" \\\n\ttDampP       %lg", idp.tDampP/fs);
		logPrintf(" \\\n\tchainLengthT %d", idp.chainLengthT);
		logPrintf(" \\\n\tchainLengthP %d", idp.chainLengthP);
		logPrintf(" \\\n\tB0           %lg", idp.B0/Bar);
	}
}
commandIonicDynamics;


struct CommandLjOverride : public Command
{
	CommandLjOverride() : Command("lj-override", "jdftx/Ionic/Optimization")
	{	format = "<rCut>";
		comments =
			"Replace electronic DFT by a Lennard-Jones only potential with cutoff <rCut> in Angstroms.\n"
			"This potential will use LJ parameters, sigma = 2^(5/6) R0 and epsilon = C6/(128 R0^6),\n"
			"where R0 and C6 are DFT-D2 parameters for each species (adjustable by command setVDW).\n"
			"This is not for real calculations, but a quick way to debug ionic-minimize,\n"
			"lattice-minimize or ionic-dynamics. Tip: set elec-cutoff to a small value to\n"
			"speed up electronic initialization (which is not bypassed for code simplicity).";
	}
	
	void process(ParamList& pl, Everything& e)
	{	pl.get(e.iInfo.ljOverride, 0., "rCut", true);
		e.iInfo.ljOverride *= Angstrom; //convert to atomic units
	}
	
	void printStatus(Everything& e, int iRep)
	{	logPrintf("%lg", e.iInfo.ljOverride/Angstrom);
	}
}
commandLjOverride;

//---- Thermostat / barostat velocities ----
struct CommandStatVelocity : public Command
{
	CommandStatVelocity(string statName) : Command(statName+"-velocity", "jdftx/Ionic/Optimization")
	{	format = "<v1> <v2> ...";
		comments =
			"Read "+statName+" internal velocities for continuing ionic dynamics.\n"
			"This command is automatically dumped with ionpos from dynamics simulations\n"
			"using Nose-Hoover chains that involve "+statName+" internal velocities.";
	}
	
	virtual diagMatrix& target(Everything& e)=0;
	
	void process(ParamList& pl, Everything& e)
	{	diagMatrix& stat = target(e);
		stat.clear();
		while(true)
		{	double v = NAN;
			pl.get(v, v, "v");
			if(std::isnan(v)) break;
			stat.push_back(v);
		}
	}
	
	void printStatus(Everything& e, int iRep)
	{	const diagMatrix& stat = target(e);
		for(const double& v: stat) logPrintf("%lg ", v);
	}
};

struct CommandThermostatVelocity : public CommandStatVelocity
{	CommandThermostatVelocity() : CommandStatVelocity("thermostat") {}
	diagMatrix& target(Everything& e) { return e.iInfo.thermostat; }
}
commandThermostatVelocity;

struct CommandBarostatVelocity : public CommandStatVelocity
{	CommandBarostatVelocity() : CommandStatVelocity("barostat")
	{	comments += "\n(The first six components are strain rate, while the rest are lattice thermostat velocities.)\n";
	}
	diagMatrix& target(Everything& e) { return e.iInfo.barostat; }
}
commandBarostatVelocity;


//----- External forces -----

EnumStringMap<IonicGaussianPotential::Geometry> igpGeometryMap
(	IonicGaussianPotential::Spherical, "Spherical",
	IonicGaussianPotential::Cylindrical, "Cylindrical",
	IonicGaussianPotential::Planar, "Planar"
);

struct CommandIonicGaussianPotential : public Command
{
	CommandIonicGaussianPotential() : Command("ionic-gaussian-potential", "jdftx/Ionic/Optimization")
	{	format = "<species> <U0> <sigma> <geometry>";
		comments = "Apply external potential and forces to ions of specified <species>.\n"
			"The potential is Gaussian with peak value <U0> Hartrees and standard\n"
			"deviation <sigma> bohrs, and is always centered at the origin. The symmetry\n"
			"of the potential is specified by <geometry>, which may be one of:\n"
			"+ Spherical: The potential is 0-D and confined near the origin.\n"
			"+ Cylindrical: The potential is 1-D and extends along the z-direction.\n"
			"+ Planar: The potential is 2-D and extends in the xy-plane.\n"
			"\n"
			"Note that the coordinates of the atoms are taken in minimum-image convention\n"
			"for the unit cell centered at the origin for the potential and force calculation.\n"
			"This command is intended primarily for applying perturbations in ionic dynamics.";
		allowMultiple = true;
		require("ion");
	}

	void process(ParamList& pl, Everything& e)
	{	IonicGaussianPotential igp;
		//Identify species:
		string spName;
		pl.get(spName, string(), "species", true);
		for(int iSpecies=0; iSpecies<int(e.iInfo.species.size()); iSpecies++)
			if(e.iInfo.species[iSpecies]->name == spName)
			{	igp.iSpecies = iSpecies;
				break;
			}
		if(igp.iSpecies < 0) //defaults to -1 in IonicGaussianPotential()
			throw string("<species> must match one of the atom-type names in the calculation");
		//Potential parameters:
		pl.get(igp.U0, 0.0, "U0", true);
		pl.get(igp.sigma, 0.0, "sigma", true);
		pl.get(igp.geometry, IonicGaussianPotential::Planar, igpGeometryMap, "geometry", true);
		vector3<> origin(0.0, 0.0, 0.0);
		string key;
		pl.get(key, string(), "origin", false);
		if(key == "origin")
		{	pl.get(origin[0], 0.0, "x", true);
			pl.get(origin[1], 0.0, "y", true);
			pl.get(origin[2], 0.0, "z", true);
		}
		igp.origin = origin;
		e.iInfo.ionicGaussianPotentials.push_back(igp);
	}

	void printStatus(Everything& e, int iRep)
	{	const IonicGaussianPotential& igp = e.iInfo.ionicGaussianPotentials[iRep];
		logPrintf("%s %lg %lg %s", e.iInfo.species[igp.iSpecies]->name.c_str(),
			igp.U0, igp.sigma, igpGeometryMap.getString(igp.geometry));
	}
}
commandIonicGaussianPotential;
