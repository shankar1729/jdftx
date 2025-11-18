/*-------------------------------------------------------------------
Copyright 2025 Andrew Diggs, Ravishankar Sundararaman

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

#include <electronic/Metadynamics.h>
#include <commands/command.h>

void MetadynamicsBond::initialize()
{	logPrintf("\n--- Initializing bond metadynamics ---\n");
	if(state_filename.length() and isReadable(state_filename))
	{	potential.clear();
		FILE* fp = fopen(state_filename.c_str(), "r");
		while(!feof(fp))
		{	double value;
			if(fscanf(fp, "%lg", &value) == 1)
				potential.push_back(value);
		}
		fclose(fp);
		logPrintf("Read bias potential on grid of length %lu from %s\n", potential.size(), state_filename.c_str());
	}
}

void MetadynamicsBond::save()
{	if(state_filename.length() and mpiWorld->isHead())
	{	logPrintf("Dumping '%s' ... ", state_filename.c_str()); logFlush();
		FILE* fp = fopen(state_filename.c_str(), "w");
		potential.print(fp, "%.15lg\n");
		fclose(fp);
		logPrintf("done.\n"); logFlush();
	}
}

inline Atom& selectAtom(std::vector<Atom>& atoms, int iSp, int at)
{	int atCount = 0;
	for(Atom& atom: atoms)
		if(atom.sp == iSp)
		{	if(atCount == at)
				return atom;
			atCount++;
		}
	die("Invalid atom specification: iSp = %d and at = %d\n", iSp, at);
}

double MetadynamicsBond::energyAndGrad(const GridInfo& gInfo, std::vector<Atom>& atoms)
{	//Compute current bond length:
	Atom& atom1 = selectAtom(atoms, iSp1, at1);
	Atom& atom2 = selectAtom(atoms, iSp2, at2);
	vector3<> dpos = atom1.pos - atom2.pos; //fractional coordinates
	dpos -= vector3<>(round(dpos)); //wrap to minimum image convention
	double bondLength = sqrt(gInfo.RTR.metric_length_squared(dpos));
	
	//Update bias potential:
	size_t min_size = size_t(std::ceil(bondLength / resolution)) + 5;
	if(potential.size() < min_size) potential.resize(min_size, 0.0); //pad bias potential as needed
	double prefactor = energy_per_step / (sqrt(2 * M_PI) * resolution);
	double r = 0.0;
	for(double& potential_i: potential)
	{	potential_i += prefactor * exp(-0.5 * std::pow((r - bondLength) / resolution, 2));
		r += resolution;
	}
	
	//Return bias potential energy and accumulate corresponding force:
	double x = bondLength / resolution; //bond length in grid units,
	int i = floor(x); //... separated into integer
	double t = x - i; //... and fractional parts
	double Ebias = potential[i] * (1.0 - t) + potential[i+1] * t;
	double Ebias_bondLength = (potential[i+1] - potential[i]) / resolution;
	vector3<> Ebias_dpos = (Ebias_bondLength / bondLength) * (gInfo.RTR * dpos);
	atom1.force -= Ebias_dpos; //negative gradient in contravariant lattice coordinates
	atom2.force += Ebias_dpos; //negative gradient in contravariant lattice coordinates
	return Ebias;
}
