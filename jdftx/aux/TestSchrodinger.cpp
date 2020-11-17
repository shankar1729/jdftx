/*-------------------------------------------------------------------
Copyright 2020 Ravishankar Sundararaman

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

#include <electronic/Everything.h>
#include <electronic/ElecMinimizer.h>
#include <electronic/ColumnBundle.h>
#include <commands/parser.h>

inline void printUsageExit(const char* errP=0)
{	if(errP) logPrintf("\nError: invalid value for <%s>\n", errP);
	logPrintf("\nUsage: TestSchrodinger <L> <Ecut> <nBands>\n\n");
	logPrintf("Solve Schrodinger equation for the nBands lowest states in a unit\n");
	logPrintf("harmonic oscillator potential within a cubic box of specified\n");
	logPrintf("size, using a plane wave basis with specified energy cutoff.\n");
	logPrintf("  <L>:      box size in bohrs\n");
	logPrintf("  <Ecut>:   kinetic energy cutoff in Hartrees\n");
	logPrintf("  <nBands>: number of bands calculated\n\n");
	exit(0);
}

inline void setHarmonicPotential(int i, vector3<> r, vector3<> r0, double* V)
{	V[i] = 0.5 * (r-r0).length_squared();
}

int main(int argc, char** argv)
{	//Get and check commandline parameters:
	if(argc != 4) printUsageExit();
	string Lstr(argv[1]), EcutStr(argv[2]), nBandsStr(argv[3]);
	double L = 0., Ecut = 0.; int nBands = 0;
	if((sscanf(Lstr.c_str(), "%lf", &L) != 1) or (L <= 0.)) printUsageExit("L");
	if((sscanf(EcutStr.c_str(), "%lf", &Ecut) != 1) or (Ecut <= 0.)) printUsageExit("Ecut");
	if((sscanf(nBandsStr.c_str(), "%d", &nBands) != 1) or (nBands <= 0)) printUsageExit("nBands");
	
	//Prepare jdftx input file in memory:
	typedef std::pair<string,string> stringPair;
	std::vector<stringPair> input;
	input.push_back(stringPair("lattice", "Cubic "+Lstr));
	input.push_back(stringPair("elec-cutoff", EcutStr));
	input.push_back(stringPair("elec-n-bands", nBandsStr));
	input.push_back(stringPair("dump", "End None"));
	input.push_back(stringPair("wavefunction", "random"));

	//Initialize:
	Everything e;
	initSystem(argc, argv);
	parse(input, e);
	e.setup();
	
	//Set the harmonic oscillator potential (r^2/2):
	nullToZero(e.eVars.Vscloc[0], e.gInfo);
	vector3<> r0 = e.gInfo.R * vector3<>(0.5,0.5,0.5); //box center
	applyFunc_r(e.gInfo, setHarmonicPotential, r0, e.eVars.Vscloc[0]->data());
	e.eVars.Vscloc[0] *= e.gInfo.dV; //internal convention for Vscloc includes grid integration factor
	logPrintf("Initialization completed successfully at t[s]: %9.2lf\n\n", clock_sec());
	
	//Perform the fixed potential calculation:
	logPrintf("\n----------- Band structure minimization -------------\n"); logFlush();
	bandMinimize(e); // Do the band-structure minimization
	logPrintf("\nConverged band energies:\n");
	e.eVars.Hsub_eigs[0].print(globalLog, "BandEig: %19.15lf\n");
	
	//Cleanup:
	finalizeSystem();
	return 0;
}
