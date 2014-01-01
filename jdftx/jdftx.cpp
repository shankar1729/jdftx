/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver

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

#include <cstdio>
#include <cmath>
#include <ctime>
#include <getopt.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/matrix.h>
#include <electronic/Dump.h>
#include <electronic/ElecMinimizer.h>
#include <electronic/LatticeMinimizer.h>
#include <electronic/InverseKohnSham.h>
#include <electronic/Vibrations.h>
#include <fluid/FluidSolver.h>
#include <core/Util.h>
#include <commands/parser.h>

//Print usage information
void printUsage(const char *name)
{	printVersionBanner();
	logPrintf("Usage: %s [options]\n",name);
	logPrintf("\n\tPerforms Joint Density Functional Theory calculations.\n\n");
	logPrintf("options:\n\n");
	logPrintf("\t-h --help               help (this output)\n");
	logPrintf("\t-v --version            version\n");
	logPrintf("\t-i --input <filename>   specify command input file, default = stdin\n");
	logPrintf("\t-o --output <filename>  specify output log file, default = stdout\n");
	logPrintf("\t-d --no-append          overwrite output file instead of appending\n");
	logPrintf("\t-t --template           print an input file template\n");
	logPrintf("\t-m --mpi-debug-log      write output from secondary MPI processes to jdftx.<proc>.mpiDebugLog (instead of /dev/null)\n");
	logPrintf("\t-n --dry-run            quit after initialization (to verify commands and other input files)\n");
	logPrintf("\t-c --cores              number of cores to use (ignored when launched using SLURM)\n");
	logPrintf("\t-s --skip-defaults      skip printing status of default commands issued automatically.\n");
	logPrintf("\n");
}

//Program entry point
int main(int argc, char** argv, char** argp)
{	mpiUtil = new MPIUtil(argc, argv);
	Everything e; //the parent data structure for, well, everything
	
	//Parse command line:
	string inputFilename, logFilename; bool appendOutput=true, dryRun=false, printDefaults=true;
	option long_options[] =
		{	{"help", no_argument, 0, 'h'},
			{"version", no_argument, 0, 'v'},
			{"input",  required_argument, 0, 'i'},
			{"output", required_argument, 0, 'o'},
			{"no-append", no_argument, 0, 'd'},
			{"template", no_argument, 0, 't'},
			{"mpi-debug-log", no_argument, 0, 'm'},
			{"dry-run", no_argument, 0, 'n'},
			{"cores", required_argument, 0, 'c'},
			{"skip-defaults", no_argument, 0, 's'},
			{0, 0, 0, 0}
		};
	while (1)
	{	int c = getopt_long(argc, argv, "hvi:o:dtmnc:s", long_options, 0);
		if (c == -1) break; //end of options
		#define RUN_HEAD(code) if(mpiUtil->isHead()) { code } delete mpiUtil;
		switch (c)
		{	case 'v': RUN_HEAD( printVersionBanner(); ) return 0;
			case 'h': RUN_HEAD( printUsage(argv[0]); ) return 0;
			case 'i': inputFilename.assign(optarg); break;
			case 'o': logFilename.assign(optarg); break;
			case 'd': appendOutput=false; break;
			case 't': RUN_HEAD( printDefaultTemplate(e); ) return 0;
			case 'm': mpiDebugLog=true; break;
			case 'n': dryRun=true; break;
			case 'c':
			{	int nCores = 0;
				if(sscanf(optarg, "%d", &nCores)==1 && nCores>0)
					nProcsAvailable=nCores;
				break;
			}
			case 's': printDefaults=false; break;
			default: RUN_HEAD( printUsage(argv[0]); ) return 1;
		}
		#undef RUN_HEAD
	}
	
	// Open the logfile (if any):
	if(logFilename.length())
	{	globalLog = fopen(logFilename.c_str(), appendOutput ? "a" : "w");
		if(!globalLog)
		{	globalLog = stdout;
			logPrintf("WARNING: Could not open log file '%s' for writing, using standard output.\n", logFilename.c_str());
		}
	}

	//Print banners, setup threads, GPUs and signal handlers
	initSystem(argc, argv);
	
	//Parse input file and setup
	ElecVars& eVars = e.eVars;
	parse(inputFilename.c_str(), e, printDefaults);
	e.setup();
	Citations::print();
	if(dryRun)
	{	logPrintf("Dry run successful: commands are valid and initialization succeeded.\n");
		finalizeSystem();
		return 0;
	}
	
	if(e.cntrl.fixed_n)
	{	//Band structure calculation - ion and fluid minimization need to be handled differently
		if(e.exCorr.needsKEdensity()) //compute (fixed) KE density for meta-GGAs
			eVars.tau = eVars.KEdensity();
		eVars.EdensityAndVscloc(e.ener);
		if(eVars.fluidSolver && eVars.fluidSolver->needsGummel())
		{	//Relies on the gummel loop, so EdensityAndVscloc would not have invoked minimize
			eVars.fluidSolver->minimizeFluid();
			eVars.EdensityAndVscloc(e.ener); //update Vscloc
		}
		e.iInfo.augmentDensityGridGrad(eVars.Vscloc); //update Vscloc atom projections for ultrasoft psp's 
		if(e.cntrl.invertKS) //Inverse Kohn-Sham problem (sequence of band structure calculations)
		{	InverseKohnSham inverseKS(e);
			inverseKS.minimize(e.inverseKSminParams);
		}
		else //Single band structure calculation
		{	logPrintf("\n----------- Band structure minimization -------------\n"); logFlush();
			elecMinimize(e); // Do the band-structure minimization
		}
	}
	else if(e.vibrations) //Bypasses ionic/lattice minimization, calls electron/fluid minimization loops at various ionic configurations
	{	e.vibrations->calculate();
	}
	else if(e.latticeMinParams.nIterations)
	{	//Lattice minimization loop (which invokes the ionic minimization loop)
		LatticeMinimizer lmin(e);
		lmin.minimize(e.latticeMinParams);
	}
	else
	{	//Ionic minimization loop (which calls electron/fluid minimization loops)
		IonicMinimizer imin(e);
		imin.minimize(e.ionicMinParams);
	}

	//Final dump:
	e.dump(DumpFreq_End, 0);
	
	finalizeSystem();
	if(globalLog && globalLog != stdout) fclose(globalLog);
	return 0;
}
