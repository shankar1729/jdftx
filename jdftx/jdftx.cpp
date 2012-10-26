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
	logPrintf("\t-i --input <filename>   specifies command input file, default = stdin\n");
	logPrintf("\t-o --output <filename>  specifies output log file, default = stdout\n");
	logPrintf("\t-d --no-append          overwrite output file instead of appending\n");
	logPrintf("\t-t --template           prints an input file template\n");
	logPrintf("\t-n --dry-run            quit after initialization (to verify commands and other input files)\n");
	logPrintf("\t-p --print-defaults     print status of default commands issued automatically.\n");
	logPrintf("\n");
}

//Program entry point
int main(int argc, char** argv, char** argp)
{
	Everything e; //the parent data structure for, well, everything
	
	//Parse command line:
	string inputFilename, logFilename; bool appendOutput=true, dryRun=false, printDefaults=false;
	option long_options[] =
		{	{"help", no_argument, 0, 'h'},
			{"version", no_argument, 0, 'v'},
			{"input",  required_argument, 0, 'i'},
			{"output", required_argument, 0, 'o'},
			{"no-append", no_argument, 0, 'd'},
			{"template", no_argument, 0, 't'},
			{"dry-run", no_argument, 0, 'n'},
			{"print-defaults", no_argument, 0, 'p'},
			{0, 0, 0, 0}
		};
	while (1)
	{	int c = getopt_long(argc, argv, "hvi:o:dtnp", long_options, 0);
		if (c == -1) break; //end of options
		switch (c)
		{	case 'v': printVersionBanner(); return 0;
			case 'h': printUsage(argv[0]); return 0;
			case 'i': inputFilename.assign(optarg); break;
			case 'o': logFilename.assign(optarg); break;
			case 'd': appendOutput=false; break;
			case 't': printDefaultTemplate(e); return 0;
			case 'n': dryRun=true; break;
			case 'p': printDefaults=true; break;
			default: printUsage(argv[0]); return 1;
		}
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
		if(e.cntrl.invertKS) //Inverse Kohn-Sham problem (sequence of band structure calculations)
		{	InverseKohnSham inverseKS(e);
			inverseKS.minimize(e.inverseKSminParams);
		}
		else //Single band structure calculation
		{	logPrintf("\n----------- Band structure minimization -------------\n"); logFlush();
			elecMinimize(e); // Do the band-structure minimization
		}
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
	if(globalLog != stdout) fclose(globalLog);
	return 0;
}
