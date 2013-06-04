/*-------------------------------------------------------------------
Copyright 2011 Deniz Gunceler

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

// Calculates the electric dipole transition elements between two slater determinants

#include <electronic/Everything.h>
#include <core/Util.h>
#include <core/Util.h>
#include <commands/parser.h>
#include <getopt.h>

#include <electronic/ColumnBundle.h>
#include <electronic/matrix.h>
#include <electronic/operators.h>
#include <electronic/Dump.h>

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
	logPrintf("\t-n --dry-run            quit after initialization (to verify commands and other input files)\n");
	logPrintf("\t-s --skip-defaults      skip printing status of default commands issued automatically.\n");
	logPrintf("\n");
}

int main(int argc, char** argv)
{

	Everything e; //the parent data structure for, well, everything

	string inputFilename, logFilename; bool appendOutput=true, printDefaults=true;
	
	option long_options[] =
	{	{"help", no_argument, 0, 'h'},
		{"version", no_argument, 0, 'v'},
		{"input",  required_argument, 0, 'i'},
		{"output", required_argument, 0, 'o'},
		{"no-append", no_argument, 0, 'd'},
		{"template", no_argument, 0, 't'},
		{"dry-run", no_argument, 0, 'n'},
		{"skip-defaults", no_argument, 0, 's'},
		{0, 0, 0, 0}
	};
	while (1)
	{	int c = getopt_long(argc, argv, "hvi:o:dtns", long_options, 0);
		if (c == -1) break; //end of options
		switch (c)
		{	case 'v': printVersionBanner(); return 0;
			case 'h': printUsage(argv[0]); return 0;
			case 'i': inputFilename.assign(optarg); break;
			case 'o': logFilename.assign(optarg); break;
			case 'd': appendOutput=false; break;
			case 't': printDefaultTemplate(e); return 0;
			case 's': printDefaults=false; break;
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
	
	initSystem(argc, argv);
	parse(inputFilename.c_str(), e, printDefaults);
	e.setup();
	
	/// ACTUAL CALCULATION HAPPERNS HERE ///
	
	int qnums = e.eInfo.qnums.size();
	
	std::vector<ColumnBundle> C1(qnums);
	std::vector<ColumnBundle> C2(qnums);
	read(C1, "C1.wfns", e);
	read(C2, "C2.wfns", e);
	
	std::vector<matrix> dipole;
	double det = 1;
	
	// Kernels for the dipole calculations
	DataRptr r0, r1, r2;
	nullToZero(r0, e.gInfo); 	nullToZero(r1, e.gInfo); 	nullToZero(r2, e.gInfo);
	applyFunc_r(e.gInfo, Moments::rn_pow_x, 0, e.gInfo.R, 1, vector3<>(0.,0.,0.), r0->data());
	applyFunc_r(e.gInfo, Moments::rn_pow_x, 1, e.gInfo.R, 1, vector3<>(0.,0.,0.), r1->data());
	applyFunc_r(e.gInfo, Moments::rn_pow_x, 2, e.gInfo.R, 1, vector3<>(0.,0.,0.), r2->data());
	
	for(size_t q =0; q<e.eInfo.qnums.size(); q++)
	{
		// If any of the blocks has zero determinant, then there's no need to compute others
		if(det == 0)
			break;
		
		int noBands = e.eInfo.nBands;
		dipole.at(q).init(noBands, noBands);
		for(int i=0; i<noBands; i++)
			for(int j=0; j<i; j++)
			{	//complex moment = dot(C1[q]);
				//dipole.at(q).set(i,j, dot());
			}
	}
	
	/// //////////////////////////////// ///
		
	finalizeSystem();
	if(globalLog != stdout) fclose(globalLog);
	return 0;
}