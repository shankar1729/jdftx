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

inline double determinant(diagMatrix& D)
{	double det = 1.;
	for(int i=0; i<D.nCols(); i++)
		det *= D[i];
	return det;
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
	int nbands = e.eInfo.nBands;
	double wInv = e.eInfo.spinType==SpinNone ? 0.5 : 1.0; //normalization factor from external to internal fillings

	// Read wavefunctions
	std::vector<ColumnBundle> C1(qnums);
	std::vector<ColumnBundle> C2(qnums);
	read(C1, "C1.wfns", e);
	read(C2, "C2.wfns", e);

	// Read fillings
	std::vector<diagMatrix> F1(qnums);
	std::vector<diagMatrix> F2(qnums);
	FILE* fp1 = fopen("F1.fillings","r");
	if(!fp1) die("Can't open file %s to read initial fillings!\n", "F1.fillings");
	FILE* fp2 = fopen("F2.fillings","r");
	if(!fp2) die("Can't open file %s to read initial fillings!\n", "F2.fillings");
	for (int q=0; q<qnums; q++)
	{	F1[q].resize(nbands); F2[q].resize(nbands);
		F1[q].scan(fp1); F2[q].scan(fp2);
		F1[q] *= wInv; F2[q] *= wInv; //NOTE: fillings are always 0 to 1 internally, but read/write 0 to 2 for SpinNone
	}

	std::vector<matrix> dipoleXMatrices;
	std::vector<matrix> dipoleYMatrices;
	std::vector<matrix> dipoleZMatrices;
	std::vector<matrix> overlapMatrices;
	complex detDipoleX = 1;
	complex detDipoleY = 1;
	complex detDipoleZ = 1;
	complex detOverlap = 1;
	
	// Real-space kernels for the dipole calculations
	DataRptr r0, r1, r2;
	nullToZero(r0, e.gInfo); 	nullToZero(r1, e.gInfo); 	nullToZero(r2, e.gInfo);
	applyFunc_r(e.gInfo, Moments::rn_pow_x, 0, e.gInfo.R, 1, vector3<>(0.,0.,0.), r0->data());
	applyFunc_r(e.gInfo, Moments::rn_pow_x, 1, e.gInfo.R, 1, vector3<>(0.,0.,0.), r1->data());
	applyFunc_r(e.gInfo, Moments::rn_pow_x, 2, e.gInfo.R, 1, vector3<>(0.,0.,0.), r2->data());
	
	for(size_t q =0; q<e.eInfo.qnums.size(); q++)
	{		
		int noBands = e.eInfo.nBands;
		dipoleXMatrices.at(q).init(noBands, noBands);
		dipoleYMatrices.at(q).init(noBands, noBands);
		dipoleZMatrices.at(q).init(noBands, noBands);
		for(int i=0; i<noBands; i++)
		{	for(int j=0; j<i; j++)
			{	
				double tol = 1e-4;
				// If either of the orbitals is empty, no need to compute anything
				if((F1[q][i] < tol) or (F2[q][j] < tol))
				{	dipoleXMatrices.at(q).set(i,j, 0.);
					dipoleXMatrices.at(q).set(j,i, 0.);
					dipoleYMatrices.at(q).set(i,j, 0.);
					dipoleYMatrices.at(q).set(j,i, 0.);
					dipoleZMatrices.at(q).set(i,j, 0.);
					dipoleZMatrices.at(q).set(j,i, 0.);
					overlapMatrices.at(q).set(i,j, 0.);
					overlapMatrices.at(q).set(j,i, 0.);
					continue;
				}
				// If both orbitals are occupied, then return the matrix element
				else if((F1[q][i] > (1.-tol)) or (F2[q][j] < (1.-tol)))
				{	
					complexDataRptr psi1 = I(C1.at(q).getColumn(i));
					complexDataRptr psi2 = I(C2.at(q).getColumn(j));
					
					complex dipoleX = (norm(detDipoleX) ? integral(psi1*r0*psi2) : 0.);
					complex dipoleY = (norm(detDipoleY) ? integral(psi1*r1*psi2) : 0.);
					complex dipoleZ = (norm(detDipoleZ) ? integral(psi1*r2*psi2) : 0.);
					complex overlap = (norm(detOverlap) ? integral(psi1*psi2) : 0.);
					
					dipoleXMatrices.at(q).set(i,j, dipoleX);
					dipoleXMatrices.at(q).set(j,i, dipoleX.conj());
					dipoleYMatrices.at(q).set(i,j, dipoleY);
					dipoleYMatrices.at(q).set(j,i, dipoleY.conj());
					dipoleZMatrices.at(q).set(i,j, dipoleZ);
					dipoleZMatrices.at(q).set(j,i, dipoleZ.conj());
					overlapMatrices.at(q).set(i,j, overlap);
					overlapMatrices.at(q).set(j,i, overlap.conj());
				}
				else
				{	die("Non-integer fillings found... Exiting.\n");
				}
			}
		}
		matrix evecs; diagMatrix eigs;
		dipoleXMatrices.at(q).diagonalize(evecs, eigs);
		detDipoleX *= determinant(eigs);
		dipoleYMatrices.at(q).diagonalize(evecs, eigs);
		detDipoleY *= determinant(eigs);
		dipoleZMatrices.at(q).diagonalize(evecs, eigs);
		detDipoleZ *= determinant(eigs);
		overlapMatrices.at(q).diagonalize(evecs, eigs);
		detOverlap *= determinant(eigs);
	}
	
	/// //////////////////////////////// ///
	
	logPrintf("\nOverlap: (%.5e, %.5e)\n", real(detOverlap), imag(detOverlap));
	logPrintf("\nDipole: <(%5.e, %.5e), (%5.e, %.5e), (%5.e, %.5e)>\n",
			  real(detDipoleX), imag(detDipoleX), real(detDipoleY), imag(detDipoleY), real(detDipoleZ), imag(detDipoleZ));
	
	
	finalizeSystem();
	if(globalLog != stdout) fclose(globalLog);
	return 0;
}