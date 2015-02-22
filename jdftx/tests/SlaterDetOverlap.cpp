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

#include <core/Util.h>
#include <commands/parser.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/matrix.h>
#include <electronic/operators.h>
#include <electronic/Dump.h>

int main(int argc, char** argv)
{	//Parse command line, initialize system and logs:
	string inputFilename; bool dryRun, printDefaults;
	initSystemCmdline(argc, argv, "Calculate electric dipole transition elements between two slater determinants.", inputFilename, dryRun, printDefaults);
	
	Everything e; //the parent data structure for, well, everything
	logSuspend(); e.elecMinParams.fpLog = nullLog;
	parse(readInputFile(inputFilename), e, printDefaults);
	e.setup();
	logResume(); e.elecMinParams.fpLog = globalLog;
	if(e.eInfo.isNoncollinear()) die("'%s' not yet implemented with noncollinear spins.\n", argv[0]);

	/// ACTUAL CALCULATION HAPPERNS HERE ///
	
	int qnums = e.eInfo.qnums.size();
	int nbands = e.eInfo.nBands;
	double wInv = e.eInfo.spinType==SpinNone ? 0.5 : 1.0; //normalization factor from external to internal fillings

	// Read wavefunctions
	std::vector<ColumnBundle> C1(qnums);
	std::vector<ColumnBundle> C2(qnums);
	init(C1, e.eInfo.nStates, e.eInfo.nBands, &(e.basis[0]), &(e.eInfo));
	init(C2, e.eInfo.nStates, e.eInfo.nBands, &(e.basis[0]), &(e.eInfo));
	read(C1, "C1.wfns", e.eInfo);
	read(C2, "C2.wfns", e.eInfo);

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
		
		// Prints fillings
		logPrintf("\n\nNormalized fillings: \n");
		logPrintf("F1 (q=%i):\n", q);
		F1.at(q).print(globalLog);
		logPrintf("\n\nF2 (q=%i):\n", q);
		F2.at(q).print(globalLog);
		
	}

	std::vector<matrix> dipoleXMatrices(qnums);
	std::vector<matrix> dipoleYMatrices(qnums);
	std::vector<matrix> dipoleZMatrices(qnums);
	std::vector<matrix> overlapMatrices(qnums);
	complex detDipoleX = 1.;
	complex detDipoleY = 1.;
	complex detDipoleZ = 1.;
	complex detOverlap = 1.;
	
	// Real-space kernels for the dipole calculations
	ScalarField r0, r1, r2;
	nullToZero(r0, e.gInfo); 	nullToZero(r1, e.gInfo); 	nullToZero(r2, e.gInfo);
	applyFunc_r(e.gInfo, Moments::rn_pow_x, 0, e.gInfo.R, 1, vector3<>(0.,0.,0.), r0->data());
	applyFunc_r(e.gInfo, Moments::rn_pow_x, 1, e.gInfo.R, 1, vector3<>(0.,0.,0.), r1->data());
	applyFunc_r(e.gInfo, Moments::rn_pow_x, 2, e.gInfo.R, 1, vector3<>(0.,0.,0.), r2->data());
	
	logPrintf("\n");
	for(size_t q =0; q<e.eInfo.qnums.size(); q++)
	{	
		logPrintf("Calculating qnum = %zu...", q);
		
		double n1 = trace(F1[q]); // No electrons in the q'th qnum for first wavefunction
		double n2 = trace(F2[q]); // No electrons in the q'th qnum for second wavefunction
		
		logPrintf("\n\t%f electrons on F1(q=%zu)\n\t%f electrons on F2(q=%zu)", n1, q, n2, q);
		
		double tol = 1e-4;
		if(abs(n1-n2) > tol)  // If different number of electrons are detected, all matrix elements are 0
		{	detDipoleX *= 0.;
			detDipoleY *= 0.;
			detDipoleZ *= 0.;
			detOverlap *= 0.;
			break;  // No need to compute other quantum numbers
		}
		
		int nfilled = round(n1);  // Get the number of filled states on each
		int ni=0, nj=0;           // Counters for the reduced matricex (just filled orbitals)
		int noBands = e.eInfo.nBands;
		dipoleXMatrices.at(q).init(nfilled, nfilled);
		dipoleYMatrices.at(q).init(nfilled, nfilled);
		dipoleZMatrices.at(q).init(nfilled, nfilled);
		overlapMatrices.at(q).init(nfilled, nfilled);
		for(int i=0; i<noBands; i++)
		{	nj = 0; // Reset nj
			if(F1[q][i] < (1.-tol))  // Check whether i'th orbital is filled
					continue;
			for(int j=0; j<noBands; j++)
			{	if(F2[q][j] < (1.-tol)) // Check whether j'th orbital is filled
					continue;
								
				// If both orbitals are occupied, then return the matrix element
				if((F1[q][i] > (1.-tol)) or (F2[q][j] > (1.-tol)))
				{	
					complexScalarField psi1 = I(C1.at(q).getColumn(i,0));
					complexScalarField psi2 = I(C2.at(q).getColumn(j,0));
					
					complex dipoleX = (norm(detDipoleX) ? integral(psi1*r0*psi2) : 0.);
					complex dipoleY = (norm(detDipoleY) ? integral(psi1*r1*psi2) : 0.);
					complex dipoleZ = (norm(detDipoleZ) ? integral(psi1*r2*psi2) : 0.);
					complex overlap = (norm(detOverlap) ? integral(psi1*psi2) : 0.);
					
					dipoleXMatrices.at(q).set(ni,nj, dipoleX);
					dipoleYMatrices.at(q).set(ni,nj, dipoleY);
					dipoleZMatrices.at(q).set(ni,nj, dipoleZ);
					overlapMatrices.at(q).set(ni,nj, overlap);
				}
				else
				{	die("Non-integer fillings found... Exiting.\n");
				}
				
				nj++;
			}
			ni++;
		}
		
		// Assert that the entire matrix have been filled
		detDipoleX *= det(dipoleXMatrices.at(q));
		detDipoleY *= det(dipoleYMatrices.at(q));
		detDipoleZ *= det(dipoleZMatrices.at(q));
		detOverlap *= det(overlapMatrices.at(q));
		
		logPrintf("\n\ndet = %.5e, %.5e\n", real(det(dipoleXMatrices.at(q))), imag(det(dipoleXMatrices.at(q))));
		dipoleXMatrices.at(q).print(globalLog,"%.2e%+.2ei\t" );
		logPrintf("\n");
		
		logPrintf("\nOverlap: (%.5e, %.5e)\n", real(detOverlap), imag(detOverlap));
		logPrintf("Dipole: <(%5.e, %.5e), (%5.e, %.5e), (%5.e, %.5e)>\n\n",
				  real(detDipoleX), imag(detDipoleX), real(detDipoleY), imag(detDipoleY), real(detDipoleZ), imag(detDipoleZ));
		logPrintf("\n\n");
	}
	
	/// //////////////////////////////// ///
	
	logPrintf("\n\n");
	logPrintf("\nOverlap: (%.5e, %.5e)", real(detOverlap), imag(detOverlap));
	logPrintf("\nDipole: <(%5.e, %.5e), (%5.e, %.5e), (%5.e, %.5e)>\n",
			  real(detDipoleX), imag(detDipoleX), real(detDipoleY), imag(detDipoleY), real(detDipoleZ), imag(detDipoleZ));
	
	
	finalizeSystem();
	return 0;
}
