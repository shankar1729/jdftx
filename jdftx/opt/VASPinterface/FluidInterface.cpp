/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman

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
#include <electronic/ColumnBundle.h>
#include <fluid/FluidSolver.h>
#include <commands/parser.h>
#include <core/Units.h>

Everything e;

//! The vast majority of fortran compilers decorate lower case names with an _
//! Edit this macro if you are using some other compiler which uses a different convention.
//! Uggggh Fortran!
//! NOTE: All C++ function names in this file must be entirely in lower case
#define DeclareFortranFunction(funcName) extern "C" void funcName##_

//! Initialize JDFTx.
//! The lattice vectors, sample counts and the scalar fields in fluidminimize.
//! are in the same order as in Fortran, and specified in Angstroms.
//! This interface takes care of giving the core of JDFTx reversed lattice directions
//! and sample counts, and for all unit conversions (energies in eV, distances in Angstrom etc.)
//! @param Rx (in, 3-vector) First lattice direction
//! @param Ry (in, 3-vector) Second lattice direction
//! @param Rz (in, 3-vector) Third lattice direction
//! @param Sx (in, integer) Number of FFT points along first lattice direction
//! @param Sy (in, integer) Number of FFT points along second lattice direction
//! @param Sz (in, integer) Number of FFT points along third lattice direction
DeclareFortranFunction(initjdftx)(double* Rx, double* Ry, double* Rz, int* Sx, int* Sy, int* Sz)
{
	//Open log file
	globalLog = fopen("FLULOG", "w");
	if(!globalLog)
	{	globalLog = stdout;
		logPrintf("WARNING: Could not open log file 'FLULOG' for writing, using standard output.\n");
	}
	
	//Initialize environment and print banner:
	const char* execName = "N/A (Running as a shared library providing fluid solvers)";
	initSystem(1, (char**)&execName);
	
	//Write a wrapper input file:
	FILE* fpWrap = fopen("FLUCAR.in", "w");
	if(!fpWrap) die("Could not write input-file wrapper. JDFTx-VASP interface needs write access to current directory.\n");
	fprintf(fpWrap,
		"lattice \\\n"
		"    %.15lf %.15lf %.15lf \\\n"
		"    %.15lf %.15lf %.15lf \\\n"
		"    %.15lf %.15lf %.15lf \n"
		"fftbox %d %d %d\n"
		"elec-cutoff 0\n"
		"symmetries none\n"
		"include FLUCAR\n",
		 Rz[0]*Angstrom, Ry[0]*Angstrom, Rx[0]*Angstrom,
		 Rz[1]*Angstrom, Ry[1]*Angstrom, Rx[1]*Angstrom,
		 Rz[2]*Angstrom, Ry[2]*Angstrom, Rx[2]*Angstrom,
		 *Sz, *Sy, *Sx);
	fclose(fpWrap);
	
	//Initialize system:
	parse("FLUCAR.in", e);
	system("rm FLUCAR.in");
	if(e.eVars.fluidParams.fluidType == FluidNone) die("No fluid model specified in FLUCAR.\n");
	if(e.iInfo.ionWidthMethod == IonInfo::IonWidthEcut)
	{	logPrintf("JDFTx interface does not have access to VASP energy cutoff:\n"
			"\tUsing FFTbox to determine nuclear width instead.\n");
		e.iInfo.ionWidthMethod = IonInfo::IonWidthFFTbox;
	}
	e.setup();
	Citations::add("JDFTx-VASP interface",
		"K. Mathew, R. Sundararaman, K. Letchworth-Weaver, T.A. Arias and R.G. Hennig (under preparation)");
	Citations::print();
}


//! Get the recommended nuclear width for the chosen fluid model
//! @param sigma (out, scalar) Recommended gaussian width
DeclareFortranFunction(getionsigma)(double* sigma)
{
	*sigma = e.iInfo.ionWidth/Angstrom;
}


//! Minimize the fluid and return the free energy and its derivatives.
//! @param Adiel (out, scalar) Fluid free energy
//! @param nCavity (in, real-space scalar field) Electron density involved in cavity determination
//! @param rhoExplicit (in, real-space scalar field) Total charge density of electronic system (valence electrons + nuclei)
//! @param Adiel_nCavity (out, real-space scalar field) Functional derivative of Adiel with respect to nCavity
//! @param Adiel_rhoExplicit (out, real-space scalar field) Functional derivative of Adiel with respect to rhoExplicit
DeclareFortranFunction(minimizefluid)(double* Adiel,
	double* nCavity, double* rhoExplicit, double* Adiel_nCavity, double* Adiel_rhoExplicit)
{
	//Convert inputs to JDFTx objects:
	ScalarField n, rho;
	nullToZero(n, e.gInfo); nullToZero(rho, e.gInfo);
	eblas_copy(n->data(), nCavity, e.gInfo.nr);
	eblas_copy(rho->data(), rhoExplicit, e.gInfo.nr);
	n *= pow(Angstrom,-3); rho *= pow(Angstrom,-3); //convert to atomic units

	//Run the fluid solver:
	logPrintf("\n---------------------- Fluid Minimization -----------------------\n");
	e.eVars.fluidSolver->set(J(rho), J(n)); n=0; rho=0;
	e.eVars.fluidSolver->minimizeFluid();
	ScalarFieldTilde A_n, A_rho; IonicGradient extraForces;
	double A = e.eVars.fluidSolver->get_Adiel_and_grad(A_rho, A_n, extraForces);
	e.dump(DumpFreq_Electronic, -1);
	
	//Convert outputs:
	*Adiel = A/eV;
	eblas_copy(Adiel_nCavity, I(A_n * (1./eV))->data(), e.gInfo.nr);
	eblas_copy(Adiel_rhoExplicit, I(A_rho * (1./eV))->data(), e.gInfo.nr);
}

