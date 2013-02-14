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

// Test fluid solver in fortran-style C++ code. 
// (Not actually written in Fortran so as to not
// introduce a fortran compiler dependency for JDFTx)

#include <cstdio>
#include <cmath>

//Forward declarations:
#define DeclareFortranFunction(funcName) extern "C" void funcName##_
DeclareFortranFunction(initjdftx)(double* Rx, double* Ry, double* Rz, int* Sx, int* Sy, int* Sz);
DeclareFortranFunction(getionsigma)(double* sigma);
DeclareFortranFunction(minimizefluid)(double* Adiel, double* nCavity, double* rhoExplicit, double* Adiel_nCavity, double* Adiel_rhoExplicit);

int main()
{
	double Rx[3] = {10., 0., 0.};
	double Ry[3] = {0., 10., 0.};
	double Rz[3] = {0., 0., 10.};
	int Sx = 14, Sy = 14, Sz = 14;
	initjdftx_(Rx, Ry, Rz, &Sx, &Sy, &Sz);
	
	double ionWidth;
	getionsigma_(&ionWidth);
	printf("ionWidth = %lg A\n", ionWidth);
	
	//Create a test "molecule" (with a slight dipole)
	int nData = Sx*Sy*Sz;
	double* n = new double[nData];
	double* rho = new double[nData];
	double* A_n = new double[nData];
	double* A_rho = new double[nData];
	double elWidth = 2*ionWidth; //width for some arbitrary electron density
	double elPrefac = 1./(elWidth*sqrt(2*M_PI));
	double ionPrefac = -1./(ionWidth*sqrt(2*M_PI));
	int i=0;
	for(int iz=-Sz/2; iz<Sz-Sz/2; iz++)
		for(int iy=-Sy/2; iy<Sy-Sy/2; iy++)
			for(int ix=-Sx/2; ix<Sx-Sx/2; ix++)
			{	double rEl = 0.25*sqrt(ix*ix + iy*iy + iz*iz);
				double rIon = 0.25*sqrt((ix-1)*(ix-1) + iy*iy + iz*iz);
				double rhoEl = elPrefac * exp(-0.5*pow(rEl/elWidth,2));
				double rhoIon = ionPrefac * exp(-0.5*pow(rIon/ionWidth,2));
				n[i] = rhoEl;
				rho[i] = rhoEl + rhoIon;
				i++;
			}
	double A;
	minimizefluid_(&A, n, rho, A_n, A_rho);
	printf("Adiel = %lg eV\n", A);
	delete[] n,
	delete[] rho;
	delete[] A_n;
	delete[] A_rho;
}
