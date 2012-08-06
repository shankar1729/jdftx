/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_BLIP_H
#define JDFTX_ELECTRONIC_BLIP_H

#include <core/Data.h>

//PW to blip conversion utility
//To use, create an object:
//	BlipConverter convert(Nx,Ny,Nz);
//and conversions can be done as:
//	vBlip = convert(vPW);
class BlipConverter
{
	int Nx, Ny, Nz;
	double *xGamma, *yGamma, *zGamma;
	double* newGamma(int N);
public:
	BlipConverter(int Nx, int Ny, int Nz);
	~BlipConverter();

	//Given a PW basis object (in real or reciprocal space) v,
	//return corresponding real-space Blip coefficient set
	// (for double or complex vectors)
	DataRptr operator()(const DataGptr& v);
	DataRptr operator()(const DataRptr& v);
	complexDataRptr operator()(const complexDataGptr& v);
	complexDataRptr operator()(const complexDataRptr& v);
};

//Compute the kinetic energy for a blip orbital phi (and set max local KE and location in unit cell)
double Tblip(const complexDataRptr& phi, double* tMax=0, int* i0max=0, int* i1max=0, int*i2max=0);

//Compute the local potential energy for blip orbital phi in blip potential V
double Vblip(const complexDataRptr& phi, const DataRptr& V);

#endif // JDFTX_ELECTRONIC_BLIP_H
