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
#include <core/vector3.h>

//PW to blip conversion utility
//To use, create an object:
//	BlipConverter convert(Nx,Ny,Nz);
//and conversions can be done as:
//	vBlip = convert(vPW);
class BlipConverter
{
	vector3<int> S;
	std::vector<double> gamma[3];
public:
	BlipConverter(const vector3<int>& S);

	//Given a PW basis object (in real or reciprocal space) v,
	//return corresponding real-space Blip coefficient set
	// (for double or complex vectors)
	DataRptr operator()(const DataGptr& v) const;
	DataRptr operator()(const DataRptr& v) const;
	complexDataRptr operator()(const complexDataGptr& v) const;
	complexDataRptr operator()(const complexDataRptr& v) const;
	
	//! Resample a scalar field on a different grid using cubic spline interpolation (in->gInfo.S must match S)
	//! Note that the input is effectively tiled to cover the entire output box (if parts of it overflow the input box)
	//! Use the wsMask to get a scalarfield that is 1 over only one of Wigner-Seitz shaped tile (so that multiplying selects a single image)
	DataRptr resample(const DataGptr& in, const GridInfo& gInfoOut) const;
	complexDataRptr resample(const complexDataGptr& in, const GridInfo& gInfoOut) const;
	
	DataRptr wsMask(const GridInfo& gInfoIn, const GridInfo& gInfoOut) const;
};

//Compute the kinetic energy for a blip orbital phi (and set max local KE and location in unit cell)
double Tblip(const complexDataRptr& phi, double* tMax=0, int* i0max=0, int* i1max=0, int*i2max=0);

//Compute the local potential energy for blip orbital phi in blip potential V
double Vblip(const complexDataRptr& phi, const DataRptr& V);

#endif // JDFTX_ELECTRONIC_BLIP_H
