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

#include <core/ScalarField.h>
#include <core/vector3.h>
#include <core/WignerSeitz.h>

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
	ScalarField operator()(const ScalarFieldTilde& v) const;
	ScalarField operator()(const ScalarField& v) const;
	complexScalarField operator()(const complexScalarFieldTilde& v) const;
	complexScalarField operator()(const complexScalarField& v) const;
};

//! Resample a scalar field from one grid to another using BLIPs
//! The periodicity of both grids are broken on their Wigner-Seitz (WS) cells centered at the origin,
//! and the values in the input WS outside the output WS are truncated, while the region in the
//! output WS outside the input WS is set to zero.
class BlipResampler
{	const GridInfo& gInfoIn;
	const GridInfo& gInfoOut;
	BlipConverter converter;
	WignerSeitz wsIn;
	WignerSeitz wsOut;
public:
	BlipResampler(const GridInfo& gInfoIn, const GridInfo& gInfoOut);
	
	ScalarField operator()(const ScalarFieldTilde& v) const;
	complexScalarField operator()(const complexScalarFieldTilde& v) const;
};

//Compute the kinetic energy for a blip orbital phi (and set max local KE and location in unit cell)
double Tblip(const complexScalarField& phi, double* tMax=0, int* i0max=0, int* i1max=0, int*i2max=0);

//Compute the local potential energy for blip orbital phi in blip potential V
double Vblip(const complexScalarField& phi, const ScalarField& V);

#endif // JDFTX_ELECTRONIC_BLIP_H
