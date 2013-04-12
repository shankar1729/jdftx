/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman, Kathleen Schwarz

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

#ifndef JDFTX_ELECTRONIC_POLARIZABILITY_H
#define JDFTX_ELECTRONIC_POLARIZABILITY_H

#include <electronic/common.h>
#include <core/vector3.h>

//! Calculate the polarizability in a convenient eigenbasis
struct Polarizability
{
	//! Response-matrix whose eigen-basis to use for final output
	enum EigenBasis
	{	NonInteracting, //!< Non-interacting response of Kohn-Sham system
		External, //!< Charge response to external electrostatic potential
		Total //!< Charge response to total electrostatic potential
	}
	eigenBasis;
	
	double Ecut; //!< energy-cutoff for occupied-valence pair densities (if zero, 4*Ecut of wavefunctions)
	int nEigs; //!< number of eigenvectors in output (if zero, output all)
	
	vector3<> dk; //!< k-point difference at which to obtain results
	
	string dkFilenamePattern; //!< if non-null, read wavefunctions and eigenvalues for offset states form here
	
	Polarizability();
	void dump(const Everything& e); //!< compute and dump polarizability eigenvectors and matrix elements
	
private:
	string dkFilename(int ik, string varName) const; //!< get the filename to read specified variable for specified k-point
	friend class PairDensityCalculator;
};

#endif // JDFTX_ELECTRONIC_POLARIZABILITY_H
