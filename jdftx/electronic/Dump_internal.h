/*-------------------------------------------------------------------
Copyright 2014 Deniz Gunceler, Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_DUMP_INTERNAL_H
#define JDFTX_ELECTRONIC_DUMP_INTERNAL_H

#include <electronic/common.h>
#include <core/ScalarFieldArray.h>
#include <vector>

//-------------------- Implemented in DumpSIC.cpp ---------------------------

class DumpSelfInteractionCorrection
{
public:
	DumpSelfInteractionCorrection(const Everything& everything);
	~DumpSelfInteractionCorrection();
	double operator()(std::vector<diagMatrix>* correctedEigenvalues); //! Evaluates the self interaction energy and (optionally) returns the corrected band eigenvalues
	double coulombExciton(int q1, int n1, int q2, int n2);  //! Approximates the coulomb excitonic contribution between two states
	void dump(const char* filenamePattern);
	bool needsTau;  //! The kinetic energy density is needed for meta-gga functionals.
private:
	const Everything* e;
	double calcSelfInteractionError(int q, int n); //! Calculates the self-interaction error of the KS orbital atthe n'th band at q'th quantum number
	std::vector<ColumnBundle> DC; //!< ColumnBundle for the derivative of the wavefunctions in each cartesian direction
};

//---------------- Implemented in DumpExcitationsMoments.cpp -----------------

void dumpExcitations(const Everything& e, const char* filename);

namespace Moments
{	void dumpMoment(const Everything& e, const char* filename, int n, vector3<> origin);
}

namespace XC_Analysis
{	ScalarFieldArray tauWeizsacker(const Everything& e);
	ScalarFieldArray spness(const Everything& e);
	ScalarFieldArray sHartree(const Everything& e);
}

//---------------- Implemented in DumpChargedDefects.cpp -----------------

struct SlabEpsilon
{	string dtotFname;
	double sigma;
	vector3<> Efield;
	
	void dump(const Everything& e, ScalarField d_tot) const;
};

struct ChargedDefect
{	struct Center
	{	vector3<> pos; //defect center in lattice coordinates
		double q, sigma; //defect electron-count and Gaussian width
	};
	std::vector<Center> center; //list of defect positions in unit cell
	
	string dtotFname; //electrostatic potential from reference neutral calculation
	
	double bulkEps; //bulk dielectric constant (Bulk mode only)
	string slabEpsFname; //slab dielectric profile (Slab mode only)
	
	double rMin, rSigma; //! Minimum distance from defect and turn-on width used for calculating alignment
	
	void dump(const Everything& e, ScalarField d_tot) const;
};


#endif // JDFTX_ELECTRONIC_DUMP_INTERNAL_H
