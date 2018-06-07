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

#include <core/ScalarFieldArray.h>
#include <core/Coulomb.h>

class Everything;
class ColumnBundle;

//! @addtogroup Output
//! @{
//! @file Dump_internal.h Implementation internals for output modules

//-------------------- Implemented in DumpSIC.cpp ---------------------------

//! Output self-interaction correction for the KS eigenvalues
class DumpSelfInteractionCorrection
{
public:
	DumpSelfInteractionCorrection(const Everything& everything);
	~DumpSelfInteractionCorrection();
	double operator()(std::vector<diagMatrix>* correctedEigenvalues); //!< Evaluates the self interaction energy and (optionally) returns the corrected band eigenvalues
	double coulombExciton(int q1, int n1, int q2, int n2);  //!< Approximates the coulomb excitonic contribution between two states
	void dump(const char* filenamePattern);
	bool needsTau;  //!< The kinetic energy density is needed for meta-gga functionals.
private:
	const Everything* e;
	double calcSelfInteractionError(int q, int n); //!< Calculates the self-interaction error of the KS orbital atthe n'th band at q'th quantum number
	std::vector<ColumnBundle> DC; //!< ColumnBundle for the derivative of the wavefunctions in each cartesian direction
};

//---------------- Implemented in DumpExcitationsMoments.cpp -----------------

//! Dump information about excitation energies and matrix elements
void dumpExcitations(const Everything& e, const char* filename);

//! Helper functions for moment calculations
namespace Moments
{	void rn_pow_x(int i, vector3<> r, int dir, matrix3<> R, double moment, vector3<> r0, double* rx); //!< set multipole spatial weights
	void dumpMoment(const Everything& e, const char* filename, int n, vector3<> origin); //!< output multipole moments
}

namespace XC_Analysis
{	ScalarFieldArray tauWeizsacker(const Everything& e); //!< Output Weizsacker KE density
	ScalarFieldArray spness(const Everything& e); //!< Output 'single-particle-ness'
	ScalarFieldArray sHartree(const Everything& e); //!< output spin Hartree potentials
}

//---------------- Implemented in DumpChargedDefects.cpp -----------------

//! Slab dielectric function calculator
struct SlabEpsilon
{	string dtotFname; //!< reference electrostatic potential filename
	double sigma; //!< smoothing width
	vector3<> Efield; //!< reference electric field
	
	void dump(const Everything& e, ScalarField d_tot) const;
};

//! Bulk dielectric function calculator
struct BulkEpsilon
{	string dtotFname; //!< reference electrostatic potential filename
	vector3<> Efield; //!< reference electric field
	
	void dump(const Everything& e, ScalarField d_tot) const;
};

//! Charged defect correction calculator
struct ChargedDefect
{	//! Model charge decsription (one for each defect in unit cell)
	struct Center
	{	vector3<> pos; //!< defect center in lattice coordinates
		double q; //!< defect electron-count
		double sigma; //!< model charge Gaussian width
	};
	std::vector<Center> center; //!< list of defect positions in unit cell
	
	CoulombParams::Geometry geometry; //!< geometry of Coulomb interaction used for correction (could differ from main calculation)
	int iDir; //!< slab trunctaion direction
	
	string dtotFname; //!< electrostatic potential from reference neutral calculation
	
	double bulkEps; //!< bulk dielectric constant (Bulk mode only)
	string slabEpsFname; //!< slab dielectric profile (Slab mode only)
	
	double rMin; //!< Minimum distance from defect used for calculating alignment
	double rSigma; //!< Turn-on width of region used for calculating alignment
	
	void dump(const Everything& e, ScalarField d_tot) const;
};



//! Parameters for BGW output
struct BGWparams
{	int nBandsDense; //!< if non-zero, use a dense ScaLAPACK solver to calculate more bands
	int blockSize; //!< block size for ScaLAPACK diagonalization
	int clusterSize; //!< maximum eigenvalue cluster size to allocate extra ScaLAPACK workspace for
	
	BGWparams() : nBandsDense(0), blockSize(32), clusterSize(10) {}
};

//! @}
#endif // JDFTX_ELECTRONIC_DUMP_INTERNAL_H
