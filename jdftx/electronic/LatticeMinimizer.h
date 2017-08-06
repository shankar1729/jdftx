/*-------------------------------------------------------------------
Copyright 2012 Deniz Gunceler, Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_LATTICEMINIMIZER_H
#define JDFTX_ELECTRONIC_LATTICEMINIMIZER_H

#include <electronic/common.h>
#include <core/Minimize.h>
#include <core/matrix3.h>

//! @addtogroup IonicSystem
//! @{
//! @file LatticeMinimizer.h Class LatticeMinimizer and related definitions

//Functions required by Minimizable<matrix3<>>
void axpy(double alpha, const matrix3<>& x, matrix3<>& y); //!< accumulate operation: Y += alpha*X
double dot(const matrix3<>& x, const matrix3<>& y); //!< inner product
matrix3<> clone(const matrix3<>& x); //!< create a copy
void randomize(matrix3<>& x); //!< initialize with random numbers

//! Lattice minimizer
class LatticeMinimizer : public Minimizable<matrix3<>>
{
public:
	LatticeMinimizer(Everything&);
	
	//Virtual functions from Minimizable:
	void step(const matrix3<>& dir, double alpha);
	double compute(matrix3<>* grad, matrix3<>* Kgrad);
	bool report(int iter);
	void constrain(matrix3<>&);
	double safeStepSize(const matrix3<>& dir) const;
	double sync(double x) const; //!< All processes minimize together; make sure scalars are in sync to round-off error

	//! Calculates the stresses along the strain directions
	std::vector<double> calculateStress();

	//!Set of independent directions in the space of all allowed strains.
	//!Their span is consistent with symmetries and truncation (if any).
	std::vector<matrix3<>> strainBasis;
private:
	Everything& e;
	matrix3<> Rorig; //!< original lattice vectors (prior to relaxation)
	matrix3<> strain; //!< minimizer state = strain relative to Rorig (i.e. R = Rorig * (1 + strain))
	bool skipWfnsDrag; //!< whether to temporarily skip wavefunction drag due to large steps (for stability)
	
	double h; //! Finite difference step size
	double centralDifference(matrix3<> direction);  //! Returns the numerical derivative along the given strain
	
	//! Updates lattice dependent quantities, but does not
	//! reconverge ionic positions or wavefunctions
	static void updateLatticeDependent(Everything& e, bool ignoreElectronic=false);
	
	friend class IonDynamics;
};

//! @}
#endif // JDFTX_ELECTRONIC_LATTICEMINIMIZER_H
