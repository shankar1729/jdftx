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

#include <electronic/IonicMinimizer.h>

//! @addtogroup IonicSystem
//! @{
//! @file LatticeMinimizer.h Class LatticeMinimizer and related definitions

//!Vector-space entry for lattice minimization (stress and forces):
struct LatticeGradient
{	matrix3<> lattice; //!< lattice component (stress)
	IonicGradient ionic; //!< ionic or internal geometry component (forces)
	
	LatticeGradient& operator*=(double scale);
};

//Functions required by Minimizable<LatticeGradient>
void axpy(double alpha, const LatticeGradient& x, LatticeGradient& y); //!< accumulate operation: Y += alpha*X
double dot(const LatticeGradient& x, const LatticeGradient& y); //!< inner product
LatticeGradient clone(const LatticeGradient& x); //!< create a copy
void randomize(LatticeGradient& x); //!< initialize with random numbers

//! Lattice minimizer
class LatticeMinimizer : public Minimizable<LatticeGradient>
{
public:
	LatticeMinimizer(Everything&);
	
	//Virtual functions from Minimizable:
	void step(const LatticeGradient& dir, double alpha);
	double compute(LatticeGradient* grad, LatticeGradient* Kgrad);
	bool report(int iter);
	void constrain(LatticeGradient&);
	double safeStepSize(const LatticeGradient& dir) const;
	double sync(double x) const; //!< All processes minimize together; make sure scalars are in sync to round-off error

	void calculateStress(); //!< calculate current stress (in Eh/a0^3 units) and store to IonInfo::stress
	double minimize(const MinimizeParams& params); //!< minor addition to Minimizable::minimize to invoke charge analysis at final positions
private:
	Everything& e;
	IonicMinimizer imin;
	matrix3<> Rorig; //!< original lattice vectors (prior to relaxation)
	matrix3<> strain; //!< minimizer state = strain relative to Rorig (i.e. R = Rorig * (1 + strain))
	vector3<> K; //!< preconditioner (scale factors on each lattice gradient direction)
	bool skipWfnsDrag; //!< whether to temporarily skip wavefunction drag due to large steps (for stability)
	
	//!Set of independent directions in the space of all allowed strains.
	//!Their span is consistent with symmetries and truncation (if any).
	std::vector<matrix3<>> strainBasis;

	double h; //! Finite difference step size
	double centralDifference(matrix3<> direction);  //! Returns the numerical derivative along the given strain
	
	//! Updates lattice dependent quantities, but does not
	//! reconverge ionic positions or wavefunctions
	static void updateLatticeDependent(Everything& e, bool ignoreElectronic=false);
	
	friend class IonDynamics;
};

//! @}
#endif // JDFTX_ELECTRONIC_LATTICEMINIMIZER_H
