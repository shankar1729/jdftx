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
#include <core/matrix.h>

//! @addtogroup IonicSystem
//! @{
//! @file LatticeMinimizer.h Class LatticeMinimizer and related definitions

//!Vector-space entry for lattice minimization (stress and forces):
struct LatticeGradient
{	matrix3<> lattice; //!< lattice component (stress)
	IonicGradient ionic; //!< ionic or internal geometry component (forces)
	diagMatrix thermostat; //!< optional extra degrees of freedom used by thermostats in IonicDynamics
	diagMatrix barostat; //!< optional extra degrees of freedom used by barostats in IonicDynamics

	void init(const class IonInfo& iInfo); //!< initialize correct sizes
	LatticeGradient& operator*=(double scale);
	LatticeGradient& operator+=(const LatticeGradient&);
	LatticeGradient operator+(const LatticeGradient&) const;
	LatticeGradient operator-(const LatticeGradient&) const;

	static double maxComponent(const LatticeGradient*); //!< MinimizeParams::maxCalculator for max force/stress component
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
	LatticeMinimizer(Everything&, bool dynamicsMode=false, bool statP=false, bool statStress=false);
	
	//Virtual functions from Minimizable:
	void step(const LatticeGradient& dir, double alpha);
	double compute(LatticeGradient* grad, LatticeGradient* Kgrad);
	bool report(int iter);
	void constrain(LatticeGradient&);
	double safeStepSize(const LatticeGradient& dir) const;
	double sync(double x) const; //!< All processes minimize together; make sure scalars are in sync to round-off error

	double minimize(const MinimizeParams& params); //!< minor addition to Minimizable::minimize to invoke charge analysis at final positions
	int nFree() { return (dynamicsMode and statP) ? 1 : int(round(trace(Pfree))); } //!< number of free lattice directions
private:
	Everything& e;
	bool dynamicsMode; //!< whether LatticeMinimizer is operating as a helper class for IonicDynamics
	bool statP; //!< in dynamicsMode, whether pressure is stat'd (hydrostatic barostat)
	bool statStress; //!<  in dynamicsMode, whether stress is stat'd (anisotropic barostat)
	IonicMinimizer imin; //!< acts as a helper to carry out force calculation, ionic steps and constraints
	matrix3<> Rorig; //!< original lattice vectors (prior to relaxation)
	matrix3<> strain; //!< minimizer state = strain relative to Rorig (i.e. R = Rorig * (1 + strain))
	bool skipWfnsDrag; //!< whether to temporarily skip wavefunction drag due to large steps (for stability)
	matrix3<> Pfree; //!< projection operator onto free directions (accounting for lattMoveScale and truncation)
	double latticeK; //!< preconditioning factor for lattice degrees of freedom
	
	//! Updates lattice dependent quantities, but does not reconverge ionic positions or wavefunctions
	static void updateLatticeDependent(Everything& e);
};

//! @}
#endif // JDFTX_ELECTRONIC_LATTICEMINIMIZER_H
