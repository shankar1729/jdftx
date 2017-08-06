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

#ifndef JDFTX_ELECTRONIC_ELECMINIMIZER_H
#define JDFTX_ELECTRONIC_ELECMINIMIZER_H

#include <core/Minimize.h>

class Everything;
class ElecInfo;
class ColumnBundle;
class matrix;

//! @addtogroup ElecSystem
//! @{
//! @file ElecMinimizer.h ElecMinimizer and helpers

//! Vector space entry for electronic minimization
struct ElecGradient
{	std::vector<ColumnBundle> C; //!< wavefunctions
	std::vector<matrix> Haux; //!< auxiliary Hamiltonian
	const ElecInfo* eInfo;
	
	void init(const Everything& e); //!< initialize C and Haux with the correct sizes for everything
	
	ElecGradient& operator*=(double alpha); //!< scalar multiply
};

//Functions required for minimize
void axpy(double alpha, const ElecGradient& x, ElecGradient& y); //!< accumulate operation: y += alpha*x
double dot(const ElecGradient& x, const ElecGradient& y, double* auxContrib=0); //!< inner product (optionally retrieve auxiliary contribution)
ElecGradient clone(const ElecGradient& x); //!< create a copy
void randomize(ElecGradient& x); //!< Initialize to random numbers

//! Variational total energy minimizer for electrons
class ElecMinimizer : public Minimizable<ElecGradient>
{
public:
	ElecMinimizer(Everything& e);
	
	//Virtual functions from Minimizable:
	void step(const ElecGradient& dir, double alpha);
	double compute(ElecGradient* grad, ElecGradient* Kgrad);
	bool report(int iter);
	void constrain(ElecGradient&);
	double sync(double x) const; //!< All processes minimize together; make sure scalars are in sync to round-off error
	
private:
	Everything& e;
	class ElecVars& eVars;
	const ElecInfo& eInfo;
	std::vector<matrix> KgradHaux; //!< latest preconditioned auxiliary gradient
	std::vector<matrix> rotPrev; //!< cumulated unitary rotations of subspace
	std::vector<matrix> rotPrevC; //!< cumulated transormation of wavefunctions (including non-unitary orthonormalization components)
	std::vector<matrix> rotPrevCinv; //!< inverse of rotPrevC (which is not just dagger, since these are not exactly unitary)
	
	bool rotExists; //!< whether rotPrev is non-trivial (not identity)
	std::shared_ptr<struct SubspaceRotationAdjust> sra; //!< Subspace rotation adjustment helper
};

void bandMinimize(Everything& e); //!< band structure minimization
void elecMinimize(Everything& e); //!< minimize electonic system
void elecFluidMinimize(Everything& e); //!< minimize electrons and fluid in a gummel loop if necessary
void convergeEmptyStates(Everything& e); //!< run bandMinimize to converge empty states (usually called from SCF / total energy calculations)

//! @}
#endif // JDFTX_ELECTRONIC_ELECMINIMIZER_H
