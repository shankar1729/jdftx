/*-------------------------------------------------------------------
Copyright 2013 Deniz Gunceler, Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_SCF_H
#define JDFTX_ELECTRONIC_SCF_H

#include <electronic/common.h>
#include <electronic/Everything.h>
#include <electronic/BandMinimizer.h>
#include <electronic/operators.h>

//! Self-Consistent Iteration for residual minimization
class SCF
{
public:
	SCF(Everything& e);
	
	//! Minimizes residual to achieve self-consistency
	void minimize();
private:
	Everything& e;
	
	matrix overlap; //! Overlap matrix of density/potential residuals
	
	//Variable that's being mixed (component names are density-like, but when mixing potential, they refer to corresponding gradient)
	struct Variable
	{	DataRptrCollection n; //electron density (or potential)
		DataRptrCollection tau; //KE density (or potential) [mGGA only]
		std::vector<matrix> rhoAtom; //atomic density matrices (or corresponding potential) [DFT+U only]
	};
	std::vector<Variable> pastVariables, pastResiduals; //!< History
	bool mixTau; //!< whether KE needs to be mixed
	
	void axpy(double alpha, const Variable& X, Variable& Y) const; // Y += alpha * X
	double dot(const Variable& X, const Variable& Y) const; //Euclidean dot product (Note metric not applied here as an optimization, as many dot products against a single residual computed)
	
	size_t variableSize() const; //!< number of bytes per variable
	void readVariable(Variable&, FILE*) const; //! read variable from file
	void writeVariable(const Variable&, FILE*) const; //! write variable to file
	Variable getVariable() const; //!< get the current variable (from ElecVars)
	void setVariable(const Variable&); //!< set the current variable (to ElecVars)
	
	RealKernel kerkerMix, diisMetric; //!< convolution kernels for kerker preconditioning and the DIIS overlap metric
	Variable applyKerker(const Variable&) const;
	Variable applyMetric(const Variable&) const;
	
	//! Pulay mixing / direct inversion in iterative subspace
	void mixDIIS();
	
	//! Updates fillings and recomputes filling energies
	void updateFillings();
	//! Updates the fillings using the maximum overlap method
	void updateMOM();
	std::vector<ColumnBundle> Cold; // Wavefunctions from the previous iteration
	
	double eigDiffRMS(const std::vector<diagMatrix>&, const std::vector<diagMatrix>&) const; //!< weigted RMS difference between two sets of eigenvalues
	void eigenShiftInit(); //!< initialize and check eigenShifts
	void eigenShiftApply(bool reverse);  //!< apply eigenshifts if false, and undo them if true

	//! Applies the single particle constraint on Vscloc
	void single_particle_constraint(double sp_constraint);
	
};

#endif