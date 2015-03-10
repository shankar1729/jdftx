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
#include <core/Pulay.h>

//! @brief Variable that is mixed during SCF
//! Component names are density-like, but when mixing potential, they refer to corresponding gradient.
struct SCFvariable
{	ScalarFieldArray n; //!< electron density (or potential)
	ScalarFieldArray tau; //!< KE density (or potential) [mGGA only]
	std::vector<matrix> rhoAtom; //!< atomic density matrices (or corresponding potential) [DFT+U only]
};

//! @brief Self-Consistent Field method for converging electronic state
class SCF : public Pulay<SCFvariable>
{
public:
	SCF(Everything& e);
	
	//! Minimizes residual to achieve self-consistency
	void minimize();
	
	static double eigDiffRMS(const std::vector<diagMatrix>&, const std::vector<diagMatrix>&, const Everything& e); //!< weigted RMS difference between two sets of eigenvalues
	
protected:
	//---- Interface to Pulay ----
	double sync(double x) const;
	double cycle(double dEprev, std::vector<double>& extraValues);
	void report(int iter);
	void axpy(double alpha, const SCFvariable& X, SCFvariable& Y) const;
	double dot(const SCFvariable& X, const SCFvariable& Y) const;
	size_t variableSize() const;
	void readVariable(SCFvariable&, FILE*) const;
	void writeVariable(const SCFvariable&, FILE*) const;
	SCFvariable getVariable() const;
	void setVariable(const SCFvariable&);
	SCFvariable precondition(const SCFvariable&) const;
	SCFvariable applyMetric(const SCFvariable&) const;

private:
	Everything& e;
	bool mixTau; //!< whether KE needs to be mixed
	RealKernel kerkerMix, diisMetric; //!< convolution kernels for kerker preconditioning and the DIIS overlap metric
	
	//! Updates fillings and recomputes filling energies
	void updateFillings();
	//! Updates the fillings using the maximum overlap method
	void updateMOM();
	std::vector<ColumnBundle> Cold; // Wavefunctions from the previous iteration
	
	double eigDiffRMS(const std::vector<diagMatrix>&, const std::vector<diagMatrix>&) const; //!< weighted RMS difference between two sets of eigenvalues
	void eigenShiftInit(); //!< initialize and check eigenShifts
	void eigenShiftApply(bool reverse);  //!< apply eigenshifts if false, and undo them if true
};

#endif