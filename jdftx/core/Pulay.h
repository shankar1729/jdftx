/*-------------------------------------------------------------------
Copyright 2015 Ravishankar Sundararaman

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

#ifndef JDFTX_CORE_PULAY_H
#define JDFTX_CORE_PULAY_H

#include <core/PulayParams.h>
#include <core/matrix.h>
#include <core/string.h>
#include <cfloat>

//! @addtogroup Algorithms
//! @{

//! @brief Pulay mixing to optimize self-consistent field optimization
//! @
template<typename Variable> class Pulay
{
public:
	Pulay(const PulayParams& pp);
	
	//! @brief Minimize energy using a self-consistent iteration
	//! @param Eprev Initial energy (optional)
	//! @param extraNames Names of extra convergence parameters to report and check every cycle
	//! @param extraThresh Thresholds for each extra convergence parameter
	//! @return Energy at final cycle
	double minimize(double Eprev=+DBL_MAX, std::vector<string> extraNames=std::vector<string>(), std::vector<double> extraThresh=std::vector<double>());
	
	void loadState(const char* filename); //!< Load the state from a single binary file
	void saveState(const char* filename) const; //!< Save the state to a single binary file
	void clearState(); //!< remove past variables and residuals
	
	//! Override to synchronize scalars over MPI processes (if the same minimization is happening in sync over many processes)
	virtual double sync(double x) const { return x; }
	
protected:
	//----- Interface specification -----

	//! @brief Single cycle of the self-consistency loop.
	//! In each subsequent cycle, Pulay will try to zero the difference
	//! between getVariable() before and after the cycle.
	//! The implementation should only do the work of computing the updated variable;
	//! debug printing and I/O, if any, should occur in report() instead.
	//! @param dEprev energy change at previous cycle, which may be used to adjust accuracy of any inner optimizations
	//! @param extraValues Specify the values of any extra convergence parameters specified in minimize()
	//! @return Energy at end of cycle
	virtual double cycle(double dEprev, std::vector<double>& extraValues)=0;
	
	virtual void report(int iter) {} //!< Override to perform optional reporting
	virtual void axpy(double alpha, const Variable& X, Variable& Y) const=0; //!< Scaled accumulate on variable
	virtual double dot(const Variable& X, const Variable& Y) const=0; //!< Euclidean dot product. Metric applied separately for efficiency.
	virtual size_t variableSize() const=0; //!< Number of bytes per variable
	virtual void readVariable(Variable&, FILE*) const=0; //!< Read variable from stream
	virtual void writeVariable(const Variable&, FILE*) const=0; //! Write variable to stream
	virtual Variable getVariable() const=0; //!< Get the current variable from state of system
	virtual void setVariable(const Variable&)=0; //!< Set the state of system to specified variable
	virtual Variable precondition(const Variable&) const=0; //!< Apply preconditioner to variable/residual
	virtual Variable applyMetric(const Variable&) const=0; //!< Apply metric to variable/residual

private:
	const PulayParams& pp; //!< Pulay parameters
	std::vector<Variable> pastVariables; //!< Previous variables
	std::vector<Variable> pastResiduals; //!< Previous residuals
	matrix overlap; //!< Overlap matrix of residuals
};

//! @}

//---------------------- Implementation ----------------------------
//!@cond

#include <core/Minimize.h>
#include <memory>

//Norm convergence check (eigenvalue-difference or residual)
//Make sure value is within tolerance for nCheck consecutive cycles
class NormCheck
{	unsigned nCheck; double threshold;
	std::deque<bool> history;
public:
	NormCheck(unsigned nCheck, double threshold) : nCheck(nCheck), threshold(fabs(threshold)) {}
	bool checkConvergence(double norm)
	{	history.push_back(fabs(norm)<threshold);
		if(history.size()==nCheck+1) history.pop_front(); //discard old unneeded elements 
		if(history.size()==nCheck)
		{	for(bool converged: history)
				if(!converged)
					return false;
			return true;
		}
		else return false;
	}
};

template<typename Variable> Pulay<Variable>::Pulay(const PulayParams& pp)
: pp(pp), overlap(pp.history, pp.history)
{
}

template<typename Variable> double Pulay<Variable>::minimize(double Eprev, std::vector<string> extraNames, std::vector<double> extraThresh)
{
	double E = sync(Eprev); Eprev = 0.;
	double dE = E-Eprev;
	assert(extraNames.size()==extraThresh.size());
	
	//Initialize convergence checkers:
	EdiffCheck ediffCheck(2, pp.energyDiffThreshold); ediffCheck.checkConvergence(E); //store the initial energy in the check's history
	NormCheck resCheck(2, pp.residualThreshold);
	std::vector<std::shared_ptr<NormCheck> > extraCheck(extraNames.size());
	for(size_t iExtra=0; iExtra<extraNames.size(); iExtra++)
		extraCheck[iExtra] = std::make_shared<NormCheck>(2, extraThresh[iExtra]);

	for(int iter=0; iter<pp.nIterations; iter++)
	{
		//If history is full, remove oldest member
		assert(pastResiduals.size() == pastVariables.size());
		if((int)pastResiduals.size() >= pp.history)
		{	size_t ndim = pastResiduals.size();
			if(ndim>1) overlap.set(0,ndim-1, 0,ndim-1, overlap(1,ndim, 1,ndim));
			pastVariables.erase(pastVariables.begin());
			pastResiduals.erase(pastResiduals.begin());
		}
		
		//Cache the old energy and variables
		Eprev = E;
		pastVariables.push_back(getVariable());

		//Perform cycle:
		std::vector<double> extraValues(extraThresh.size());
		E = sync(cycle(dE, extraValues));
		dE = E - Eprev;
		for(auto& v: extraValues) v = sync(v);
			
		//Calculate and cache residual:
		double residualNorm = 0.;
		{	Variable residual = getVariable(); axpy(-1., pastVariables.back(), residual);
			pastResiduals.push_back(residual);
			residualNorm = sync(sqrt(dot(residual,residual)));
		}
		
		//Print energy and convergence parameters:
		fprintf(pp.fpLog, "%sCycle: %2i   %s: ", pp.linePrefix, iter, pp.energyLabel);
		fprintf(pp.fpLog, pp.energyFormat, E);
		fprintf(pp.fpLog, "   d%s: %+.3e", pp.energyLabel, dE);
		fprintf(pp.fpLog, "   |Residual|: %.3e", residualNorm);
		for(size_t iExtra=0; iExtra<extraNames.size(); iExtra++)
			fprintf(pp.fpLog, "   |%s|: %.3e", extraNames[iExtra].c_str(), extraValues[iExtra]);
		fprintf(pp.fpLog, "  t[s]: %9.2lf", clock_sec());
		fprintf(pp.fpLog, "\n"); fflush(pp.fpLog);
		
		//Optional reporting:
		report(iter);
		
		//Check for convergence and update variable:
		if(std::isnan(E))
		{	fprintf(pp.fpLog, "%sE=%le. Stopping ...\n\n", pp.linePrefix, E);
			return E;
		}
		bool converged = false;
		if(!converged && ediffCheck.checkConvergence(E))
		{	fprintf(pp.fpLog, "%sConverged (|Delta E|<%le for 2 iters).\n\n", pp.linePrefix, pp.energyDiffThreshold);
			converged = true;
		}
		if(!converged && resCheck.checkConvergence(residualNorm))
		{	fprintf(pp.fpLog, "%sConverged (|Residual|<%le for 2 iters).\n\n", pp.linePrefix, pp.residualThreshold);
			converged = true;
		}
		if(!converged && extraNames.size())
			for(size_t iExtra=0; iExtra<extraNames.size(); iExtra++)
				if(extraCheck[iExtra]->checkConvergence(extraValues[iExtra]))
				{	fprintf(pp.fpLog, "%sConverged (|%s|<%le for 2 iters).\n\n", pp.linePrefix, extraNames[iExtra].c_str(), extraThresh[iExtra]);
					converged = true;
					break;
				}
		fflush(pp.fpLog);
		if(converged || killFlag) break; //converged or manually interrupted
		
		//---- DIIS/Pulay mixing -----
			
		//Update the overlap matrix
		size_t ndim = pastResiduals.size();
		Variable MlastResidual = applyMetric(pastResiduals.back());
		for(size_t j=0; j<ndim; j++)
		{	double thisOverlap = dot(pastResiduals[j], MlastResidual);
			overlap.set(j, ndim-1, thisOverlap);
			overlap.set(ndim-1, j, thisOverlap);
		}
		
		//Invert the residual overlap matrix to get the minimum of residual
		matrix cOverlap(ndim+1, ndim+1); //Add row and column to enforce normalization constraint
		cOverlap.set(0, ndim, 0, ndim, overlap(0, ndim, 0, ndim));
		for(size_t j=0; j<ndim; j++)
		{	cOverlap.set(j, ndim, 1);
			cOverlap.set(ndim, j, 1);
		}
		cOverlap.set(ndim, ndim, 0);
		matrix cOverlap_inv = inv(cOverlap);
		
		//Update variable:
		Variable v;
		for(size_t j=0; j<ndim; j++)
		{	double alpha = cOverlap_inv.data()[cOverlap_inv.index(j, ndim)].real();
			axpy(alpha, pastVariables[j], v);
			axpy(alpha, precondition(pastResiduals[j]), v);
		}
		setVariable(v);
	}
	return E;
}

template<typename Variable> void Pulay<Variable>::loadState(const char* filename)
{
	size_t nBytesCycle = 2 * variableSize(); //number of bytes per history entry
	size_t nBytesFile = fileSize(filename);
	size_t ndim = nBytesFile / nBytesCycle;
	size_t dimOffset = 0;
	if(int(ndim) > pp.history)
	{	dimOffset = ndim - pp.history;
		ndim = pp.history;
	}
	if(nBytesFile % nBytesCycle != 0)
		die("Pulay history file '%s' does not contain an integral multiple of the mixed variables and residuals.\n", filename);
	fprintf(pp.fpLog, "%sReading %lu past variables and residuals from '%s' ... ", pp.linePrefix, ndim, filename); logFlush();
	pastVariables.resize(ndim);
	pastResiduals.resize(ndim);
	FILE* fp = fopen(filename, "r");
	if(dimOffset) fseek(fp, dimOffset*nBytesCycle, SEEK_SET);
	for(size_t idim=0; idim<ndim; idim++)
	{	readVariable(pastVariables[idim], fp);
		readVariable(pastResiduals[idim], fp);
	}
	fclose(fp);
	fprintf(pp.fpLog, "done.\n"); fflush(pp.fpLog);
	//Compute overlaps of loaded history:
	for(size_t i=0; i<ndim; i++)
	{	Variable Mresidual_i = applyMetric(pastResiduals[i]);
		for(size_t j=0; j<=i; j++)
		{	double thisOverlap = dot(pastResiduals[j], Mresidual_i);
			overlap.set(i,j, thisOverlap);
			overlap.set(j,i, thisOverlap);
		}
	}
}

template<typename Variable> void Pulay<Variable>::saveState(const char* filename) const
{
	if(mpiWorld->isHead())
	{	FILE* fp = fopen(filename, "w");
		for(size_t idim=0; idim<pastVariables.size(); idim++)
		{	writeVariable(pastVariables[idim], fp);
			writeVariable(pastResiduals[idim], fp);
		}
		fclose(fp);
	}
}

template<typename Variable> void Pulay<Variable>::clearState()
{	pastVariables.clear();
	pastResiduals.clear();
}

//!@endcond

#endif //JDFTX_CORE_PULAY_H
