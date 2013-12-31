/*-------------------------------------------------------------------
Copyright 2013 Deniz Gunceler

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
		
		//! Updates fillings and recomputes filling energies
		void updateFillings();
		
		//! Mixes densities or Vscloc and Vtau
		void mixPlain(DataRptrCollection& variable_n, DataRptrCollection& variable_tau, 
							DataRptrCollection& prevVariable_n, DataRptrCollection& prevVariable_tau, double mixFraction=0.5);
		
		//! Uses direct inversion in the iterative subspace to extrapolate to a new density/potential
		void mixDIIS(DataRptrCollection& variable_n, DataRptrCollection& variable_tau, 
				  std::vector<DataRptrCollection>& pastVariables_n, std::vector<DataRptrCollection>& pastVariables_tau,
				  std::vector<DataRptrCollection>& pastResiduals, std::vector< DataRptrCollection >& pastResiduals_tau);
};

#endif