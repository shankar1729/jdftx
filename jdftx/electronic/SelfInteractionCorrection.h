/*-------------------------------------------------------------------#include <electronic/common.h>
Copyright 2012 Deniz Gunceler

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

#ifndef JDFTX_ELECTRONIC_SELFINTERACTIONCORRECTION
#define JDFTX_ELECTRONIC_SELFINTERACTIONCORRECTION

#include <electronic/common.h>
#include <core/Data.h>
#include <vector>

class SelfInteractionCorrection
{
	public:
		SelfInteractionCorrection(const Everything& everything);
		~SelfInteractionCorrection();
		
		double operator()(std::vector<diagMatrix>* correctedEigenvalues); //! Evaluates the self interaction energy and (optionally) returns the corrected band eigenvalues
		
		double coulombExciton(int q1, int n1, int q2, int n2);  //! Approximates the coulomb excitonic contribution between two states
		
		void dump(const char* filenamePattern);
		
		bool needsTau;  //! The kinetic energy density is needed for meta-gga functionals.
		
	private:
		const Everything* e;
		
		double calcSelfInteractionError(int q, int n); //! Calculates the self-interaction error of the KS orbital atthe n'th band at q'th quantum number
		
		//! ColumnBundle for the derivative of the wavefunctions in each cartesian direction
		std::vector<ColumnBundle> DC;
	
};

#endif // JDFTX_ELECTRONIC_SELFINTERACTIONCORRECTION
