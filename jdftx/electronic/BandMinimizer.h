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

#ifndef JDFTX_ELECTRONIC_BANDMINIMIZER_H
#define JDFTX_ELECTRONIC_BANDMINIMIZER_H

#include <electronic/common.h>
#include <core/Minimize.h>
#include <electronic/ColumnBundle.h>

class BandMinimizer : public Minimizable<ColumnBundle>
{	
	public:
		BandMinimizer(Everything& e, int qActive, bool precond=true);
	
		double compute(ColumnBundle* grad);
		void step(const ColumnBundle& dir, double alpha);
		ColumnBundle precondition(const ColumnBundle& grad);
		bool report(int iter);
		void constrain(ColumnBundle&);
		
		int qActive;  //! Quantum number of the subspace that is being minimized

		ColumnBundle Kgrad; //! Preconditioned gradient
	private:
		Everything& e;
		ElecVars& eVars;
		ElecInfo& eInfo;
		
		bool precond;	
};

#endif