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

#include <core/Minimize.h>

class ColumnBundle;
class Everything;

//! @addtogroup ElecSystem
//! @{

//! Conjugate gradients eigen-solver
class BandMinimizer : public Minimizable<ColumnBundle>
{
public:
	BandMinimizer(Everything& e, int q); //!< Construct band-structure minimizer for quantum number q

	//Interface for Minimizable:
	double compute(ColumnBundle* grad, ColumnBundle* Kgrad);
	void step(const ColumnBundle& dir, double alpha);
	void constrain(ColumnBundle&);

private:
	Everything& e;
	class ElecVars& eVars;
	const class ElecInfo& eInfo;
	int q; //!< Current quantum number
};

//! @}
#endif
