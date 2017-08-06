/*-------------------------------------------------------------------
Copyright 2014 Ravishankar Sundararaman

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


#ifndef JDFTX_ELECTRONIC_BANDDAVIDSON_H
#define JDFTX_ELECTRONIC_BANDDAVIDSON_H

#include <core/Minimize.h>

class Everything;

//! @addtogroup ElecSystem
//! @{

//! Davidson eigensolver
class BandDavidson
{
public:
	BandDavidson(Everything& e, int q); //!< Construct Davidson eigenvalue solver for quantum number q
	void minimize(); //!< Converge eigenproblem with tolerance set by e.elecMinParams
	
private:
	Everything& e;
	class ElecVars& eVars;
	const class ElecInfo& eInfo;
	int q;  //!< Current quantum number
};

//! @}
#endif // JDFTX_ELECTRONIC_BANDDAVIDSON_H