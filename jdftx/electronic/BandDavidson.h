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

#include <electronic/common.h>
#include <core/Minimize.h>
#include <electronic/ColumnBundle.h>

class BandDavidson
{	
	public:
		BandDavidson(Everything& e, int qActive);
		void minimize(); //!< converge eigenvalues till trace(Hsub); tolerance picked up from e.elecMinParams
		
		int qActive;  //!< Quantum number of the subspace that is being minimized
	private:
		Everything& e;
		ElecVars& eVars;
		ElecInfo& eInfo;
};

#endif // JDFTX_ELECTRONIC_BANDDAVIDSON_H