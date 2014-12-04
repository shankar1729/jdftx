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

#ifndef JDFTX_PHONON_PHONON_H
#define JDFTX_PHONON_PHONON_H

#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>

class Phonon
{
public:
	Everything e, eSup; //data for original unit cell and supercell
	vector3<int> sup; //phonon supercell 
	
	void setup();
};

#endif //JDFTX_PHONON_PHONON_H