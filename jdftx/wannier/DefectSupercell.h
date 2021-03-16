/*-------------------------------------------------------------------
Copyright 2021 Ravishankar Sundararaman

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

#ifndef JDFTX_WANNIER_DEFECT_SUPERCELL_H
#define JDFTX_WANNIER_DEFECT_SUPERCELL_H

#include <core/vector3.h>
#include <core/string.h>
#include <memory>

class DefectSupercell
{
public:
	string name; //defect name
	string inFile; //input file of supercell calculation
	vector3<int> sup; //supercell of defect relative to unit cell
	vector3<> xCenter; //lattice coordinates of supercell where defect is centered
	double q; //net electron count of defect
	
	void initialize(const class Wannier*);
	
private:
	const class Wannier* wannier;
	const class Everything* e;
	std::shared_ptr<class Everything> eSup;
};

#endif //JDFTX_WANNIER_DEFECT_SUPERCELL_H
