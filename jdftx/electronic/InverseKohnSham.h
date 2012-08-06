/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_INVERSEKOHNSHAM_H
#define JDFTX_ELECTRONIC_INVERSEKOHNSHAM_H

#include <core/Minimize.h>
#include <core/DataCollection.h>
#include <electronic/common.h>

class InverseKohnSham : public Minimizable<DataRptrCollection>
{
public:
	InverseKohnSham(Everything&);
	
	//Interface for minimize:
	void step(const DataRptrCollection& dir, double alpha);
	double compute(DataRptrCollection* grad);
	DataRptrCollection precondition(const DataRptrCollection& grad);
	void constrain(DataRptrCollection& dir);
	bool report(int iter);
private:
	Everything& e;
	DataRptrCollection n; //band structure density (should equal ElecVars::n at end)
	RealKernel gaussCutoff; //gaussian cutoff kernel (of width Control::invertKS_sigma)
	//Chi guess:
	std::shared_ptr<class InvertChi> invertChi; //Compute inverse chi using linear CG
};

#endif // JDFTX_ELECTRONIC_INVERSEKOHNSHAM_H
