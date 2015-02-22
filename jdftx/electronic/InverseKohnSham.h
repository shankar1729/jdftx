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
#include <core/ScalarFieldArray.h>
#include <electronic/common.h>

class InverseKohnSham : public Minimizable<ScalarFieldArray>
{
public:
	InverseKohnSham(Everything&);
	
	//Interface for minimize:
	void step(const ScalarFieldArray& dir, double alpha);
	double compute(ScalarFieldArray* grad);
	ScalarFieldArray precondition(const ScalarFieldArray& grad);
	void constrain(ScalarFieldArray& dir);
	bool report(int iter);
	double sync(double x) const; //!< All processes minimize together; make sure scalars are in sync to round-off error
private:
	Everything& e;
	ScalarFieldArray n; //band structure density (should equal ElecVars::n at end)
	RealKernel gaussCutoff; //gaussian cutoff kernel (of width Control::invertKS_sigma)
	//Chi guess:
	std::shared_ptr<class InvertChi> invertChi; //Compute inverse chi using linear CG
};

#endif // JDFTX_ELECTRONIC_INVERSEKOHNSHAM_H
