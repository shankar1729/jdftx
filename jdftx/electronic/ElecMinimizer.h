/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

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

#ifndef JDFTX_ELECTRONIC_ELECMINIMIZER_H
#define JDFTX_ELECTRONIC_ELECMINIMIZER_H

#include <electronic/common.h>
#include <electronic/ColumnBundle.h>
#include <electronic/matrix.h>
#include <core/Minimize.h>

struct ElecGradient
{
	std::vector<ColumnBundle> Y; 
	std::vector<matrix> B;
	
	void init(Everything& e); //!< initialize Y and B with the correct sizes for everything
	void init(int nStates); //!< initialize Y and B arrays, but with empty matrices and ColumnBundles
	
	ElecGradient& operator*=(double alpha); //!< scalar multiply
};

//Functions required for minimize
void axpy(double alpha, const ElecGradient& x, ElecGradient& y); //!< accumulate operation: Y += alpha*X
double dot(const ElecGradient& x, const ElecGradient& y); //!< inner product
ElecGradient clone(const ElecGradient& x); //!< create a copy
void randomize(ElecGradient& x); //!< Initialize to random numbers

class ElecMinimizer : public Minimizable<ElecGradient>
{
	Everything& e;
	ElecVars& eVars;
	ElecInfo& eInfo;
	bool precond;
	ElecGradient Kgrad;
public:
	ElecMinimizer(Everything& e, bool precond=true);
	
	//Virtual functions from Minimizable:
	void step(const ElecGradient& dir, double alpha);
	double compute(ElecGradient* grad);
	ElecGradient precondition(const ElecGradient& grad);
	bool report(int iter);
	void constrain(ElecGradient&);
};

void elecMinimize(Everything& e); //!< minimize electonic system
void elecFluidMinimize(Everything& e); //!< minimize electrons and fluid in a gummel loop if necessary

#endif // JDFTX_ELECTRONIC_ELECMINIMIZER_H
