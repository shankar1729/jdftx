/*-------------------------------------------------------------------
Copyright 2012 Deniz Gunceler, Ravishankar Sundararaman

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

#include <electronic/LatticeMinimizer.h>
#include <core/Random.h>

//Functions required by Minimizable<matrix3<>>
void axpy(double alpha, const matrix3<>& x, matrix3<>& y) { y += alpha * x; }
double dot(const matrix3<>& x, const matrix3<>& y) { return trace(x*(~y)); }
matrix3<> clone(const matrix3<>& x) { return x; }
void randomize(matrix3<>& x) { for(int i=0; i<3; i++) for(int j=0; j<3; j++) x(i,j) = Random::normal(); }

//-------------  class LatticeMinimizer -----------------


LatticeMinimizer::LatticeMinimizer(Everything& e) : e(e)
{
}

void LatticeMinimizer::step(const matrix3<>& dir, double alpha)
{
}

double LatticeMinimizer::compute(matrix3<>* grad)
{
	return 0.;
}


bool LatticeMinimizer::report(int iter)
{	return false;
}

void LatticeMinimizer::constrain(matrix3<>&)
{
}
