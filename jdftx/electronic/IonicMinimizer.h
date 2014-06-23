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

#ifndef JDFTX_ELECTRONIC_IONICMINIMIZER_H
#define JDFTX_ELECTRONIC_IONICMINIMIZER_H

#include <electronic/common.h>
#include <electronic/RadialFunction.h>
#include <core/Minimize.h>
#include <core/vector3.h>
#include <core/matrix3.h>
#include <vector>

//! Object to hold all the forces
struct IonicGradient : std::vector< std::vector< vector3<> > >
{
	void init(const IonInfo&); //!< initialize to zeroes with the correct species and atom numbers for iInfo
	void print(const Everything&, FILE*, const char* prefix="force") const;
	
	IonicGradient& operator*=(double);
	IonicGradient& operator+=(const IonicGradient&);
	
	IonicGradient operator*(double) const;
	IonicGradient operator+(const IonicGradient&) const;
	IonicGradient operator-(const IonicGradient&) const;
};

void axpy(double alpha, const IonicGradient& x, IonicGradient& y); //!< accumulate operation: Y += alpha*X
double dot(const IonicGradient& x, const IonicGradient& y); //!< inner product
IonicGradient clone(const IonicGradient& x); //!< create a copy
void randomize(IonicGradient& x); //!< initialize with random numbers

IonicGradient operator*(const matrix3<>&, const IonicGradient&); //!< coordinate transformations


class IonicMinimizer : public Minimizable<IonicGradient>
{	Everything& e;
	
public:
	IonicMinimizer(Everything& e);
	//Virtual functions from Minimizable:
	void step(const IonicGradient& dir, double alpha);
	double compute(IonicGradient* grad);
	IonicGradient precondition(const IonicGradient& grad);
	bool report(int iter);
	void constrain(IonicGradient&);
	double sync(double x) const; //!< All processes minimize together; make sure scalars are in sync to round-off error
	
	double minimize(const MinimizeParams& params); //!< minor addition to Minimizable::minimize to invoke charge analysis at final positions
private:
	bool populationAnalysisPending; //!< report() has requested a charge analysis output that is yet to be done
};

#endif // JDFTX_ELECTRONIC_IONICMINIMIZER_H
