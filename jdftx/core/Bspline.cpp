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

#include <core/Bspline.h>

Bspline::Bspline(double xStart, double xEnd, int nIntervals):nCoeff(nIntervals),xStart(xStart)
{	h = (xEnd-xStart)/nCoeff;
	invh = 1.0/h;
	coeff = new double[4*nCoeff];
}

Bspline::~Bspline()
{	delete[] coeff;
}

void Bspline::init()
{	double x = xStart;
	double f = getValue(x);
	double fp = h*getDeriv(x);
	double* a = coeff;

	for(int i=0; i<nCoeff; i++)
	{	double xNext = x + h;
		double fNext = getValue(xNext);
		double fpNext = h*getDeriv(xNext);

		a[0] = f;
		a[1] = fp;
		a[2] = 3*(fNext-f) - fpNext - 2*fp;
		a[3] = 2*(f-fNext) + fpNext + fp;

		x = xNext;
		f = fNext;
		fp = fpNext;
		a += 4;
	}
}

void Bspline::maxErrorEstimate(FILE* fp, const char* plotFileName)
{	double maxRelErr=0.0, maxAbsErr=0.0; int iMaxRelErr=0, iMaxAbsErr=0;
	double maxRelErrP=0.0, maxAbsErrP=0.0; int iMaxRelErrP=0, iMaxAbsErrP=0;

	for(int i=0; i<nCoeff; i++)
	{	double x = xStart+(i+0.5)*h;

		double err = fabs(getValue(x)-value(x));
		double relerr = fabs(err/getValue(x));
		if(err > maxAbsErr) { maxAbsErr = err; iMaxAbsErr=i; }
		if(relerr > maxRelErr) { maxRelErr = relerr; iMaxRelErr=i; }

		double errP = fabs(getDeriv(x)-deriv(x));
		double relerrP = fabs(err/getDeriv(x));
		if(errP > maxAbsErrP) { maxAbsErrP = errP; iMaxAbsErrP=i; }
		if(relerrP > maxRelErrP) { maxRelErrP = relerrP; iMaxRelErrP=i; }
	}
	fprintf(fp, "Max value error: abs=%le at %d, rel=%le at %d\n", maxAbsErr, iMaxAbsErr, maxRelErr, iMaxRelErr);
	fprintf(fp, "Max deriv error: abs=%le at %d, rel=%le at %d\n", maxAbsErrP, iMaxAbsErrP, maxRelErrP, iMaxRelErrP);

	if(plotFileName)
	{	FILE* fpPlot = fopen(plotFileName, "w");
		for(double x=xStart; x<xStart+h*nCoeff; x+=0.1*h)
			fprintf(fpPlot, "%le\t%le\t%le\n", x, getValue(x), value(x));
		fclose(fpPlot);
	}
}
