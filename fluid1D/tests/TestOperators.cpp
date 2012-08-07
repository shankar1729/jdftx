/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

This file is part of Fluid1D.

Fluid1D is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Fluid1D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Fluid1D.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#include <core1D/Operators.h>
#include <core/Util.h>
#include <cmath>

int main(int argc, char** argv)
{	initSystem(argc, argv);

	GridInfo gInfo(GridInfo::Spherical, 256, 0.25);
	
	{	puts("\nTest 1: Norm, dot product");
		ScalarField r1(&gInfo);
		initRandom(r1);
		printf("\tNorm of random vector = %le (expectation value: %le)\n", nrm2(r1), sqrt(gInfo.S));
		printf("\tSelf dot product = %le (expectation value: %le)\n", dot(r1, r1), double(gInfo.S));

		puts("\nTest 2: Linear combines");
		ScalarField r2 = 3.0*r1 - r1;
		ScalarField r3 = r1 + r2;
		r2 -= 2.0*((r3 - 0.5*r1*2.0) - r1);
		printf("\tLinear combination norm = %le (should be 0 within roundoff)\n",  nrm2(r2));

		puts("\nTest 3: Transform tests");
		initRandom(r2);
		ScalarFieldTilde g1 = J(r1);
		printf("\tTransform inverse relative error: %le\n", nrm2(I(g1)-r1) / nrm2(r1));
		printf("\tRelative error between JA.O(JA) and A.JdagOJ(A): %le\n", dot(g1,O(g1))/dot(r1,JdagOJ(r1))-1.);
		printf("\tRelative error between A.Idag(B) and I(A).B: %le\n", dot(g1,Idag(r2))/dot(r2,I(g1))-1.);
		printf("\tRelative error between A.IDdag(B) and ID(A).B: %le\n", dot(g1,IDdag(r2))/dot(r2,ID(g1))-1.);
		printf("\tRelative error between A.IDDdag(B) and IDD(A).B: %le\n", dot(g1,IDDdag(r2))/dot(r2,IDD(g1))-1.);
		printf("\tRelative error between A.Jdag(B) and J(A).B: %le\n", dot(g1,J(r2))/dot(r2,Jdag(g1))-1.);
	}

	{	puts("\nTest 4: Accuracy of the ID and IDD operators:");
		ScalarField x(&gInfo), Dx(&gInfo), DDx(&gInfo);
		double *xData = x.data(), *DxData = Dx.data(), *DDxData = DDx.data();
		for(int i=0; i<gInfo.S; i++)
		{	double t = gInfo.r[i] / gInfo.rMax;
			xData[i] = t*t*(3 - 2*t); //a test function
			DxData[i] = 6*t*(1-t) / gInfo.rMax; //it's analytic vector derivative
			//it's analytic 'special' tensor second derivative:
			switch(gInfo.coord)
			{	case GridInfo::Spherical:   DDxData[i] = (-6*t) / pow(gInfo.rMax,2); break;
				case GridInfo::Cylindrical: DDxData[i] = (3-9*t) / pow(gInfo.rMax,2); break;
				case GridInfo::Planar:      DDxData[i] = (6-12*t) / pow(gInfo.rMax,2); break;
			}
		}
		ScalarFieldTilde xTilde = J(x);
		ScalarField numDx = ID(xTilde), numDDx = IDD(xTilde);
		ScalarField errDx = numDx - Dx, errDDx = numDDx - DDx;
		printf("\tRelative error between ID and analytic derivative of cubic spline: %le\n", sqrt(dot(errDx,JdagOJ(errDx))/dot(Dx,JdagOJ(Dx))));
		printf("\tRelative error between IDD and analytic second derivative of cubic spline: %le\n", sqrt(dot(errDDx,JdagOJ(errDDx))/dot(DDx,JdagOJ(DDx))));
		double *numDxData = numDx.data(), *numDDxData = numDDx.data();
		FILE* fp = fopen("test-derivative-transform", "w");
		for(int i=0; i<gInfo.S; i++)
			fprintf(fp, "%lf\t%le\t%le\t%le\t%le\t%le\n", gInfo.r[i], xData[i], DxData[i], numDxData[i], DDxData[i], numDDxData[i]);
		fclose(fp);
	}
}
