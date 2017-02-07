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

#include <core/GridInfo.h>
#include <core/Data.h>
#include <cmath>
#include <gsl/gsl_sf.h>

GridInfo::GridInfo(GridInfo::CoordinateSystem coord, int S, double hMean)
: coord(coord), S(S), rMax(S*hMean), r(S), G(S), w(S), wTilde(S)
{
	switch(coord)
	{
		case Spherical:
		{
			//Compute roots of spherical Bessel functions
			std::vector<double> x(S+1), y(S+1);
			for(int i=0; i<=S; i++)
			{	if(i) x[i] = gsl_sf_bessel_zero_Jnu(0.5, i); //j0(x) = sqrt(pi/2x) * J_0.5(x)
				y[i] = gsl_sf_bessel_zero_Jnu(1.5, i); //j0'(x) = -j1(x) = -sqrt(pi/2x) * J_1.5(x)
			}
			//Initialize basis and quadrature grid
			for(int i=0; i<S; i++)
			{	r[i] = rMax*x[i+1]/y[S];
				G[i] = y[i]/rMax;
				w[i] =  4 * pow(M_PI/gsl_sf_bessel_j1(x[i+1]),2) * pow(rMax/y[S],3);
				wTilde[i] = 1. / ((i ? 2 : 4.0/3) * M_PI * pow(rMax,3) * pow(gsl_sf_bessel_j0(y[i]),2));
			}
			//Setup transform matrices
			matI.resize(S*S); auto elemI = matI.begin();
			matID.resize(S*S); auto elemID = matID.begin();
			matIDD.resize(S*S); auto elemIDD = matIDD.begin();
			for(int j=0; j<S; j++) for(int i=0; i<S; i++)
			{	double Gr = r[j] * G[i];
				*(elemI++) = gsl_sf_bessel_j0(Gr);
				*(elemID++) = -G[i] * gsl_sf_bessel_j1(Gr);
				*(elemIDD++) = pow(G[i],2) * gsl_sf_bessel_j2(Gr); // j0''-j0'/r
			}
			break;
		}
		
		case Cylindrical:
		{
			//Compute roots of Bessel functions
			std::vector<double> X(S+1), Y(S+1);
			for(int i=0; i<=S; i++)
			{	if(i) X[i] = gsl_sf_bessel_zero_J0(i);
				Y[i] = gsl_sf_bessel_zero_J1(i); //J0'(x) = -J1(x)
			}
			//Initialize basis and quadrature grid
			for(int i=0; i<S; i++)
			{	r[i] = rMax*X[i+1]/Y[S];
				G[i] = Y[i]/rMax;
				w[i] =  4*M_PI * pow(rMax/(Y[S]*gsl_sf_bessel_J1(X[i+1])), 2);
				wTilde[i] = 1. / (M_PI * pow(rMax*gsl_sf_bessel_J0(Y[i]), 2));
			}
			//Setup transform matrices
			matI.resize(S*S); auto elemI = matI.begin();
			matID.resize(S*S); auto elemID = matID.begin();
			matIDD.resize(S*S); auto elemIDD = matIDD.begin();
			for(int j=0; j<S; j++) for(int i=0; i<S; i++)
			{	double Gr = r[j] * G[i];
				*(elemI++) = gsl_sf_bessel_J0(Gr);
				*(elemID++) = -G[i] * gsl_sf_bessel_J1(Gr);
				*(elemIDD++) = pow(G[i],2) * (3*gsl_sf_bessel_Jn(2,Gr)-gsl_sf_bessel_J0(Gr))/4; //J0''-J0'/2r
			}
			break;
		}
		
		case Planar:
		{
			//Initialize basis and quadrature grid
			for(int i=0; i<S; i++)
			{	r[i] = (i+0.5)*hMean;
				G[i] = i * (M_PI/rMax);
				w[i] = hMean;
				wTilde[i] = 2./((i ? 1 : 2) * rMax);
			}
			//Setup transform plans
			ScalarField temp(this);
			ScalarFieldTilde tempTilde(this);
			planPlanarIdag  = fftw_plan_r2r_1d(S, temp.data(), tempTilde.data(), FFTW_REDFT10, FFTW_MEASURE); //DCT-II  (DCT)
			planPlanarI     = fftw_plan_r2r_1d(S, tempTilde.data(), temp.data(), FFTW_REDFT01, FFTW_MEASURE); //DCT-III (IDCT)
			planPlanarIDdag = fftw_plan_r2r_1d(S, temp.data(), tempTilde.data(), FFTW_RODFT10, FFTW_MEASURE); //DST-II  (DST)
			planPlanarID    = fftw_plan_r2r_1d(S, tempTilde.data(), temp.data(), FFTW_RODFT01, FFTW_MEASURE); //DST-III (IDST)
			break;
		}
	}
}

GridInfo::~GridInfo()
{
	switch(coord)
	{
		case Spherical:
		case Cylindrical:
			break;
			
		case Planar:
			fftw_destroy_plan(planPlanarI);
			fftw_destroy_plan(planPlanarIdag);
			fftw_destroy_plan(planPlanarID);
			fftw_destroy_plan(planPlanarIDdag);
			break;
	}
}

double GridInfo::Volume() const
{	switch(coord)
	{	case Spherical: return (4*M_PI/3) * pow(rMax,3);
		case Cylindrical: return M_PI * pow(rMax,2);
		case Planar: return rMax;
	}
	return 0.;
}
