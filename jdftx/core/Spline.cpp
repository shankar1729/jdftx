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

#include <core/Spline.h>

namespace QuinticSpline
{
	std::vector<double> getCoeff(const std::vector<double>& samples, bool oddExtension)
	{	//Setup a penta-diagonal system of equations for quintic-Blip coefficeints:
		std::vector<double> x = samples; //destructible copy
		int N = x.size();
		std::vector<double> tmpBand[5]; std::vector<double>* band = tmpBand+2;
		for(int m=-2; m<=2; m++) band[m].resize(N-abs(m));
		const double off1 = 13./33;
		const double off2 = 1./66;
		for(int i=0; i<N; i++)   band[0][i] = 1.;
		for(int i=0; i<N-1; i++) band[1][i] = band[-1][i] = off1;
		for(int i=0; i<N-2; i++) band[2][i] = band[-2][i] = off2;
		//Mirror BC at rho=0:
		if(oddExtension) x[0] = 0.;
		double sign = oddExtension ? -1. : 1.;
		band[0][1] += sign * off2;
		band[1][0] += sign * off1;
		band[2][0] += sign * off2;
		//Natural BC at rho=rhoMax: (third and fourth derivatives vanish)
		const double extrap[2][3] = {{+1,-3,+3}, {+3,-8,+6}};
		band[-1][N-3] += off2 * extrap[0][0];
		band[ 0][N-2] += off2 * extrap[0][1];
		band[ 1][N-2] += off2 * extrap[0][2];
		band[-2][N-3] += off1 * extrap[0][0] + off2 * extrap[1][0];
		band[-1][N-2] += off1 * extrap[0][1] + off2 * extrap[1][1];
		band[ 0][N-1] += off1 * extrap[0][2] + off2 * extrap[1][2];
		//Pentadiagonal solve: Forward sweep:
		for(int i=0; i<N; i++)
			for(int m=1; m<=2; m++)
				if(i+m<N)
				{	register double rot = band[-m][i]/band[0][i]; //gauss elimination factor
					for(int n=1; n<m; n++)
						band[n-m][i+n] -= rot*band[n][i];
					for(int n=m; n<=2; n++)
						if(i+n<N)
							band[n-m][i+m] -= rot*band[n][i];
					x[i+m] -= rot*x[i];
				}
		//Pentadiagonal solve: Reverse sweep:
		for(int i=N-1; i>=0; i--)
		{	register double invDiag = 1.0/band[0][i];
			for(int m=1; m<=2; m++)
				if(i+m<N) x[i] -= band[m][i]*x[i+m];
			x[i] *= invDiag;
		}
		//Mirror two coefficients at beginning:
		std::vector<double> coeff;
		coeff.push_back(sign * x[2]);
		coeff.push_back(sign * x[1]);
		coeff.insert(coeff.end(), x.begin(), x.end());
		//Extrapolate two coefficients for natural BC at rhoMax:
		coeff.push_back(extrap[0][0]*x[N-3] + extrap[0][1]*x[N-2] + extrap[0][2]*x[N-1]);
		coeff.push_back(extrap[1][0]*x[N-3] + extrap[1][1]*x[N-2] + extrap[1][2]*x[N-1]);
		return coeff;
	}
}
