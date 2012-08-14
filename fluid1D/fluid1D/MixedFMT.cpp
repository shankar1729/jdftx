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

#include <fluid1D/MixedFMT.h>
#include <core1D/Operators.h>
#include <cmath>

//White-Bear mark II FMT scale functions (and derivatives) [replace with f2(x)=f3(x)=1 for standard Tarazona FMT]:
double WB_f2(double x, double& f2_x)
{	if(x<0.002)
	{	f2_x = x*((2./9) + x*(3./18 + x*(4./30)));
		return 1 + x*x*((1./9) + x*(1./18 + x*(1./30)));
	}
	else
	{	f2_x = (-1./3)*(x*(2+x) + 2*log(1-x)) / (x*x);
		return 1 + (1./3)*(2-x + 2*(1-x)*log(1-x)/x);
	}
}
double WB_f3(double x, double& f3_x)
{	if(x<0.005)
	{	f3_x = -4./9 + x*(2./18 + x*(3./45 + x*(4./90)));
		return 1 + x*(-4./9 + x*(1./18 + x*(1./45 + x*(1./90))));
	}
	else
	{	f3_x = 2*(1-x) * (x*(2+x) + 2*log(1-x)) / (3*x*x*x);
		return 1 - (x*(2+x*(-3+x*2)) + 2*(1-x)*(1-x)*log(1-x)) / (3*x*x);
	}
}

double phiFMT_calc(int i, const double* wArr,
	const double *n0arr, const double *n1arr, const double *n2arr, const double *n3arr,
	const double* n1vArr, const double* n2vArr, const double* n2mArr,
	double *grad_n0arr, double *grad_n1arr, double *grad_n2arr, double *grad_n3arr,
	double* grad_n1vArr, double* grad_n2vArr, double* grad_n2mArr)
{
	double n0 = n0arr[i];
	double n1 = n1arr[i];
	double n2 = n2arr[i];
	double n3 = n3arr[i];
	if(n0<0. || n1<0. || n2<0. || n3<0.) return 0.;
	double n1v = n1vArr[i];
	double n2v = n2vArr[i];
	double n2m = n2mArr[i];
	double pole = 1./(1-n3); //The following is derived easily using: d(pole^N)/d(n3) = N pole^(N+1)

	double f2prime, f2 = WB_f2(n3, f2prime), comb2 = (n1*n2-n1v*n2v);
	double f3prime, f3 = WB_f3(n3, f3prime), comb3 = (n2*(n2*n2-3*n2v*n2v) + n2m*(6*n2v*n2v-n2m*n2m));
	double w = wArr[i];
	double wphi    = w * (n0*log(pole) + pole*(f2*comb2 + pole*((1./((24*M_PI)))*f3*comb3)));
	grad_n0arr[i] += w * log(pole);
	grad_n1arr[i] += w * pole*(f2*n2);
	grad_n2arr[i] += w * pole*(f2*n1 + pole*((1./(8*M_PI))*f3*(n2*n2-n2v*n2v)));
	grad_n3arr[i] += w * pole*(
		n0 + pole*(f2*comb2 + pole*((1./(12*M_PI))*f3*comb3))
		+ f2prime*comb2 + pole*((1./(24*M_PI))*f3prime*comb3));
	grad_n1vArr[i] += w * (-pole*f2)*n2v;
	grad_n2vArr[i] += w * (-pole)*(f2*n1v + (pole*f3*(n2-2*n2m)/(4*M_PI))*n2v );
	grad_n2mArr[i] += w * pole*pole*f3*(2*n2v*n2v-n2m*n2m)*(1./(8*M_PI));
	return wphi;
}

double PhiFMT(const ScalarField& n0, const ScalarField& n1, const ScalarField& n2,
	const ScalarFieldTilde& n3tilde, const ScalarFieldTilde& n1vTilde, const ScalarFieldTilde& n2mTilde,
	ScalarField& grad_n0, ScalarField& grad_n1, ScalarField& grad_n2,
	ScalarFieldTilde& grad_n3tilde, ScalarFieldTilde& grad_n1vTilde, ScalarFieldTilde& grad_n2mTilde)
{
	const GridInfo& gInfo = *(n0.gInfo);
	ScalarField n3 = I(n3tilde);
	ScalarField n1v = ID(n1vTilde);
	ScalarField n2v = ID(-n3tilde);
	ScalarField n2m = IDD(n2mTilde);

	ScalarField grad_n3, grad_n1v, grad_n2v, grad_n2m;
	nullToZero(grad_n0, gInfo); nullToZero(grad_n1, gInfo); nullToZero(grad_n2, gInfo); nullToZero(grad_n3, gInfo);
	nullToZero(grad_n1v, gInfo); nullToZero(grad_n2v, gInfo); nullToZero(grad_n2m, gInfo);

	double Phi = serialAccumulate(phiFMT_calc, gInfo.S, gInfo.w.data(),
			n0.data(), n1.data(), n2.data(), n3.data(), n1v.data(), n2v.data(), n2m.data(),
			grad_n0.data(), grad_n1.data(), grad_n2.data(), grad_n3.data(),
			grad_n1v.data(), grad_n2v.data(), grad_n2m.data());
	n3=0; n1v=0; n2v=0; n2m=0; //no longer need these weighted densities (clean up)

	grad_n2mTilde += IDDdag(grad_n2m); grad_n2m=0;
	grad_n1vTilde += IDdag(grad_n1v); grad_n1v=0;
	grad_n3tilde += ( Idag(grad_n3) - IDdag(grad_n2v) ); grad_n3=0; grad_n2v=0;
	return Phi;
}

double phiFMTuniform(double n0, double n1, double n2, double n3,
	double& grad_n0, double& grad_n1, double& grad_n2, double& grad_n3)
{
	double one = 1., zero=0., dump=0.;
	return phiFMT_calc(0, &one, &n0, &n1, &n2, &n3, &zero, &zero, &zero,
		&grad_n0, &grad_n1, &grad_n2, &grad_n3, &dump, &dump, &dump);
}

//-------------- Bonding corrections ------------

double phiBond_calc(int i, double Rhm, double scale, const double* wArr,
	const double *n0arr, const double *n2arr, const double *n3arr, const double* n2vArr,
	double *grad_n0arr, double *grad_n2arr, double *grad_n3arr, double* grad_n2vArr)
{
	double n0 = n0arr[i];
	double n2 = n2arr[i];
	double n3 = n3arr[i];
	if(n0<0. || n2<0. || n3<0.) return 0.;
	double n2v = n2vArr[i];
	double pole = 1.0/(1-n3); //The following is derived easily using: d(pole^N)/d(n3) = N pole^(N+1)
	//Vector correction factor and derivatives:
	double n2vsq = n2v*n2v;
	double zeta = (n2<=0 || n2vsq>n2*n2) ? 0. : (1-n2vsq/(n2*n2));
	double zeta_n2 = zeta ? 2*n2vsq/(n2*n2*n2) : 0.;
	double zeta_n2v = zeta ? (-2/(n2*n2))*n2v : 0.;
	//Compute the contact correlation function and its derivatives:
	double gContact      = pole*(1 + zeta * pole*(Rhm*n2 + pole*(2.0/9)*pow(Rhm*n2,2) ));
	double gContact_zeta = pow(pole,2) * (Rhm*n2 + pole*(2.0/9)*pow(Rhm*n2,2) );
	double gContact_n2   = zeta * pow(pole,2)*(Rhm + pole*(4.0/9)*pow(Rhm,2)*n2 ) + gContact_zeta * zeta_n2;
	double gContact_n3   = pow(pole,2)*(1 + zeta * pole*(2*Rhm*n2 + pole*(6.0/9)*pow(Rhm*n2,2) ));
	double gContact_n2v  = gContact_zeta * zeta_n2v;
	//Compute the bonding corrections and its derivatives:
	double wscale = scale * wArr[i];
	double wphi     = -wscale*n0*log(gContact);
	grad_n0arr[i]  += -wscale*log(gContact);
	grad_n2arr[i]  += -wscale*n0*gContact_n2/gContact;
	grad_n3arr[i]  += -wscale*n0*gContact_n3/gContact;
	grad_n2vArr[i] += -wscale*n0*gContact_n2v/gContact;
	return wphi;
}

double PhiBond(double Rhm, double scale, const ScalarField& n0mol, const ScalarField& n2, const ScalarFieldTilde& n3tilde,
	ScalarField& grad_n0mol, ScalarField& grad_n2, ScalarFieldTilde& grad_n3tilde)
{
	const GridInfo& gInfo = *(n0mol.gInfo);
	//Compute n3 and n2v in real space from n3tilde:
	ScalarField n3 = I(n3tilde);
	ScalarField n2v = ID(-n3tilde);
	//Bonding correction and gradient:
	ScalarField grad_n3, grad_n2v;
	nullToZero(grad_n0mol, gInfo);
	nullToZero(grad_n2, gInfo);
	nullToZero(grad_n3, gInfo);
	nullToZero(grad_n2v, gInfo);
	
	double Phi = serialAccumulate(phiBond_calc, gInfo.S, Rhm, scale, gInfo.w.data(),
			n0mol.data(), n2.data(), n3.data(), n2v.data(),
			grad_n0mol.data(), grad_n2.data(), grad_n3.data(), grad_n2v.data());
	n3=0; n2v=0; //no longer need these weighted densities (clean up)
	//Propagate grad_n2v and grad_n3 to grad_n3tilde:
	grad_n3tilde += ( Idag(grad_n3) - IDdag(grad_n2v) );
	return Phi;
}

double phiBondUniform(double Rhm, double scale, double n0mol, double n2, double n3,
	double& grad_n0mol, double& grad_n2, double& grad_n3)
{	double one = 1., zero = 0., dump = 0.;
	return phiBond_calc(0, Rhm, scale, &one,
		&n0mol, &n2, &n3, &zero, &grad_n0mol, &grad_n2, &grad_n3, &dump);
}
