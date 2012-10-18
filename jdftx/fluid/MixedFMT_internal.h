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

#ifndef JDFTX_FLUID_MIXEDFMT_INTERNAL_H
#define JDFTX_FLUID_MIXEDFMT_INTERNAL_H

#include <vector>
#include <core/matrix3.h>
#include <core/tensor3.h>

__hostanddev__ void tensorKernel_calc(int i, const vector3<int> iG, bool nyq, const matrix3<> G,
	const complex* nTilde, tensor3<complex*> mTilde)
{	complex minus_nTilde = nyq ? complex(0,0) : -nTilde[i];
	vector3<> Gvec = iG*G;
	double Gsq = Gvec.length_squared();
	mTilde.xy()[i] = minus_nTilde*Gvec.x()*Gvec.y();
	mTilde.yz()[i] = minus_nTilde*Gvec.y()*Gvec.z();
	mTilde.zx()[i] = minus_nTilde*Gvec.z()*Gvec.x();
	mTilde.xxr()[i] = minus_nTilde*(Gvec.x()*Gvec.x() - (1.0/3)*Gsq);
	mTilde.yyr()[i] = minus_nTilde*(Gvec.y()*Gvec.y() - (1.0/3)*Gsq);
}

__hostanddev__ void tensorKernel_grad_calc(int i, const vector3<int> iG, bool nyq, const matrix3<> G,
	tensor3<const complex*> grad_mTilde, complex* grad_nTilde)
{	complex temp = complex(0,0);
	if(!nyq)
	{	vector3<> Gvec = iG*G;
		double Gsq = Gvec.length_squared();
		temp += grad_mTilde.xy()[i]*Gvec.x()*Gvec.y();
		temp += grad_mTilde.yz()[i]*Gvec.y()*Gvec.z();
		temp += grad_mTilde.zx()[i]*Gvec.z()*Gvec.x();
		temp += grad_mTilde.xxr()[i]*(Gvec.x()*Gvec.x() - (1.0/3)*Gsq);
		temp += grad_mTilde.yyr()[i]*(Gvec.y()*Gvec.y() - (1.0/3)*Gsq);
	}
	grad_nTilde[i] = -temp;
}

//Compute vT*m*v for a vector v and a symmetric traceless tensor m
__hostanddev__ double mul_vTmv(const tensor3<>& m, const vector3<>& v)
{	return 2*(m.xy()*v.x()*v.y() + m.yz()*v.y()*v.z() + m.zx()*v.z()*v.x())
		+ m.xxr()*v.x()*v.x() + m.yyr()*v.y()*v.y() - (m.xxr()+m.yyr())*v.z()*v.z();
}
//Accumulate gradient of above function
__hostanddev__ void mul_vTmv_grad(const double grad_mul, const tensor3<>& m, const vector3<>& v,
	tensor3<>& grad_m, vector3<>& grad_v)
{	grad_m.xy() += (2*grad_mul)*v.x()*v.y();
	grad_m.yz() += (2*grad_mul)*v.y()*v.z();
	grad_m.zx() += (2*grad_mul)*v.z()*v.x();
	grad_m.xxr() += grad_mul*(v.x()*v.x()-v.z()*v.z());
	grad_m.yyr() += grad_mul*(v.y()*v.y()-v.z()*v.z());
	grad_v.x() += (2*grad_mul)*(m.xy()*v.y() + m.zx()*v.z() + m.xxr()*v.x());
	grad_v.y() += (2*grad_mul)*(m.xy()*v.x() + m.yz()*v.z() + m.yyr()*v.y());
	grad_v.z() += (2*grad_mul)*(m.yz()*v.y() + m.zx()*v.x()- (m.xxr()+m.yyr())*v.z());
}


//Compute tr(m^3) for a symmetric traceless tensor m (See ~/Water1D/FMT_tensorWeights.m for expressions)
__hostanddev__ double trace_cubed(const tensor3<>& m)
{	return -3.0 * (-2.0*m.xy()*m.yz()*m.zx() + pow(m.yz(),2)*m.xxr() - pow(m.xy(),2)*(m.xxr()+m.yyr()) + m.yyr()*(pow(m.zx(),2)+m.xxr()*(m.xxr()+m.yyr())));
}
//Accumulate gradient of above function
__hostanddev__ void trace_cubed_grad(const double grad_trace, const tensor3<>& m, tensor3<>& grad_m)
{	grad_m.xy() += (6*grad_trace)* (m.yz()*m.zx()+m.xy()*(m.xxr()+m.yyr()));
	grad_m.yz() += (6*grad_trace)* (m.xy()*m.zx()-m.yz()*m.xxr());
	grad_m.zx() += (6*grad_trace)* (m.xy()*m.yz()-m.zx()*m.yyr());
	grad_m.xxr() += (-3*grad_trace)* (-pow(m.xy(),2)+pow(m.yz(),2)+m.yyr()*(2*m.xxr()+m.yyr()));
	grad_m.yyr() += (-3*grad_trace)* (-pow(m.xy(),2)+pow(m.zx(),2)+m.xxr()*(m.xxr()+2*m.yyr()));
}

//White-Bear mark II FMT scale functions (and derivatives) [replace with f2(x)=f3(x)=1 for standard Tarazona FMT]:
__hostanddev__ double WB_f2(double x, double& f2_x)
{	if(x<0.002)
	{	f2_x = x*((2./9) + x*(3./18 + x*(4./30)));
		return 1 + x*x*((1./9) + x*(1./18 + x*(1./30)));
	}
	else
	{	f2_x = (-1./3)*(x*(2+x) + 2*log(1-x)) / (x*x);
		return 1 + (1./3)*(2-x + 2*(1-x)*log(1-x)/x);
	}
}
__hostanddev__ double WB_f3(double x, double& f3_x)
{	if(x<0.005)
	{	f3_x = -4./9 + x*(2./18 + x*(3./45 + x*(4./90)));
		return 1 + x*(-4./9 + x*(1./18 + x*(1./45 + x*(1./90))));
	}
	else
	{	f3_x = 2*(1-x) * (x*(2+x) + 2*log(1-x)) / (3*x*x*x);
		return 1 - (x*(2+x*(-3+x*2)) + 2*(1-x)*(1-x)*log(1-x)) / (3*x*x);
	}
}

__hostanddev__ double phiFMT_calc(int i,
	const double *n0arr, const double *n1arr, const double *n2arr, const double *n3arr,
	vector3<const double*> n1vArr, vector3<const double*> n2vArr, tensor3<const double*> n2mArr,
	double *grad_n0arr, double *grad_n1arr, double *grad_n2arr, double *grad_n3arr,
	vector3<double*> grad_n1vArr, vector3<double*> grad_n2vArr, tensor3<double*> grad_n2mArr)
{
	double n0 = n0arr[i];
	double n1 = n1arr[i];
	double n2 = n2arr[i];
	double n3 = n3arr[i];
	if(n0<0. || n1<0. || n2<0. || n3<0.) return 0.;
	if(n3>=1.) return NAN;
	vector3<> n1v = loadVector(n1vArr, i);
	vector3<> n2v = loadVector(n2vArr, i);
	tensor3<> n2m = loadTensor(n2mArr, i);
	double tensorPart = mul_vTmv(n2m, n2v) - 0.5*trace_cubed(n2m);

	double n1v_n2v = dot(n1v, n2v);
	double n2vsq = n2v.length_squared();
	double pole = 1./(1-n3); //The following is derived easily using: d(pole^N)/d(n3) = N pole^(N+1)

	double f2prime, f2 = WB_f2(n3, f2prime), comb2 = (n1*n2-n1v_n2v);
	double f3prime, f3 = WB_f3(n3, f3prime), comb3 = (n2*(n2*n2-3*n2vsq) + 9.*tensorPart);
	double phi    = n0*log(pole) + pole*(f2*comb2 + pole*((1./((24*M_PI)))*f3*comb3));
	double phi_n0 = log(pole);
	double phi_n1 = pole*(f2*n2);
	double phi_n2 = pole*(f2*n1 + pole*((1./(8*M_PI))*f3*(n2*n2-n2vsq)));
	double phi_n3 = pole*(n0 + pole*(f2*comb2 + pole*((1./(12*M_PI))*f3*comb3)))
		+ pole*(f2prime*comb2 + pole*((1./(24*M_PI))*f3prime*comb3));
	vector3<> phi_n1v = (-pole*f2)*n2v;
	vector3<> phi_n2v = (-pole)*(f2*n1v + (pole*f3*n2/(4*M_PI))*n2v );
	double phi_tensorPart = pole*pole*f3*(9./(24*M_PI));

	tensor3<> phi_n2m;
	mul_vTmv_grad(phi_tensorPart, n2m, n2v, phi_n2m, phi_n2v);
	trace_cubed_grad(-0.5*phi_tensorPart, n2m, phi_n2m);
	accumTensor(phi_n2m, grad_n2mArr, i);
	accumVector(phi_n2v, grad_n2vArr, i);
	accumVector(phi_n1v, grad_n1vArr, i);
	grad_n3arr[i] += phi_n3;
	grad_n2arr[i] += phi_n2;
	grad_n1arr[i] += phi_n1;
	grad_n0arr[i] += phi_n0;
	return phi;
}

__hostanddev__ double phiBond_calc(int i, double Rhm, double scale,
	const double *n0arr, const double *n2arr, const double *n3arr, vector3<const double*> n2vArr,
	double *grad_n0arr, double *grad_n2arr, double *grad_n3arr, vector3<double*> grad_n2vArr)
{
	double n0 = n0arr[i];
	double n2 = n2arr[i];
	double n3 = n3arr[i];
	if(n0<0. || n2<0. || n3<0.) return 0.;
	double pole = 1.0/(1-n3); //The following is derived easily using: d(pole^N)/d(n3) = N pole^(N+1)
	//Vector correction factor and derivatives:
	vector3<> n2v = loadVector(n2vArr, i);
	double n2vsq = n2v.length_squared();
	double zeta = (n2<=0 || n2vsq>n2*n2) ? 0.0 : (1-n2vsq/(n2*n2));
	double zeta_n2 = zeta ? 2*n2vsq/(n2*n2*n2) : 0.0;
	vector3<> zeta_n2v = zeta ? (-2/(n2*n2))*n2v : vector3<>();
	//Compute the contact correlation function and its derivatives:
	double gContact      = pole*(1 + zeta * pole*(Rhm*n2 + pole*(2.0/9)*pow(Rhm*n2,2) ));
	double gContact_n2   = zeta * pow(pole,2)*(Rhm + pole*(4.0/9)*pow(Rhm,2)*n2 );
	double gContact_n3   = pow(pole,2)*(1 + zeta * pole*(2*Rhm*n2 + pole*(6.0/9)*pow(Rhm*n2,2) ));
	double gContact_zeta = pow(pole,2) * (Rhm*n2 + pole*(2.0/9)*pow(Rhm*n2,2) );
	//Compute the bonding corrections and its derivatives:
	double phi     = -scale*n0*log(gContact);
	double phi_n0   = -scale*log(gContact);
	double phi_n2   = -scale*n0*gContact_n2/gContact;
	double phi_n3   = -scale*n0*gContact_n3/gContact;
	double phi_zeta = -scale*n0*gContact_zeta/gContact;
	//Accumulate the gradients and return the answer:
	grad_n0arr[i] += phi_n0;
	grad_n2arr[i] += phi_n2 + phi_zeta * zeta_n2;
	grad_n3arr[i] += phi_n3;
	accumVector(phi_zeta * zeta_n2v, grad_n2vArr, i);
	return phi;
}

#endif // JDFTX_FLUID_MIXEDFMT_INTERNAL_H
