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

#ifndef JDFTX_ELECTRONIC_NONLOCALPCM_INTERNAL_H
#define JDFTX_ELECTRONIC_NONLOCALPCM_INTERNAL_H

#include <core/matrix3.h>

//! Symmetric matrix in 3D
template<typename scalar=double> class symmMatrix3
{
	scalar v[6];
public:
	//Accessors
	__hostanddev__ scalar& operator[](int k) { return v[k]; }
	__hostanddev__ const scalar& operator[](int k) const { return v[k]; }
	__hostanddev__ scalar& xx() { return v[0]; }
	__hostanddev__ scalar& yy() { return v[1]; }
	__hostanddev__ scalar& zz() { return v[2]; }
	__hostanddev__ scalar& yz() { return v[3]; }
	__hostanddev__ scalar& zx() { return v[4]; }
	__hostanddev__ scalar& xy() { return v[5]; }
	__hostanddev__ const scalar& xx() const { return v[0]; }
	__hostanddev__ const scalar& yy() const { return v[1]; }
	__hostanddev__ const scalar& zz() const { return v[2]; }
	__hostanddev__ const scalar& yz() const { return v[3]; }
	__hostanddev__ const scalar& zx() const { return v[4]; }
	__hostanddev__ const scalar& xy() const { return v[5]; }
	
	symmMatrix3(std::vector<scalar> a) { for(int k=0; k<6; k++) v[k]=a[k]; }
};

//! Load matrix from a const symmetric matrix field
template<typename scalar> __hostanddev__ matrix3<scalar> loadSymmMatrix(const symmMatrix3<const scalar*>& mArr, int i)
{	return matrix3<scalar>(
		mArr.xx()[i], mArr.xy()[i], mArr.zx()[i],
		mArr.xy()[i], mArr.yy()[i], mArr.yz()[i],
		mArr.zx()[i], mArr.yz()[i], mArr.zz()[i]);
}
//! Load matrix from a symmetric matrix field
template<typename scalar> __hostanddev__ matrix3<scalar> loadSymmMatrix(const symmMatrix3<scalar*>& mArr, int i)
{	return matrix3<scalar>(
		mArr.xx()[i], mArr.xy()[i], mArr.zx()[i],
		mArr.xy()[i], mArr.yy()[i], mArr.yz()[i],
		mArr.zx()[i], mArr.yz()[i], mArr.zz()[i]);
}
//! Store matrix to a symmetric matrix field
template<typename scalar> __hostanddev__ void storeSymmMatrix(const matrix3<scalar>& m, symmMatrix3<scalar*>& mArr, int i)
{	mArr.xx()[i] = m(0,0);
	mArr.yy()[i] = m(1,1);
	mArr.zz()[i] = m(2,2);
	mArr.yz()[i] = m(1,2);
	mArr.zx()[i] = m(2,0);
	mArr.xy()[i] = m(0,1);
}
//! Accumulate tensor onto a symmetric matrix field
template<typename scalar> __hostanddev__ void accumSymmMatrix(const matrix3<scalar>& m, symmMatrix3<scalar*>& mArr, int i)
{	mArr.xx()[i] += m(0,0);
	mArr.yy()[i] += m(1,1);
	mArr.zz()[i] += m(2,2);
	mArr.yz()[i] += m(1,2);
	mArr.zx()[i] += m(2,0);
	mArr.xy()[i] += m(0,1);
}

//! Compute surface tension reduction factor f(H1,H2):
__hostanddev__ double cavitation_f(double H1, double H2, double& f_H1, double& f_H2)
{	const double X = 11.2;
	const double A = 25, B = 7.6, C = 28, D = 39;
	const double P = 0.047, Q = 2800, R = 20, S = 25;
	//Symmetric part rational function:
	double sym, sym_H1, sym_H2;
	{	double num = 1 + A*H2 - B*H1*H1;
		if(num < 0.) { f_H1=f_H2=0.; return 0.; }
		double denInv = 1./(1 + H2*(C + H2*D)), den_H2 = C + H2*(2*D);
		sym = sqrt(num*denInv);
		sym_H1 = (0.5/sym)*(-2*B)*H1*denInv;
		sym_H2 = (0.5/sym)*denInv*(A - num*denInv*den_H2);
	}
	//Assymetric part rational function
	double asym, asym_H2;
	{	double num = 1. + Q*H2;
		double denInv = 1./(1 + H2*(R + H2*S)), den_H2 = R + H2*(2*S);
		asym = sqrt(num*denInv);
		asym_H2 = (0.5/asym)*denInv*(Q - num*denInv*den_H2);
	}
	//put together reduction factor:
	double switchFac = P*tanh(X*H1), switchFac_H1 = (P*X)/pow(cosh(H1*X),2);
	f_H1 = sym_H1 - asym * switchFac_H1;
	f_H2 = sym_H2 - asym_H2 * switchFac;
	return sym - asym * switchFac;
}

//!Return the cavitation energy density given the first and second cartesian
//!derivatives of shape, and accumulate the derivative w.r.t those
__hostanddev__ double cavitationEnergy_calc(int i,
	vector3<const double*> Dshape, symmMatrix3<const double*> DDshape,
	vector3<double*> grad_Dshape, symmMatrix3<double*> grad_DDshape)
{	
	const double sigma0 = 4.62e-5;
	vector3<> Ds = loadVector(Dshape, i);
	matrix3<> DDs = loadSymmMatrix(DDshape, i);
	double ds = Ds.length();
	if(ds < 1e-8) return 0.;
	double dsInv = 1./ds;
	double dsInvSq = dsInv*dsInv;
	
	//Compute relevant scalar contractions:
	double Ds_DDs_Ds = DDs.metric_length_squared(Ds), trDDs = trace(DDs);
	matrix3<> DDsSq = DDs * DDs;
	double Ds_DDsSq_Ds = DDsSq.metric_length_squared(Ds), trDDsSq = trace(DDsSq);
	
	//Compute curvatures:
	double H1 = 0.5*dsInv*(trDDs - dsInvSq*Ds_DDs_Ds);
	double H2 = 0.5*dsInvSq*(trDDsSq + dsInvSq*(-2*Ds_DDsSq_Ds + dsInvSq*Ds_DDs_Ds*Ds_DDs_Ds));
	
	//Compute cavitation energy density:
	double f_H1, f_H2, f = cavitation_f(H1, H2, f_H1, f_H2);
	double Ecav = sigma0 * ds * f;
	double Ecav_ds = sigma0 * f;
	double Ecav_H1 = sigma0 * ds * f_H1;
	double Ecav_H2 = sigma0 * ds * f_H2;
	
	//Propagate curvature derivatives to scalar contractions:
	Ecav_ds += Ecav_H1 * (-0.5*dsInvSq*(trDDs - 3*dsInvSq*Ds_DDs_Ds));
	Ecav_ds += Ecav_H2 * (-dsInv*dsInvSq*(trDDsSq + dsInvSq*(-4*Ds_DDsSq_Ds + 3*dsInvSq*Ds_DDs_Ds*Ds_DDs_Ds)));
	double Ecav_trDDs = Ecav_H1 * (0.5*dsInv);
	double Ecav_Ds_DDs_Ds = dsInvSq*( Ecav_H1*(-0.5*dsInv) + Ecav_H2*(dsInvSq*dsInvSq*Ds_DDs_Ds) );
	double Ecav_trDDsSq = Ecav_H2 * (0.5*dsInvSq);
	double Ecav_Ds_DDsSq_Ds = Ecav_H2 * (-dsInvSq*dsInvSq);

	//Propagate derivatives wrt scalar contractions to Ds and DDs:
	vector3<> Ecav_Ds
		= (Ecav_ds * dsInv) * Ds
		+ (Ecav_Ds_DDs_Ds * 2) * (DDs * Ds)
		+ (Ecav_Ds_DDsSq_Ds * 2) * (DDsSq * Ds);
	matrix3<> DsDsDag = outer(Ds,Ds);
	matrix3<> Ecav_DDs
		= Ecav_trDDs * matrix3<>(1,1,1)
		+ (Ecav_trDDsSq * 2) * DDs;
		+ Ecav_Ds_DDs_Ds * DsDsDag
		+ Ecav_Ds_DDsSq_Ds * (DDs*DsDsDag + DsDsDag*DDs);
	
	//Store results:
	accumVector(Ecav_Ds, grad_Dshape, i);
	accumSymmMatrix(Ecav_DDs, grad_DDshape, i);
	return Ecav;
}

#endif // JDFTX_ELECTRONIC_NONLOCALPCM_INTERNAL_H
