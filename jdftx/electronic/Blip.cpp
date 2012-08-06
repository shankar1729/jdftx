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

#include <electronic/common.h>
#include <electronic/Blip.h>
#include <core/GridInfo.h>
#include <core/vector3.h>
#include <core/Thread.h>
#include <core/Operators.h>

//-----------------------------------------------
//---------- PW to Blip conversion --------------
//-----------------------------------------------

double* BlipConverter::newGamma(int N)
{
	double* gamma = new double[N];
	gamma[0] = 2.0/3;
	for(int i=1; i<N; i++)
	{
		double k = (2*M_PI/N) * (i<N/2 ? i : N-i);
		gamma[i] = (1.0/6) * pow(k*k/(1-cos(k)),2);
	}
	return gamma;
}


BlipConverter::BlipConverter(int Nx, int Ny, int Nz) : Nx(Nx), Ny(Ny), Nz(Nz)
{
	xGamma = newGamma(Nx);
	yGamma = newGamma(Ny);
	zGamma = newGamma(Nz);
}


BlipConverter::~BlipConverter()
{
	delete[] xGamma;
	delete[] yGamma;
	delete[] zGamma;
}


//Given a complex PW basis object, return corresponding real-space Blip coefficient set
complexDataRptr BlipConverter::operator()(const complexDataGptr& vTilde)
{	complex* vTildeData = vTilde->data();
	int index=0;
	for(int ix=0; ix<Nx; ix++)
		for(int iy=0; iy<Ny; iy++)
			for(int iz=0; iz<Nz; iz++)
				vTildeData[index++] *= (xGamma[ix]*yGamma[iy]*zGamma[iz]);
	return I(vTilde);
}
complexDataRptr BlipConverter::operator()(const complexDataRptr& v)
{	return (*this)(J(v));
}

//Given a real PW basis object v, return corresponding real-space Blip coefficient set
DataRptr BlipConverter::operator()(const DataGptr& vTilde)
{	complex* vTildeData = vTilde->data();
	int index=0;
	for(int ix=0; ix<Nx; ix++)
		for(int iy=0; iy<Ny; iy++)
			for(int iz=0; iz<=Nz/2; iz++)
				vTildeData[index++] *= (xGamma[ix]*yGamma[iy]*zGamma[iz]);
	return I(vTilde);
}
DataRptr BlipConverter::operator()(const DataRptr& v)
{	return (*this)(J(v));
}


//---------------------------------------------------------------------------
//---------- Analytical integration via orthogonal polynomials --------------
//---------------------------------------------------------------------------

//coefficients of the tricubic polynomial in each cell expressed in Legendre polynomials
template<typename scalar> struct TriCubic
{	scalar c[4][4][4];

	//convert 1d blip coeffs to legendre coefficients (in-place)
	static inline void convert1d(scalar& coeff0, scalar& coeff1, scalar& coeff2, scalar& coeff3)
	{
		register scalar a=coeff0, b=coeff1, c=coeff2, d=coeff3;
		coeff0 = (1.0/16)*(a+d) + (11.0/16)*(b+c);
		coeff1 = (9.0/80)*(d-a) + (33.0/80)*(c-b);
		coeff2 = (1.0/16)*(a+d-b-c);
		coeff3 = (1.0/80)*(d-a) + (3.0/80)*(b-c);
	}

	TriCubic(const scalar* v, const vector3<int>& S, int i0, int i1, int i2)
	{	//Load coefficients:
		int j0[4], j1[4], j2[4];
		#define WrapIndex(d) j##d[k]=i##d+k; if(j##d[k]>=S[d]) j##d[k] -= S[d];
		for(int k=0; k<4; k++)
		{	WrapIndex(0);
			WrapIndex(1);
			WrapIndex(2);
		}
		#undef WrapIndex
		for(int k0=0; k0<4; k0++)
		{
			register int offs0 = j0[k0]*S[1];
			for(int k1=0; k1<4; k1++)
			{
				register int offs01 = (offs0+j1[k1])*S[2];
				for(int k2=0; k2<4; k2++)
					c[k0][k1][k2] = v[offs01+j2[k2]];
			}
		}
		//Recombine into Legendre coefficients
		for(int i0=0; i0<4; i0++) for(int i1=0; i1<4; i1++) convert1d(c[i0][i1][0], c[i0][i1][1], c[i0][i1][2], c[i0][i1][3]);
		for(int i0=0; i0<4; i0++) for(int i2=0; i2<4; i2++) convert1d(c[i0][0][i2], c[i0][1][i2], c[i0][2][i2], c[i0][3][i2]);
		for(int i1=0; i1<4; i1++) for(int i2=0; i2<4; i2++) convert1d(c[0][i1][i2], c[1][i1][i2], c[2][i1][i2], c[3][i1][i2]);
	}
};

inline double Tcell(const TriCubic<complex>& phi, const matrix3<>* Tmat)
{
	register double T00=0, T11=0, T22=0; //diagonal terms in oblique laplacian
	register complex T01=0, T12=0, T20=0; //off-diagonal terms in oblique laplacian

	static const double w[4] = {1.0, 1.0/3, 1.0/5, 1.0/7}; //weights when no derivatives

	//00 term:
	for(int i1=0; i1<4; i1++)
		for(int i2=0; i2<4; i2++)
			T00 += w[i1]*w[i2]*
				( 4.0*norm(phi.c[1][i1][i2]+phi.c[3][i1][i2])
				+ 12.0*norm(phi.c[2][i1][i2])
				+ 20.0*norm(phi.c[3][i1][i2]) );
	//11 term:
	for(int i0=0; i0<4; i0++)
		for(int i2=0; i2<4; i2++)
			T11 += w[i0]*w[i2]*
				( 4.0*norm(phi.c[i0][1][i2]+phi.c[i0][3][i2])
				+ 12.0*norm(phi.c[i0][2][i2])
				+ 20.0*norm(phi.c[i0][3][i2]) );
	//22 term:
	for(int i0=0; i0<4; i0++)
		for(int i1=0; i1<4; i1++)
			T11 += w[i0]*w[i1]*
				( 4.0*norm(phi.c[i0][i1][1]+phi.c[i0][i1][3])
				+ 12.0*norm(phi.c[i0][i1][2])
				+ 20.0*norm(phi.c[i0][i1][3]) );

	static const int ja[4] = {1,2,3,3};
	static const int jb[4] = {0,1,2,0};
	for(int i0=0; i0<4; i0++)
		for(int i1=0; i1<4; i1++)
			for(int i2=0; i2<4; i2++)
			{
				T12 += 4.0*w[i0]*conj(phi.c[i0][ja[i1]][jb[i2]])*phi.c[i0][jb[i1]][ja[i2]];
				T20 += 4.0*w[i1]*conj(phi.c[ja[i0]][i1][jb[i2]])*phi.c[jb[i0]][i1][ja[i2]];
				T01 += 4.0*w[i2]*conj(phi.c[ja[i0]][jb[i1]][i2])*phi.c[jb[i0]][ja[i1]][i2];
			}

	return (*Tmat)(0,0)*T00 + (*Tmat)(1,1)*T11 + (*Tmat)(2,2)*T22
		+ 2*( (*Tmat)(1,2)*T12.real() + (*Tmat)(2,0)*T20.real() + (*Tmat)(0,1)*T01.real() );
}


struct SparseTensorElem
{
	int a,b,c;
	double w;

	SparseTensorElem(int a1, int b1, int c1, double w)
	:a(a1-1),b(b1-1),c(c1-1),w(w) //shifting from 1-based indices (from Mathematica)
	{}
};

inline double Vcell(const TriCubic<complex>& phi, const TriCubic<double>& V)
{
	static const int nSparse=15;
	static const SparseTensorElem s[nSparse] =
	{
		{1,1,1, 1.0},
		{1,2,2, 1.0/3},
		{1,3,3, 1.0/5},
		{1,4,4, 1.0/7},
		{2,2,1, 2.0/3},
		{2,3,2, 4.0/15},
		{2,4,3, 6.0/35},
		{3,2,2, 2.0/15},
		{3,3,1, 2.0/5},
		{3,3,3, 2.0/35},
		{3,4,2, 6.0/35},
		{3,4,4, 4.0/105},
		{4,3,2, 6.0/35},
		{4,4,1, 2.0/7},
		{4,4,3, 8.0/105}
	};
	complex result=0.0;
	for(int i0=0; i0<nSparse; i0++)
		for(int i1=0; i1<nSparse; i1++)
			for(int i2=0; i2<nSparse; i2++)
				result += s[i0].w * s[i1].w * s[i2].w
					* V.c[s[i0].a][s[i1].a][s[i2].a]
					* phi.c[s[i0].b][s[i1].b][s[i2].b]
					* conj(phi.c[s[i0].c][s[i1].c][s[i2].c]);
	return result.real();
}


//------------------------------------------------------------------------
//---------- Outer threaded loop over cells (both methods)  --------------
//------------------------------------------------------------------------

//Local potential energy for one slice (to be threaded over slices)
double Vblip_sub(int i0, const vector3<int> S, const complex* phi, const double* V)
{
	double res=0.0;
	for(int i1=0; i1<S[1]; i1++)
		for(int i2=0; i2<S[2]; i2++)
			res += Vcell(TriCubic<complex>(phi, S, i0, i1, i2), TriCubic<double>(V, S, i0, i1, i2));
	return res;
}


//Compute the local potential energy for blip orbital phi in blip potential V
double Vblip(const complexDataRptr& phi, const DataRptr& V)
{	assert(phi->gInfo.S == V->gInfo.S);
	return phi->gInfo.dV * threadedAccumulate(Vblip_sub, phi->gInfo.S[0], phi->gInfo.S, phi->data(), V->data());
}


//Kinetic energy for one slice (to be threaded over slices)
double Tblip_sub(int i0, const vector3<int> S, const complex* phi, const matrix3<>* Tmat,
	double* tMaxPtr, int* i0maxPtr, int* i1maxPtr, int* i2maxPtr, std::mutex* m)
{
	double res=0.0;
	double tMax=0.0; int i0max=0, i1max=0, i2max=0;
	for(int i1=0; i1<S[1]; i1++)
		for(int i2=0; i2<S[2]; i2++)
		{	double t = Tcell(TriCubic<complex>(phi, S, i0, i1, i2), Tmat);
			res += t;
			if(t>tMax) { tMax=t; i0max=i0; i1max=i1; i2max=i2; }
		}
	m->lock();
	if(tMaxPtr && tMax>*tMaxPtr)
	{	*tMaxPtr=tMax;
		if(i0maxPtr) *i0maxPtr=i0max;
		if(i1maxPtr) *i1maxPtr=i1max;
		if(i2maxPtr) *i2maxPtr=i2max;
	}
	m->unlock();
	return res;
}


//Compute the kinetic energy for a blip orbital phi
double Tblip(const complexDataRptr& phi, double* tMax, int* i0max, int* i1max, int* i2max)
{
	matrix3<> diagSinv(1.0/phi->gInfo.S[0], 1.0/phi->gInfo.S[1], 1.0/phi->gInfo.S[2]);
	matrix3<> H = phi->gInfo.R * diagSinv;
	matrix3<> invHT = inv(~H);
	matrix3<> Tmat = 0.5*(~invHT)*invHT;
	std::mutex m; if(tMax) *tMax=0.0;
	return phi->gInfo.dV * threadedAccumulate(Tblip_sub, phi->gInfo.S[0], phi->gInfo.S,
		phi->data(), &Tmat, tMax, i0max, i1max, i2max, &m);
}

