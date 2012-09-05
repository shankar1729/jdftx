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

#ifndef JDFTX_CORE_SIMPLEX_H
#define JDFTX_CORE_SIMPLEX_H

//! @file Simplex.h Compute the fourier transform of simplex-shaped theta functions

#include <array>
#include <algorithm>
#include <core/scalar.h>

//! n-dimensional simplex with one vertex at origin (+ its image under inversion)
template<int n> struct Simplex
{
	typedef std::array<double,n> Point;
	std::array<Point,n> v; //!< vertices of simplex (excpet for the implicit one at the origin)
	double V; //!< parallelopiped volume prefactor
	
	void init(); //!< Call after setting vertices to initialize volume prefactor
	static double dot(const Point& p1, const Point& p2); //!< dot product
	
	//! Get fourier transform of simplex (+ its image under inversion) at a specified G-vector
	double getTilde(const Point& G) const;

};

//---------------------- Implementation ----------------------------
//!@cond

//Hard-coded determinants for relevant sizes
inline double det(const std::array<std::array<double,2>,2>& M)
{	return M[0][0]*M[1][1] - M[0][1]*M[1][0];
}
inline double det(const std::array<std::array<double,3>,3>& M)
{	return M[0][0] * (M[1][1]*M[2][2] - M[1][2]*M[2][1])
		+  M[0][1] * (M[1][2]*M[2][0] - M[1][0]*M[2][2])
		+  M[0][2] * (M[1][0]*M[2][1] - M[1][1]*M[2][0]);
}
//Set parallelopiped-volume prefactor for simplex:
template<int n> void Simplex<n>::init()
{	V = fabs(det(v));
}


//Dot product
template<int n> double Simplex<n>::dot(const Simplex<n>::Point& p1, const Simplex<n>::Point& p2)
{	double ret = 0.;
	for(int i=0; i<n; i++)
		ret += p1[i] * p2[i];
	return ret;
}


//Homogeneous polynomial of degree d in {x} (Note: template recursion - no runtime recursion overhead)
template<int nx> double homPoly(int d, const double* x)
{	double ret = 0., xi = 1.;
	for(int i=0; i<d; i++)
	{	ret += xi * homPoly<nx-1>(d-i, x+1);
		xi *= x[0];
	}
	ret += xi; //homPoly(0, ...) = 1
	return ret;
}
template<> inline double homPoly<0>(int d, const double* x) { return d ? 0. : 1.; } //template-recursion base case
template<> inline double homPoly<1>(int d, const double* x) { return pow(x[0],d); } //template-recursion base case

//Get contribution from a single cluster in a simplex (called by Simplex::getTilde via getSimplexTilde_mSwitch)
#define accumClusterArgs const std::array<double,n+1>& d, \
		const typename std::array<double,n+1>::const_iterator& iStart, \
		const typename std::array<double,n+1>::const_iterator& iEnd, \
		double dMean, double& result
template<int n, int m> void getSimplexTilde_accumCluster(accumClusterArgs)
{	double den = 1.;
	std::array<double,n+1-m> x; //array of individual denominators
	auto xIter = x.begin();
	for(auto j=d.begin(); j!=iStart; j++) den *= -( *(xIter++) = 1./(*j - dMean) );
	for(auto j=iEnd; j!=d.end(); j++) den *= -( *(xIter++) = 1./(*j - dMean) );
	//Add contributions from homogeneous polynomials in x of various order:
	complex phase = cis(n*(0.5*M_PI) - dMean);
	for(int k=0; k<m; k++)
	{	if(k) den /= k; //Extra factor of k! relative to original denominator
		result += (den * homPoly<n+1-m>(m-1-k, x.begin())) * phase.real();
		phase *= complex(0,-1); //so that phase.real() = cos((n-k)*(0.5*M_PI) - dMean)
	}
}
//Called by Simplex::getTilde for each cluster. This struct switches m to a
//compile time constant so as to call getSimplexTilde_accumCluster<n,m>(...)
template<int n, int M> struct getSimplexTilde_mSwitch
{	void operator()(int m, accumClusterArgs)
	{	if(m==M) getSimplexTilde_accumCluster<n,M>(d, iStart, iEnd, dMean, result);
		else getSimplexTilde_mSwitch<n,M-1>()(m, d, iStart, iEnd, dMean, result);
	}
};
template<int n> struct getSimplexTilde_mSwitch<n,0>
{	void operator()(int m, accumClusterArgs) {} //template recursion end - will never be called
};
#undef accumClusterArgs

//Get fourier transform of simplex (+ its image under inversion) at a specified G-vector
template<int n> double Simplex<n>::getTilde(const Simplex<n>::Point& G) const
{	//Get list of dot products:
	std::array<double,n+1> d;
	for(int i=0; i<n; i++)
		d[i] = dot(v[i], G);
	d[n] = 0.; //corresponding to origin
	//Sort list:
	std::sort(d.begin(), d.end());
	//Detect clusters and accumulate contributions:
	const double dTol = 1e-6; //cluster tolerance
	double ret = 0.;
	for(auto iStart=d.begin(); iStart!=d.end();)
	{	//Detect cluster:
		double dPrev = *iStart, dSum = 0.;
		auto iEnd = iStart;
		while(iEnd!=d.end() && *iEnd<dPrev+dTol)
			dSum += (dPrev = *(iEnd++));
		int m = int(iEnd-iStart); //multiplicity
		double dMean = dSum/m;
		//Get contribution from cluster:
		//(use templates to agressively optimize code for each cluster size)
		getSimplexTilde_mSwitch<n,n+1>()(m, d, iStart, iEnd, dMean, ret);
		iStart = iEnd; //advance to next cluster
	}
	//Return result with prefactor:
	return (2*V) * ret;
}

//!@endcond
#endif // JDFTX_CORE_SIMPLEX_H
