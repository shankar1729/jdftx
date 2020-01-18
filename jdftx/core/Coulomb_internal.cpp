/*-------------------------------------------------------------------
Copyright 2020 Ravishankar Sundararaman

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

#include <core/Coulomb_internal.h>
#include <core/LoopMacros.h>
#include <core/Thread.h>

template<typename Coulomb_calc>
void coulombAnalytic_thread(size_t iStart, size_t iStop, vector3<int> S, const matrix3<>& GGT, const Coulomb_calc& calc, complex* data)
{	THREAD_halfGspaceLoop
	(	data[i] *= calc(iG, GGT);
	)
}
#define DECLARE_coulombAnalytic(Type) \
	void coulombAnalytic(vector3<int> S, const matrix3<>& GGT, const Coulomb##Type##_calc& calc, complex* data) \
	{	threadLaunch(coulombAnalytic_thread<Coulomb##Type##_calc>, S[0]*S[1]*(1+S[2]/2), S, GGT, calc, data); \
	}
DECLARE_coulombAnalytic(Periodic)
DECLARE_coulombAnalytic(Slab)
DECLARE_coulombAnalytic(Spherical)
#undef DECLARE_coulombAnalytic


template<typename Coulomb_calc>
void coulombAnalyticStress_thread(size_t iStart, size_t iStop, vector3<int> S, const matrix3<>& GGT, const Coulomb_calc& calc,
	const complex* X, const complex* Y, symmetricMatrix3<>* grad_RTR)
{
	THREAD_halfGspaceLoop
	(	double weight = ((iG[2]==0) or (2*iG[2]==S[2])) ? 1 : 2; //weight factor for points in reduced reciprocal space of real scalar fields
		grad_RTR[i] = (weight * real(X[i].conj() * Y[i])) * calc.latticeGradient(iG, GGT);
	)
}
#define DECLARE_coulombAnalyticStress(Type) \
	void coulombAnalyticStress(vector3<int> S, const matrix3<>& GGT, const Coulomb##Type##_calc& calc, \
		const complex* X, const complex* Y, symmetricMatrix3<>* grad_RTR) \
	{	\
		threadLaunch(coulombAnalyticStress_thread<Coulomb##Type##_calc>, S[0]*S[1]*(1+S[2]/2), S, GGT, calc, X, Y, grad_RTR); \
	}
DECLARE_coulombAnalyticStress(Periodic)
DECLARE_coulombAnalyticStress(Slab)
#undef DECLARE_coulombAnalyticStress


template<typename Exchange_calc>
void exchangeAnalytic_thread(size_t iStart, size_t iStop, vector3<int> S, const matrix3<>& GGT, const Exchange_calc& calc,
	complex* data, const vector3<>& kDiff, double Vzero, double thresholdSq)
{	THREAD_fullGspaceLoop
	(	double kplusGsq = GGT.metric_length_squared(iG + kDiff);
		data[i] *= kplusGsq<thresholdSq ? Vzero : calc(kplusGsq);
	)
}
//Specialization for slab mode:
template<> void exchangeAnalytic_thread<ExchangeSlab_calc>(size_t iStart, size_t iStop,
	vector3<int> S, const matrix3<>& GGT, const ExchangeSlab_calc& calc,
	complex* data, const vector3<>& kDiff, double Vzero, double thresholdSq)
{	THREAD_fullGspaceLoop
	(	data[i] *= calc(iG, GGT, kDiff, Vzero, thresholdSq);
	)
}
#define DECLARE_exchangeAnalytic(Type) \
	void exchangeAnalytic(vector3<int> S, const matrix3<>& GGT, const Exchange##Type##_calc& calc, \
		complex* data, const vector3<>& kDiff, double Vzero, double thresholdSq) \
	{	\
		threadLaunch(exchangeAnalytic_thread<Exchange##Type##_calc>, \
			S[0]*S[1]*S[2], S, GGT, calc, data, kDiff, Vzero, thresholdSq); \
	}
DECLARE_exchangeAnalytic(Periodic)
DECLARE_exchangeAnalytic(PeriodicScreened)
DECLARE_exchangeAnalytic(Spherical)
DECLARE_exchangeAnalytic(SphericalScreened)
DECLARE_exchangeAnalytic(Slab)
#undef DECLARE_exchangeAnalytic


void multRealKernel_thread(size_t iStart, size_t iStop,
	vector3<int> S, const double* kernel, complex* data)
{	THREAD_fullGspaceLoop( multRealKernel_calc(i, iG, S, kernel, data); )
}
void multRealKernel(vector3<int> S, const double* kernel, complex* data)
{	threadLaunch(multRealKernel_thread, S[0]*S[1]*S[2], S, kernel, data);
}

//Multiply a complexScalarFieldTilde by a kernel sampled with offset and rotation by rot
void multTransformedKernel_thread(size_t iStart, size_t iStop,
	vector3<int> S, const double* kernel, complex* data, const vector3<int>& offset)
{	THREAD_fullGspaceLoop( multTransformedKernel_calc(i, iG, S, kernel, data, offset); )
}
void multTransformedKernel(vector3<int> S, const double* kernel, complex* data, const vector3<int>& offset)
{	threadLaunch(multTransformedKernel_thread, S[0]*S[1]*S[2], S, kernel, data, offset);
}
