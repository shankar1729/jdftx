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

#include <core/CoulombPeriodic.h>
#include <core/CoulombSlab.h>
#include <core/CoulombWire.h>
#include <core/CoulombIsolated.h>
#include <core/Coulomb_internal.h>
#include <core/Coulomb_ExchangeEval.h>
#include <core/LoopMacros.h>
#include <core/BlasExtra.h>
#include <core/Thread.h>

std::shared_ptr<Coulomb> CoulombParams::createCoulomb(const GridInfo& gInfo) const
{	if(geometry != Periodic)
		logPrintf("\n---------- Setting up coulomb interaction ----------\n");
	switch(geometry)
	{	case Periodic:    return std::make_shared<CoulombPeriodic>(gInfo, *this);
		case Slab:        return std::make_shared<CoulombSlab>(gInfo, *this);
		case Wire:        return std::make_shared<CoulombWire>(gInfo, *this);
		case Cylindrical: return std::make_shared<CoulombCylindrical>(gInfo, *this);
		case Isolated:    return std::make_shared<CoulombIsolated>(gInfo, *this);
		case Spherical:   return std::make_shared<CoulombSpherical>(gInfo, *this);
		default: return 0; //never encountered (to suppress warning)
	}
}

vector3<bool> CoulombParams::isTruncated() const
{	switch(geometry)
	{	case Isolated:
		case Spherical:
			return vector3<bool>(true, true, true);
		case Cylindrical:
		case Wire:
		{	vector3<bool> result(true, true, true);
			result[iDir] = false;
			return result;
		}
		case Slab:
		{	vector3<bool> result(false, false, false);
			result[iDir] = true;
			return result;
		}
		case Periodic:
		default:
			return vector3<bool>(false, false, false);
	}
}


//--------------- class Coulomb ----------------

DataGptr Coulomb::operator()(const DataGptr& in) const
{	DataGptr out(in->clone()); //create destructible copy
	return (*this)((DataGptr&&)out);
}

complexDataGptr Coulomb::operator()(complexDataGptr&& in, vector3<> kDiff, double omega) const
{	auto exEvalOmega = exchangeEval.find(omega);
	assert(exEvalOmega != exchangeEval.end());
	return (*exEvalOmega->second)((complexDataGptr&&)in, kDiff);
}

complexDataGptr Coulomb::operator()(const complexDataGptr& in, vector3<> kDiff, double omega) const
{	complexDataGptr out(in->clone()); //create destructible copy
	return (*this)((complexDataGptr&&)out, kDiff, omega);
}

double Coulomb::energyAndGrad(std::vector<Atom>& atoms) const
{	if(!ewald) ((Coulomb*)this)->ewald = createEwald(gInfo.R, atoms.size());
	return ewald->energyAndGrad(atoms);
}

Coulomb::Coulomb(const GridInfo& gInfo, const CoulombParams& params)
: gInfo(gInfo), params(params)
{
}

void Coulomb::initExchangeEval()
{	//Initialize Exchange evaluators if required
	for(double omega: params.omegaSet)
		exchangeEval[omega]  = std::make_shared<ExchangeEval>(gInfo, params, *this, omega);
}


//-------- CPU implementation of Coulomb_internal.h --------

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
		threadLaunch(threadOperators ? 0 : 1, exchangeAnalytic_thread<Exchange##Type##_calc>, \
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
{	threadLaunch(threadOperators ? 0 : 1, multRealKernel_thread, S[0]*S[1]*S[2], S, kernel, data);
}

//Multiply a complexDataGptr by a kernel sampled with offset and rotation by rot
void multTransformedKernel_thread(size_t iStart, size_t iStop,
	vector3<int> S, const double* kernel, complex* data,
	const vector3<int>& offset, const matrix3<int>& rot)
{	THREAD_fullGspaceLoop( multTransformedKernel_calc(i, iG, S, kernel, data, offset, rot); )
}
void multTransformedKernel(vector3<int> S, const double* kernel, complex* data,
	const vector3<int>& offset, const matrix3<int>& rot)
{	threadLaunch(threadOperators ? 0 : 1, multTransformedKernel_thread, S[0]*S[1]*S[2],
		S, kernel, data, offset, rot);
}
