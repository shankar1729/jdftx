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
#include <core/LoopMacros.h>
#include <core/Thread.h>

std::shared_ptr<Coulomb> CoulombTruncationParams::createCoulomb(const GridInfo& gInfo) const
{	if(type != Periodic)
		logPrintf("\n---------- Setting up coulomb interaction ----------\n");
	switch(type)
	{	case Periodic:  return std::make_shared<CoulombPeriodic>(gInfo, *this);
		case Slab:      return std::make_shared<CoulombSlab>(gInfo, *this);
		case Wire:      return std::make_shared<CoulombWire>(gInfo, *this);
		case Isolated:  return std::make_shared<CoulombIsolated>(gInfo, *this);
		case Spherical: return std::make_shared<CoulombSpherical>(gInfo, *this);
		default: return 0; //never encountered (to suppress warning)
	}
}

Coulomb::Coulomb(const GridInfo& gInfo, const CoulombTruncationParams& params)
: gInfo(gInfo), params(params)
{
}

DataGptr Coulomb::operator()(const DataGptr& in) const
{	DataGptr out(in->clone()); //create destructible copy
	return (*this)((DataGptr&&)out);
}


//-------- CPU implementation of Coulomb_internal.h --------

template<typename Coulomb_calc>
void coulombAnalytic_thread(int iStart, int iStop, vector3<int> S, const matrix3<>& GGT, const Coulomb_calc& calc, complex* data)
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
