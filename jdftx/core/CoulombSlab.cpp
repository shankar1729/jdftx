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

#include <core/CoulombSlab.h>
#include <core/Coulomb_internal.h>
#include <core/BlasExtra.h>

CoulombSlab::CoulombSlab(const GridInfo& gInfo, const CoulombTruncationParams& params)
: Coulomb(gInfo, params)
{	//TODO: Check if iDir is orthogonal to the rest
	die("Not yet implemented.\n");
}

DataGptr CoulombSlab::operator()(DataGptr&& in) const
{	int iDir = params.iDir;
	double hlfL = 0.5*gInfo.R(iDir,iDir);
	callPref(coulombAnalytic)(gInfo.S, gInfo.GGT, CoulombSlab_calc(iDir, hlfL), in->dataPref(false));
	return in;
}

double CoulombSlab::energyAndGrad(std::vector<Coulomb::PointCharge>& pointCharges) const
{	//TODO: 2D periodic Ewald sum
	return 0.;
}
