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

#include <core/CoulombWire.h>
#include <core/Operators.h>
#include <core/Util.h>

CoulombWire::CoulombWire(const GridInfo& gInfo, const CoulombTruncationParams& params)
: Coulomb(gInfo, params), ws(gInfo.R), Vc(gInfo)
{	//TODO: Lots!
	die("Not yet implemented.\n");
}

DataGptr CoulombWire::operator()(DataGptr&& in) const
{	return Vc * in;
}

double CoulombWire::energyAndGrad(std::vector<Coulomb::PointCharge>& pointCharges) const
{	//TODO: 1D periodic Ewald sum
	return 0.;
}
