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

#ifndef JDFTX_CORE_COULOMBISOLATED_H
#define JDFTX_CORE_COULOMBISOLATED_H

#include <core/GridInfo.h>
#include <core/WignerSeitz.h>
#include <core/string.h>

//! Coulomb interaction for an isolated system (no periodicity)
class CoulombIsolated
{
public:
	CoulombIsolated(const GridInfo& gInfo, double borderWidth, string filename=string());
	
	DataGptr operator()(const DataGptr&); //!< apply isolated Coulomb kernel
	
private:
	const GridInfo& gInfo;
	WignerSeitz ws;
	RealKernel Vc;
};

#endif // JDFTX_CORE_COULOMBISOLATED_H
