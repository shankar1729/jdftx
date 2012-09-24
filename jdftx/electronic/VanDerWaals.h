/*-------------------------------------------------------------------
Copyright 2012 Deniz Gunceler

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

#include <utility>

#include <core/vector3.h>
#include <core/Util.h>
#include <core/Coulomb.h>

#include <electronic/common.h>
#include <string.h>


class VanDerWaals
{
	public:
		
		void vanDerWaals();
		void setup(const Everything &everything);
		
		//! Returns the Van der Waals energy and gradient (in cartesian coordinates)
		//! s6 is the scaling parameter that depends on the Ex-Corr functional used
		double VDWEnergyAndGrad(std::vector<Atom>& atoms, string EXCorr);
		
	private:
		
		const Everything* e;
		
		//! C6 and R0 parameters for the VDW interactions
		struct VDWParameters
		{	double C6;
			double R0;
		};
		//! Returns the C6 coefficient and R0
		std::vector<VDWParameters> vdwParams;
		std::pair<double, double> getC6AndR0(int atomicNumber);
		
		//! Scaling factor that depends on the EXCorr functional
		std::map<string, double> scalingFactor;
};