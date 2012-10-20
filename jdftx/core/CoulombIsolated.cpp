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

#include <core/CoulombIsolated.h>
#include <core/Coulomb_internal.h>
#include <core/CoulombKernel.h>
#include <core/Operators.h>
#include <core/Util.h>
#include <core/Thread.h>
#include <core/LoopMacros.h>
#include <core/BlasExtra.h>

//! Analog of ewald sum for isolated systems
//! (no Ewald trick required; prefixed 'Ewald' for consistent naming)
struct EwaldIsolated
{	const GridInfo& gInfo;
	const WignerSeitz& ws; //!< Wigner-Seitz cell
	double Rc; //!< cutoff radius for spherical mode (used for ion overlap checks only)
	
	EwaldIsolated(const GridInfo& gInfo, const WignerSeitz& ws, double Rc=0.)
	: gInfo(gInfo), ws(ws), Rc(Rc)
	{
	}
	
	double energyAndGrad(std::vector<Atom>& atoms) const
	{	if(!atoms.size()) return 0.;
		double E = 0.;
		//Shift all points into a Wigner-Seitz cell centered on one of the atoms; choice of this atom
		//is irrelevant if every atom lies in the WS cell of the other with a consistent translation:
		vector3<> pos0 = atoms[0].pos;
		for(Atom& a: atoms)
			a.pos = pos0 + ws.restrict(a.pos - pos0);
		//Loop over all pairs of atoms:
		for(unsigned i=0; i<atoms.size(); i++)
		{	Atom& a1 = atoms[i];
			for(unsigned j=0; j<i; j++)
			{	Atom& a2 = atoms[j];
				vector3<> x = a1.pos - a2.pos; //lattice coords
				double rSq = gInfo.RTR.metric_length_squared(x), r = sqrt(rSq);
				if((!(ws.restrict(x)==x)) || (Rc && (r>=Rc)))
					die("Separation between atoms %d and %d lies outside extent of coulomb truncation.\n", i, j);
				double dE = (a1.Z * a2.Z) / r;
				vector3<> dF = (gInfo.RTR * x) * (dE/rSq);
				E += dE;
				a1.force += dF;
				a2.force -= dF;
			}
		}
		return E;
	}
};

//----------------- class CoulombIsolated ---------------------

CoulombIsolated::CoulombIsolated(const GridInfo& gInfo, const CoulombParams& params)
: Coulomb(gInfo, params), ws(gInfo.R), Vc(gInfo)
{	//Compute kernel:
	CoulombKernel(gInfo.R, gInfo.S, params.isTruncated()).compute(Vc.data, ws);
	Vc.set();
	initExchangeEval();
}

DataGptr CoulombIsolated::operator()(DataGptr&& in) const
{	return Vc * in;
}

double CoulombIsolated::energyAndGrad(std::vector<Atom>& atoms) const
{	return EwaldIsolated(gInfo, ws).energyAndGrad(atoms);
}


//----------------- class CoulombSpherical ---------------------

CoulombSpherical::CoulombSpherical(const GridInfo& gInfo, const CoulombParams& params)
: Coulomb(gInfo, params), ws(gInfo.R), Rc(params.Rc)
{	double RcMax = ws.inRadius();
	if(Rc > RcMax)
		die("Spherical truncation radius %lg exceeds Wigner-Seitz cell in-radius of %lg bohrs.\n", Rc, RcMax);
	if(!Rc) Rc = RcMax;
	logPrintf("Initialized spherical truncation of radius %lg bohrs\n", Rc);
	initExchangeEval();
}

DataGptr CoulombSpherical::operator()(DataGptr&& in) const
{	callPref(coulombAnalytic)(gInfo.S, gInfo.GGT, CoulombSpherical_calc(Rc), in->dataPref(false));
	return in;
}

double CoulombSpherical::energyAndGrad(std::vector<Atom>& atoms) const
{	return EwaldIsolated(gInfo, ws, Rc).energyAndGrad(atoms);
}
