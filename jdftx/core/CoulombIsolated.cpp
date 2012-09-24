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
//! (no Ewald trick required, just for consistent naming)
struct EwaldIsolated
{	const GridInfo& gInfo;
	const WignerSeitz& ws;
	bool wsTruncated; //true => Wigner-Seitz truncation, false => spherical
	double criticalDist; //borderWidth for Wigner-Seitz, Rc for spherical
	
	EwaldIsolated(const GridInfo& gInfo, const WignerSeitz& ws, bool wsTruncated, double criticalDist)
	: gInfo(gInfo), ws(ws), wsTruncated(wsTruncated), criticalDist(criticalDist)
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
		//Loop over all pairs of pointcharges:
		for(unsigned i=0; i<atoms.size(); i++)
		{	Atom& a1 = atoms[i];
			for(unsigned j=0; j<i; j++)
			{	Atom& a2 = atoms[j];
				vector3<> x = a1.pos - a2.pos; //lattice coords
				double rSq = gInfo.RTR.metric_length_squared(x), r = sqrt(rSq);
				if(wsTruncated)
				{	if(ws.boundaryDistance(x) <= criticalDist)
						die("Separation between atoms %d and %d lies in the truncation border + margin.\n", i, j);
				}
				else
				{	if(r >= criticalDist)
						die("Atoms %d and %d are separated by r = %lg >= Rc-ionMargin = %lg bohrs.\n", i, j, r, criticalDist);
				}
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
{
	//Select gauss-smoothing parameter:
	double maxBorderWidth = sqrt(0.5) * ws.inRadius();
	if(params.borderWidth > maxBorderWidth)
		die("Border width %lg bohrs must be less than %lg bohrs (Wigner-Seitz cell in-radius/sqrt(2)).\n",
			params.borderWidth, maxBorderWidth);
	double sigmaBorder = params.borderWidth / CoulombKernelDesc::nSigmasPerWidth;
	logPrintf("Selecting gaussian width %lg bohrs (for border width %lg bohrs).\n", sigmaBorder, params.borderWidth);

	//Create kernel description:
	vector3<bool> isTruncated(true, true, true);
	vector3<> sigmaBorders(sigmaBorder, sigmaBorder, sigmaBorder);
	CoulombKernelDesc kernelDesc(gInfo.R, gInfo.S, isTruncated, sigmaBorders);
	
	if(!kernelDesc.loadKernel(Vc.data, params.filename)) //Try reading the kernel
	{	kernelDesc.computeKernel(Vc.data, ws); //Compute the kernel
		kernelDesc.saveKernel(Vc.data, params.filename); //Save kernel if requested
	}
	Vc.set();
	initExchangeEval();
}

DataGptr CoulombIsolated::operator()(DataGptr&& in) const
{	return Vc * in;
}

double CoulombIsolated::energyAndGrad(std::vector<Atom>& atoms) const
{	return EwaldIsolated(gInfo, ws, true, params.borderWidth + params.ionMargin).energyAndGrad(atoms);
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
{	return EwaldIsolated(gInfo, ws, false, Rc - params.ionMargin).energyAndGrad(atoms);
}
