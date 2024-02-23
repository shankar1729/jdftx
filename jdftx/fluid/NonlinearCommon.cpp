/*-------------------------------------------------------------------
Copyright 2024 Ravishankar Sundararaman

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

#include <fluid/NonlinearCommon.h>
#include <fluid/PCM_internal.h>


NonlinearCommon::NonlinearCommon(const FluidSolverParams& fsp, double epsBulk)
: pMol(0.), ionNbulk(0.), ionZ(0.), screeningEval(0), dielectricEval(0)
{
	const auto& solvent = fsp.solvents[0];
	pMol = solvent->pMol ? solvent->pMol : solvent->molecule.getDipole().length();
	
	//Initialize dielectric evaluation class:
	dielectricEval = new NonlinearPCMeval::Dielectric(fsp.linearDielectric,
		fsp.T, solvent->Nbulk, pMol, epsBulk, solvent->epsInf);

	//Check and setup ionic screening:
	if(fsp.cations.size() > 1) die("NonlinearPCM currently only supports a single cationic component.\n");
	if(fsp.anions.size() > 1) die("NonlinearPCM currently only supports a single anionic component.\n");
	assert(fsp.anions.size() == fsp.cations.size()); //this should be ensured by charge neutrality check in FluidSolver constructor
	if(fsp.cations.size())
	{	//Ensure charge balanced:
		if(fabs(fsp.cations[0]->molecule.getCharge() + fsp.anions[0]->molecule.getCharge())>1e-12)
			die("NonlinearPCM currently only supports charge-balanced (Z:Z) electrolytes.\n");
		ionNbulk = fsp.cations[0]->Nbulk;
		ionZ = fsp.anions[0]->molecule.getCharge();
		double VhsCation = fsp.cations[0]->molecule.getVhs(); if(!VhsCation) VhsCation = (4*M_PI/3)*pow(fsp.cations[0]->Rvdw,3);
		double VhsAnion = fsp.anions[0]->molecule.getVhs(); if(!VhsAnion) VhsAnion = (4*M_PI/3)*pow(fsp.anions[0]->Rvdw,3);
		screeningEval = new NonlinearPCMeval::Screening(fsp.linearScreening, fsp.T, ionNbulk, ionZ, VhsCation, VhsAnion, epsBulk);
	}
	else
	{	ionNbulk = 0.;
		screeningEval = 0;
	}
	
	//Dielectric lookup table
	double dxMapped = 1./512;
	std::vector<double> samples;
	for(double xMapped=0.; xMapped<=1.; xMapped+=dxMapped)
	{	double x;
		if(xMapped==0.) x = 1e-12;
		else if(xMapped==1.) x = 1e+12;
		else x = xMapped/(1.-xMapped); //inverse of xMapped = x / (1 + x)
		double eps = dielectricEval->eps_from_x(x), frac, logsinch;
		dielectricEval->calcFunctions(eps, frac, logsinch);
		double energy = dielectricEval->NT * (
			logsinch - 0.5 * dielectricEval->alpha * std::pow(eps * frac, 2)
			+ 0.5 * dielectricEval->X * (x * x)
		);
		samples.push_back(energy/(x*x));
	}
	dielEnergyLookup.init(0, samples, dxMapped);
	
	//Screening lookup table
	if(screeningEval)
	{
		double dVmapped = 1./512;
		std::vector<double> samples;
		for(double Vmapped=-1.; Vmapped<=1.; Vmapped+=dVmapped)
		{	if(fabs(Vmapped)==1)
				samples.push_back(0.);
			else
			{	double V = Vmapped / (1 - Vmapped*Vmapped); //inverse of Vmapped = 2V/(1 + sqrt(1 + 4V^2))
				double x = screeningEval->x_from_V(V);
				
				//Corresponding energy for the nonlinear screening:
				double f_x, f = screeningEval->fHS(x, f_x);
				double energy =
					( exp(+V - f_x * screeningEval->x0minus)
					+ exp(-V - f_x * screeningEval->x0plus)
					+ f_x * x - f - 2.0); //upto factor of NT; zero energy at V=0 and energy -> infty at infinite |V|
				samples.push_back(1./(1. + energy)); //map from [0, infty) -> [1, 0), with zero at Vmapped = +/-1
			}
		}
		ionEnergyLookup.init(1, samples, dVmapped);
	}
}

NonlinearCommon::~NonlinearCommon()
{	delete dielectricEval;
	dielEnergyLookup.free();
	if(screeningEval)
	{	delete screeningEval;
		ionEnergyLookup.free();
	}
}
