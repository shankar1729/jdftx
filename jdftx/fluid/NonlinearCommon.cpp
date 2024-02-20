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
	std::vector<double> gSamples, energySamples;
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
		gSamples.push_back(eps/x);
		energySamples.push_back(energy/(x*x));
	}
	gLookup.init(0, gSamples, dxMapped);
	dielEnergyLookup.init(0, energySamples, dxMapped);
	
// 	for(double xMapped=0.001; xMapped<1.0; xMapped+=0.01)
// 	{	double x = xMapped/(1.-xMapped);
// 		double E = x / dielectricEval->pByT;
// 		vector3<> Evec(0.0, 0.0, E);
// 		double s = 1.0;
// 		
// 		//Epsilon from phiToState:
// 		double epsilon = 0.0;
// 		vector3<const double*> Dphi_const(&Evec[0], &Evec[1], &Evec[2]);
// 		dielectricEval->phiToState_calc(0, Dphi_const, &s, gLookup, false, vector3<double*>(NULL, NULL, NULL), &epsilon);
// 		
// 		//Chi from apply:
// 		double A = 0.0;
// 		vector3<double*> Dphi(&Evec[0], &Evec[1], &Evec[2]);
// 		dielectricEval->apply_calc(0, dielEnergyLookup, &s, Dphi, &A);
// 		const double epsilon2 = (4*M_PI) * Dphi[2][0] / E;
// 		
// 		logPrintf("%.6lf %12.9lf %12.9lf\n", x, epsilon, epsilon2);
// 	}
// 	die("Testing.\n");

// 	FILE* fp = fopen("energyTest.dat", "w");
// 	for(double x=0.001; x<1.0; x+=0.001)
// 		fprintf(fp, "%lg %le\n", x, dielEnergyLookup(x));
// 	fclose(fp);
// 	die("Testing.\n");
	
	//Screening lookup table
	if(screeningEval)
	{
		double dVmapped = 1./1024;
		std::vector<double> samples, energySamples;
		for(double Vmapped=-1.; Vmapped<=1.; Vmapped+=dVmapped)
		{	if(fabs(Vmapped)==1)
			{	samples.push_back(0.);
				energySamples.push_back(0.);
			}
			else
			{	//double V = std::pow(Vmapped/(1.-Vmapped*Vmapped), 3.); //inverse of Vmapped = copysign(2cbrt(V) / (1 + sqrt(1 + (2cbrt(V))^2)), V)
				double V = Vmapped / (1 - Vmapped*Vmapped); //inverse of Vmapped = 2V/(1 + sqrt(1 + 4V^2))
				double x = screeningEval->x_from_V(V);
				samples.push_back(1.0 - x); //makes result -> 0 at edges, because x -> 1 for infinite |V|
				
				//Corresponding energy for the nonlinear screening:
				double f_x, f = screeningEval->fHS(x, f_x);
				double energy = screeningEval->NT * 
					( exp(+V - f_x * screeningEval->x0minus)
					+ exp(-V - f_x * screeningEval->x0plus)
					+ f_x * x - f);
				energySamples.push_back(1.0 / energy); //makes result -> 0 at edges, because energy -> infty at infinite |V|
			}
		}
		xLookup.init(1, samples, dVmapped);
		ionEnergyLookup.init(1, energySamples, dVmapped);
		
		/*
		for(double Vmapped=-0.992; Vmapped<1.0; Vmapped+=0.005)
		{	double V = Vmapped / (1 - Vmapped*Vmapped);
			const double phi = V / screeningEval->ZbyT;
			const double s = 1.0;
			
			//kappa^2 from phiToState:
			double kappaSq = 0.0;
			screeningEval->phiToState_calc(0, &phi, &s, xLookup, false, NULL, NULL, &kappaSq);
			
			//kappa^2 from apply:
			double A=0.0, A_phi=0.0;
			screeningEval->apply_calc(0, ionEnergyLookup, &s, &phi, &A, &A_phi);
			const double kappaSq2 = (4*M_PI) * A_phi / phi;
			
			logPrintf("%.6lf %12.9lf %12.9lf\n", V, kappaSq, kappaSq2);
		}
			
		FILE* fp = fopen("energyTest.dat", "w");
		for(double Vmapped=-1.0; Vmapped<=1.0; Vmapped+=0.0001)
			fprintf(fp, "%lg %le %le\n", Vmapped,
				xLookup(1.0 + Vmapped),
				ionEnergyLookup(1.0 + Vmapped));
		fclose(fp);
		die("Testing.\n");
		*/
	}
}

NonlinearCommon::~NonlinearCommon()
{	delete dielectricEval;
	gLookup.free();
	dielEnergyLookup.free();
	if(screeningEval)
	{	delete screeningEval;
		xLookup.free();
		ionEnergyLookup.free();
	}
}
