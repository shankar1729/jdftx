/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

This file is part of Fluid1D.

Fluid1D is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Fluid1D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Fluid1D.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#include <fluid/IdealGasPomega.h>
#include <fluid/FluidMixture.h>
#include <fluid/Euler.h>
#include <core/BlasExtra.h>

IdealGasPomega::IdealGasPomega(Fex* fex, double xBulk, const SO3quad& quad, const TranslationOperator& trans)
: IdealGas(quad.nOrientations(),fex,xBulk), quad(quad), trans(trans)
{
	site0mult = 0;
	while(molecule->site[site0mult].index==0) site0mult++;
}

void IdealGasPomega::initState(const ScalarField* Vex, ScalarField* logPomega, double scale, double Elo, double Ehi) const
{	ScalarFieldCollection Veff(molecule->nIndices); nullToZero(Veff, gInfo);
	for(int k=0; k<molecule->nIndices; k++)
	{	Veff[k] += V[k];
		Veff[k] += Vex[k];
	}
	double Emin=+DBL_MAX, Emax=-DBL_MAX, Emean=0.0;
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		ScalarField Emolecule;
		//Sum the potentials collected over sites for each orientation:
		for(int i=0; i<molecule->nSites; i++)
		{	const Site& s = molecule->site[i];
			trans.S_axpy(rot*s.pos, 1.0, Veff[s.index], Emolecule);
		}
		//Accumulate stats and cap:
		Emean += quad.weight(o) * integral(Emolecule)/gInfo.Volume();
		double Emin_o, Emax_o;
		eblas_capMinMax(gInfo.S, Emolecule.data(), Emin_o, Emax_o, Elo, Ehi);
		if(Emin_o<Emin) Emin=Emin_o;
		if(Emax_o>Emax) Emax=Emax_o;
		//Set contributions to the state (with appropriate scale factor):
		logPomega[o] = (-scale/T) * Emolecule;
	}
	//Print stats:
	logPrintf("\tIdealGasPomega[%s] single molecule energy: min = %le, max = %le, mean = %le\n",
		   molecule->name.c_str(), Emin, Emax, Emean);
}

void IdealGasPomega::getDensities(const ScalarField* logPomega, ScalarField* N, double& P) const
{	ScalarFieldCollection wN; nullToZero(wN, gInfo, molecule->nIndices);
	P = 0.;
	double& S = ((IdealGasPomega*)this)->S;
	S=0.0;
	//Loop over orientations:
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		ScalarField N_o = (quad.weight(o) * Nbulk) * exp(logPomega[o]); //contribution form this orientation
		ScalarField wN_o = DiagJdagOJ1(N_o);
		//Accumulate wN_o to each site density with appropriate translations:
		for(int i=0; i<molecule->nSites; i++)
			trans.Sdag_axpy(rot*molecule->site[i].pos, 1.0, wN_o, wN[molecule->site[i].index]);
		//Accumulate contributions to the entropy:
		S += dot(wN_o, logPomega[o]);
		//Accumulate the cell dipole moments:
		double pHat = rot(2,2);
		P += pHat * integral(N_o);
	}
	for(int k=0; k<molecule->nIndices; k++)
		N[k] = DiagJdagOJ1inv(wN[k]);
}

double IdealGasPomega::compute(const ScalarField* logPomega, const ScalarField* N, ScalarField* grad_N,
	const double& P, double& grad_P, const double Nscale, double& grad_Nscale) const
{	double PhiNI = 0.0;
	//Add contributions due to external potentials:
	for(int j=0; j<molecule->nIndices; j++)
		if(V[j])
		{	ScalarField JdagOJV = JdagOJ(V[j]);
			grad_N[j] += JdagOJV;
			PhiNI += dot(N[j], JdagOJV);
		}
	//Contributions due to uniform electric field:
	grad_P -= fluidMixture.Eexternal * molecule->get_dipole();
	PhiNI -= fluidMixture.Eexternal * P * molecule->get_dipole();
	//KE and mu:
	grad_N[0] -= DiagJdagOJ1(mu/site0mult, gInfo);
	PhiNI -= (T+mu)*integral(N[0])/site0mult;
	//Entropy (this part deals with Nscale explicitly, so need to increment grad_Nscale):
	grad_Nscale += T*S;
	PhiNI += Nscale*T*S;
	return PhiNI;
}

void IdealGasPomega::convertGradients(const ScalarField* logPomega, const ScalarField* N,
	const ScalarField* grad_N, double grad_P, ScalarField* grad_logPomega, const double Nscale) const
{
	ScalarFieldCollection grad_wN(molecule->nIndices);
	for(int k=0; k<molecule->nIndices; k++)
		grad_wN[k] = DiagJdagOJ1inv(grad_N[k]);
	//Loop over orientations:
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		ScalarField grad_wN_o; //gradient w.r.t wN_o (as calculated in getDensities)
		//Collect the contributions from each grad_wN in grad_wN_o
		for(int i=0; i<molecule->nSites; i++)
			trans.S_axpy(rot*molecule->site[i].pos, 1.0, grad_wN[molecule->site[i].index], grad_wN_o);
		//Collect the contributions the entropy:
		grad_wN_o += T * logPomega[o];
		//Collect the contribution from grad_P:
		double pHat = rot(2,2);
		grad_wN_o += grad_P * pHat;
		//Propagate grad_wN_o to grad_logPomega[o]:
		grad_logPomega[o] = Diag(DiagJdagOJ1((quad.weight(o) * Nbulk * Nscale) * exp(logPomega[o]))) * grad_wN_o;
	}
}

