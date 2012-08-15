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

#include <fluid1D/IdealGasPsiAlpha.h>
#include <fluid/Euler.h>
#include <core/BlasExtra.h>
#include <core1D/Operators.h>

//! count the number of independent variables
int countIndep(const Molecule* molecule)
{	int nIndep = 0;
	for(int i=0; i<molecule->nSites; i++)
	{	if(((!i) || (molecule->site[i].index != molecule->site[i-1].index)) //this is a distinct site
			&& molecule->site[i].prop->indepSite) //should contribute to indep variables
			nIndep++;
	}
	return nIndep;
}

IdealGasPsiAlpha::IdealGasPsiAlpha(Fex* fex, double xBulk, const SO3quad& quad, const TranslationOperator& trans)
: IdealGas(countIndep(fex->getMolecule()),fex,xBulk), quad(quad), trans(trans)
{
	//Set the indep indices and count the site multiplicities:
	siteToIndep.assign(molecule->nSites, -1);
	densityToIndep.assign(molecule->nIndices, -1);
	indepMult.assign(nIndep, 0);
	indepMultTot=0;
	int cur_nIndep = 0;
	for(int i=0; i<molecule->nSites; i++)
	{	if(((!i) || (molecule->site[i].index != molecule->site[i-1].index)) //this is a distinct site
			&& molecule->site[i].prop->indepSite) //should contribute to indep variables
			cur_nIndep++;
		if(molecule->site[i].prop->indepSite)
		{	siteToIndep[i] = cur_nIndep-1;
			densityToIndep[molecule->site[i].index] = cur_nIndep-1;
			indepMult[siteToIndep[i]]++;
			indepMultTot++;
		}
	}
}

void IdealGasPsiAlpha::initState(const ScalarField* Vex, ScalarField* psi, double scale, double Elo, double Ehi) const
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
	}
	//Print stats:
	logPrintf("\tIdealGasPsiAlpha[%s] single molecule energy: min = %le, max = %le, mean = %le\n",
		   molecule->name.c_str(), Emin, Emax, Emean);
	//Initialize the state (simply a constant factor times the potential):
	for(int k=0; k<nIndep; k++) psi[k] *= 0.0;
	for(int i=0; i<molecule->nIndices; i++) psi[densityToIndep[i]] = (-scale/T)*Veff[i];
}

void IdealGasPsiAlpha::getDensities(const ScalarField* psi, ScalarField* N, double& P) const
{	ScalarFieldCollection wN; nullToZero(wN, gInfo, molecule->nIndices);
	//Loop over orientations:
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		ScalarField sum_psi; //the exponent in the boltzmann factor
		//Collect psi's from each site in this orientation:
		for(int i=0; i<molecule->nSites; i++)
			if(siteToIndep[i]>=0)
				trans.S_axpy(rot*molecule->site[i].pos, 1.0, psi[siteToIndep[i]], sum_psi);
		ScalarField wN_o = DiagJdagOJ1(quad.weight(o) * Nbulk * exp(sum_psi)); //contribution from this orientation
		//Accumulate wN_o to each site density with appropriate translations:
		for(int i=0; i<molecule->nSites; i++)
			trans.Sdag_axpy(rot*molecule->site[i].pos, 1.0, wN_o, wN[molecule->site[i].index]);
	}
	for(int k=0; k<molecule->nIndices; k++)
		N[k] = DiagJdagOJ1inv(wN[k]);
}

double IdealGasPsiAlpha::compute(const ScalarField* psi, const ScalarField* N, ScalarField* grad_N,
	const double& P, double& grad_P, const double Nscale, double& grad_Nscale) const
{	double PhiNI = 0.0;
	bool muAdded = false;
	for(int j=0; j<molecule->nIndices; j++)
	{	ScalarField PhiNI_Nj;
		int i = densityToIndep[j];
		if(V[j]) PhiNI_Nj += JdagOJ(V[j]); //all site densities may have an external potential
		if(i>=0)
		{	PhiNI_Nj += DiagJdagOJ1(T*psi[i] - (muAdded ? 0. : mu/indepMult[i])); //entropy and mu
			if(!muAdded)
			{	PhiNI -= T*integral(N[j])/indepMult[i]; //KE (Note: its d/dpsi cancels the del/delpsi of the entropy term)
				muAdded = true;
			}
		}
		if(PhiNI_Nj)
		{	grad_N[j] += PhiNI_Nj;
			PhiNI += dot(N[j], PhiNI_Nj);
		}
	}
	return PhiNI;
}

void IdealGasPsiAlpha::convertGradients(const ScalarField* psi, const ScalarField* N,
	const ScalarField* grad_N, double grad_P, ScalarField* grad_psi, const double Nscale) const
{
	ScalarFieldCollection grad_wN(molecule->nIndices);
	for(int k=0; k<molecule->nIndices; k++)
		grad_wN[k] = DiagJdagOJ1inv(grad_N[k]);
	for(int k=0; k<nIndep; k++) initZero(grad_psi[k], gInfo); 
	//Loop over orientations:
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		ScalarField grad_wN_o; //gradient w.r.t wN_o (as calculated in getDensities)
		//Collect the contributions from each grad_wN in grad_wN_o
		for(int i=0; i<molecule->nSites; i++)
			trans.S_axpy(rot*molecule->site[i].pos, 1.0, grad_wN[molecule->site[i].index], grad_wN_o);
		//Calculate wN_o again (with Nscale this time):
		ScalarField sum_psi;
		for(int i=0; i<molecule->nSites; i++)
			if(siteToIndep[i]>=0)
				trans.S_axpy(rot*molecule->site[i].pos, 1.0, psi[siteToIndep[i]], sum_psi);
		ScalarField wN_o = DiagJdagOJ1((quad.weight(o) * Nbulk * Nscale) * exp(sum_psi)); //contribution from this orientation
		//Accumulate wN_o * grad_wN_o into each grad_psi with appropriate translations:
		ScalarField grad_psi_term = Diag(wN_o) * grad_wN_o;
		for(int i=0; i<molecule->nSites; i++)
			if(siteToIndep[i]>=0)
				trans.Sdag_axpy(rot*molecule->site[i].pos, 1.0, grad_psi_term, grad_psi[siteToIndep[i]]);
	}
}

