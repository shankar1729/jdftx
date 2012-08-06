/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

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

#include <fluid/IdealGasPsiAlpha.h>
#include <fluid/Euler.h>

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

void IdealGasPsiAlpha::initState(const DataRptr* Vex, DataRptr* psi, double scale, double Elo, double Ehi) const
{	DataRptrCollection Veff(molecule->nIndices); nullToZero(Veff, gInfo);
	for(int k=0; k<molecule->nIndices; k++)
	{	Veff[k] += V[k];
		Veff[k] += Vex[k];
	}
	double Emin=+DBL_MAX, Emax=-DBL_MAX, Emean=0.0;
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		DataRptr Emolecule;
		//Sum the potentials collected over sites for each orientation:
		for(int i=0; i<molecule->nSites; i++)
		{	const Site& s = molecule->site[i];
			trans.taxpy(-(rot*s.pos), 1.0, Veff[s.index], Emolecule);
		}
		//Accumulate stats and cap:
		Emean += quad.weight(o) * sum(Emolecule)/gInfo.nr;
		double Emin_o, Emax_o;
		callPref(eblas_capMinMax)(gInfo.nr, Emolecule->dataPref(), Emin_o, Emax_o, Elo, Ehi);
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

void IdealGasPsiAlpha::getDensities(const DataRptr* psi, DataRptr* N, vector3<>& P) const
{	for(int k=0; k<molecule->nIndices; k++) N[k]=0;
	//Loop over orientations:
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		DataRptr sum_psi; //the exponent in the boltzmann factor
		//Collect psi's from each site in this orientation:
		for(int i=0; i<molecule->nSites; i++)
			if(siteToIndep[i]>=0)
				trans.taxpy(-rot*molecule->site[i].pos, 1.0, psi[siteToIndep[i]], sum_psi);
		DataRptr N_o = quad.weight(o) * Nbulk * exp(sum_psi); //contribution from this orientation
		//Accumulate N_o to each site density with appropriate translations:
		for(int i=0; i<molecule->nSites; i++)
			trans.taxpy(rot*molecule->site[i].pos, 1.0, N_o, N[molecule->site[i].index]);
	}
}

double IdealGasPsiAlpha::compute(const DataRptr* psi, const DataRptr* N, DataRptr* grad_N,
	const vector3<>& P, vector3<>& grad_P, const double Nscale, double& grad_Nscale) const
{	double PhiNI = 0.0;
	bool muAdded = false;
	for(int j=0; j<molecule->nIndices; j++)
	{	DataRptr PhiNI_Nj;
		int i = densityToIndep[j];
		PhiNI_Nj += V[j]; //all site densities may have an external potential
		if(i>=0)
		{	PhiNI_Nj += T*psi[i]; //entropy part
			if(!muAdded)
			{	//Add terms associated with whole molecule (mu and KE) to first indep
				PhiNI_Nj -= mu/indepMult[i];
				PhiNI -= T*integral(N[j])/indepMult[i];
				muAdded = true;
			}
		}
		if(PhiNI_Nj)
		{	grad_N[j] += PhiNI_Nj;
			PhiNI += gInfo.dV*dot(N[j], PhiNI_Nj);
		}
	}
	return PhiNI;
}

void IdealGasPsiAlpha::convertGradients(const DataRptr* psi, const DataRptr* N,
	const DataRptr* grad_N, vector3<> grad_P, DataRptr* grad_psi, const double Nscale) const
{
	for(int k=0; k<nIndep; k++) grad_psi[k]=0;
	//Loop over orientations:
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		DataRptr grad_N_o; //gradient w.r.t N_o (as calculated in getDensities)
		//Collect the contributions from each grad_N in grad_N_o
		for(int i=0; i<molecule->nSites; i++)
			trans.taxpy(-rot*molecule->site[i].pos, 1.0, grad_N[molecule->site[i].index], grad_N_o);
		//Calculate N_o again (with Nscale this time):
		DataRptr sum_psi;
		for(int i=0; i<molecule->nSites; i++)
			if(siteToIndep[i]>=0)
				trans.taxpy(-rot*molecule->site[i].pos, 1.0, psi[siteToIndep[i]], sum_psi);
		DataRptr N_o = (quad.weight(o) * Nbulk * Nscale) * exp(sum_psi); //contribution from this orientation
		//Accumulate N_o * grad_N_o into each grad_psi with appropriate translations:
		DataRptr grad_psi_term = N_o * grad_N_o;
		for(int i=0; i<molecule->nSites; i++)
			if(siteToIndep[i]>=0)
				trans.taxpy(rot*molecule->site[i].pos, 1.0, grad_psi_term, grad_psi[siteToIndep[i]]);
	}
}

