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

#include <fluid/IdealGasMuEps.h>
#include <fluid/Euler.h>

IdealGasMuEps::IdealGasMuEps(Fex* fex, double xBulk, const SO3quad& quad, const TranslationOperator& trans)
: IdealGas(4,fex,xBulk), quad(quad), trans(trans)
{
	site0mult = 0;
	while(molecule->site[site0mult].index==0) site0mult++;
}

void IdealGasMuEps::initState(const DataRptr* Vex, DataRptr* mueps, double scale, double Elo, double Ehi) const
{	DataRptrCollection Veff(molecule->nIndices); nullToZero(Veff, gInfo);
	for(int i=0; i<4; i++) mueps[i]=0;
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
		//Add contributions to the state:
		vector3<> pHat = rot*vector3<>(0,0,1); //reference dipole is along z
		mueps[0] += quad.weight(o) * Emolecule;
		for(int k=0; k<3; k++)
			mueps[k+1] += (pHat[k]*quad.weight(o)) * Emolecule;
	}
	//Set the scale factor:
	for(int k=0; k<4; k++) mueps[k] *= (-scale/T);
	//Print stats:
	logPrintf("\tIdealGasMuEps[%s] single molecule energy: min = %le, max = %le, mean = %le\n",
		   molecule->name.c_str(), Emin, Emax, Emean);
}

void IdealGasMuEps::getDensities(const DataRptr* mueps, DataRptr* N, vector3<>& P) const
{	for(int k=0; k<molecule->nIndices; k++) N[k]=0;
	P = vector3<>(0,0,0);
	double& S = ((IdealGasMuEps*)this)->S;
	S=0.0;
	//Loop over orientations:
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		vector3<> pHat = rot*vector3<>(0,0,1);
		DataRptr sum_mueps; //the exponent in the boltzmann factor
		sum_mueps += mueps[0];
		for(int k=0; k<3; k++)
			sum_mueps += pHat[k]*mueps[k+1];
		DataRptr N_o = quad.weight(o) * Nbulk * exp(sum_mueps); //contribution from this orientation
		//Accumulate N_o to each site density with appropriate translations:
		for(int i=0; i<molecule->nSites; i++)
			trans.taxpy(rot*molecule->site[i].pos, 1.0, N_o, N[molecule->site[i].index]);
		//Accumulate contributions to the entropy:
		S += gInfo.dV*dot(N_o, mueps[0]);
		for(int k=0; k<3; k++)
			S += pHat[k] * gInfo.dV*dot(N_o, mueps[k+1]);
		//Accumulate the cell dipole moments:
		P += pHat * integral(N_o);
	}
}

double IdealGasMuEps::compute(const DataRptr* mueps, const DataRptr* N, DataRptr* grad_N,
	const vector3<>& P, vector3<>& grad_P, const double Nscale, double& grad_Nscale) const
{	double PhiNI = 0.0;
	//Add contributions due to external potentials:
	for(int j=0; j<molecule->nIndices; j++)
		if(V[j])
		{	grad_N[j] += V[j];
			PhiNI += gInfo.dV*dot(N[j], V[j]);
		}
	//Contributions due to uniform electric field:
	grad_P -= Eexternal * molecule->get_dipole();
	PhiNI -= dot(Eexternal, P) * molecule->get_dipole();
	//KE and mu:
	grad_N[0] -= mu/site0mult;
	PhiNI -= (T+mu)*integral(N[0])/site0mult;
	//Entropy (this part deals with Nscale explicitly, so need to increment grad_Nscale):
	grad_Nscale += T*S;
	PhiNI += Nscale*T*S;
	return PhiNI;
}

void IdealGasMuEps::convertGradients(const DataRptr* mueps, const DataRptr* N,
	const DataRptr* grad_N, vector3<> grad_P, DataRptr* grad_mueps, const double Nscale) const
{	for(int k=0; k<4; k++) grad_mueps[k]=0;
	//Loop over orientations:
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		vector3<> pHat = rot*vector3<>(0,0,1);
		DataRptr grad_N_o; //gradient w.r.t N_o (as calculated in getDensities)
		//Collect the contributions from each grad_N in grad_N_o
		for(int i=0; i<molecule->nSites; i++)
			trans.taxpy(-rot*molecule->site[i].pos, 1.0, grad_N[molecule->site[i].index], grad_N_o);
		//Collect the contributions the entropy:
		grad_N_o += T*mueps[0];
		for(int k=0; k<3; k++)
			grad_N_o += (T*pHat[k])*mueps[k+1];
		//Collect the contribution from grad_P:
		grad_N_o += dot(grad_P, pHat);
		//Calculate N_o again:
		DataRptr sum_mueps;
		sum_mueps += mueps[0];
		for(int k=0; k<3; k++)
			sum_mueps += pHat[k]*mueps[k+1];
		DataRptr N_o = (quad.weight(o) * Nbulk * Nscale) * exp(sum_mueps);
		//Accumulate N_o * grad_N_o into each component of grad_mueps with appropriate weights:
		DataRptr grad_mueps_term = N_o * grad_N_o;
		grad_mueps[0] += grad_mueps_term;
		for(int k=0; k<3; k++)
			grad_mueps[k+1] += pHat[k] * grad_mueps_term;
	}
}
