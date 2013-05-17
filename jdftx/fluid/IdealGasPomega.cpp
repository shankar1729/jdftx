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

#include <fluid/IdealGasPomega.h>
#include <fluid/Euler.h>
#include <electronic/operators.h>

IdealGasPomega::IdealGasPomega(const FluidMixture* fluidMixture, const FluidComponent* comp, const SO3quad& quad, const TranslationOperator& trans, unsigned nIndepOverride)
: IdealGas(nIndepOverride ? nIndepOverride : quad.nOrientations(), fluidMixture, comp), quad(quad), trans(trans), pMol(molecule.getDipole())
{
}

string IdealGasPomega::representationName() const
{	return "Pomega";
}

void IdealGasPomega::initState_o(int o, const matrix3<>& rot, double scale, const DataRptr& Eo, DataRptr* logPomega) const
{	logPomega[o] += (-scale/T) * Eo;
}

void IdealGasPomega::getDensities_o(int o, const matrix3<>& rot, const DataRptr* logPomega, DataRptr& logPomega_o) const
{	logPomega_o += logPomega[o];
}

void IdealGasPomega::convertGradients_o(int o, const matrix3<>& rot, const DataRptr& Phi_logPomega_o, DataRptr* Phi_logPomega) const
{	Phi_logPomega[o] += Phi_logPomega_o;
}


void IdealGasPomega::initState(const DataRptr* Vex, DataRptr* indep, double scale, double Elo, double Ehi) const
{	for(int k=0; k<nIndep; k++)indep[k]=0;
	DataRptrCollection Veff(molecule.sites.size()); nullToZero(Veff, gInfo);
	for(unsigned i=0; i<molecule.sites.size(); i++)
	{	Veff[i] += V[i];
		Veff[i] += Vex[i];
	}
	double Emin=+DBL_MAX, Emax=-DBL_MAX, Emean=0.0;
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		DataRptr Emolecule;
		//Sum the potentials collected over sites for each orientation:
		for(unsigned i=0; i<molecule.sites.size(); i++)
			for(vector3<> pos: molecule.sites[i]->positions)
				trans.taxpy(-(rot*pos), 1., Veff[i], Emolecule);
		//Accumulate stats and cap:
		Emean += quad.weight(o) * sum(Emolecule)/gInfo.nr;
		double Emin_o, Emax_o;
		callPref(eblas_capMinMax)(gInfo.nr, Emolecule->dataPref(), Emin_o, Emax_o, Elo, Ehi);
		if(Emin_o<Emin) Emin=Emin_o;
		if(Emax_o>Emax) Emax=Emax_o;
		//Set contributions to the state (with appropriate scale factor):
		initState_o(o, rot, scale, Emolecule, indep);
	}
	//Print stats:
	logPrintf("\tIdealGas%s[%s] single molecule energy: min = %le, max = %le, mean = %le\n",
		   representationName().c_str(), molecule.name.c_str(), Emin, Emax, Emean);
}

void IdealGasPomega::getDensities(const DataRptr* indep, DataRptr* N, vector3<>& P0) const
{	for(unsigned i=0; i<molecule.sites.size(); i++) N[i]=0;
	double& S = ((IdealGasPomega*)this)->S;
	S=0.0;
	DataRptrVec P;
	//Loop over orientations:
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		DataRptr logPomega_o; getDensities_o(o, rot, indep,logPomega_o);
		DataRptr N_o = (quad.weight(o) * Nbulk) * exp(logPomega_o); //contribution form this orientation
		//Accumulate N_o to each site density with appropriate translations:
		for(unsigned i=0; i<molecule.sites.size(); i++)
			for(vector3<> pos: molecule.sites[i]->positions)
				trans.taxpy(rot*pos, 1., N_o, N[i]);
		//Accumulate contributions to the entropy:
		S += gInfo.dV*dot(N_o, logPomega_o);
		//Accumulate the polarization density:
		if(pMol.length_squared()) P += (rot * pMol) * N_o;
	}
	
	//Compute and cache dipole correlation correction:
	IdealGasPomega* cache = ((IdealGasPomega*)this);
	if(pMol.length_squared())
	{	P0 = sumComponents(P) / gInfo.nr;
		cache->Ecorr_P = I(molecule.mfKernel*(molecule.mfKernel*(corrPrefac*J(P))));
		cache->Ecorr = 0.5*gInfo.dV*dot(cache->Ecorr_P, P);
	}
	else
	{	P0 = vector3<>();
		cache->Ecorr = 0;
		cache->Ecorr_P = 0;
	}
}

double IdealGasPomega::compute(const DataRptr* indep, const DataRptr* N, DataRptr* Phi_N, const double Nscale, double& Phi_Nscale) const
{	double PhiNI = 0.0;
	//Add contributions due to external potentials:
	for(unsigned i=0; i<molecule.sites.size(); i++)
		if(V[i])
		{	Phi_N[i] += V[i];
			PhiNI += gInfo.dV*dot(N[i], V[i]);
		}
	//KE and mu:
	double invSite0mult = 1./molecule.sites[0]->positions.size();
	Phi_N[0] -= mu * invSite0mult;
	PhiNI -= (T+mu)*integral(N[0])* invSite0mult;
	//Entropy and correlation correction (this part deals with Nscale explicitly, so need to increment Phi_Nscale):
	Phi_Nscale += (T*S + Ecorr);
	PhiNI += Nscale*(T*S + Ecorr);
	return PhiNI;
}

void IdealGasPomega::convertGradients(const DataRptr* indep, const DataRptr* N, const DataRptr* Phi_N, const vector3<>& Phi_P0, DataRptr* Phi_indep, const double Nscale) const
{	for(int k=0; k<nIndep; k++) Phi_indep[k]=0;
	//Loop over orientations:
	for(int o=0; o<quad.nOrientations(); o++)
	{	matrix3<> rot = matrixFromEuler(quad.euler(o));
		DataRptr logPomega_o; getDensities_o(o, rot, indep, logPomega_o);
		DataRptr N_o = (quad.weight(o) * Nbulk * Nscale) * exp(logPomega_o);
		DataRptr Phi_N_o; //gradient w.r.t N_o (as calculated in getDensities)
		//Collect the contributions from each Phi_N in Phi_N_o
		for(unsigned i=0; i<molecule.sites.size(); i++)
			for(vector3<> pos: molecule.sites[i]->positions)
				trans.taxpy(-rot*pos, 1., Phi_N[i], Phi_N_o);
		//Collect the contributions from the entropy:
		Phi_N_o += T*logPomega_o;
		//Collect the contribution from Phi_P0 and Ecorr_P:
		if(pMol.length_squared()) Phi_N_o += dot(rot * pMol, Nscale*Ecorr_P) + dot(rot * pMol, Phi_P0);
		//Propagate Phi_N_o to Phi_logPomega_o and then to Phi_indep:
		convertGradients_o(o, rot, N_o*Phi_N_o, Phi_indep);
	}
}

