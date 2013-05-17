/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver

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

#include <fluid/FluidMixture.h>
#include <fluid/IdealGas.h>
#include <electronic/operators.h>

extern string rigidMoleculeCDFT_ScalarEOSpaper;

FluidMixture::FluidMixture(const GridInfo& gInfo, const double T)
: gInfo(gInfo), T(T), verboseLog(false), Qtol(1e-12), nIndepIdgas(0), nDensities(0), polarizable(false)
{
	logPrintf("Initializing fluid mixture at T=%lf K ...\n", T/Kelvin);
	Citations::add("Rigid-molecule density functional theory framework", rigidMoleculeCDFT_ScalarEOSpaper);
}

FluidMixture::~FluidMixture()
{	
}

void FluidMixture::initialize(double p, double epsBulkOverride, double epsInfOverride)
{	logPrintf("Adjusting fluid pressure to p=%lf bar\n", p/Bar);
	//Compute the maximum possible density (core packed limit)
	double Nguess=0., n3=0.;
	assert(component.size());
	for(const FluidComponent* c: component)
	{	Nguess += c->Nbulk;
		n3 += c->Nbulk * c->molecule.getVhs();
	}
	const double mulStep = 0.99;
	double Nstart = Nguess;
	if(n3 > mulStep) Nstart *= mulStep/n3; //ensure that mixtyure doesn;t exceed packing limit
	double pTest = compute_p(Nstart);
	//Find an interval of N that brackets P:
	double Nlo, Nhi;
	if(pTest > p)
	{	Nlo = Nstart;
		do
		{	Nhi = Nlo;
			Nlo = mulStep*Nhi;
			pTest = compute_p(Nlo);
		}
		while(pTest>p);
	}
	else
	{	Nhi = Nstart;
		do
		{	Nlo = Nhi;
			Nhi = Nlo/mulStep;
			pTest = compute_p(Nhi);
		}
		while(pTest<p);
	}
	//Bisect on N to get the pressure:
	double Ntot;
	do
	{	Ntot = 0.5*(Nlo + Nhi);
		pTest = compute_p(Ntot);
		if(pTest<p) Nlo = Ntot;
		else Nhi = Ntot;
	}
	while(Nhi-Nlo>1e-12*Nhi);
	//Set the chemical potentials and bulk densities for each component:
	std::vector<double> Nmol(component.size()), Phi_Nmol(component.size());
	for(unsigned ic=0; ic<component.size(); ic++) Nmol[ic] = (Ntot/Nguess)*component[ic]->Nbulk;
	computeUniformEx(Nmol, Phi_Nmol);
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const FluidComponent& c = *component[ic];
		c.idealGas->Nbulk = Nmol[ic];
		c.idealGas->mu = Phi_Nmol[ic];
		logPrintf("   Component '%s' at bulk density %le bohr^-3\n", c.molecule.name.c_str(), Nmol[ic]);
	}
	this->p = p;
	
	//Determine dipole correlation factor:
	double epsBulk = epsBulkOverride, epsInf = epsInfOverride;
	if(!epsBulk)
	{	epsBulk = 1.;
		for(const auto& c: component)
			epsBulk += (c->idealGas->get_Nbulk()/c->pureNbulk(T)) * (c->epsBulk - 1.);
	}
	if(!epsInf)
	{	epsInf = 1.;
		for(const auto& c: component)
			epsInf += (c->idealGas->get_Nbulk()/c->pureNbulk(T)) * (c->epsInf - 1.);
	}
	double chiRot = 0., chiPol = 0.;
	for(const FluidComponent* c: component)
	{	chiRot += c->idealGas->get_Nbulk() * c->molecule.getDipole().length_squared()/(3.*T);
		chiPol += c->idealGas->get_Nbulk() * c->molecule.getAlphaTot();
	}
	Crot = (epsBulk>epsInf && chiRot) ? (epsBulk-epsInf)/(4.*M_PI*chiRot) : 1.;
	Cpol = (epsInf>1. && chiPol) ? (epsInf-1.)/(4.*M_PI*chiPol) : 1.;
	logPrintf("   Local polarization-density correlation factors, Crot: %lg  Cpol: %lg\n", Crot, Cpol);
	for(const FluidComponent* c: component)
	{	double pMolSq = c->molecule.getDipole().length_squared();
		if(pMolSq) c->idealGas->corrPrefac = (1./Crot-1.)*3*T/(pMolSq*c->idealGas->get_Nbulk());
	}
	
	//Initialize preconditioners:
	Kindep.resize(component.size());
	for(unsigned ic=0; ic<component.size(); ic++)
	{	//Determine second derivative of total excess functional:
		double Nbulk = Nmol[ic];
		const double dNfac = 1e-4;
		Nmol[ic]=Nbulk*(1.+dNfac); computeUniformEx(Nmol, Phi_Nmol); double Phi_Np = Phi_Nmol[ic];
		Nmol[ic]=Nbulk*(1.-dNfac); computeUniformEx(Nmol, Phi_Nmol); double Phi_Nm = Phi_Nmol[ic];
		Nmol[ic]=Nbulk;
		double NNPhi_NN = Nbulk * (Phi_Np - Phi_Nm) / (2.*dNfac);
		if(NNPhi_NN < -0.5*Nbulk*T) NNPhi_NN = -0.5*Nbulk*T;
		//Set preconditioner:
		Kindep[ic] = 1./(gInfo.dV * (Nbulk*T + NNPhi_NN));
	}
	if(polarizable)
	{	double chiPol = 0.;
		for(const FluidComponent* c: component)
			chiPol += c->idealGas->get_Nbulk() * c->molecule.getAlphaTot();
		Keps = 1./(gInfo.dV * chiPol);
	}
}

const std::vector<const FluidComponent*>& FluidMixture::getComponents() const
{	return component;
}


void FluidMixture::addComponent(FluidComponent* comp)
{	component.push_back(comp);
	//Set the offsets for this component:
	comp->offsetIndep = nIndepIdgas;
	comp->offsetDensity = nDensities;
	//Update the totals, which become the offset for the next component
	nIndepIdgas += comp->idealGas->nIndep;
	nDensities += comp->molecule.sites.size();
	//Update the polarizable flag:
	polarizable |= bool(comp->molecule.getAlphaTot());
}

void FluidMixture::addFmix(const Fmix* fmix)
{	fmixArr.push_back(fmix);
}


void FluidMixture::initState(double scale, double Elo, double Ehi)
{	//Compute the effective nonlinear coupling potential for the uniform fluid:
	DataRptrCollection Vex(nDensities);
	{	//Get the uniform fluid site densities:
		std::vector<double> Nbulk(nDensities);
		for(const FluidComponent* c: component)
			for(unsigned i=0; i<c->molecule.sites.size(); i++)
				Nbulk[c->offsetDensity + i] = c->idealGas->get_Nbulk() * c->molecule.sites[i]->positions.size();
		DataGptrCollection Ntilde(nDensities);
		nullToZero(Ntilde, gInfo);
		for(unsigned i=0; i<nDensities; i++) Ntilde[i]->setGzero(Nbulk[i]);
		//Set Vex to the difference between the potentials returned by Fmix::compute() and Fmix::computeUniform():
		DataGptrCollection Phi_Ntilde(nDensities);
		std::vector<double> Phi_Nbulk(nDensities);
		for(const Fmix* fmix: fmixArr)
		{	fmix->compute(Ntilde, Phi_Ntilde);
			fmix->computeUniform(Nbulk, Phi_Nbulk);
		}
		for(unsigned i=0; i<nDensities; i++)
			if(Phi_Ntilde[i]) Vex[i] = Jdag(Phi_Ntilde[i]) - Phi_Nbulk[i];
		//Add electrostatic coupling:
		if(rhoExternal)
		{	DataGptr dExternal = (-4*M_PI)*Linv(O(rhoExternal));
			for(const FluidComponent* c: component)
				for(unsigned i=0; i<c->molecule.sites.size(); i++)
				{	const Molecule::Site& s = *(c->molecule.sites[i]);
					if(s.chargeKernel) Vex[c->offsetDensity+i] += I(s.chargeKernel * dExternal);
				}
		}
	}
	state.assign(get_nIndep(), 0);
	//Call initState for each component
	logPrintf("\n----- FluidMixture::initState() -----\n");
	for(const FluidComponent* c: component)
		c->idealGas->initState(&Vex[c->offsetDensity], &state[c->offsetIndep], scale, Elo, Ehi);
	
	//Guess polarizability independent variables:
	if(polarizable)
	{	DataGptr rhoTot; vector3<> Prot0; double epsInf=1.;
		for(unsigned ic=0; ic<component.size(); ic++)
		{	const FluidComponent& c = *(component[ic]);
			DataRptrCollection N(c.molecule.sites.size()); vector3<> P0c;
			c.idealGas->getDensities(&state[c.offsetIndep], &N[0], P0c);
			Prot0 += P0c;
			for(unsigned i=0; i<c.molecule.sites.size(); i++)
			{	const Molecule::Site& s = *(c.molecule.sites[i]);
				if(s.chargeKernel)
					rhoTot += s.chargeKernel(0) * (c.molecule.mfKernel * J(N[i]));
				if(s.polKernel && rhoExternal)
					rhoTot += (Cpol * s.alpha) * divergence(c.molecule.mfKernel * J(N[i] * I(gradient(s.polKernel*coulomb(rhoExternal)))));
				epsInf += 4*M_PI * (integral(N[i])/gInfo.detR) * Cpol * s.alpha;
			}
		}
		DataGptrVec epsMF = gradient((-1./epsInf)*coulomb(rhoTot));
		vector3<> epsMF0 = (Eexternal - 4*M_PI*Prot0)/epsInf;
		for(int k=0; k<3; k++) state[nIndepIdgas+k] = I(epsMF[k]) + epsMF0[k];
	}
	logPrintf("\n");
}

void FluidMixture::loadState(const char* filename)
{	nullToZero(state, gInfo, get_nIndep());
	loadFromFile(state, filename);
}

void FluidMixture::saveState(const char* filename) const
{	saveToFile(state, filename);
}

FluidMixture::Outputs::Outputs(DataRptrCollection* N, vector3<>* electricP,
	DataGptr* Phi_rhoExternal, DataRptrCollection* psiEff, EnergyComponents* Phi)
:N(N),electricP(electricP),Phi_rhoExternal(Phi_rhoExternal),psiEff(psiEff),Phi(Phi)
{
}

double FluidMixture::getFreeEnergy(Outputs outputs) const
{	DataRptrCollection Phi_state;
	return (*this)(state, Phi_state, outputs);
	//calls operator() and provides whichever outputs are requested
	//before calling must create an Outputs structure with the required quantities allocated.
}

void FluidMixture::step(const DataRptrCollection& dir, double alpha)
{	axpy(alpha, dir, state);
}

double FluidMixture::compute(DataRptrCollection* grad)
{	DataRptrCollection tempGrad;
	return (*this)(state, grad ? *grad : tempGrad, Outputs());
}

DataRptrCollection FluidMixture::precondition(const DataRptrCollection& grad)
{	DataRptrCollection Kgrad(get_nIndep());
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const FluidComponent& c = *component[ic];
		for(unsigned k=c.offsetIndep; k<c.offsetIndep+c.idealGas->nIndep; k++)
			Kgrad[k] = Kindep[ic]*grad[k];
	}
	for(unsigned k=nIndepIdgas; k<get_nIndep(); k++)
		Kgrad[k] = Keps * grad[k];
	return Kgrad;
}

double FluidMixture::compute_p(double Ntot) const
{	std::vector<double> Nmol(component.size()), Phi_Nmol(component.size());
	double Nguess=0.;
	for(const FluidComponent* c: component) Nguess += c->Nbulk;
	for(unsigned ic=0; ic<component.size(); ic++) Nmol[ic] = (Ntot/Nguess)*component[ic]->Nbulk;
	//p = N T + Sum_ic Nic d(aEx)/dNic - aEx where aEx is Helmholtz energy density
	double p = Ntot*T - computeUniformEx(Nmol, Phi_Nmol);
	for(unsigned ic=0; ic<component.size(); ic++)
		p += Nmol[ic]*Phi_Nmol[ic];
	return p;
}
