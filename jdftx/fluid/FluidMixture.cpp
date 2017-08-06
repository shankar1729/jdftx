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
#include <gsl/gsl_multiroots.h>

extern string rigidMoleculeCDFT_ScalarEOSpaper;

FluidMixture::FluidMixture(const GridInfo& gInfo, const double T)
: gInfo(gInfo), T(T), verboseLog(false), useMFKernel(false), Qtol(1e-12), nIndepIdgas(0), nDensities(0), polarizable(false)
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
	std::vector<double> Nmol0(component.size(), 0.), Phi_Nmol0(component.size());
	computeUniformEx(Nmol0, Phi_Nmol0);
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const FluidComponent& c = *component[ic];
		c.idealGas->Nbulk = Nmol[ic];
		c.idealGas->mu = Phi_Nmol[ic];
		logPrintf("   Component '%s' at bulk density %le bohr^-3 (with vapor pressure ~ %.2lg KPa and chemical potential %.2lg H)\n",
			c.molecule.name.c_str(), Nmol[ic], c.idealGas->Nbulk*T*exp((Phi_Nmol[ic]-Phi_Nmol0[ic])/T)/KPascal, Phi_Nmol[ic]-Phi_Nmol0[ic]);
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
	Cpol = chiPol ? (epsInf-1.)/(4.*M_PI*chiPol) : 1.;
	if(!Cpol) polarizable = false;
	logPrintf("   Local polarization-density correlation factors, Crot: %lg  Cpol: ", Crot);
	if(polarizable) logPrintf("%lg\n", Cpol); else logPrintf("none/disabled\n");
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
	ScalarFieldArray Vex(nDensities);
	{	//Get the uniform fluid site densities:
		std::vector<double> Nbulk(nDensities);
		for(const FluidComponent* c: component)
			for(unsigned i=0; i<c->molecule.sites.size(); i++)
				Nbulk[c->offsetDensity + i] = c->idealGas->get_Nbulk() * c->molecule.sites[i]->positions.size();
		ScalarFieldTildeArray Ntilde(nDensities);
		nullToZero(Ntilde, gInfo);
		for(unsigned i=0; i<nDensities; i++) Ntilde[i]->setGzero(Nbulk[i]);
		//Set Vex to the difference between the potentials returned by Fmix::compute() and Fmix::computeUniform():
		ScalarFieldTildeArray Phi_Ntilde(nDensities);
		std::vector<double> Phi_Nbulk(nDensities);
		for(const Fmix* fmix: fmixArr)
		{	fmix->compute(Ntilde, Phi_Ntilde);
			fmix->computeUniform(Nbulk, Phi_Nbulk);
		}
		for(unsigned i=0; i<nDensities; i++)
			if(Phi_Ntilde[i]) Vex[i] = Jdag(Phi_Ntilde[i]) - Phi_Nbulk[i];
		//Add electrostatic coupling:
		if(rhoExternal)
		{	ScalarFieldTilde dExternal = coulomb(rhoExternal);
			for(const FluidComponent* c: component)
				for(unsigned i=0; i<c->molecule.sites.size(); i++)
				{	const Molecule::Site& s = *(c->molecule.sites[i]);
				  //initialization from electrostatic potential doesn't work for charged species
				  if(s.chargeKernel && (!c->molecule.getCharge())) Vex[c->offsetDensity+i] += I(s.chargeKernel * dExternal);
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
	{	ScalarFieldTilde rhoTot; vector3<> Prot0; double epsInf=1.;
		for(unsigned ic=0; ic<component.size(); ic++)
		{	const FluidComponent& c = *(component[ic]);
			ScalarFieldArray N(c.molecule.sites.size()); vector3<> P0c;
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
		VectorFieldTilde epsMF = gradient((-1./epsInf)*coulomb(rhoTot));
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
{	if(mpiUtil->isHead()) saveToFile(state, filename);
}

FluidMixture::Outputs::Outputs(ScalarFieldArray* N, vector3<>* electricP,
	ScalarFieldTilde* Phi_rhoExternal, ScalarFieldArray* psiEff, EnergyComponents* Phi)
:N(N),electricP(electricP),Phi_rhoExternal(Phi_rhoExternal),psiEff(psiEff),Phi(Phi)
{
}

double FluidMixture::getFreeEnergy(Outputs outputs) const
{	ScalarFieldArray Phi_state;
	return (*this)(state, Phi_state, outputs);
	//calls operator() and provides whichever outputs are requested
	//before calling must create an Outputs structure with the required quantities allocated.
}

void FluidMixture::step(const ScalarFieldArray& dir, double alpha)
{	axpy(alpha, dir, state);
}

double FluidMixture::compute(ScalarFieldArray* grad, ScalarFieldArray* Kgrad)
{	ScalarFieldArray tempGrad;
	double E = (*this)(state, grad ? *grad : tempGrad, Outputs());
	//Compute preconditioned gradient:
	if(Kgrad)
	{	*Kgrad = clone(grad ? *grad : tempGrad);
		for(unsigned ic=0; ic<component.size(); ic++)
		{	const FluidComponent& c = *component[ic];
			for(unsigned k=c.offsetIndep; k<c.offsetIndep+c.idealGas->nIndep; k++)
				Kgrad->at(k) *= Kindep[ic];
		}
		for(unsigned k=nIndepIdgas; k<get_nIndep(); k++)
			Kgrad->at(k) *= Keps;
	}
	return E;
}

double FluidMixture::sync(double x) const
{	mpiUtil->bcast(x);
	return x;
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

//-------- Helper utilities for adjusting fluid mixture to vapor-liquid equilibrium

struct BoilingPressureSolver
{	const FluidMixture& fm;
	int nComponents;
	std::vector<double> xLiq;
	std::vector<double> Nliq, Nvap, aPrimeLiq, aPrimeVap; double Pliq, Pvap, NliqTot; //Dependent quantities set by compute:
	gsl_multiroot_fsolver* fsolver;
	gsl_multiroot_function solverFunc;
	
	BoilingPressureSolver(const FluidMixture& fm, std::vector<double> xLiq, double NliqGuess, double NvapGuess)
	: fm(fm), nComponents(fm.component.size()), xLiq(xLiq),
	Nliq(nComponents), Nvap(nComponents), aPrimeLiq(nComponents), aPrimeVap(nComponents)
	{
		//Initialize state:
		gsl_vector* params = gsl_vector_alloc(nComponents+1);
		for(int i=0; i<nComponents; i++)
			gsl_vector_set(params, i, log(NvapGuess * xLiq[i])); //Use liquid mole fractions as guess for vapor
		gsl_vector_set(params, nComponents, log(NliqGuess));
		//Create solver:
		fsolver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, nComponents+1);
		solverFunc.f = &BoilingPressureSolver::errFunc;
		solverFunc.n = nComponents+1;
		solverFunc.params = this;
		gsl_multiroot_fsolver_set(fsolver, &solverFunc, params);
		gsl_vector_free(params);
	}
	
	~BoilingPressureSolver()
	{	gsl_multiroot_fsolver_free(fsolver);
	}
	
	void solve(double tol)
	{	printState(0);
		for(int iter=1; iter<100; iter++)
		{	int status = gsl_multiroot_fsolver_iterate(fsolver);
			printState(iter);
			if(status) die("Boiling Pressure solver stuck - try different guesses for the densities.\n");
			if(gsl_multiroot_test_residual(fsolver->f, tol) != GSL_CONTINUE) break;
		}
		if(gsl_multiroot_test_residual(fsolver->f, tol) != GSL_SUCCESS)
			die("Boiling Pressure solver failed to converge - try different guesses for the densities.\n");
	}
	
	static int errFunc(const gsl_vector *params, void *bpSolver, gsl_vector* err)
	{	((BoilingPressureSolver*)bpSolver)->compute(params, err);
		return 0;
	}
	
	void compute(const gsl_vector* params, gsl_vector* err)
	{	//Set the densities:
		for(int i=0; i<nComponents; i++) Nvap[i] = exp(gsl_vector_get(params, i));
		NliqTot = exp(gsl_vector_get(params, nComponents));
		for(int i=0; i<nComponents; i++) Nliq[i] = NliqTot * xLiq[i];
		//Compute the bulk excess free energy density and gradient
		double aVap = fm.computeUniformEx(Nvap, aPrimeVap);
		double aLiq = fm.computeUniformEx(Nliq, aPrimeLiq);
		//Compute the differences in pressures and chemical potentials:
		Pvap = -aVap;
		Pliq = -aLiq;
		for(int i=0; i<nComponents; i++)
		{	Pvap += Nvap[i] * (fm.T + aPrimeVap[i]);
			Pliq += Nliq[i] * (fm.T + aPrimeLiq[i]);
			gsl_vector_set(err, i, log(Nliq[i]/Nvap[i]) + (aPrimeLiq[i]-aPrimeVap[i])/fm.T); //chemical potential difference / T (dimensionless)
		}
		gsl_vector_set(err, nComponents, (Pliq-Pvap)/(NliqTot*fm.T)); //dimensionless pressure difference
	}
	
	void printState(size_t iter)
	{	logPrintf("\tBPsolve: Iter: %lu  NliqTot: %.3le  Nvap: (", iter, exp(gsl_vector_get(fsolver->x, nComponents)));
		for(int i=0; i<nComponents; i++) logPrintf(" %.3le", exp(gsl_vector_get(fsolver->x, i)));
		logPrintf(")  DeltaP/NliqT: %.3le  DeltaMu/T: (", gsl_vector_get(fsolver->f, nComponents));
		for(int i=0; i<nComponents; i++) logPrintf(" %.3le", gsl_vector_get(fsolver->f, i));
		logPrintf(")\n");
	}
};

double FluidMixture::getBoilingPressure(double NliqGuess, double NvapGuess, std::vector<double>* Nvap) const
{	logPrintf("Finding vapor-liquid equilibrium state points:\n"); logFlush();
	//Collect and normalize mole fractions:
	double xTot = 0.; std::vector<double> xLiq(component.size());
	for(size_t ic=0; ic<component.size(); ic++)
	{	xLiq[ic] = component[ic]->Nbulk;
		xTot += xLiq[ic];
	}
	for(double& x: xLiq) x /= xTot;
	
	BoilingPressureSolver bpSolver(*this, xLiq, NliqGuess, NvapGuess);
	bpSolver.solve(1e-8);
	
	logPrintf("At equilibrium:\n\tPliq = %le bar, Pvap = %le bar\n", bpSolver.Pliq/Bar, bpSolver.Pvap/Bar);
	for(size_t ic=0; ic<component.size(); ic++)
	{	const FluidComponent& c = *component[ic];
		logPrintf("\tComponent '%s': Nliq = %le bohr^-3, Nvap = %le bohr^-3\n",
			c.molecule.name.c_str(), bpSolver.Nliq[ic], bpSolver.Nvap[ic]);
	}
	if(Nvap) *Nvap = bpSolver.Nvap;
	return bpSolver.Pliq;
}
