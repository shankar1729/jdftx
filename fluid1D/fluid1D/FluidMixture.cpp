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

#include <fluid1D/FluidMixture.h>
#include <fluid1D/MixedFMT.h>
#include <core/EnergyComponents.h>
#include <core/BlasExtra.h>
#include <core/Random.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

FluidMixture::FluidMixture(const GridInfo& gInfo, const double T)
:gInfo(gInfo),T(T),verboseLog(false),nDensities(0),nIndepIdgas(0),polarizable(false)
{	logPrintf("Initializing fluid mixture at T=%lf K ...\n", T/Kelvin);
}

void FluidMixture::setPressure(double p, double Nguess)
{	logPrintf("Adjusting fluid pressure to p=%lf bar\n", p/Bar);
	//Normalize mole fractions:
	double xBulkTot=0.0;
	for(std::vector<Component>::iterator c=component.begin(); c!=component.end(); c++) xBulkTot += c->idealGas->xBulk;
	for(std::vector<Component>::iterator c=component.begin(); c!=component.end(); c++) c->idealGas->xBulk /= xBulkTot;
	//Compute the maximum possible density (core packed limit)
	double zi3=0.0; //n3/N for the bulk mixture
	for(std::vector<Component>::const_iterator c=component.begin(); c!=component.end(); c++)
	{	for(int j=0; j<c->molecule->nSites; j++)
		{	const SiteProperties& s = *(c->molecule->site[j].prop);
			if(s.sphereRadius) zi3 += s.w3->at(0) * c->idealGas->xBulk;
		}
	}
	double Nmax = 1.0/zi3;
	const double mulStep = 0.99;
	if(Nguess > mulStep*Nmax) Nguess = mulStep*Nmax;
	double pTest = compute_p(Nguess);
	//Find an interval of N that brackets P:
	double Nlo, Nhi;
	if(pTest > p)
	{	Nlo = Nguess;
		do
		{	Nhi = Nlo;
			Nlo = mulStep*Nhi;
			pTest = compute_p(Nlo);
		}
		while(pTest>p);
	}
	else
	{	Nhi = Nguess;
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
	std::vector<double> Nmol(component.size()), grad_Nmol(component.size());
	for(unsigned ic=0; ic<component.size(); ic++) Nmol[ic] = Ntot*component[ic].idealGas->xBulk;
	computeUniformEx(Nmol, grad_Nmol);
	for(unsigned ic=0; ic<component.size(); ic++)
	{	Component& c = component[ic];
		c.idealGas->Nbulk = Nmol[ic];
		c.idealGas->mu = grad_Nmol[ic];
		logPrintf("\tComponent '%s' at bulk density %le bohr^-3\n", c.molecule->name.c_str(), Nmol[ic]);
	}
	this->p = p;
}

//-------- Helper utilities for adjusting fluid mixture to vapor-liquid equilibrium

struct BoilingPressureSolver
{	const FluidMixture& fm;
	int nComponents;
	std::vector<double> xLiq;
	std::vector<double> Nliq, Nvap, aPrimeLiq, aPrimeVap; double Pliq, Pvap, NliqTot; //Dependent quantities set by compute:
	gsl_multiroot_fsolver* fsolver;
		
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
		gsl_multiroot_function solverFunc = {&BoilingPressureSolver::errFunc, nComponents+1, this};
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

double FluidMixture::setBoilingPressure(std::vector<double>* Nvap, double tol, double NliqGuess, double NvapGuess)
{	logPrintf("Finding vapor-liquid equilibrium state points:\n"); logFlush();
	//Collect and normalize mole fractions:
	double xTot = 0.; std::vector<double> xLiq(component.size());
	for(size_t i=0; i<component.size(); i++)
	{	xLiq[i] = component[i].idealGas->xBulk;
		xTot += xLiq[i];
	}
	for(size_t i=0; i<component.size(); i++) xLiq[i] /= xTot;
	
	BoilingPressureSolver bpSolver(*this, xLiq, NliqGuess, NvapGuess);
	bpSolver.solve(tol);
	
	logPrintf("At equilibrium:\n\tPliq = %le bar, Pvap = %le bar\n", bpSolver.Pliq/Bar, bpSolver.Pvap/Bar);
	for(size_t ic=0; ic<component.size(); ic++)
	{	Component& c = component[ic];
		logPrintf("\tComponent '%s': Nliq = %le bohr^-3, Nvap = %le bohr^-3\n",
			c.molecule->name.c_str(), bpSolver.Nliq[ic], bpSolver.Nvap[ic]);
	}
	setPressure(bpSolver.Pliq, bpSolver.NliqTot); //set the pressure
	if(Nvap) *Nvap = bpSolver.Nvap;
	return bpSolver.Pliq;
}


unsigned FluidMixture::get_offsetDensity(const Fex* fex) const
{	for(std::vector<Component>::const_iterator c=component.begin(); c!=component.end(); c++)
		if(fex==c->fex)
			return c->offsetDensity;
	assert(!"Could not find excess functional in fluid mixture."); return 0;
}

unsigned int FluidMixture::get_nComponents() const
{	return component.size();
}

const FluidMixture::Component& FluidMixture::get_component(unsigned int c) const
{	return component[c];
}


void FluidMixture::addComponent(IdealGas* idealGas, const Fex* fex)
{	Component comp = {idealGas, fex, fex->getMolecule(), nIndepIdgas, nDensities};
	comp.indexedSite.resize(comp.molecule->nIndices);
	comp.indexedSiteMultiplicity.resize(comp.molecule->nIndices);
	for(int i=0; i<comp.molecule->nSites; i++)
	{	comp.indexedSite[comp.molecule->site[i].index] = comp.molecule->site[i].prop;
		comp.indexedSiteMultiplicity[comp.molecule->site[i].index]++;
		polarizable |= (comp.molecule->site[i].prop->alpha && comp.molecule->site[i].prop->alphaKernel);
	}
	component.push_back(comp);
	//Update the totals, which become the offset for the next component
	nIndepIdgas += idealGas->nIndep;
	nDensities += comp.molecule->nIndices;
}

void FluidMixture::addFmix(const Fmix* fmix)
{	fmixArr.push_back(fmix);
}


void FluidMixture::initState(double scale, double Elo, double Ehi)
{	//Compute the effective nonlinear coupling potential for the uniform fluid:
	ScalarFieldCollection Vex(nDensities); //TODO: Interface with the electronic side to get the coupling potential here
	//Call initState for each component
	nullToZero(state, gInfo, get_nIndep());
	logPrintf("\n----- FluidMixture::initState() -----\n");
	for(std::vector<Component>::iterator c=component.begin(); c!=component.end(); c++)
		c->idealGas->initState(&Vex[c->offsetDensity], &state[c->offsetIndep], scale, Elo, Ehi);
	logPrintf("\n");
}

void FluidMixture::loadState(const char* filename)
{	nullToZero(state, gInfo, get_nIndep());
	loadFromFile(state, filename);
}

void FluidMixture::saveState(const char* filename) const
{	saveToFile(state, filename);
}

FluidMixture::Outputs::Outputs(ScalarFieldCollection* N, double* electricP, ScalarFieldCollection* psiEff)
:N(N),electricP(electricP),psiEff(psiEff)
{
}

double FluidMixture::operator()(const ScalarFieldCollection& indep, ScalarFieldCollection& grad_indep, Outputs outputs) const
{	
	//logPrintf("indep.size: %d nIndep: %d\n",indep.size(),nIndep);
	assert(indep.size()==get_nIndep());

	//---------- Compute site densities from the independent variables ---------
	std::vector<double> P(component.size()); //total dipole moment per component
	ScalarFieldTildeCollection Ntilde(nDensities); //site densities (in reciprocal space)
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const Component& c = component[ic];
		ScalarFieldCollection N(c.molecule->nIndices);
		c.idealGas->getDensities(&indep[c.offsetIndep], &N[0], P[ic]);
		for(int j=0; j<c.molecule->nIndices; j++)
		{	//Replace negative densities with 0:
			double Nmin, Nmax;
			eblas_capMinMax(gInfo.S, N[j].data(), Nmin, Nmax, 0.);
			//store site densities in fourier space
			Ntilde[c.offsetDensity+j] = J(N[j]);
			N[j] = 0; //Really skimp on memory!
		}
	}

	//----------- Handle density constraints ------------
	std::vector<double> Nscale(component.size(), 1.0); //density scale factor that satisfies the constraint
	std::vector<std::vector<double> > NscalePrime(component.size(),
		std::vector<double>(component.size(),0.0)); //jacobian of Nscale w.r.t the uncorrected molecule counts
	std::vector<string> names; //list of molecule names
	
	//Find fixed N and charged species:
	for(unsigned ic=0; ic<component.size(); ic++)
	{
		const Component& c = component[ic];
		double Qmolecule = c.molecule->get_charge();
		if(c.idealGas->Nnorm>0 || Qmolecule)
		{
			double N0 = integral(Ntilde[c.offsetDensity])/c.indexedSiteMultiplicity[0];
			if(c.idealGas->Nnorm>0)
			{
				Nscale[ic] = c.idealGas->Nnorm/N0;
				NscalePrime[ic][ic] = -c.idealGas->Nnorm/pow(N0,2);
			}
		}
	}
	std::vector<double> grad_Nscale(component.size(), 0.0); //accumulate explicit derivatives w.r.t Nscale here

	//Apply the scale factors to the site densities
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const Component& c = component[ic];
		for(int j=0; j<c.molecule->nIndices; j++)
			Ntilde[c.offsetDensity+j] *= Nscale[ic];
		P[ic] *= Nscale[ic];
	}

	EnergyComponents Phi; //the grand free energy (with component information)
	ScalarFieldTildeCollection grad_Ntilde(nDensities); //gradients (partial derivative) w.r.t reciprocal space site densities
	std::vector<double> grad_P(component.size()); //gradient (partial derivative) w.r.t cell dipole

	//--------- Compute the (scaled) mean field coulomb interaction --------
	{	ScalarFieldTilde Od; //total electrostatic potential (with a factor of O)
		ScalarFieldTilde rho; //total charge density
		
		double electricPtot = 0.; //total electric dipole moment in cell
		double geomFacPtot = (gInfo.coord==GridInfo::Planar ? 1./gInfo.rMax : 0.); //geometrical factor in dipole-sum contribution
		
		for(unsigned ic=0; ic<component.size(); ic++)
		{	const Component& c = component[ic];
			ScalarFieldTilde rho_c; //total charge from this fluid component
			for(int j=0; j<c.molecule->nIndices; j++)
			{	const SiteProperties& s = *c.indexedSite[j];
				if(s.chargeZ && s.chargeKernel)
					rho_c += s.chargeZ * (*s.chargeKernel * Ntilde[c.offsetDensity+j]);
			}
			if(!rho_c) continue;
			ScalarFieldTilde Od_c = O(-4*M_PI*Linv(O(rho_c))); //total electrostatic potential from this fluid component
			Od += Od_c;
			rho += rho_c;
			electricPtot += c.molecule->get_dipole()*P[ic];

			//Add the extra correlation contribution within component
			double corrFac = c.fex->get_aDiel() - 1;
			ScalarFieldTilde Od_cCorr = corrFac * Od_c;
			Phi["Coulomb"] += 0.5*dot(rho_c, Od_cCorr);
			for(int j=0; j<c.molecule->nIndices; j++)
			{	const SiteProperties& s = *c.indexedSite[j];
				if(s.chargeZ && s.chargeKernel)
					grad_Ntilde[c.offsetDensity+j] += s.chargeZ * (*s.chargeKernel * Od_cCorr);
			}
			Phi["PsqCell"] += 0.5 * 4*M_PI*corrFac * pow(P[ic],2) * pow(c.molecule->get_dipole(),2) * geomFacPtot;
			grad_P[ic] = 4*M_PI*corrFac * P[ic] * pow(c.molecule->get_dipole(),2) * geomFacPtot;
		}
		//Now add the true mean-field contribution from Od to all the site densities:
		if(rho)
		{	Phi["Coulomb"] += 0.5*dot(rho,Od);
			Phi["PsqCell"] += 0.5 * 4*M_PI * pow(electricPtot,2) * geomFacPtot;
			for(unsigned ic=0; ic<component.size(); ic++)
			{	const Component& c = component[ic];
				for(int j=0; j<c.molecule->nIndices; j++)
				{	const SiteProperties& s = *c.indexedSite[j];
					if(s.chargeZ && s.chargeKernel)
						grad_Ntilde[c.offsetDensity+j] += s.chargeZ * (*s.chargeKernel * Od);
				}
				grad_P[ic] += 4*M_PI * electricPtot * c.molecule->get_dipole() * geomFacPtot;
			}
		}
		if(outputs.electricP) *outputs.electricP = electricPtot;
	}

	//--------- Hard sphere mixture and bonding -------------
	{	//Compute the FMT weighted densities:
		ScalarFieldTilde n0tilde, n1tilde, n2tilde, n3tilde, n1vTilde, n2mTilde;
		std::vector<ScalarField> n0mol(component.size(), 0); //partial n0 for molecules that need bonding corrections
		std::vector<int> n0mult(component.size(), 0); //number of sites which contribute to n0 for each molecule
		std::vector<std::map<double,int> > bond(component.size()); //sets of bonds for each molecule
		bool bondsPresent = false; //whether bonds are present for any molecule
		for(unsigned ic=0; ic<component.size(); ic++)
		{	const Component& c = component[ic];
			bond[ic] = c.molecule->getBonds();
			ScalarFieldTilde n0molTilde;
			for(int j=0; j<c.molecule->nIndices; j++)
			{	const SiteProperties& s = *(c.indexedSite[j]);
				if(s.sphereRadius)
				{	const ScalarFieldTilde& Nsite = Ntilde[c.offsetDensity+j];
					n0mult[ic] += c.indexedSiteMultiplicity[j];
					n0molTilde += *s.w0  * Nsite;
					n1tilde    += *s.w1  * Nsite;
					n2tilde    += *s.w2  * Nsite;
					n3tilde    += *s.w3  * Nsite;
					n1vTilde   += *s.w1v * Nsite;
					n2mTilde   += *s.w2m * Nsite;
				}
			}
			if(n0molTilde) n0tilde += n0molTilde;
			if(bond[ic].size())
			{	n0mol[ic] = I(n0molTilde);
				bondsPresent = true;
			}
		}
		if(n0tilde) //at least one sphere in the mixture
		{	ScalarField n0 = I(n0tilde); n0tilde=0;
			ScalarField n1 = I(n1tilde); n1tilde=0;
			ScalarField n2 = I(n2tilde); n2tilde=0;
			ScalarField grad_n0, grad_n1, grad_n2; ScalarFieldTilde grad_n3tilde, grad_n1vTilde, grad_n2mTilde;
			//Compute the sphere mixture free energy:
			Phi["MixedFMT"] += T * PhiFMT(n0, n1, n2, n3tilde, n1vTilde, n2mTilde,
				grad_n0, grad_n1, grad_n2, grad_n3tilde, grad_n1vTilde, grad_n2mTilde);
			//Bonding corrections if required
			if(bondsPresent)
			{	for(unsigned ic=0; ic<component.size(); ic++)
				{	const Component& c = component[ic];
					ScalarField grad_n0mol;
					for(std::map<double,int>::iterator b=bond[ic].begin(); b!=bond[ic].end(); b++)
						Phi["Bonding"] += T * PhiBond(b->first, b->second*1.0/n0mult[ic],
							n0mol[ic], n2, n3tilde, grad_n0mol, grad_n2, grad_n3tilde);
					if(grad_n0mol)
					{	//Propagate gradient w.r.t n0mol[ic] to the site densities:
						ScalarFieldTilde grad_n0molTilde = Idag(grad_n0mol);
						for(int j=0; j<c.molecule->nIndices; j++)
						{	const SiteProperties& s = *(c.indexedSite[j]);
							if(s.sphereRadius)
								grad_Ntilde[c.offsetDensity+j] += T * (*s.w0  * grad_n0molTilde);
						}
					}
				}
			}
			//Accumulate gradients w.r.t weighted densities to site densities:
			ScalarFieldTilde grad_n0tilde = Idag(grad_n0); grad_n0=0;
			ScalarFieldTilde grad_n1tilde = Idag(grad_n1); grad_n1=0;
			ScalarFieldTilde grad_n2tilde = Idag(grad_n2); grad_n2=0;
			for(std::vector<Component>::const_iterator c=component.begin(); c!=component.end(); c++)
			{	for(int j=0; j<c->molecule->nIndices; j++)
				{	const SiteProperties& s = *(c->indexedSite[j]);
					if(s.sphereRadius)
					{	ScalarFieldTilde& grad_Nsite = grad_Ntilde[c->offsetDensity+j];
						grad_Nsite += T * (*s.w0  * grad_n0tilde);
						grad_Nsite += T * (*s.w1  * grad_n1tilde);
						grad_Nsite += T * (*s.w2  * grad_n2tilde);
						grad_Nsite += T * (*s.w3  * grad_n3tilde);
						grad_Nsite += T * (*s.w1v * grad_n1vTilde);
						grad_Nsite += T * (*s.w2m * grad_n2mTilde);
					}
				}
			}
		}
	}

	//---------- Excess functionals --------------
	for(std::vector<Component>::const_iterator c=component.begin(); c!=component.end(); c++)
		Phi["Fex("+c->molecule->name+")"] += c->fex->compute(&Ntilde[c->offsetDensity], &grad_Ntilde[c->offsetDensity]);

	//--------- Mixing functionals --------------
	for(std::vector<const Fmix*>::const_iterator fmix=fmixArr.begin(); fmix!=fmixArr.end(); fmix++)
		Phi["Fmix("+(*fmix)->getName()+")"] += (*fmix)->compute(Ntilde, grad_Ntilde);

	//--------- PhiNI ---------
	nullToZero(grad_Ntilde, gInfo);
	if(outputs.N) outputs.N->resize(nDensities);
	//Put the site densities and gradients back in real space
	ScalarFieldCollection N(nDensities);
	ScalarFieldCollection grad_N(nDensities);
	for(unsigned i=0; i<nDensities; i++)
	{	N[i] = I(Ntilde[i]); Ntilde[i]=0;
		grad_N[i] = Jdag(grad_Ntilde[i]); grad_Ntilde[i] = 0;
		if(outputs.N) (*outputs.N)[i] = N[i]; //Link site-density to return pointer if necessary
	}
	//Estimate psiEff based on gradients, if requested
	if(outputs.psiEff)
	{	outputs.psiEff->clear();
		for(const Component& c: component)
		{	bool muAdded = false;
			int curSiteMult=0; //multiplicity of current site
			for(int i=0; i<c.molecule->nSites; i++)
			{	curSiteMult++;
				if(!(c.molecule->site[i].prop->indepSite)) continue;
				if((i+1!=c.molecule->nSites)
					&& (c.molecule->site[i].index == c.molecule->site[i+1].index))
					continue; //this is not the last of a chain of identical sites
				int iDensity = c.molecule->site[i].index;
				ScalarField psiCur = grad_N[c.offsetDensity+iDensity] + c.idealGas->V[iDensity];
				if(!muAdded)
				{	psiCur -= c.idealGas->mu/curSiteMult;
					muAdded = true;
				}
				psiCur *= (-1./T);
				outputs.psiEff->push_back(psiCur);
				curSiteMult=0; //will be starting a new site after this
			}
		}
	}
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const Component& c = component[ic];
		Phi["PhiNI("+c.molecule->name+")"] +=
			c.idealGas->compute(&indep[c.offsetIndep], &N[c.offsetDensity], &grad_N[c.offsetDensity],
				P[ic], grad_P[ic], Nscale[ic], grad_Nscale[ic]);

		//Fixed N correction to entropy:
		if(Nscale[ic]!=1.0)
		{	double deltaTs = T*log(Nscale[ic])/c.indexedSiteMultiplicity[0];
			grad_N[c.offsetDensity] += deltaTs;
			Phi["PhiNI("+c.molecule->name+")"] += integral(N[c.offsetDensity])*deltaTs;
		}
	}
	//Add in the implicit contributions to grad_Nscale
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const Component& ci = component[ic];
		bool anyNonzero=false;
		for(unsigned jc=0; jc<component.size(); jc++)
			if(NscalePrime[ic][jc])
				anyNonzero=true;
		if(anyNonzero)
		{	grad_Nscale[ic] += P[ic] * grad_P[ic];
			for(int ki=0; ki<ci.molecule->nIndices; ki++)
				grad_Nscale[ic] += dot(N[ci.offsetDensity+ki], grad_N[ci.offsetDensity+ki]) / Nscale[ic];
		}
	}
	//Propagate gradients from Nscale to N:
	for(unsigned jc=0; jc<component.size(); jc++)
	{	const Component& cj = component[jc];
		double grad_Ncontrib = 0.0;
		for(unsigned ic=0; ic<component.size(); ic++)
			if(NscalePrime[ic][jc])
				grad_Ncontrib += grad_Nscale[ic] * NscalePrime[ic][jc];
		if(grad_Ncontrib)
			grad_N[cj.offsetDensity] += DiagJdagOJ1(grad_Ncontrib / (Nscale[jc] * cj.indexedSiteMultiplicity[0]), gInfo);
	}

	//Propagate gradients from grad_N to grad_indep
	grad_indep.resize(get_nIndep());
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const Component& c = component[ic];
		c.idealGas->convertGradients(&indep[c.offsetIndep], &N[c.offsetDensity],
			&grad_N[c.offsetDensity], grad_P[ic], &grad_indep[c.offsetIndep], Nscale[ic]);
	}
	if(polarizable) nullToZero(grad_indep[nIndepIdgas], gInfo);
	
	Phi["+pV"] += p * gInfo.Volume(); //background correction

	if(verboseLog) Phi.print(globalLog, true, "\t\t\t\t%15s = %25.16lf\n");
	return Phi;
}

double FluidMixture::getFreeEnergy(Outputs outputs) const
{	ScalarFieldCollection grad_state;
	return (*this)(state, grad_state, outputs);
	//calls operator() and provides whichever outputs are requested
	//before calling must create an Outputs structure with the required quantities allocated.
}

void FluidMixture::step(const ScalarFieldCollection& dir, double alpha)
{	axpy(alpha, dir, state);
}

double FluidMixture::compute(ScalarFieldCollection* grad)
{	ScalarFieldCollection tempGrad;
	return (*this)(state, grad ? *grad : tempGrad, Outputs());
}

ScalarFieldCollection FluidMixture::precondition(const ScalarFieldCollection& grad)
{	ScalarFieldCollection Kgrad(grad.size());
	for(unsigned i=0; i<grad.size(); i++)
		Kgrad[i] = DiagJdagOJ1inv(grad[i]);
    return Kgrad;
}

double FluidMixture::computeUniformEx(const std::vector<double>& Nmol, std::vector<double>& grad_Nmol) const
{	//--------- Compute the site densities ----------
	std::vector<double> N(nDensities);
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const Component& c = component[ic];
		for(int j=0; j<c.molecule->nSites; j++)
			N[c.offsetDensity + c.molecule->site[j].index] += Nmol[ic];
	}

	EnergyComponents phi; std::vector<double> grad_N(nDensities);

	//-------- Mean field coulomb -------- (No contribution in uniform fluid)

	//-------- Hard sphere/bonding -----------
	double n0=0.0, n1=0.0, n2=0.0, n3=0.0;
	std::vector<double> n0mol(component.size(), 0.0); //partial n0 for molecules that need bonding corrections
	std::vector<int> n0mult(component.size(), 0); //number of sites which contribute to n0 for each molecule
	std::vector<std::map<double,int> > bond(component.size()); //sets of bonds for each molecule
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const Component& c = component[ic];
		bond[ic] = c.molecule->getBonds();
		for(int j=0; j<c.molecule->nIndices; j++)
		{	const SiteProperties& s = *(c.indexedSite[j]);
			if(s.sphereRadius)
			{	double Nsite = N[c.offsetDensity+j];
				n0mult[ic] += c.indexedSiteMultiplicity[j];
				n0mol[ic]  += s.w0->data()[0]  * Nsite;
				n1         += s.w1->data()[0]  * Nsite;
				n2         += s.w2->data()[0]  * Nsite;
				n3         += s.w3->data()[0]  * Nsite;
			}
		}
		n0 += n0mol[ic];
	}
	if(n0)
	{	double grad_n0=0.0, grad_n1=0.0, grad_n2=0.0, grad_n3=0.0;
		phi["MixedFMT"] = T * phiFMTuniform(n0, n1, n2, n3, grad_n0, grad_n1, grad_n2, grad_n3);
		//Bonding corrections:
		for(unsigned ic=0; ic<component.size(); ic++)
		{	const Component& c = component[ic];
			double grad_n0mol=0.0;
			for(std::map<double,int>::iterator b=bond[ic].begin(); b!=bond[ic].end(); b++)
				phi["Bonding"] += T * phiBondUniform(b->first, b->second*1.0/n0mult[ic],
					n0mol[ic], n2, n3, grad_n0mol, grad_n2, grad_n3);
			if(grad_n0mol)
			{	//Propagate gradient w.r.t n0mol[ic] to the site densities:
				for(int j=0; j<c.molecule->nIndices; j++)
				{	const SiteProperties& s = *(c.indexedSite[j]);
					if(s.sphereRadius)
						grad_N[c.offsetDensity+j] += T * (s.w0->data()[0]  * grad_n0mol);
				}
			}
		}
		//Convert FMT weighted gradients to site density gradients:
		for(std::vector<Component>::const_iterator c=component.begin(); c!=component.end(); c++)
		{	for(int j=0; j<c->molecule->nIndices; j++)
			{	const SiteProperties& s = *(c->indexedSite[j]);
				if(s.sphereRadius)
				{	double& grad_Nsite = grad_N[c->offsetDensity+j];
					grad_Nsite += T * (s.w0->data()[0]  * grad_n0);
					grad_Nsite += T * (s.w1->data()[0]  * grad_n1);
					grad_Nsite += T * (s.w2->data()[0]  * grad_n2);
					grad_Nsite += T * (s.w3->data()[0]  * grad_n3);
				}
			}
		}
	}

	//---------- Excess functionals --------------
	for(std::vector<Component>::const_iterator c=component.begin(); c!=component.end(); c++)
		phi["Fex("+c->molecule->name+")"] += c->fex->computeUniform(&N[c->offsetDensity], &grad_N[c->offsetDensity]);

	//--------- Mixing functionals --------------
	for(std::vector<const Fmix*>::const_iterator fmix=fmixArr.begin(); fmix!=fmixArr.end(); fmix++)
		phi["Fmix("+(*fmix)->getName()+")"] += (*fmix)->computeUniform(N, grad_N);

	//--------- Convert site density gradients to molecular density gradients ---------
	for(unsigned ic=0; ic<component.size(); ic++)
	{	const Component& c = component[ic];
		grad_Nmol[ic] = 0.0;
		for(int j=0; j<c.molecule->nSites; j++)
			grad_Nmol[ic] += grad_N[c.offsetDensity + c.molecule->site[j].index];
	}

	//phi.print(gInfo.fpLog, true, "\t\t\t\t%15s = %12.4le\n"); //Uncomment when debugging to get contributions

	return phi;
}

double FluidMixture::compute_p(double Ntot) const
{	std::vector<double> Nmol(component.size()), grad_Nmol(component.size());
	for(unsigned ic=0; ic<component.size(); ic++) Nmol[ic] = Ntot*component[ic].idealGas->xBulk;
	//p = N T + Sum_ic Nic d(aEx)/dNic - aEx where aEx is Helmholtz energy density
	double p = Ntot*T - computeUniformEx(Nmol, grad_Nmol);
	for(unsigned ic=0; ic<component.size(); ic++)
		p += Nmol[ic]*grad_Nmol[ic];
	return p;
}
