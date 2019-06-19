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

#include <electronic/Everything.h>
#include <electronic/VanDerWaals.h>
#include <core/ScalarFieldIO.h>
#include <core/VectorField.h>
#include <core/SphericalHarmonics.h>
#include <fluid/SaLSA.h>
#include <fluid/PCM_internal.h>
#include <gsl/gsl_linalg.h>
#include <cstring>

struct MultipoleResponse
{	int l; //!< angular momentum
	int iSite; //!< site index (-1 => molecule center, -2 => ion center)
	int siteMultiplicity; //!< number of sites sharing index iSite
	
	RadialFunctionG V; //!< radial part of response eigenfunctions scaled by square root of eigenvalue and 1/G^l
	
	MultipoleResponse(int l, int iSite, int siteMultiplicity, const std::vector<double>& Vsamples, double dG)
	: l(l), iSite(iSite), siteMultiplicity(siteMultiplicity)
	{	V.init(0, Vsamples, dG);
	}
	
	~MultipoleResponse()
	{	V.free();
	}
	
	//index into shape or site arrays, as appropriate for iSite
	template<typename T> T& selectSite(std::vector<T>& shape, std::vector<T>& site) const
	{	return  iSite<0 ? shape[std::min(-iSite, int(shape.size())) -1] : site[iSite];
	}
	template<typename T> const T& selectSite(const std::vector<T>& shape, const std::vector<T>& site) const
	{	return  iSite<0 ? shape[std::min(-iSite, int(shape.size())) -1] : site[iSite];
	}
};


SaLSA::SaLSA(const Everything& e, const FluidSolverParams& fsp)
: PCM(e, fsp), siteShape(fsp.solvents[0]->molecule.sites.size())
{	
	assert(!useGummel()); //Non-variational energy: cannot use Gummel loop!
	logPrintf("   Initializing non-local response weight functions:\n");
	const double dG = gInfo.dGradial, Gmax = gInfo.GmaxGrid;
	unsigned nGradial = unsigned(ceil(Gmax/dG))+5;

	//Initialize fluid molecule's spherically-averaged electron density kernel:
	const auto& solvent = fsp.solvents[0];
	std::vector<double> nFluidSamples(nGradial);
	for(unsigned i=0; i<nGradial; i++)
	{	double G = i*dG;
		nFluidSamples[i] = 0.;
		for(const auto& site: solvent->molecule.sites)
		{	double nTilde = site->elecKernel(G);
			for(const vector3<>& r: site->positions)
				nFluidSamples[i] += nTilde * bessel_jl(0, G*r.length());
		}
	}
	nFluid.init(0, nFluidSamples, dG);
	
	//Determine dipole correlation factors:
	double chiRot = 0., chiPol = 0.;
	for(const auto& c: fsp.components)
	{	chiRot += c->Nbulk * c->molecule.getDipole().length_squared()/(3.*fsp.T);
		chiPol += c->Nbulk * c->molecule.getAlphaTot();
	}
	double sqrtCrot = (epsBulk>epsInf && chiRot) ? sqrt((epsBulk-epsInf)/(4.*M_PI*chiRot)) : 1.;
	double epsInfEff = chiRot ? epsInf : epsBulk; //constrain to epsBulk for molecules with no rotational susceptibility
	double sqrtCpol = (epsInfEff>1. && chiPol) ? sqrt((epsInfEff-1.)/(4.*M_PI*chiPol)) : 1.;
	
	//Rotational and translational response (includes ionic response):
	const double bessel_jl_by_Gl_zero[4] = {1., 1./3, 1./15, 1./105}; //G->0 limit of j_l(G)/G^l
	for(unsigned iComp=0; iComp<fsp.components.size(); iComp++)
	{	const std::shared_ptr<FluidComponent>& c = fsp.components[iComp];
		for(int l=0; l<=fsp.lMax; l++)
		{	//Calculate radial densities for all m:
			gsl_matrix* V = gsl_matrix_calloc(nGradial, 2*l+1); //allocate and set to zero
			double prefac = sqrt(4.*M_PI*c->Nbulk/fsp.T);
			for(unsigned iG=0; iG<nGradial; iG++)
			{	double G = iG*dG;
				for(const auto& site: c->molecule.sites)
				{	double Vsite = prefac * site->chargeKernel(G);
					for(const vector3<>& r: site->positions)
					{	double rLength = r.length();
						double bessel_jl_by_Gl = G ? bessel_jl(l,G*rLength)/pow(G,l) : bessel_jl_by_Gl_zero[l]*pow(rLength,l);
						vector3<> rHat = (rLength ? 1./rLength : 0.) * r;
						for(int m=-l; m<=+l; m++)
							*gsl_matrix_ptr(V,iG,l+m) += Vsite * bessel_jl_by_Gl * Ylm(l,m, rHat);
					}
				}
			}
			//Scale dipole active modes:
			for(int lm=0; lm<2l+1; lm++)
				if(l==1 && fabs(gsl_matrix_get(V,0,lm))>1e-6)
					for(unsigned iG=0; iG<nGradial; iG++)
						*gsl_matrix_ptr(V,iG,lm) *= sqrtCrot;
			//Get linearly-independent non-zero modes by performing an SVD:
			gsl_vector* S = gsl_vector_alloc(2*l+1);
			gsl_matrix* U = gsl_matrix_alloc(2*l+1, 2*l+1);
			gsl_matrix* tmpMat = gsl_matrix_alloc(2*l+1, 2*l+1);
			gsl_vector* tmpVec = gsl_vector_alloc(2*l+1);
			gsl_linalg_SV_decomp_mod(V, tmpMat, U, S, tmpVec);
			gsl_vector_free(tmpVec);
			gsl_matrix_free(tmpMat);
			gsl_matrix_free(U);
			//Add response functions for non-singular modes:
			for(int mode=0; mode<2*l+1; mode++)
			{	double Smode = gsl_vector_get(S, mode);
				if(Smode*Smode < 1e-3) break;
				std::vector<double> Vsamples(nGradial);
				for(unsigned iG=0; iG<nGradial; iG++)
					Vsamples[iG] = Smode * gsl_matrix_get(V, iG, mode);
				response.push_back(std::make_shared<MultipoleResponse>(l, -(iComp+1), 1, Vsamples, dG));
			}
			gsl_vector_free(S);
			gsl_matrix_free(V);
		}
	}
	
	//Polarizability response:
	for(unsigned iSite=0; iSite<solvent->molecule.sites.size(); iSite++)
	{	const Molecule::Site& site = *(solvent->molecule.sites[iSite]);
		if(site.polKernel)
		{	std::vector<double> Vsamples(nGradial);
			double prefac = sqrtCpol * sqrt(solvent->Nbulk * site.alpha);
			for(unsigned iG=0; iG<nGradial; iG++)
				Vsamples[iG] = prefac * site.polKernel(iG*dG);
			response.push_back(std::make_shared<MultipoleResponse>(1, iSite, site.positions.size(), Vsamples, dG));
		}
	}
	
	const double GzeroTol = 1e-12;
	
	//Compute bulk properties and print summary:
	double epsBulk = 1.; double k2factor = 0.; std::map<int,int> lCount;
	for(const std::shared_ptr<MultipoleResponse>& resp: response)
	{	lCount[resp->l]++;
		double respGzero = (4*M_PI) * pow(resp->V(0), 2) * resp->siteMultiplicity;
		if(resp->l==0) k2factor += respGzero;
		if(resp->l==1) epsBulk += respGzero;
	}
	for(auto lInfo: lCount)
		logPrintf("      l: %d  #weight-functions: %d\n", lInfo.first, lInfo.second);
	logPrintf("   Bulk dielectric-constant: %lg", epsBulk);
	if(k2factor > GzeroTol) logPrintf("   screening-length: %lg bohrs.\n", sqrt(epsBulk/k2factor));
	else logPrintf("\n");
	if(fsp.lMax >= 1) assert(fabs(epsBulk-this->epsBulk) < 1e-3); //verify consistency of correlation factors
	assert(fabs(k2factor-this->k2factor) < 1e-3); //verify consistency of site charges
	
	//Initialize preconditioner kernel:
	std::vector<double> KkernelSamples(nGradial);
	for(unsigned i=0; i<nGradial; i++)
	{	double G = i*dG, G2=G*G;
		//Compute diagonal part of the hessian ( 4pi(Vc^-1 + chi) ):
		double diagH = G2;
		for(const auto& resp: response)
			diagH += pow(G2,resp->l) * pow(resp->V(G), 2);
		//Set its inverse square-root as the preconditioner:
		KkernelSamples[i] = (diagH>GzeroTol) ? 1./sqrt(diagH) : 0.;
	}
	Kkernel.init(0, KkernelSamples, dG);
	
	//MPI division:
	TaskDivision(response.size(), mpiWorld).myRange(rStart, rStop);
}

SaLSA::~SaLSA()
{	nFluid.free();
	Kkernel.free();
}


ScalarFieldTilde SaLSA::chi(const ScalarFieldTilde& phiTilde) const
{	ScalarFieldTilde rhoTilde;
	for(int r=rStart; r<rStop; r++)
	{	const MultipoleResponse& resp = *response[r];
		const ScalarField& s = resp.selectSite(shape, siteShape);
		if(resp.l>6) die("Angular momenta l > 6 not supported.\n");
		double prefac = pow(-1,resp.l) * 4*M_PI/(2*resp.l+1);
		rhoTilde -= prefac * (resp.V * lDivergence(J(s * I(lGradient(resp.V * phiTilde, resp.l))), resp.l));
	}
	nullToZero(rhoTilde, gInfo); rhoTilde->allReduceData(mpiWorld, MPIUtil::ReduceSum);
	return rhoTilde;
}

ScalarFieldTilde SaLSA::hessian(const ScalarFieldTilde& phiTilde) const
{	return (-1./(4*M_PI*gInfo.detR)) * L(phiTilde) - chi(phiTilde);
}

ScalarFieldTilde SaLSA::precondition(const ScalarFieldTilde& rTilde) const
{	return Kkernel*(J(epsInv*I(Kkernel*rTilde)));
}

double SaLSA::sync(double x) const
{	mpiWorld->bcast(x);
	return x;
}

void SaLSA::set_internal(const ScalarFieldTilde& rhoExplicitTilde, const ScalarFieldTilde& nCavityTilde)
{
	this->rhoExplicitTilde = rhoExplicitTilde; zeroNyquist(this->rhoExplicitTilde);
	
	//Compute cavity shape function (0 to 1)
	nCavity = I(nFluid * (nCavityTilde + getFullCore()));
	updateCavity();

	//Compute site shape functions with the spherical ansatz:
	const auto& solvent = fsp.solvents[0];
	for(unsigned iSite=0; iSite<solvent->molecule.sites.size(); iSite++)
		siteShape[iSite] = I(Sf[iSite] * J(shape[0]));
	
	//Update the inhomogeneity factor of the preconditioner
	epsInv = inv(1. + (epsBulk-1.)*shape[0]);
	
	//Initialize the state if it hasn't been loaded:
	if(!state) nullToZero(state, gInfo);
}


void SaLSA::minimizeFluid()
{
	logPrintf("\tSaLSA fluid occupying %lf of unit cell:", integral(shape[0])/gInfo.detR); logFlush();
	fprintf(e.fluidMinParams.fpLog, "\n\tWill stop at %d iterations, or sqrt(|r.z|)<%le\n",
		e.fluidMinParams.nIterations, e.fluidMinParams.knormThreshold);
	int nIter = solve(rhoExplicitTilde, e.fluidMinParams);
	logPrintf("\tCompleted after %d iterations at t[s]: %9.2lf\n", nIter, clock_sec());
}


double SaLSA::get_Adiel_and_grad_internal(ScalarFieldTilde& Adiel_rhoExplicitTilde, ScalarFieldTilde& Adiel_nCavityTilde, IonicGradient* extraForces) const
{
	EnergyComponents& Adiel = ((SaLSA*)this)->Adiel;
	const ScalarFieldTilde& phi = state; // that's what we solved for in minimize

	//First-order correct estimate of electrostatic energy:
	ScalarFieldTilde phiExt = coulomb(rhoExplicitTilde);
	Adiel["Electrostatic"] = -0.5*dot(phi, O(hessian(phi))) + dot(phi - 0.5*phiExt, O(rhoExplicitTilde));
	
	//Gradient w.r.t rhoExplicitTilde:
	Adiel_rhoExplicitTilde = phi - phiExt;

	//The "cavity" gradient is computed by chain rule via the gradient w.r.t to the shape function:
	const auto& solvent = fsp.solvents[0];
	ScalarFieldArray Adiel_shape(shape.size()); ScalarFieldArray Adiel_siteShape(solvent->molecule.sites.size());
	for(int r=rStart; r<rStop; r++)
	{	const MultipoleResponse& resp = *response[r];
		ScalarField& Adiel_s = resp.selectSite(Adiel_shape, Adiel_siteShape);
		if(resp.l>6) die("Angular momenta l > 6 not supported.\n");
		double prefac = 0.5 * 4*M_PI/(2*resp.l+1);
		ScalarFieldArray IlGradVphi = I(lGradient(resp.V * phi, resp.l));
		for(int lpm=0; lpm<(2*resp.l+1); lpm++)
			Adiel_s -= prefac * (IlGradVphi[lpm]*IlGradVphi[lpm]);
	}
	for(unsigned iSite=0; iSite<solvent->molecule.sites.size(); iSite++)
		if(Adiel_siteShape[iSite])
			Adiel_shape[0] += I(Sf[iSite] * J(Adiel_siteShape[iSite]));
	for(ScalarField& A_s : Adiel_shape)
	{	nullToZero(A_s, gInfo);
		A_s->allReduceData(mpiWorld, MPIUtil::ReduceSum);
	}
	
	//Propagate shape gradients to A_nCavity:
	ScalarField Adiel_nCavity;
	propagateCavityGradients(Adiel_shape, Adiel_nCavity, Adiel_rhoExplicitTilde, extraForces);
	Adiel_nCavityTilde = nFluid * J(Adiel_nCavity);
	
	accumExtraForces(extraForces, Adiel_nCavityTilde);
	return Adiel;
}

void SaLSA::getSusceptibility_internal(const std::vector<complex>& omega, std::vector<SusceptibilityTerm>& susceptibility, ScalarFieldArray& sArr, bool elecOnly) const
{	const FluidComponent& solvent = *(fsp.solvents[0]);
	susceptibility.clear();
	sArr.clear();
	for(const std::shared_ptr<MultipoleResponse>& resp: response)
	{	SusceptibilityTerm st;
		//Add shape function:
		const ScalarField& shapeCur = resp->selectSite(shape, siteShape);
		st.iSite = sArr.size();
		sArr.push_back(shapeCur);
		//Response parameters:
		st.l = resp->l;
		st.w = &(resp->V);
		//Frequency dependence:
		double nucPrefac = 0., elPrefac = 0.;
		if(resp->iSite >= 0)
			elPrefac = 1.;
		else
			nucPrefac = elecOnly ? 0. : 1.;
		st.prefactor = solvent.getChiPrefactor(omega, nucPrefac, elPrefac);
		susceptibility.push_back(st);
	}
}

void SaLSA::loadState(const char* filename)
{	ScalarField Istate(ScalarFieldData::alloc(gInfo));
	loadRawBinary(Istate, filename); //saved data is in real space
	state = J(Istate);
}

void SaLSA::saveState(const char* filename) const
{	if(mpiWorld->isHead()) saveRawBinary(I(state), filename); //saved data is in real space
}

void SaLSA::dumpDensities(const char* filenamePattern) const
{	PCM::dumpDensities(filenamePattern);
	
	//Dump effective site densities:
	double chiRot = 0., chiPol = 0.;
	for(const auto& c: fsp.components)
	{	chiRot += c->Nbulk * c->molecule.getDipole().length_squared()/(3.*fsp.T);
		chiPol += c->Nbulk * c->molecule.getAlphaTot();
	}
	double sqrtCrot = (epsBulk>epsInf && chiRot) ? sqrt((epsBulk-epsInf)/(4.*M_PI*chiRot)) : 1.;
	const double bessel_jl_by_Gl_zero[4] = {1., 1./3, 1./15, 1./105}; //G->0 limit of j_l(G)/G^l
	
	for(unsigned iComp=0; iComp<fsp.components.size(); iComp++)
	{	const std::shared_ptr<FluidComponent>& c = fsp.components[iComp];
		int iShape = std::min(iComp, unsigned(shape.size()-1)); //index to relevant cavity
		ScalarFieldTildeArray Ntilde(c->molecule.sites.size());
		for(int l=0; l<=fsp.lMax; l++)
		{	double prefac = sqrt(4.*M_PI*c->Nbulk/fsp.T) * (l==1 ? sqrtCrot : 1.);
			for(int m=-l; m<=+l; m++)
			{	const double dG = 0.02, Gmax = gInfo.GmaxGrid;
				unsigned nGradial = unsigned(ceil(Gmax/dG))+5;
				std::vector<double> VtotSamples(nGradial);
				std::vector< std::vector<double> > VsiteSamples(c->molecule.sites.size(), std::vector<double>(nGradial));
				for(unsigned iG=0; iG<nGradial; iG++)
				{	double G = iG*dG;
					for(unsigned iSite=0; iSite<c->molecule.sites.size(); iSite++)
					{	const auto& site = c->molecule.sites[iSite];
						for(const vector3<>& r: site->positions)
						{	double rLength = r.length();
							double bessel_jl_by_Gl = G ? bessel_jl(l,G*rLength)/pow(G,l) : bessel_jl_by_Gl_zero[l]*pow(rLength,l);
							vector3<> rHat = (rLength ? 1./rLength : 0.) * r;
							double Ylm_rHat = Ylm(l, m, rHat);
							VtotSamples[iG] += prefac * site->chargeKernel(G) * bessel_jl_by_Gl * Ylm_rHat;
							VsiteSamples[iSite][iG] += prefac * bessel_jl_by_Gl * Ylm_rHat;
						}
					}
				}
				
				RadialFunctionG Vtot; Vtot.init(0, VtotSamples, dG);
				ScalarFieldTilde temp = lDivergence(J(shape[iShape] * I(lGradient(Vtot * state, l))), l);
				Vtot.free();
				for(unsigned iSite=0; iSite<c->molecule.sites.size(); iSite++)
				{	RadialFunctionG Vsite; Vsite.init(0, VsiteSamples[iSite], dG);
					Ntilde[iSite] -= (pow(-1,l) * 4*M_PI/(2*l+1)) * (Vsite * temp);
					Vsite.free();
				}
			}
		}
		char filename[256]; 
		for(unsigned j=0; j<c->molecule.sites.size(); j++)
		{	const Molecule::Site& s = *(c->molecule.sites[j]);
			ScalarField N = I(Ntilde[j]) + c->Nbulk*(iComp ? shape[iShape] : siteShape[j]);
			ostringstream oss; oss << "N_" << c->molecule.name;
			if(c->molecule.sites.size()>1) oss << "_" << s.name;
			sprintf(filename, filenamePattern, oss.str().c_str());
			logPrintf("Dumping '%s' ... ", filename); logFlush();
			if(mpiWorld->isHead()) saveRawBinary(N, filename);
			
			{	//debug sphericalized site densities
				ostringstream oss; oss << "Nspherical_" << c->molecule.name;
				if(c->molecule.sites.size()>1) oss << "_" << s.name;
				sprintf(filename, filenamePattern, oss.str().c_str());
				saveSphericalized(&N,1,filename);
			}
			logPrintf("Done.\n"); logFlush();
		}
	}
}
