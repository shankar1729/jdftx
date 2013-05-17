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

#include <core/DataIO.h>
#include <core/DataMultiplet.h>
#include <electronic/NonlocalPCM.h>
#include <electronic/PCM_internal.h>
#include <electronic/Everything.h>
#include <electronic/SphericalHarmonics.h>
#include <electronic/VanDerWaals.h>
#include <electronic/operators.h>
#include <gsl/gsl_linalg.h>

struct MultipoleResponse
{	int l; //!< angular momentum
	int iSite; //!< site index (-1 => molecule center)
	int siteMultiplicity; //!< number of sites sharing index iSite
	
	RadialFunctionG V; //!< radial part of response eigenfunctions scaled by square root of eigenvalue and 1/G^l
	
	MultipoleResponse(int l, int iSite, int siteMultiplicity, const std::vector<double>& Vsamples, double dG)
	: l(l), iSite(iSite), siteMultiplicity(siteMultiplicity)
	{	V.init(0, Vsamples, dG);
	}
	
	~MultipoleResponse()
	{	V.free();
	}
};

template<int l, int m> void set_Ylm(const vector3<> qHat, double& result) { result = Ylm<l,m>(qHat); }

NonlocalPCM::NonlocalPCM(const Everything& e, const FluidSolverParams& fsp)
: PCM(e, fsp), siteShape(fsp.solvents[0]->molecule.sites.size())
{	
	logPrintf("   Initializing non-local response weight functions:\n");
	const double dG = 0.02, Gmax = e.gInfo.GmaxGrid;
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
	const double bessel_jl_by_Gl_zero[3] = {1., 1./3, 1./15}; //G->0 limit of j_l(G)/G^l
	for(const auto& c: fsp.components)
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
						{	double Ylm_rHat=0.; SwitchTemplate_lm(l, m, set_Ylm, (rHat, Ylm_rHat))
							*gsl_matrix_ptr(V,iG,l+m) += Vsite * bessel_jl_by_Gl * Ylm_rHat;
						}
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
				response.push_back(std::make_shared<MultipoleResponse>(l, -1, 1, Vsamples, dG));
			}
			gsl_vector_free(S);
			gsl_matrix_free(V);
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
	assert(fabs(epsBulk-this->epsBulk) < 1e-3); //verify consistency of correlation factors
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
}

NonlocalPCM::~NonlocalPCM()
{	nFluid.free();
	Kkernel.free();
}


DataGptr NonlocalPCM::chi(const DataGptr& phiTilde) const
{	DataGptr rhoTilde;
	for(const auto& resp: response)
	{	const DataRptr& s = resp->iSite<0 ? shape : siteShape[resp->iSite];
		switch(resp->l)
		{	case 0: rhoTilde -= resp->V * J(s*I(resp->V * phiTilde)); break;
			case 1: rhoTilde += resp->V * divergence(J(s*I(gradient(resp->V * phiTilde)))); break;
			case 2: rhoTilde -= 1.5 * (resp->V * tensorDivergence(J(s*I(tensorGradient(resp->V * phiTilde))))); break;
			default: die("NonlocalPCM: Angular momentum l=%d not implemented.\n", resp->l);
		}
	}
	return rhoTilde;
}


DataGptr NonlocalPCM::hessian(const DataGptr& phiTilde) const
{	return (-1./(4*M_PI*e.gInfo.detR)) * L(phiTilde) - chi(phiTilde);
}

DataGptr NonlocalPCM::precondition(const DataGptr& rTilde) const
{	return Kkernel*(J(epsInv*I(Kkernel*rTilde)));
}


void NonlocalPCM::set(const DataGptr& rhoExplicitTilde, const DataGptr& nCavityTilde)
{
	this->rhoExplicitTilde = clone(rhoExplicitTilde); zeroNyquist(this->rhoExplicitTilde);
	
	//Compute cavity shape function (0 to 1)
	nCavity = I(nFluid * nCavityTilde);
	updateCavity();

	//Compute site shape functions with the spherical ansatz:
	const auto& solvent = fsp.solvents[0];
	for(unsigned iSite=0; iSite<solvent->molecule.sites.size(); iSite++)
		siteShape[iSite] = I(Sf[iSite] * J(shape));
	
	logPrintf("\tNonlocalPCM fluid occupying %lf of unit cell:", integral(shape)/e.gInfo.detR);
	logFlush();

	//Update the inhomogeneity factor of the preconditioner
	epsInv = inv(1. + (epsBulk-1.)*shape);
	
	//Initialize the state if it hasn't been loaded:
	if(!state) nullToZero(state, e.gInfo);
}


void NonlocalPCM::minimizeFluid()
{
	fprintf(e.fluidMinParams.fpLog, "\n\tWill stop at %d iterations, or sqrt(|r.z|)<%le\n",
		e.fluidMinParams.nIterations, e.fluidMinParams.knormThreshold);
	int nIter = solve(rhoExplicitTilde, e.fluidMinParams);
	logPrintf("\tCompleted after %d iterations.\n", nIter);
}

double NonlocalPCM::get_Adiel_and_grad(DataGptr& Adiel_rhoExplicitTilde, DataGptr& Adiel_nCavityTilde, IonicGradient& extraForces) const
{
	EnergyComponents& Adiel = ((NonlocalPCM*)this)->Adiel;
	const DataGptr& phi = state; // that's what we solved for in minimize

	//First-order correct estimate of electrostatic energy:
	DataGptr phiExt = coulomb(rhoExplicitTilde);
	Adiel["Electrostatic"] = -0.5*dot(phi, O(hessian(phi))) + dot(phi - 0.5*phiExt, O(rhoExplicitTilde));
	
	//Gradient w.r.t rhoExplicitTilde:
	Adiel_rhoExplicitTilde = phi - phiExt;

	//The "cavity" gradient is computed by chain rule via the gradient w.r.t to the shape function:
	const auto& solvent = fsp.solvents[0];
	DataRptr Adiel_shape; DataRptrCollection Adiel_siteShape(solvent->molecule.sites.size());
	for(const std::shared_ptr<MultipoleResponse>& resp: response)
	{	DataRptr& Adiel_s = resp->iSite<0 ? Adiel_shape : Adiel_siteShape[resp->iSite];
		switch(resp->l)
		{	case 0:
			{	DataRptr IVphi = I(resp->V * phi);
				Adiel_s -= 0.5 * (IVphi * IVphi);
				break;
			}
			case 1:
			{	DataRptrVec IgradVphi = I(gradient(resp->V * phi));
				Adiel_s -= 0.5 * lengthSquared(IgradVphi);
				break;
			}
			case 2:
			{	DataRptrTensor ItgradVphi = I(tensorGradient(resp->V * phi));
				Adiel_s -= 0.5 * 1.5
					* 2*( ItgradVphi[0]*ItgradVphi[0] + ItgradVphi[1]*ItgradVphi[1] + ItgradVphi[2]*ItgradVphi[2]
						+ ItgradVphi[3]*ItgradVphi[3] + ItgradVphi[4]*ItgradVphi[4] + ItgradVphi[3]*ItgradVphi[4]);
				break;
			}
			default:
				die("NonlocalPCM: Angular momentum l=%d not implemented.\n", resp->l);
		}
	}
	for(unsigned iSite=0; iSite<solvent->molecule.sites.size(); iSite++)
		if(Adiel_siteShape[iSite])
			Adiel_shape += I(Sf[iSite] * J(Adiel_siteShape[iSite]));
	
	//Propagate shape gradients to A_nCavity:
	DataRptr Adiel_nCavity;
	propagateCavityGradients(Adiel_shape, Adiel_nCavity);
	Adiel_nCavityTilde = nFluid * J(Adiel_nCavity);
	
	if(vdwForces) extraForces = *vdwForces;
	return Adiel;
}

void NonlocalPCM::loadState(const char* filename)
{	DataRptr Istate(DataR::alloc(e.gInfo));
	loadRawBinary(Istate, filename); //saved data is in real space
	state = J(Istate);
}

void NonlocalPCM::saveState(const char* filename) const
{	saveRawBinary(I(state), filename); //saved data is in real space
}
