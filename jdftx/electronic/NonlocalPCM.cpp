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
#include <electronic/operators.h>
#include <gsl/gsl_linalg.h>

struct MultipoleResponse
{	int l; //!< angular momentum
	int iSite; //!< site index (-1 => molecule center)
	int siteMultiplicity; //!< number of sites sharing index iSite
	
	RadialFunctionG V; //!< radial part of response eigenfunction scaled by square root of eigenvalue and 1/G^l
	std::vector<double> Vsamples; //!< radial G-samples corresponding to V
	
	template<typename Func, typename... Args>
	MultipoleResponse(int l, int iSite, int siteMultiplicity, double dG, double Gmax, const Func& func, Args... args)
	: l(l), iSite(iSite), siteMultiplicity(siteMultiplicity), Vsamples(unsigned(ceil(Gmax/dG))+5)
	{
		for(unsigned i=0; i<Vsamples.size(); i++)
			Vsamples[i] = func(i*dG, args...);
		V.init(0, Vsamples, dG); //note that the G=0 behaviour of V is effectively l=0 since the G^l has been factored out
	}
	
	//version without initializing radial function
	MultipoleResponse(int l, int iSite, int siteMultiplicity, double dG, double Gmax)
	: l(l), iSite(iSite), siteMultiplicity(siteMultiplicity), Vsamples(unsigned(ceil(Gmax/dG))+5)
	{
	}
	
	~MultipoleResponse()
	{	V.free();
	}
};

//Fourier transform of cuspless exponential
inline double cusplessExpTilde(double G, double norm, double a)
{	double aG = a*G;
	double den = 1./(1.+aG*aG);
	return norm * den*den*den;
}

//Fourier transform of gaussian
inline double gaussTilde(double G, double norm, double sigma)
{	double sigmaG = sigma*G;
	return norm * exp(-0.5*sigmaG*sigmaG);
}

//Radial kernel for spherically-averaged electron density of solvent molecule (l=0 only)
inline double nFluid_calc(double G, const std::vector<FluidSolverParams::PcmSite>* pcmSites)
{	double nFluid = 0.;
	for(const auto& site: *pcmSites)
	{	double nTilde = cusplessExpTilde(G, site.Zel, site.aEl);
		for(const vector3<>& r: site.pos)
			nFluid += nTilde * bessel_jl(0, G*r.length());
	}
	return nFluid;
}

template<int l, int m> void set_Ylm(const vector3<> qHat, double& result) { result = Ylm<l,m>(qHat); }

NonlocalPCM::NonlocalPCM(const Everything& e, const FluidSolverParams& fsp)
: PCM(e, fsp), siteShape(fsp.pcmSite.size())
{	
	logPrintf("   Initializing non-local response weight functions:\n");
	const double dG = 0.02, Gmax = e.iInfo.GmaxLoc;
	unsigned nGradial = unsigned(ceil(Gmax/dG))+5;

	//Initialize fluid molecule's spherically-averaged electron density kernel:
	nFluid.init(0, dG, e.iInfo.GmaxLoc, nFluid_calc, &(fsp.pcmSite));
	
	//Determine dipole correlation factors:
	vector3<> pMol; double alphaTot = 0.;
	for(const auto& site: fsp.pcmSite)
	{	for(const vector3<>& r: site.pos) pMol += r * (site.Zel-site.Znuc);
		alphaTot += site.alpha * site.pos.size();
	}
	double l1prefac = pMol.length_squared()>1e-6 ? sqrt((fsp.epsBulk-fsp.epsInf)*3./pMol.length_squared()) : 0.;
	double polPrefac = alphaTot ? sqrt(((l1prefac ? fsp.epsInf : fsp.epsBulk) - 1.)/(4.*M_PI*alphaTot)) : 0.;
	
	//Rotational response:
	const double bessel_jl_by_Gl_zero[3] = {1., 1./3, 1./15}; //G->0 limit of j_l(G)/G^l
	for(int l=0; l<=fsp.npcmParams.lMax; l++)
	{	//Calculate radial densities for all m:
		gsl_matrix* V = gsl_matrix_calloc(nGradial, 2*l+1); //allocate and set to zero
		double prefac = (l==1 && l1prefac) ? l1prefac : sqrt(4.*M_PI*fsp.Nbulk/fsp.T);
		for(unsigned iG=0; iG<nGradial; iG++)
		{	double G = iG*dG;
			for(const auto& site: fsp.pcmSite)
			{	double VsiteTilde = prefac * (cusplessExpTilde(G, site.Zel, site.aEl) - gaussTilde(G, site.Znuc, site.sigmaNuc));
				for(const vector3<>& r: site.pos)
				{	double rLength = r.length();
					double bessel_jl_by_Gl = G ? bessel_jl(l,G*rLength)/pow(G,l) : bessel_jl_by_Gl_zero[l]*pow(rLength,l);
					vector3<> rHat = (rLength ? 1./rLength : 0.) * r;
					for(int m=-l; m<=+l; m++)
					{	double Ylm_rHat=0.; SwitchTemplate_lm(l, m, set_Ylm, (rHat, Ylm_rHat))
						*gsl_matrix_ptr(V,iG,l+m) += VsiteTilde * bessel_jl_by_Gl * Ylm_rHat;
					}
				}
			}
		}
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
			auto resp = std::make_shared<MultipoleResponse>(l, -1, 1, dG, Gmax);
			for(unsigned iG=0; iG<resp->Vsamples.size(); iG++)
				resp->Vsamples[iG] = Smode * gsl_matrix_get(V, iG, mode);
			resp->V.init(0, resp->Vsamples, dG);
			response.push_back(resp);
		}
		gsl_vector_free(S);
		gsl_matrix_free(V);
	}
	
	//Polarizability response:
	for(unsigned iSite=0; iSite<fsp.pcmSite.size(); iSite++)
	{	const auto& site = fsp.pcmSite[iSite];
		if(site.alpha)
			response.push_back(std::make_shared<MultipoleResponse>(1, iSite, site.pos.size(), dG, Gmax, cusplessExpTilde, polPrefac*sqrt(site.alpha), site.aPol));
	}
	
	//Ionic (monopole response)
	if(fsp.ionicConcentration)
	{	//Currently add gaussian ions with sigma = (1./4) ionic radius
		//TODO: get true ion charge density profiles
		double norm = fsp.ionicZelectrolyte * sqrt(fsp.ionicConcentration / fsp.T);
		response.push_back(std::make_shared<MultipoleResponse>(0, -1, 1, dG, Gmax, gaussTilde, norm, 0.25*fsp.ionicRadiusPlus));
		response.push_back(std::make_shared<MultipoleResponse>(0, -1, 1, dG, Gmax, gaussTilde, norm, 0.25*fsp.ionicRadiusMinus));
	}
	
	const double GzeroTol = 1e-12;
	
	//Compute bulk properties and print summary:
	double epsBulk = 1.; double k2factor = 0.; std::map<int,int> lCount;
	for(const std::shared_ptr<MultipoleResponse>& resp: response)
	{	lCount[resp->l]++;
		double respGzero = (4*M_PI) * pow(resp->Vsamples[0], 2) * resp->siteMultiplicity;
		if(resp->l==0) k2factor += respGzero;
		if(resp->l==1) epsBulk += respGzero;
	}
	for(auto lInfo: lCount)
		logPrintf("      l: %d  #weight-functions: %d\n", lInfo.first, lInfo.second);
	logPrintf("   Bulk dielectric-constant: %lg", epsBulk);
	if(k2factor > GzeroTol) logPrintf("   screening-length: %lg bohrs.\n", sqrt(epsBulk/k2factor));
	else logPrintf("\n");
	assert(fabs(epsBulk-fsp.epsBulk) < 1e-3); //verify consistency of correlation factors
	
	//Initialize preconditioner kernel:
	std::vector<double> KkernelSamples(nGradial);
	for(unsigned i=0; i<KkernelSamples.size(); i++)
	{	double G = i*dG, G2=G*G;
		//Compute diagonal part of the hessian ( 4pi(Vc^-1 + chi) ):
		double diagH = G2;
		for(const auto& resp: response)
			diagH += pow(G2,resp->l) * resp->Vsamples[i]*resp->Vsamples[i];
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
	for(const std::shared_ptr<MultipoleResponse>& resp: response)
	{	const DataRptr& s = resp->iSite<0 ? shape : siteShape[resp->iSite];
		switch(resp->l)
		{	case 0:
			{	rhoTilde -= resp->V * J(s*I(resp->V * phiTilde));
				break;
			}
			case 1:
			{	rhoTilde += resp->V * divergence(J(s*I(gradient(resp->V * phiTilde))));
				break;
			}
			case 2:
			{	rhoTilde -= 1.5 * (resp->V * tensorDivergence(J(s*I(tensorGradient(resp->V * phiTilde)))));
				break;
			}
			default:
				die("NonlocalPCM: Angular momentum l=%d not implemented.\n", resp->l);
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
	for(unsigned iSite=0; iSite<fsp.pcmSite.size(); iSite++)
		siteShape[iSite] = I(Sf[iSite] * J(shape));
	
	logPrintf("\tNonlocalPCM fluid occupying %lf of unit cell:", integral(shape)/e.gInfo.detR);
	logFlush();

	//Update the inhomogeneity factor of the preconditioner
	epsInv = inv(1. + (fsp.epsBulk-1.)*shape);
	
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
	const DataGptr& phiTilde = state; // that's what we solved for in minimize

	//The "electrostatic" gradient is the potential due to the bound charge alone:
	Adiel_rhoExplicitTilde = phiTilde - (-4*M_PI)*Linv(O(rhoExplicitTilde));
	Adiel["Electrostatic"] = 0.5*dot(Adiel_rhoExplicitTilde, O(rhoExplicitTilde)) //True energy if phi was an exact solution
		+ 0.5*dot(O(phiTilde), rhoExplicitTilde - hessian(phiTilde)); //First order residual correction (remaining error is second order)
	

	//The "cavity" gradient is computed by chain rule via the gradient w.r.t to the shape function:
	DataRptr Adiel_shape; DataRptrCollection Adiel_siteShape(fsp.pcmSite.size());
	for(const std::shared_ptr<MultipoleResponse>& resp: response)
	{	DataRptr& Adiel_s = resp->iSite<0 ? Adiel_shape : Adiel_siteShape[resp->iSite];
		switch(resp->l)
		{	case 0:
			{	DataRptr IVphi = I(resp->V * phiTilde);
				Adiel_s -= 0.5 * (IVphi * IVphi);
				break;
			}
			case 1:
			{	DataRptrVec IgradVphi = I(gradient(resp->V * phiTilde));
				Adiel_s -= 0.5 * lengthSquared(IgradVphi);
				break;
			}
			case 2:
			{	DataRptrTensor ItgradVphi = I(tensorGradient(resp->V * phiTilde));
				Adiel_s -= 0.5 * 1.5
					* 2*( ItgradVphi[0]*ItgradVphi[0] + ItgradVphi[1]*ItgradVphi[1] + ItgradVphi[2]*ItgradVphi[2]
						+ ItgradVphi[3]*ItgradVphi[3] + ItgradVphi[4]*ItgradVphi[4] + ItgradVphi[3]*ItgradVphi[4]);
				break;
			}
			default:
				die("NonlocalPCM: Angular momentum l=%d not implemented.\n", resp->l);
		}
	}
	for(unsigned iSite=0; iSite<fsp.pcmSite.size(); iSite++)
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

void NonlocalPCM::dumpDensities(const char* filenamePattern) const
{	string filename(filenamePattern);
	filename.replace(filename.find("%s"), 2, "Shape");
	logPrintf("Dumping '%s'... ", filename.c_str());  logFlush();
	saveRawBinary(shape, filename.c_str());
	logPrintf("done.\n"); logFlush();
}
