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

#include <core/Coulomb_ExchangeEval.h>
#include <core/CoulombIsolated.h>
#include <core/CoulombWire.h>
#include <core/CoulombKernel.h>
#include <core/Coulomb_internal.h>
#include <core/LatticeUtils.h>
#include <core/LoopMacros.h>
#include <core/BlasExtra.h>
#include <core/Util.h>
#include <gsl/gsl_sf.h>

//-------------- Auxiliary Function method --------------------------
// Based on [ P. Carrier, S. Rohra and A. Gorling, Phys. Rev. B 75, 205126 (2007) ]
// for 3D, and generalized to other geometries by R. Sundararaman et al, PRB 87, 165122 (2013)

// Reciprocal space function periodic on Brillouin zone which ~ G^2 near G=0
// and is strictly > 0 everywhere else in the zone for any crystal system
// g is a dimensionless coordinate (coefficients to reciprocal lattice vectors)
// GGT is the reciprocal space metric for this dimensionless coordinate
inline double effGsq(const vector3<>& g, const matrix3<>& GGT)
{	vector3<> sinPi; for(int i=0; i<3; i++) sinPi[i] = sin(M_PI*g[i]);
	vector3<> sin2Pi; for(int i=0; i<3; i++) sin2Pi[i] = sin(2*M_PI*g[i]);
	return (1./(M_PI*M_PI)) *
		( sinPi[0]*sinPi[0]*GGT(0,0) + sinPi[1]*sinPi[1]*GGT(1,1) + sinPi[2]*sinPi[2]*GGT(2,2)
		+ 0.5*(sin2Pi[0]*sin2Pi[1]*GGT(0,1) + sin2Pi[1]*sin2Pi[2]*GGT(1,2) + sin2Pi[2]*sin2Pi[0]*GGT(2,0)));
}

//Singular function for 3D periodic systems (from Carrier et al., but generalized
//to include screened exchange, which although not singular, benefits from regularization
//when the effective supercell is smaller than the length scale of screening)
inline double fSingular3D(const vector3<>& g, const matrix3<>& GGT, double omegaSq)
{	double Gsq = effGsq(g, GGT);
	if(omegaSq) return Gsq ? (1.-exp(-0.25*Gsq/omegaSq))*(4*M_PI)/Gsq : M_PI/omegaSq;
	else return Gsq ? (4*M_PI)/Gsq : 0.;
}

//Singular function for 2D periodic systems (slab geometry)
//Assumes g=0 along the truncated direction
inline double fSingular2D(const vector3<>& g, const matrix3<>& GGT, double omegaSq)
{	double G = sqrt(effGsq(g, GGT));
	if(omegaSq)
	{	double hlfOmegaInv=0.5/sqrt(omegaSq);
		return (2*M_PI) * hlfOmegaInv * erf_by_x(hlfOmegaInv*G);
	}
	else return G ? 2*M_PI/G : 0.;
}

//Singular function for 1D periodic systems (wire geometry)
//Assumes g=0 along the truncated directions
inline double fSingular1D(const vector3<>& g, const matrix3<>& GGT, double omegaSq)
{	double Gsq = effGsq(g, GGT);
	if(omegaSq) return Gsq ? gsl_sf_expint_Ei(-0.25*Gsq/omegaSq)-log(Gsq) : M_EULER-log(4.*omegaSq);
	else return Gsq ? log(4)-2*M_EULER-log(Gsq) : 0.;
}

//Integrate fSingular on a box of size scale and center gCenter (in reciprocal lattice coordinates)]
template<typename FSingular> double fSingularIntegralBox(FSingular fSingular, vector3<bool> isTruncated,
	double scale, const vector3<> gCenter, const matrix3<>& GGT, double omegaSq)
{	//Weights and abscissae of the 15-point gauss quadrature:
	const int N = 7;
	static const double w[N+1] =
	{	0.030753241996117268354628393577204, 0.070366047488108124709267416450667,
		0.107159220467171935011869546685869, 0.139570677926154314447804794511028,
		0.166269205816993933553200860481209, 0.186161000015562211026800561866423,
		0.198431485327111576456118326443839, 0.202578241925561272880620199967519
	};
	static const double x[N+1] =
	{	0.987992518020485428489565718586613, 0.937273392400705904307758947710209,
		0.848206583410427216200648320774217, 0.724417731360170047416186054613938,
		0.570972172608538847537226737253911, 0.394151347077563369897207370981045,
		0.201194093997434522300628303394596, 0.000000000000000000000000000000000
	};
	int N0 = isTruncated[0] ? 0 : N;
	int N1 = isTruncated[1] ? 0 : N;
	int N2 = isTruncated[2] ? 0 : N;
	double ret = 0.0;
	double h = 0.5 * scale; 
	vector3<> g;
	for(int i0=-N0; i0<=N0; i0++)
	{	g[0] = gCenter[0] + h*x[N-abs(i0)]*(i0>0?1:-1);
		double w0 = (N0 ? h*w[N-abs(i0)] : 1.);
		for(int i1=-N1; i1<=N1; i1++)
		{	g[1] = gCenter[1] + h*x[N-abs(i1)]*(i1>0?1:-1);
			double w01 = w0 * (N1 ? h*w[N-abs(i1)] : 1.);
			for(int i2=-N2; i2<=N2; i2++)
			{	g[2] = gCenter[2] + h*x[N-abs(i2)]*(i2>0?1:-1);
				double w012 = w01 * (N2 ? h*w[N-abs(i2)] : 1.);
				ret += w012 * fSingular(g, GGT, omegaSq);
			}
		}
	}
	return ret;
}

//Integrate fSingular between the box of size scale and scale/3 centered at the origin (in reciprocal lattice coordinates)
template<typename FSingular> double fSingularIntegralBoxDiff(FSingular fSingular, vector3<bool> isTruncated,
	double scale, const matrix3<>& GGT, double omegaSq)
{	double scaleBy3 = scale/3.;
	int N0 = isTruncated[0] ? 0 : 1;
	int N1 = isTruncated[1] ? 0 : 1;
	int N2 = isTruncated[2] ? 0 : 1;
	double ret = 0.0;
	vector3<int> ig;
	for(ig[0]=-N0; ig[0]<=N0; ig[0]++)
	for(ig[1]=-N1; ig[1]<=N1; ig[1]++)
	for(ig[2]=-N2; ig[2]<=N2; ig[2]++)
		if(ig.length_squared()) //except the center box
			ret += fSingularIntegralBox<>(fSingular, isTruncated, scaleBy3, scaleBy3*ig, GGT, omegaSq);
	return ret;
}

//Compute the difference between the integral of fSingular (scaled appropriately) and its sum over kmesh
template<typename FSingular> double fSingularIntegralMinusSum(FSingular fSingular, vector3<bool> isTruncated,
	const matrix3<>& GGT, double omegaSq, const std::vector<vector3<>>& kmesh)
{	//Compute the integral of fSingular over the brillouin zone
	double fSingularIntegral = 0.;
	for(double scale=1.0; scale>1e-16; scale/=3.)
		fSingularIntegral += fSingularIntegralBoxDiff(fSingular, isTruncated, scale, GGT, omegaSq);
	//Compute the sum over the k-point mesh:
	double fSingularSum = 0.;
	for(unsigned i=0; i<kmesh.size(); i++)
		fSingularSum += fSingular(kmesh[i]-kmesh[0], GGT, omegaSq);
	//Return difference:
	return fSingularIntegral * kmesh.size() - fSingularSum;
}

//Return the Vzero correction for specified lattice vectors, kmesh, truncation mode etc.
inline double getAuxiliaryFunctionVzero(matrix3<> R, vector3<bool> isTruncated, int iDir, const std::vector<vector3<>>& kmesh, double omegaSq)
{	//Compute lattice vector-dependent quantities (recomputed here for lattice derivative):
	double detR = fabs(det(R));
	matrix3<> G = 2*M_PI*inv(R);
	matrix3<> GGT = G*(~G);
	//Switch based on dimension:
	int nPeriodic=0;
	for(int k=0; k<3; k++)
		if(!isTruncated[k])
			nPeriodic++;
	switch(nPeriodic)
	{	case 0:
			assert(!"Auxiliary function method meaningless for isolated geometry.\n");
			return 0.;
		case 1:
			return (detR/R.column(iDir).length()) //transverse area to untruncated axis
				* fSingularIntegralMinusSum(fSingular1D, isTruncated, GGT, omegaSq, kmesh);
		case 2:
			return R.column(iDir).length() //truncated axis length
				* fSingularIntegralMinusSum(fSingular2D, isTruncated, GGT, omegaSq, kmesh);
		case 3:
			return fSingularIntegralMinusSum(fSingular3D, isTruncated, GGT, omegaSq, kmesh);
		default:
			return 0.; //Should never reach here (just to avoid compiler warnings)
	}
}

//---------------- Spherical truncation for screened exchange --------------
double erfcIntegrand(double r, void* params)
{	double omega = *((double*)params);
	return erfc(omega * r);
}
double truncatedErfcTilde(double G, double omega, double Rc)
{	if(G==0)
	{	double RcSq = Rc * Rc;
		double omegaSq = omega * omega;
		return M_PI * (2*Rc*(Rc - exp(-RcSq*omegaSq)/(omega*sqrt(M_PI)))
			+ (1./omegaSq - 2*RcSq) * erf(omega*Rc));
	}
	assert(G > 0.);
	const int nBisections = 20;
	gsl_integration_workspace* iWS = gsl_integration_workspace_alloc(nBisections);
	gsl_integration_qawo_table* qawoTable = gsl_integration_qawo_table_alloc(G, Rc, GSL_INTEG_SINE, nBisections);
	gsl_function f = { &erfcIntegrand, &omega };
	double result, err;
	gsl_integration_qawo(&f, 0., 0., 1e-12*std::max(0.1,G*Rc), nBisections, iWS, qawoTable, &result, &err);
	gsl_integration_qawo_table_free(qawoTable);
	gsl_integration_workspace_free(iWS);
	return (4*M_PI/G) * result;
}

//---------------- Wigner-Seitz truncated exchange ---------------

//Extract unit cell exchange kernel for a particualr k-point difference dk
//from supercell kernel in dataSuper, and store it in data
void extractExchangeKernel_thread(size_t iStart, size_t iStop, const vector3<>& dk,
	const vector3<int>& S, const vector3<int>& Ssuper, const matrix3<int>& super,
	const double* dataSuper, double* data, const symmetricMatrix3<>* dataSuper_RRT, symmetricMatrix3<>* data_RRT)
{	//Get the integer vector corresponding to dk in the supercell:
	double err;
	vector3<int> dkSuper = round(dk * super, &err); //note: transformation is by transpose of super
	assert(err < symmThreshold);
	//Collect data from supercell to unit-cell with k-point offset:
	THREAD_fullGspaceLoop
	(	vector3<int> iGsuper = iG * super + dkSuper;
		//Reduce to inversion-symmetry reduced index in supercell data:
		if(iGsuper[2]<0) iGsuper = -iGsuper;
		if(iGsuper[1]<0) iGsuper[1] += Ssuper[1];
		if(iGsuper[0]<0) iGsuper[0] += Ssuper[0];
		size_t iSuper = iGsuper[2] + (1+Ssuper[2]/2) * size_t(iGsuper[1] + Ssuper[1]*iGsuper[0]);
		data[i] = dataSuper[iSuper];
		if(data_RRT) data_RRT[i] = dataSuper_RRT[iSuper];
	)
}


//-------------------- class ExchangeEval -----------------------

ExchangeEval::ExchangeEval(const GridInfo& gInfo, const CoulombParams& params, const Coulomb& coulomb, double omega)
: gInfo(gInfo), omega(omega), VcGamma(0), VcGamma_RRT(0), kernelData(0)
{
	if(!omega) logPrintf("\n-------- Setting up exchange kernel --------\n");
	else logPrintf("\n--- Setting up screened exchange kernel (omega = %lg) ---\n", omega);
	
	//Obtain supercell parameters, and adjust for mesh embedding where necessary:
	assert(params.supercell);
	const matrix3<int>& super = params.supercell->super;
	const std::vector< vector3<> >& kmesh = params.supercell->kmesh;
	matrix3<> Rsuper = gInfo.R * super; //this could differ from supercell->Rsuper, because the embedding gInfo.R is scaled up from the original gInfo.R
	
	//Check supercell:
	if(params.geometry != CoulombParams::Periodic)
	{	vector3<bool> isTruncated = params.isTruncated();
		//Make sure that the supercell and unit cell are identical in truncated directions:
		for(int k=0; k<3; k++)
			if(isTruncated[k])
			{	vector3<int> comb = super.column(k); //linear combinations of lattice-vectors in k^th supercell vector
				if((comb.length_squared()!=1) || (abs(comb[k])!=1))
				{	string dirName(3, '0'); dirName[k] = '1';
					die("More than one k-point along truncated lattice direction %s\n", dirName.c_str());
				}
			}
	}
	
	//Select kernel mode:
	switch(params.exchangeRegularization)
	{	//G=0 correction based modes:
		case CoulombParams::None:
		case CoulombParams::AuxiliaryFunction:
		case CoulombParams::ProbeChargeEwald:
			switch(params.geometry)
			{	case CoulombParams::Periodic:
					kernelMode = PeriodicKernel; break;
				case CoulombParams::Slab:
					kernelMode = SlabKernel; break;
				case CoulombParams::Wire:
				case CoulombParams::Cylindrical:
					kernelMode = NumericalKernel; break;
				case CoulombParams::Spherical:
					kernelMode = SphericalKernel; break;
				case CoulombParams::Isolated:
					kernelMode = WignerSeitzGammaKernel; break;
			}
			if(params.geometry==CoulombParams::Periodic)
			{	if(params.exchangeRegularization==CoulombParams::AuxiliaryFunction)
					Citations::add("Auxiliary function method for exact exchange in arbitrary lattice systems", "P. Carrier et al., Phys. Rev. B 75, 205126 (2007)");
				if(params.exchangeRegularization==CoulombParams::ProbeChargeEwald)
					Citations::add("Probe-charge Ewald method for exact exchange", "J. Paier et al, J. Chem. Phys. 122, 234102 (2005)");
			}
			else if(params.exchangeRegularization!=CoulombParams::None)
				Citations::add("Exact exchange in reduced-dimensionality systems", wsTruncationPaper);
			break;
		//Truncation based modes:
		case CoulombParams::SphericalTruncated:
			Citations::add("Spherical truncated method for exact exchange", "J. Spencer et al, Phys. Rev. B 77, 193110 (2008)");
			kernelMode = SphericalKernel; break;
		case CoulombParams::WignerSeitzTruncated:
			Citations::add("Wigner-Seitz truncated method for exact exchange", wsTruncationPaper);
			kernelMode = NumericalKernel; break;
	}
	
	//Perform G=0 handling if required
	double VzeroCorrection = 0.;
	matrix3<> VzeroCorrection_RRT; //optional lattice derivative
	double detRsuper = fabs(det(Rsuper));
	
	if(params.exchangeRegularization==CoulombParams::AuxiliaryFunction)
	{	double omegaSq = omega*omega;
		vector3<bool> isTruncated = params.isTruncated();
		VzeroCorrection = getAuxiliaryFunctionVzero(gInfo.R, isTruncated, params.iDir, kmesh, omegaSq);
		//Compute lattice derivative if requied:
		if(params.computeStress)
		{	const double h = 1e-5; //finite difference size
			matrix3<> id(1,1,1); //identity
			std::vector<matrix3<>> strainBasis(6);
			for(int iDir=0; iDir<3; iDir++)
			{	//Tensions:
				strainBasis[iDir](iDir,iDir) = 1.;
				int jDir = (iDir+1)%3;
				int kDir = (iDir+2)%3;
				//Shears:
				strainBasis[iDir+3](jDir,kDir) = sqrt(0.5);
				strainBasis[iDir+3](kDir,jDir) = sqrt(0.5);
			}
			for(matrix3<> strainDir: strainBasis)
			{	double VzeroCorrection_dir = (0.5/h) * 
					( getAuxiliaryFunctionVzero((id+h*strainDir)*gInfo.R, isTruncated, params.iDir, kmesh, omegaSq)
					- getAuxiliaryFunctionVzero((id-h*strainDir)*gInfo.R, isTruncated, params.iDir, kmesh, omegaSq) );
				VzeroCorrection_RRT += VzeroCorrection_dir * strainDir;
			}
		}
	}
	
	if(params.exchangeRegularization==CoulombParams::ProbeChargeEwald)
	{	double Eperiodic = 0.; //Periodic interaction of a point charge in supercell geometry
		matrix3<> Eperiodic_RRT;
		if(omega) //Directly compute the periodic interaction in real space
		{	matrix3<> invRsuper = inv(Rsuper);
			vector3<bool> isTruncated = params.isTruncated();
			vector3<int> Nreal(0,0,0);
			double rMax = CoulombKernel::nSigmasPerWidth * sqrt(0.5)/omega;
			for(int k=0; k<3; k++)
				if(!isTruncated[k]) //no sampling along truncated directions
					Nreal[k] = 1+ceil(rMax * invRsuper.row(k).length());
			//Loop over neighbouring cells in real space:
			matrix3<> RsuperTRsuper = (~Rsuper)*Rsuper;
			vector3<int> iR; //integer cell number
			double omegaSq = omega*omega;
			for(iR[0]=-Nreal[0]; iR[0]<=Nreal[0]; iR[0]++)
				for(iR[1]=-Nreal[1]; iR[1]<=Nreal[1]; iR[1]++)
					for(iR[2]=-Nreal[2]; iR[2]<=Nreal[2]; iR[2]++)
					{	double rSq = RsuperTRsuper.metric_length_squared(iR);
						if(!rSq) continue; //exclude self-interaction
						double r = sqrt(rSq);
						Eperiodic += 0.5*erfc(omega*r)/r;
						if(params.computeStress)
						{	vector3<> rVec = Rsuper * iR;
							double minus_E_r_by_r = 0.5*(erfc(omega*r)/r + (2./sqrt(M_PI))*omega*exp(-omegaSq*rSq))/rSq;
							Eperiodic_RRT -= minus_E_r_by_r * outer(rVec,rVec);
						}
					}
		}
		else //Use the appropriate Ewald method
		{	std::vector<Atom> atoms(1, Atom(1., vector3<>())); //single unit point charge
			logSuspend();
			Eperiodic = coulomb.createEwald(Rsuper, 1)->energyAndGrad(atoms, params.computeStress ? &Eperiodic_RRT : 0);
			logResume();
			//Correction for G=0 difference between cylinder and wire truncation modes:
			if(params.geometry == CoulombParams::Cylindrical)
			{	double rho0 = ((CoulombCylindrical&)coulomb).Rc; //cylinder mode uses this as reference rho in logarithmic singularity
				double L = Rsuper.column(params.iDir).length();
				Eperiodic -= log(rho0) / L;
				if(params.computeStress)
				{	vector3<> zHat = (1./L) * Rsuper.column(params.iDir);
					Eperiodic_RRT += (log(rho0) / L) * outer(zHat,zHat);
				}
			}
		}
		VzeroCorrection = (-2.*detRsuper) * Eperiodic;
		if(params.computeStress)
			VzeroCorrection_RRT = (-2.*detRsuper) * (Eperiodic_RRT + matrix3<>(1,1,1)*Eperiodic);
	}
	
	if(VzeroCorrection) logPrintf("Vxx(G=0) correction = %le\n", VzeroCorrection/detRsuper);
	
	//Setup the kernel depending on the kernel mode:
	switch(kernelMode)
	{	case PeriodicKernel:
		{	Vzero = VzeroCorrection + (omega ? M_PI/(omega*omega) : 0.);
			Vzero_RRT = VzeroCorrection_RRT;
			logPrintf("3D periodic kernel with %s G=0.\n", VzeroCorrection ? "modified" : "unmodified");
			break;
		}
		case SphericalKernel:
		{	//Select the truncation radius:
			if(params.exchangeRegularization==CoulombParams::None)
				Rc = ((CoulombSpherical&)coulomb).Rc; //same radius as truncation for Hartree/Vloc
			else
				Rc = pow((3./(4*M_PI)) * detRsuper, 1./3); //from supercell volume
			logPrintf("Truncated on a sphere of radius %lg bohrs.\n", Rc);
			if(omega)
			{	//Initialize look-up table for erfc fourier transform
				//--- Determine end of range for look-up table:
				double Gmax = 0.0;
				vector3<int> c;
				for(c[0]=-1; c[0]<=1; c[0]+=2) for(c[1]=-1; c[1]<=1; c[1]+=2) for(c[2]=-1; c[2]<=1; c[2]+=2)
				{	vector3<> f; for(int k=0; k<3; k++) f[k] = c[k]*(gInfo.S[k]/2 + 1);
					double G = sqrt(gInfo.GGT.metric_length_squared(f));
					if(G>Gmax) Gmax=G;
				}
				//--- Compute samples and spline coefficients:
				const double dG = 0.01 * (2*M_PI)/Rc;
				int nSamples = int(ceil(Gmax/dG)) + 10;
				std::vector<double> samples(nSamples);
				for(int i=0; i<nSamples; i++)
					samples[i] = truncatedErfcTilde(dG*i, omega, Rc);
				Vzero = samples[0]; //Note: this mode will always have VzeroCorrection = 0.
				sphericalScreenedCoeff = ManagedArray<double>(QuinticSpline::getCoeff(samples));
				sphericalScreenedCalc.coeff = sphericalScreenedCoeff.dataPref();
				sphericalScreenedCalc.dGinv = 1.0/dG;
				sphericalScreenedCalc.nSamples = nSamples;
				sphericalScreenedCalc.Rc = Rc;
				sphericalScreenedCalc.erfcOmegaRc_4piBy3 = (4*M_PI/3) * erfc(omega*Rc);
			}
			else Vzero = (2*M_PI) * Rc*Rc;  //Note: this mode will always have VzeroCorrection = 0.
			//Optionally contribute G=0 contribution to stress:
			if(params.computeStress)
				Vzero_RRT = matrix3<>(1.,1.,1.) * (4*M_PI/3) * Rc*Rc * erfc(Rc * omega); //since Rc ~ detR ^ (1./3)
			break;
		}
		case SlabKernel:
		{	slabCalc.iDir = params.iDir;
			slabCalc.hlfL = 0.5 * gInfo.R.column(params.iDir).length();
			logPrintf("Truncated on a slab of thickness %lg bohrs.\n", 2*slabCalc.hlfL);
			if(omega)
			{	//Initialize look-up tables (quintic splines) of the difference between the screened
				//and the analytic unscreened kernels as a function of Gplane for each iG[iDir]:
				//--- Determine end of range for look-up table:
				double GplaneMax = 0.0;
				vector3<int> c;
				for(c[0]=-1; c[0]<=1; c[0]+=2) for(c[1]=-1; c[1]<=1; c[1]+=2) for(c[2]=-1; c[2]<=1; c[2]+=2)
				{	vector3<> f; for(int k=0; k<3; k++) f[k] = c[k]*(gInfo.S[k]/2 + 1);
					f[params.iDir] = 0.; //project out truncated direction
					double Gplane = sqrt(gInfo.GGT.metric_length_squared(f));
					if(Gplane>GplaneMax) GplaneMax=Gplane;
				}
				//--- Compute samples for each Gplane:
				const double dG = 0.01 * (2*M_PI)/slabCalc.hlfL;
				int nSamples = int(ceil(GplaneMax/dG)) + 10;
				int nAxis = gInfo.S[params.iDir]/2+1;
				double hAxis = 2.*slabCalc.hlfL/gInfo.S[params.iDir];
				std::vector<std::vector<double>> samples(nAxis, std::vector<double>(nSamples));
				ManagedArray<double> fftMem; fftMem.init(nAxis);
				double* fftArr = fftMem.data();
				fftw_plan planDCT = fftw_plan_r2r_1d(nAxis, fftArr, fftArr, FFTW_REDFT00, FFTW_ESTIMATE);
				for(int iSample=0; iSample<nSamples; iSample++)
				{	double Gplane = iSample*dG;
					//Initialize the partial fourier transform of the smooth part (-erf(omega r)/r before truncation):
					for(int iAxis=0; iAxis<nAxis; iAxis++)
					{	double z = iAxis*hAxis;
						if(Gplane==0.)
							fftArr[iAxis] = (1./(omega*sqrt(M_PI))) * exp(-omega*omega*z*z) + z*erf(omega*z);
						else
						{	fftArr[iAxis] = (-0.5/Gplane) * ( (Gplane*z > 100.) ? 0.
								: ( exp(Gplane*z) * erfc(0.5*Gplane/omega + omega*z)
								 + exp(-Gplane*z) * erfc(0.5*Gplane/omega - omega*z) ) );
						}
					}
					//Transform the truncated direction to reciprocal space:
					fftw_execute(planDCT);
					for(int iAxis=0; iAxis<nAxis; iAxis++)
						samples[iAxis][iSample] = (2*M_PI * hAxis) * fftArr[iAxis];
				}
				fftw_destroy_plan(planDCT);
				//--- Initialize splines for each iAxis:
				Vzero = samples[0][0] + (-2.*M_PI)*pow(slabCalc.hlfL,2) + VzeroCorrection;
				for(int iSample=0; iSample<nSamples; iSample++) samples[0][iSample] *= (iSample * dG); //handle 1/Gplane sinbularity in the Gz=0 plane
				samples[0][0] = (-4.*M_PI)*slabCalc.hlfL; //replace regularized version with limit of singular one multiplied out as above
				std::vector<double> coeff, coeffSub;
				for(int iAxis=0; iAxis<nAxis; iAxis++)
				{	coeffSub = QuinticSpline::getCoeff(samples[iAxis]); //coefficients for a given iAxis
					coeff.insert(coeff.end(), coeffSub.begin(), coeffSub.end());
				}
				slabCoeff = ManagedArray<double>(coeff);
				slabCalc.coeff = slabCoeff.dataPref();
				slabCalc.dGinv = 1.0/dG;
				slabCalc.nSamples = nSamples;
				slabCalc.nCoeff = coeffSub.size();
			}
			else Vzero = (-2.*M_PI)*pow(slabCalc.hlfL,2) + VzeroCorrection; //Everything else handled analytically
			break;
		}
		case WignerSeitzGammaKernel:
		{	if(abs(det(super)) != 1)
				die("Exact-exchange in Isolated geometry should be used only with a single k-point.\n");
			if(omega) //Create an omega-screened version (but gamma-point only):
			{	VcGamma = new RealKernel(gInfo);
				symmetricMatrix3<>* VcGamma_RRTdata = 0;
				if(params.computeStress)
				{	VcGamma_RRT = new ManagedArray<symmetricMatrix3<>>();
					VcGamma_RRT->init(gInfo.nG);
					VcGamma_RRTdata = VcGamma_RRT->data();
				}
				CoulombKernel(gInfo.R, gInfo.S, params.isTruncated(), omega).compute(VcGamma->data(), ((CoulombIsolated&)coulomb).ws, VcGamma_RRTdata);
			}
			else //use the same kernel as hartree/Vloc
			{	VcGamma = &((CoulombIsolated&)coulomb).Vc; 
				if(params.computeStress)
					VcGamma_RRT = &((CoulombIsolated&)coulomb).Vc_RRT; 
				logPrintf("Using previously initialized isolated coulomb kernel.\n");
			}
		}
		case NumericalKernel:
		{	//Create the kernel on the k-point supercell:
			vector3<bool> isTruncated = params.exchangeRegularization==CoulombParams::WignerSeitzTruncated
				? vector3<bool>(true, true, true) //All directions truncated for Wigner-Seitz truncated method
				: params.isTruncated(); //Same truncation geometry as Hartree/Vloc for G=0 based methods
			//--- set up supercell sample count:
			vector3<int> Ssuper(0,0,0), s; //loop over vertices of parallelopiped:
			for(s[0]=-1; s[0]<=1; s[0]+=2)
			for(s[1]=-1; s[1]<=1; s[1]+=2)
			for(s[2]=-1; s[2]<=1; s[2]+=2)
			{	vector3<> iG;
				for(int k=0; k<3; k++)
					iG[k] = s[k] * (gInfo.S[k]/2 + 1); //include margin for k-point
				vector3<> iGsuper = (~super) * iG;
				for(int k=0; k<3; k++)
				{	int Ssuper_k = 2*std::abs(iGsuper[k]);
					if(Ssuper_k > Ssuper[k])
						Ssuper[k] = Ssuper_k;
				}
			}
			logPrintf("Creating Wigner-Seitz truncated kernel on k-point supercell with sample count ");
			Ssuper.print(globalLog, " %d");
			//Note: No FFTs of dimensions Ssuper are required (so no need to make it fftSuitable())
			//--- create kernel on supercell:
			size_t nGsuper = Ssuper[0]*(Ssuper[1]*size_t(1+Ssuper[2]/2));
			double* dataSuper = new double[nGsuper];
			if(!dataSuper) die_alone("Out of memory. (need %.1lfGB for supercell exchange kernel)\n", nGsuper*1e-9*sizeof(double));
			symmetricMatrix3<>* dataSuper_RRT = 0;
			if(params.computeStress)
			{	dataSuper_RRT = new symmetricMatrix3<>[nGsuper];
				if(!dataSuper_RRT) die_alone("Out of memory. (need %.1lfGB for lattice derivatives of supercell exchange kernel)\n", 6*nGsuper*1e-9*sizeof(double));
			}
			WignerSeitz wsSuper(Rsuper);
			CoulombKernel(Rsuper, Ssuper, isTruncated, omega).compute(dataSuper, wsSuper, dataSuper_RRT);
			dataSuper[0] += VzeroCorrection; //For slab/wire geometry kernels in AuxiliaryFunction/ProbeChargeEwald methods
			if(params.computeStress) dataSuper_RRT[0] += VzeroCorrection_RRT;
			
			//Construct k-point difference mesh:
			for(const vector3<>& kpoint: kmesh)
			{	vector3<> dk = kpoint - kmesh.front();
				for(int k=0; k<3; k++) dk[k] -= floor(dk[k] + 0.5); //reduce to fundamental zone:
				 dkArr.push_back(dk);
			}
			
			//Split supercell kernel into one for each k-point difference:
			logPrintf("Splitting supercell kernel to unit-cell with k-points ... "); logFlush();
			size_t nKernelData = dkArr.size() * gInfo.nr;
			kernelData.init(nKernelData);
			if(params.computeStress)
				kernelData_RRT.init(nKernelData);
			for(size_t i=0; i<dkArr.size(); i++)
				threadLaunch(extractExchangeKernel_thread, gInfo.nr, dkArr[i],
					gInfo.S, Ssuper, super, dataSuper, kernelData.data() + i*gInfo.nr,
					dataSuper_RRT, params.computeStress ? kernelData_RRT.data() + i*gInfo.nr : 0);
			delete[] dataSuper;
			if(params.computeStress)
				delete[] dataSuper_RRT;
			logPrintf("Done.\n");
			break;
		}
	}
}

ExchangeEval::~ExchangeEval()
{	
	if(VcGamma && omega) delete VcGamma;
	if(VcGamma_RRT && omega) delete VcGamma_RRT;
}


void multTransformedKernel(complexScalarFieldTilde& X, const double* kernel, const vector3<int>& offset)
{	assert(X);
	if(!offset.length_squared())
		callPref(eblas_zmuld)(X->gInfo.nr, kernel, 1, X->dataPref(false), 1);
	else
		callPref(multTransformedKernel)(X->gInfo.S, kernel, X->dataPref(false), offset);
}


complexScalarFieldTilde ExchangeEval::operator()(complexScalarFieldTilde&& in, vector3<> kDiff) const
{
	#define CALL_exchangeAnalytic(calc) callPref(exchangeAnalytic)(gInfo.S, gInfo.GGT, calc, in->dataPref(false), kDiff, Vzero, symmThresholdSq)
	switch(kernelMode)
	{	case PeriodicKernel:
		{	if(omega) CALL_exchangeAnalytic(ExchangePeriodicScreened_calc(omega));
			else CALL_exchangeAnalytic(ExchangePeriodic_calc());
			break;
		}
		case SphericalKernel:
		{	if(omega) CALL_exchangeAnalytic(sphericalScreenedCalc);
			else CALL_exchangeAnalytic(ExchangeSpherical_calc(Rc));
			break;
		}
		case SlabKernel:
		{	CALL_exchangeAnalytic(slabCalc);
			break;
		}
		case WignerSeitzGammaKernel:
		{	assert(kDiff.length_squared() < symmThresholdSq); //gamma-point only
			callPref(multRealKernel)(gInfo.S, VcGamma->dataPref(), in->dataPref(false));
			break;
		}
		case NumericalKernel:
		{	//Find the appropriate kDiff:
			bool kDiffFound = false;
			for(unsigned ik=0; ik<dkArr.size(); ik++)
				if(circDistanceSquared(dkArr[ik], kDiff) < symmThresholdSq)
				{	//Find the integer offset, if any:
					double err;
					vector3<int> offset = round(dkArr[ik] - kDiff, &err);
					assert(err < symmThreshold);
					//Multiply kernel:
					multTransformedKernel(in, kernelData.dataPref() + gInfo.nr * ik, offset);
					kDiffFound = true;
					break;
				}
			assert(kDiffFound);
			break;
		}
	}
	#undef CALL_exchangeAnalytic
	return in;
}

matrix3<> ExchangeEval::latticeGradient(const complexScalarFieldTilde& X, vector3<> kDiff) const
{
	#define RETURN_exchangeAnalyticStress(calc) \
	{	ManagedArray<symmetricMatrix3<>> result; result.init(gInfo.nr, isGpuEnabled()); \
		callPref(exchangeAnalyticStress)(gInfo.S, gInfo.G, calc, X->dataPref(), result.dataPref(), kDiff, symmThresholdSq); \
		matrix3<> resultSum = callPref(eblas_sum)(gInfo.nr, result.dataPref()); \
		if(kDiff.length_squared() < symmThresholdSq) resultSum += Vzero_RRT * X->getGzero().norm(); \
		return gInfo.detR * resultSum; \
	}
	switch(kernelMode)
	{	case PeriodicKernel:
		{	if(omega) RETURN_exchangeAnalyticStress(ExchangePeriodicScreened_calc(omega))
			else RETURN_exchangeAnalyticStress(ExchangePeriodic_calc())
		}
		case SphericalKernel:
		{	if(omega) RETURN_exchangeAnalyticStress(sphericalScreenedCalc)
			else RETURN_exchangeAnalyticStress(ExchangeSpherical_calc(Rc))
		}
		//case SlabKernel:
			//RETURN_exchangeAnalyticStress(slabCalc)
		case WignerSeitzGammaKernel:
		{	assert(kDiff.length_squared() < symmThresholdSq); //gamma-point only
			ManagedArray<symmetricMatrix3<>> result; result.init(gInfo.nr, isGpuEnabled());
			callPref(realKernelStress)(gInfo.S, VcGamma_RRT->dataPref(), X->dataPref(), result.dataPref());
			matrix3<> resultSum = callPref(eblas_sum)(gInfo.nr, result.dataPref());
			return gInfo.detR * resultSum;
		}
		case NumericalKernel:
		{	//Find kernel with the appropriate kDiff:
			for(unsigned ik=0; ik<dkArr.size(); ik++)
				if(circDistanceSquared(dkArr[ik], kDiff) < symmThresholdSq)
				{	//Find the integer offset, if any:
					double err;
					vector3<int> offset = round(dkArr[ik] - kDiff, &err);
					assert(err < symmThreshold);
					//Compute stress:
					ManagedArray<symmetricMatrix3<>> result; result.init(gInfo.nr, isGpuEnabled());
					callPref(transformedKernelStress)(gInfo.S, kernelData_RRT.dataPref() + gInfo.nr * ik, X->dataPref(), result.dataPref(), offset);
					matrix3<> resultSum = callPref(eblas_sum)(gInfo.nr, result.dataPref());
					return gInfo.detR * resultSum;
				}
			assert(!"Encountered invalid kDiff");
			return matrix3<>();
		}
		default:
			die_alone("Lattice gradient not yet implemented for this kernel mode.");
	}
	#undef RETURN_exchangeAnalyticStress
	return matrix3<>(); //to avoid compiler warnings
}
