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

#include <core/CoulombWire.h>
#include <core/Operators.h>
#include <core/Util.h>
#include <core/Bspline.h>
#include <core/LoopMacros.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>

//!Check orthogonality and return lattice direction name (declared in CoulombSlab.cpp)
extern string checkOrthogonality(const GridInfo& gInfo, int iDir);

//! Compute Cbar_k^sigma - the gaussian convolved cylindrical coulomb kernel - by numerical quadrature
struct Cbar
{
	Cbar() { iWS = gsl_integration_workspace_alloc(maxIntervals); }
	~Cbar() { gsl_integration_workspace_free(iWS); }
	
	//! Compute Cbar_k^sigma(rho)
	double operator()(double k, double sigma, double rho)
	{	assert(k >= 0.);
		assert(sigma > 0.);
		assert(rho >= 0.);
		if(k == 0) //Use closed form in terms of the exponential integral function:
		{	const double xMax = 700.; //threshold (with some margin) to underflow in expint_E1
			double hlfSigmaInvSq = 0.5/(sigma*sigma);
			double x = hlfSigmaInvSq*rho*rho;
			if(x < 3.5e-3) return (M_EULER + log(hlfSigmaInvSq)) - x*(1. - x*(1./4 - x*(1./18 - x*(1./96))));
			else return -2.*log(rho) - (x>xMax ? 0. : gsl_sf_expint_E1(x));
		}
		else
		{	double R = rho/sigma;
			double K = k*sigma;
			if(R*(R-2*K) > 100.)
				return 2. * gsl_sf_bessel_K0_scaled(k*rho) * exp(-k*rho);
			std::pair<double,double> params;
			gsl_function f; f.params = &params;
			if(R < 1.)
			{	params.first  = R;
				params.second = K;
				f.function = &integrandSmallRho;
			}
			else
			{	params.first  = R*R;
				params.second = K*R;
				f.function = &integrandLargeRho;
			}
			double result, err;
			gsl_integration_qagiu(&f, 0., 0., 1e-13, maxIntervals, iWS, &result, &err);
			return 2 * exp(-0.5*(K*K + R*R)) * result;
		}
	}
	
private:
	double absTol, relTol; //!< Absolute and relative tolerance
	static const size_t maxIntervals = 1000; //!< Size of integration workspace
	gsl_integration_workspace* iWS; //!< Integration workspace
	
	//! Integrand for rho < sigma
	static double integrandSmallRho(double t,  void* params)
	{	const std::pair<double,double>& p = *((std::pair<double,double>*)params);
		const double& R = p.first; //R = (rho/sigma)
		const double& K = p.second; //K = (k*sigma)
		return t * exp(-0.5*t*t + t*(R-K)) * gsl_sf_bessel_I0_scaled(R*t) * gsl_sf_bessel_K0_scaled(K*t);
	}

	//! Integrand for rho > sigma
	static double integrandLargeRho(double t, void* params)
	{	const std::pair<double,double>& p = *((std::pair<double,double>*)params);
		const double& Rsq = p.first; //Rsq = (rho/sigma)^2
		const double& KR  = p.second; //KR = (k*sigma)*(rho/sigma) = k*rho
		return t * Rsq * exp(-0.5*Rsq*t*t + t*(Rsq-KR)) * gsl_sf_bessel_I0_scaled(Rsq*t) * gsl_sf_bessel_K0_scaled(KR*t);
	}
};

//! Look-up table for Cbar_k^sigma(rho) for specific values of k and sigma
struct Cbar_k_sigma
{
	Cbar_k_sigma(double k, double sigma, double rhoMax)
	{	assert(rhoMax > 0.);
		//Pick grid and initialize sample values:
		drho = 0.03*sigma; //With 5th order splines, this guarantees rel error ~ 1e-14 typical, 1e-12 max
		drhoInv = 1./drho;
		isLog = (k != 0.); //When k!=0, samples are positive and interpolate on the logarithm
		std::vector<double> x(size_t(drhoInv*rhoMax)+10);
		Cbar cbar;
		for(size_t i=0; i<x.size(); i++)
		{	double c = cbar(k, sigma, i*drho);
			if(isLog) x[i] = (c>0 ? log(c) : (i ? x[i-1] : log(DBL_MIN)));
			else x[i] = c;
		}
		coeff = QuinticSpline::getCoeff(x);
	}
	
	//! Get value:
	double value(double rho) const
	{	double f = QuinticSpline::value(coeff.data(), drhoInv * rho);
		return isLog ? exp(f) : f;
	}
	
	//! Get derivative:
	double deriv(double rho) const
	{	double fp = QuinticSpline::deriv(coeff.data(), drhoInv * rho) * drhoInv;
		return isLog ? fp * value(rho) : fp;
	}
	
private:
	friend int main();
	double drho, drhoInv; bool isLog;
	std::vector<double> coeff;
};



//! 1D Ewald sum
struct EwaldWire
{
	const GridInfo& gInfo;
	int iDir; //!< truncated direction
	const WignerSeitz& ws; //!< Wigner-Seitz cell
	bool wsTruncated; //!< true => Wigner-Seitz truncation, false => cylindrical
	double criticalDist; //!< borderWidth+ionMargin for Wigner-Seitz, Rc-ionMargin for cylindrical

	double sigma; //!< gaussian width for Ewald sums
	vector3<int> Nreal; //!< max unit cell indices for real-space sum
	vector3<int> Nrecip; //!< max unit cell indices for reciprocal-space sum
	
	std::vector<std::shared_ptr<Cbar_k_sigma>> cbar_k_sigma;
	
	EwaldWire(const GridInfo& gInfo, int iDir, const WignerSeitz& ws, bool wsTruncated, double criticalDist)
	: gInfo(gInfo), iDir(iDir), ws(ws), wsTruncated(wsTruncated), criticalDist(criticalDist)
	{	logPrintf("\n---------- Setting up 1D ewald sum ----------\n");
		//Determine optimum gaussian width for 1D Ewald sums:
		// From below, the number of reciprocal cells ~ |R.column[iDir]|
		//    and number of real space cells ~ |G.row[iDir]|
		// including the fact that a term in the reciprocal space sum
		// costs roughly 10 times as much as that in the real space sum
		sigma = sqrt(10.*gInfo.R.column(iDir).length() / gInfo.G.row(iDir).length());
		logPrintf("Optimum gaussian width for ewald sums = %lf bohr.\n", sigma);
		
		//Carry real space sums to Rmax = 10 sigma and Gmax = 10/sigma
		//This leads to relative errors ~ 1e-22 in both sums, well within double precision limits
		for(int k=0; k<3; k++)
		{	Nreal[k]  = (k!=iDir) ? 0 : 1+ceil(10. * gInfo.G.row(k).length() * sigma / (2*M_PI));
			Nrecip[k] = (k!=iDir) ? 0 : 1+ceil(10. * gInfo.R.column(k).length() / (2*M_PI*sigma));
		}
		logPrintf("Real space sums over %d unit cells with max indices ", 2*Nreal[iDir]+1);
		Nreal.print(globalLog, " %d ");
		logPrintf("Reciprocal space sums over %d terms with max indices ", Nrecip[iDir]+1);
		Nrecip.print(globalLog, " %d ");
		
		//Initialize Cbar_k^sigma look-up tables:
		cbar_k_sigma.resize(Nrecip[iDir]+1);
		vector3<int> iG(0,0,0);
		double rhoMax = ws.circumRadius(iDir);
		for(iG[iDir]=0; iG[iDir]<=Nrecip[iDir]; iG[iDir]++)
		{	double k = sqrt(gInfo.GGT.metric_length_squared(iG));
			cbar_k_sigma[iG[iDir]] = std::make_shared<Cbar_k_sigma>(k, sigma, rhoMax);
		}
	}
	
	double energyAndGrad(std::vector<Atom>& atoms) const
	{	if(!atoms.size()) return 0.;
		double eta = sqrt(0.5)/sigma, etaSq=eta*eta;
		//Position independent terms: (Self-energy correction)
		double ZsqTot = 0.;
		for(const Atom& a: atoms)
			ZsqTot += a.Z * a.Z;
		double E = -0.5 * ZsqTot * eta * (2./sqrt(M_PI));
		//Reduce positions to first unit cell:
		//Shift all points in the truncated directions into the 2D Wigner-Seitz cell
		//centered on one of the atoms; choice of this atom is irrelevant if every atom
		//lies in the WS cell of the other with a consistent translation:
		vector3<> pos0 = atoms[0].pos;
		for(Atom& a: atoms)
			a.pos = pos0 + ws.restrict(a.pos - pos0);
		//Real space sum:
		vector3<int> iR(0,0,0); //integer cell number
		for(const Atom& a2: atoms)
			for(Atom& a1: atoms)
				for(iR[iDir]=-Nreal[iDir]; iR[iDir]<=Nreal[iDir]; iR[iDir]++)
				{	vector3<> x = iR + (a1.pos - a2.pos);
					double rSq = gInfo.RTR.metric_length_squared(x);
					if(!rSq) continue; //exclude self-interaction
					double r = sqrt(rSq);
					E += 0.5 * a1.Z * a2.Z * erfc(eta*r)/r;
					a1.force += (gInfo.RTR * x) *
						(a1.Z * a2.Z * (erfc(eta*r)/r + (2./sqrt(M_PI))*eta*exp(-etaSq*rSq))/rSq);
				}
		//Reciprocal space sum:
		double volPrefac = 0.5 / sqrt(gInfo.RTR(iDir,iDir));
		for(unsigned i1=0; i1<atoms.size(); i1++)
		{	Atom& a1 = atoms[i1];
			for(unsigned i2=0; i2<=i1; i2++)
			{	Atom& a2 = atoms[i2];
				double prefac = volPrefac * a1.Z * a2.Z * (i1==i2 ? 1 : 2);
				vector3<> r12 = a1.pos - a2.pos;
				vector3<> rho12vec = r12; rho12vec[iDir] = 0.; //projected to truncation plane
				double rho12 = sqrt(gInfo.RTR.metric_length_squared(rho12vec));
				if(wsTruncated)
				{	if(ws.boundaryDistance(rho12vec, iDir) <= criticalDist)
						die("Separation between atoms %d and %d lies in the truncation border + margin.\n", i1, i2);
				}
				else
				{	if(rho12 >= criticalDist)
						die("Atoms %d and %d are separated by rho = %lg >= Rc-ionMargin = %lg bohrs.\n", i1, i2, rho12, criticalDist);
				}
				double E12 = 0.; vector3<> E12_r12(0.,0.,0.); //energy and gradient from this pair
				vector3<int> iG(0,0,0); //integer reciprocal cell number (only iG[iDir] will be non-zero)
				for(iG[iDir]=0; iG[iDir]<=Nrecip[iDir]; iG[iDir]++)
				{	//1D structure factor term and derivative
					double c, s; sincos((2*M_PI)*dot(iG,r12), &s, &c);
					if(iG[iDir]) { c *= 2.; s *= 2.; } //include contribution from -iG[iDir]
					//Contribution from truncated directions:
					double rhoTerm = cbar_k_sigma[iG[iDir]]->value(rho12);
					double rhoTermPrime = cbar_k_sigma[iG[iDir]]->deriv(rho12);
					//Update energy and forces:
					E12 += prefac * c * rhoTerm;
					E12_r12 += (prefac * -s * rhoTerm * (2*M_PI)) * iG;
						+ (prefac * c * rhoTermPrime * (rho12 ? 1./rho12 : 0.)) * rho12vec;
				}
				E += E12;
				a1.force -= E12_r12;
				a2.force += E12_r12;
			}
		}
		return E;
	}
};


//----------------- class CoulombWire ---------------------

//Threaded initialization for CoulombWire
struct CoulombWire_init
{	//Data arrays:
	complex* padArr; double* padRealArr; //fft array on padded grid
	complex* denseArr; double* denseRealArr; //fft array on dense (unpadded) grid
	double* Vc; //output kernel (3D fftw c2r layout)
	fftw_plan fftPlanC2R, fftPlanR2C;
	//Geometry:
	int iDir, jDir, kDir;
	vector3<int> S, Spad, Sdense;
	matrix3<> Rplanar, RpadPlanar, GGT;
	std::vector<Simplex<2>>* simplexArr; //2D Wigner-Seitz cell simplicial tesselation
	WignerSeitz *wsDense, *wsPad; //2D wigner-Seitz cells
	double sigma, sigmaBorder;
	
	void computePlane(int iPlane)
	{
		double kz = iPlane * (2*M_PI/Rplanar.column(iDir).length());
		//Short-cut numerical truncation if K_0(kz rho) becomes zero well inside WS cell:
		vector3<int> pitch;
		pitch[2] = 1;
		pitch[1] = pitch[2] * (1 + S[2]/2);
		pitch[0] = pitch[1] * S[1];
		if(kz * (wsDense->inRadius(iDir) - 10.*sigmaBorder) > 50.) //=> K_0 < 1e-22
		{	vector3<int> iG; iG[iDir] = iPlane;
			for(iG[jDir]=1-S[jDir]/2; iG[jDir]<=S[jDir]/2; iG[jDir]++)
				for(iG[kDir]=1-S[kDir]/2; iG[kDir]<=S[kDir]/2; iG[kDir]++)
				{	double curV = (4*M_PI) / GGT.metric_length_squared(iG);
					//Save to the appropriate locations in Vc (2 sign combinations):
					for(int si=0; si<2; si++)
					{	vector3<int> iv = iG;
						if(si && iv[iDir]) { iv[iDir] = S[iDir] - iv[iDir]; }
						if(iv[jDir] < 0) iv[jDir] += S[jDir];
						if(iv[kDir] < 0) iv[kDir] += S[kDir];
						if(iv[2] <= S[2]/2)
							Vc[dot(iv, pitch)] = curV;
					}
				}
			return;
		}
		
		//Compute smoothed WignerSeitz-shaped theta function in Fourier space:
		const double invApad = RpadPlanar.column(iDir).length() / fabs(det(RpadPlanar));
		vector3<int> iGpad(0,0,0);
		matrix3<> GTpadPlanar = (2*M_PI) * ~inv(RpadPlanar);
		complex* theta = padArr;
		for(iGpad[jDir]=0;;)
		{	for(iGpad[kDir]=0; iGpad[kDir]<=Spad[kDir]/2; iGpad[kDir]++)
			{	vector3<> G = GTpadPlanar * iGpad; //reciprocal lattice vector in planar cartesian coords
				Simplex<2>::Point Gpoint({{ G[0], G[1] }}); //convert to Simplex<2>::Point (third direction is truncated)
				double curTheta = 0.;
				for(const Simplex<2>& simplex: *simplexArr)
					curTheta += simplex.getTilde(Gpoint);
				*(theta++) = invApad * curTheta * exp(-0.5*G.length_squared()*sigmaBorder*sigmaBorder);
			}
			iGpad[jDir]++;
			if(2*iGpad[jDir]>Spad[jDir]) iGpad[jDir] -= Spad[jDir];
			if(iGpad[jDir]==0) break;
		}
		fftw_execute_dft_c2r(fftPlanC2R, (fftw_complex*)padArr, padRealArr);
		
		//Multiply by Cbar_k_sigma(rho) in real-space and fold into dense array:
		int pitchPad = 2*(1+Spad[kDir]/2);
		int pitchDense = 2*(1+Sdense[kDir]/2);
		vector3<> invSpad, invSdense;
		for(int k=0; k<3; k++)
		{	invSpad[k] = 1./Spad[k];
			invSdense[k] = 1./Sdense[k];
		}
		matrix3<> hPad; //mesh offset vectors
		for(int k=0; k<3; k++)
			hPad.set_col(k, RpadPlanar.column(k) / Spad[k]);
		matrix3<> hThPad = (~hPad) * hPad; //metric in mesh coordinates
		double dA = fabs(det(Rplanar)) / (Rplanar.column(iDir).length() * Sdense[jDir] * Sdense[kDir]);
		Cbar_k_sigma cbar_k_sigma(kz, sigma, wsPad->circumRadius(iDir));
		vector3<int> iv(0,0,0);
		memset(denseArr, 0, sizeof(complex)*Sdense[jDir]*(1+Sdense[kDir]/2));
		for(iv[jDir]=0; iv[jDir]<Spad[jDir]; iv[jDir]++)
			for(iv[kDir]=0; iv[kDir]<Spad[kDir]; iv[kDir]++)
			{	//Compute index mappings:
				vector3<int> ivPad = wsPad->restrict(iv, Spad, invSpad); //position in mesh coordinates within padded WignerSeitz cell
				vector3<int> ivDense = wsDense->restrict(ivPad, Sdense, invSdense); //position in mesh coordinates within orig WS cell
				for(int k=0; k<3; k++)
				{	ivDense[k] = ivDense[k] % Sdense[k];
					if(ivDense[k]<0) ivDense[k] += Sdense[k];
				}
				int iPad = iv[kDir] + pitchPad * iv[jDir];
				int iDense = ivDense[kDir] + pitchDense * ivDense[jDir];
				//Accumulate theta multiplied by erf/r to the mapped position:
				double rho = sqrt(hThPad.metric_length_squared(ivPad)); //distance of minmal periodic image from origin
				denseRealArr[iDense] += dA * cbar_k_sigma.value(rho) * padRealArr[iPad];
			}
		fftw_execute_dft_r2c(fftPlanR2C, denseRealArr, (fftw_complex*)denseArr);
		
		//Add analytic short-ranged parts in fourier space (and down-sample to final resolution):
		pitchDense = 1+Sdense[kDir]/2;
		vector3<int> iG; iG[iDir] = iPlane;
		for(iG[jDir]=1-S[jDir]/2; iG[jDir]<=S[jDir]/2; iG[jDir]++)
			for(iG[kDir]=0; iG[kDir]<=S[kDir]/2; iG[kDir]++)
			{	//Collect the data from the dense grid transform:
				double curV = denseArr[iG[kDir] + pitchDense*(iG[jDir]<0 ? iG[jDir]+Sdense[jDir] : iG[jDir])].real();
				//Add the analytical short-ranged part:
				double Gsq = GGT.metric_length_squared(iG);
				curV += (4*M_PI) * (Gsq ? (1.-exp(-0.5*sigma*sigma*Gsq))/Gsq : 0.5*sigma*sigma);
				//Save to the appropriate locations in Vc (4 sign combinations):
				for(int si=0; si<2; si++)
				{	for(int sk=0; sk<2; sk++)
					{	vector3<int> iv = iG;
						if(si && iv[iDir]) { iv[iDir] = S[iDir] - iv[iDir]; }
						if(sk && iv[kDir]) { iv[kDir] = S[kDir] - iv[kDir]; iv[jDir] = -iv[jDir]; }
						if(iv[jDir] < 0) iv[jDir] += S[jDir];
						if(iv[2] <= S[2]/2)
							Vc[dot(iv, pitch)] = curV;
					}
				}
			}
	}
	
	static void thread(int iThread, int nThreads, CoulombWire_init* cwInitArr,
		int nPlanes, int* nPlanesDone, std::mutex* m)
	{
		while(true)
		{	//Get next available job:
			m->lock();
			int iPlane = (*nPlanesDone)++;
			m->unlock();
			if(iPlane >= nPlanes)
				break; //job queue empty
			//Perform job:
			cwInitArr[iThread].computePlane(iPlane);
		}
	}
};

CoulombWire::CoulombWire(const GridInfo& gInfo, const CoulombParams& params)
: Coulomb(gInfo, params), ws(gInfo.R), Vc(gInfo)
{	//Check orthogonality
	string dirName = checkOrthogonality(gInfo, params.iDir);
	
	//Read precomputed kernel from file if supplied
	if(params.filename.length())
	{	FILE* fp = fopen(params.filename.c_str(), "rb");
		if(fp)
		{	matrix3<> R; int iDir; vector3<int> S; double bw;
			fread(&R, sizeof(matrix3<>), 1, fp);
			fread(&iDir, sizeof(int), 1, fp);
			fread(&S, sizeof(vector3<int>), 1, fp);
			fread(&bw, sizeof(double), 1, fp);
			#define CHECKerror(errorCondition, paramName) \
				else if(errorCondition) \
					logPrintf("Precomputed coulomb kernel file '%s' has different " paramName " (recomputing it now)\n", params.filename.c_str());
			if(false) {} //dummy condition to start a chain of else-if's
			CHECKerror(R != gInfo.R, "lattice vectors")
			CHECKerror(iDir != params.iDir, "truncation direction")
			CHECKerror(!(S == gInfo.S), "sample count")
			CHECKerror(bw != params.borderWidth, "border width")
			else if(fread(Vc.data, sizeof(double), gInfo.nG, fp) != unsigned(gInfo.nG))
				logPrintf("Error reading precomputed coulomb kernel from '%s' (computing it now)\n", params.filename.c_str());
			else
			{	logPrintf("Successfully read precomputed coulomb kernel from '%s'\n", params.filename.c_str());
				Vc.set();
				initExchangeEval();
				return;
			}
			#undef CHECKparam
		}
		else logPrintf("Could not open precomputed coulomb kernel file '%s' (computing it now)\n", params.filename.c_str());
	}
	
	//Direction indices:
	int iDir = params.iDir;
	int jDir = (iDir + 1) % 3;
	int kDir = (iDir + 2) % 3;
	
	//Select gauss-smoothing parameter:
	double maxBorderWidth = 0.5 * ws.inRadius(iDir);
	if(params.borderWidth > maxBorderWidth)
		die("Border width %lg bohrs must be less than %lg bohrs (half the Wigner-Seitz cell in-radius).\n",
			params.borderWidth, maxBorderWidth);
	double sigmaBorder = 0.1 * params.borderWidth;
	double sigma = 0.1*ws.inRadius(iDir) - sigmaBorder; //so that 10(sigma+sigmaBorder) < inRadius
	logPrintf("Selecting gaussian width %lg bohrs (for border width %lg bohrs).\n", sigmaBorder, params.borderWidth);

	//Set up dense integration grids:
	logPrintf("Setting up FFT grids: ");
	double Gnyq = 10./sigmaBorder; //lower bound on Nyquist frequency
	vector3<int> Sdense; //dense fft sample count
	matrix3<> Rpad; vector3<int> Spad; //padded lattice vectors and sample count
	for(int k=0; k<3; k++)
		if(k == iDir)
		{	Sdense[k] = gInfo.S[k];
			Spad[k] = gInfo.S[k];
			Rpad.set_col(k, gInfo.R.column(k));
		}
		else
		{	Sdense[k] = std::max(gInfo.S[k], 2*int(ceil(Gnyq * gInfo.R.column(k).length() / (2*M_PI))));
			while(!fftSuitable(Sdense[k])) Sdense[k]+=2; //pick the next even number suitable for FFT
			//Pad the super cell by border width:
			Spad[k] = Sdense[k] + 2*ceil(Sdense[k]*params.borderWidth/gInfo.R.column(k).length());
			while(!fftSuitable(Spad[k])) Spad[k]+=2;
			Rpad.set_col(k, gInfo.R.column(k) * (Spad[k]*1./Sdense[k]));
		}
	logPrintf("%d x %d, and %d x %d padded.\n", Sdense[jDir], Sdense[kDir], Spad[jDir], Spad[kDir]);
	int nGdense = Sdense[jDir] * (1+Sdense[kDir]/2);
	int nGpad = Spad[jDir] * (1+Spad[kDir]/2); //number of symmetry reduced planar reciprocal lattice vectors
	int nrPad = Spad[jDir] * Spad[kDir]; //number of planar real space points
	
	//Construct padded Wigner-Seitz cell, and check border:
	logPrintf("For padded lattice, "); WignerSeitz wsPad(Rpad);
	std::vector<vector3<>> vArr = ws.getVertices();
	matrix3<> invRpad = inv(Rpad);
	for(vector3<> v: vArr)
		if(wsPad.boundaryDistance(wsPad.restrict(invRpad*v), iDir) < 0.9*params.borderWidth)
			die("\nPadded Wigner-Seitz cell does not fit inside Wigner-Seitz cell of padded lattice.\n"
				"This can happen for reducible lattice vectors; HINT: the reduced lattice vectors,\n"
				"if different from the input lattice vectors, are printed during Symmetry setup.\n");
	
	//Plan Fourier transforms:
	assert(nrPad > 0); //overflow check
	complex* tempArr = (complex*)fftw_malloc(sizeof(complex)*nGpad);
	#define INSUFFICIENT_MEMORY_ERROR \
		die("Insufficient memory (need %.1fGB). Hint: try increasing border width.\n", \
			(nGpad+nGdense)*1e-9*nProcsAvailable*sizeof(complex));
	if(!tempArr) INSUFFICIENT_MEMORY_ERROR
	logPrintf("Planning fourier transforms ... "); logFlush();
	fftw_plan_with_nthreads(1); //Multiple simultaneous single threaded fourier transforms
	fftw_plan fftPlanC2R = fftw_plan_dft_c2r_2d(Spad[jDir], Spad[kDir], (fftw_complex*)tempArr, (double*)tempArr, FFTW_ESTIMATE);
	fftw_plan fftPlanR2C = fftw_plan_dft_r2c_2d(Sdense[jDir], Sdense[kDir], (double*)tempArr, (fftw_complex*)tempArr, FFTW_ESTIMATE);
	fftw_free(tempArr);
	logPrintf("Done.\n");
	
	//Setup threads for initializing each plane perpendicular to truncated direction:
	logPrintf("Computing wire[%s]-truncated coulomb kernel ... ", dirName.c_str()); logFlush();
	std::vector<CoulombWire_init> cwInitArr(nProcsAvailable);
	std::vector<Simplex<2>> simplexArr = ws.getSimplices(iDir);
	for(CoulombWire_init& c: cwInitArr)
	{	//Allocate data arrays:
		c.padArr = (complex*)fftw_malloc(sizeof(complex)*nGpad);
		c.denseArr = (complex*)fftw_malloc(sizeof(complex)*nGdense);
		if(!c.padArr || !c.denseArr) INSUFFICIENT_MEMORY_ERROR
		#undef INSUFFICIENT_MEMORY_ERROR
		c.padRealArr = (double*)c.padArr;
		c.denseRealArr = (double*)c.denseArr;
		c.Vc = Vc.data;
		c.fftPlanC2R = fftPlanC2R;
		c.fftPlanR2C = fftPlanR2C;
		//Copy geometry definitions:
		c.iDir = iDir; c.jDir = jDir; c.kDir = kDir;
		c.S = gInfo.S; c.Sdense = Sdense; c.Spad = Spad;
		c.Rplanar = ws.getRplanar(iDir);
		c.RpadPlanar = wsPad.getRplanar(iDir);
		c.GGT = gInfo.GGT;
		c.simplexArr = &simplexArr;
		c.wsDense = &ws; c.wsPad = &wsPad;
		c.sigma = sigma; c.sigmaBorder = sigmaBorder;
	}
	
	//Launch threads:
	std::mutex mJobCount; int nPlanesDone = 0; //for job management
	threadLaunch(CoulombWire_init::thread, 0, cwInitArr.data(),
		1+gInfo.S[iDir]/2, &nPlanesDone, &mJobCount);
	Vc.set();
	fftw_destroy_plan(fftPlanC2R);
	fftw_destroy_plan(fftPlanR2C);
	
	//Cleanup threads:
	for(CoulombWire_init& c: cwInitArr)
	{	fftw_free(c.padArr);
		fftw_free(c.denseArr);
	}
	logPrintf("Done.\n");
	
	//Save kernel if requested:
	if(params.filename.length())
	{	logPrintf("Saving wire[%s]-truncated coulomb kernel to '%s' ... ", dirName.c_str(), params.filename.c_str()); logFlush();
		FILE* fp = fopen(params.filename.c_str(), "wb");
		if(!fp) die("could not open file for writing.\n");
		fwrite(&gInfo.R, sizeof(matrix3<>), 1, fp);
		fwrite(&params.iDir, sizeof(int), 1, fp);
		fwrite(&gInfo.S, sizeof(vector3<int>), 1, fp);
		fwrite(&params.borderWidth, sizeof(double), 1, fp);
		fwrite(Vc.data, sizeof(double), gInfo.nG, fp);
		fclose(fp);
		logPrintf("Done.\n");
	}
	initExchangeEval();
}

DataGptr CoulombWire::operator()(DataGptr&& in) const
{	return Vc * in;
}

double CoulombWire::energyAndGrad(std::vector<Atom>& atoms) const
{	if(!ewald)
		((CoulombWire*)this)->ewald = std::make_shared<EwaldWire>(gInfo, params.iDir, ws, true, params.borderWidth + params.ionMargin);
	return ewald->energyAndGrad(atoms);
}

//----------------- class CoulombCylindrical ---------------------

void setVcylindrical(int iStart, int iStop, vector3<int> S, const matrix3<> GGT, int iDir, double Rc, double* Vc)
{	THREAD_halfGspaceLoop
	(	double Gsq = GGT.metric_length_squared(iG);
		double GaxisSq = GGT(iDir,iDir) * iG[iDir] * iG[iDir];
		double GplaneSq = Gsq - GaxisSq;
		double RGaxis = GaxisSq>0. ? Rc*sqrt(GaxisSq) : 0.; //safe sqrt to prevent NaN from roundoff errors
		double RGplane = GplaneSq>0. ? Rc*sqrt(GplaneSq) : 0.; //safe sqrt to prevent NaN from roundoff errors
		if(iG[iDir])
			Vc[i] = (4*M_PI/Gsq) * (1.
				+ RGplane * gsl_sf_bessel_J1(RGplane) * gsl_sf_bessel_K0(RGaxis)
				- RGaxis * gsl_sf_bessel_J0(RGplane) * gsl_sf_bessel_K1(RGaxis) );
		else
		{	Vc[i] = GplaneSq
				? (4*M_PI/GplaneSq) * (1. - gsl_sf_bessel_J0(RGplane)
					- (RGplane) * gsl_sf_bessel_J1(RGplane) * log(Rc) )
				: M_PI*Rc*Rc * (1. - 2*log(Rc));
		}
	)
}

CoulombCylindrical::CoulombCylindrical(const GridInfo& gInfo, const CoulombParams& params)
: Coulomb(gInfo, params), ws(gInfo.R), Rc(params.Rc), Vc(gInfo)
{	//Check orthogonality:
	string dirName = checkOrthogonality(gInfo, params.iDir);
	//Check the truncation radius:
	double RcMax = ws.inRadius(params.iDir);
	if(Rc > RcMax)
		die("Cylindrical truncation radius %lg exceeds 2D Wigner-Seitz cell in-radius of %lg bohrs.\n", Rc, RcMax);
	if(!Rc) Rc = RcMax;
	//Set the kernel:
	threadLaunch(setVcylindrical, gInfo.nG, gInfo.S, gInfo.GGT, params.iDir, Rc, Vc.data);
	Vc.set();
	logPrintf("Initialized cylindrical truncation of radius %lg bohrs with axis along lattice direction %s\n", Rc, dirName.c_str());
	initExchangeEval();
}

DataGptr CoulombCylindrical::operator()(DataGptr&& in) const
{	return Vc * in;
}

double CoulombCylindrical::energyAndGrad(std::vector<Atom>& atoms) const
{	if(!ewald)
		((CoulombCylindrical*)this)->ewald = std::make_shared<EwaldWire>(gInfo, params.iDir, ws, false, Rc - params.ionMargin);
	return ewald->energyAndGrad(atoms);
}
