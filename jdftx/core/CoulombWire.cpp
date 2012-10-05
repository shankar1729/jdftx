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
#include <core/CoulombKernel.h>
#include <core/Coulomb_internal.h>
#include <core/Operators.h>
#include <core/Util.h>
#include <core/Bspline.h>
#include <core/LoopMacros.h>
#include <gsl/gsl_sf.h>

//Check orthogonality and return lattice direction name
string checkOrthogonality(const GridInfo& gInfo, int iDir)
{	string dirName(3, '0'); dirName[iDir] = '1';
	if((!WignerSeitz::isOrthogonal(gInfo.R.column(iDir),gInfo.R.column((iDir+1)%3)))
	|| (!WignerSeitz::isOrthogonal(gInfo.R.column(iDir),gInfo.R.column((iDir+2)%3))) )
		die("Lattice direction %s is not perpendicular to the other two basis vectors.\n", dirName.c_str());
	return dirName;
}



//--------------- class Cbar ----------

Cbar::Cbar() { iWS = gsl_integration_workspace_alloc(maxIntervals); }
Cbar::~Cbar() { gsl_integration_workspace_free(iWS); }
	
//Compute Cbar_k^sigma(rho)
double Cbar::operator()(double k, double sigma, double rho)
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

//Integrand for rho < sigma
double Cbar::integrandSmallRho(double t,  void* params)
{	const std::pair<double,double>& p = *((std::pair<double,double>*)params);
	const double& R = p.first; //R = (rho/sigma)
	const double& K = p.second; //K = (k*sigma)
	return t * exp(-0.5*t*t + t*(R-K)) * gsl_sf_bessel_I0_scaled(R*t) * gsl_sf_bessel_K0_scaled(K*t);
}

//Integrand for rho > sigma
double Cbar::integrandLargeRho(double t, void* params)
{	const std::pair<double,double>& p = *((std::pair<double,double>*)params);
	const double& Rsq = p.first; //Rsq = (rho/sigma)^2
	const double& KR  = p.second; //KR = (k*sigma)*(rho/sigma) = k*rho
	return t * Rsq * exp(-0.5*Rsq*t*t + t*(Rsq-KR)) * gsl_sf_bessel_I0_scaled(Rsq*t) * gsl_sf_bessel_K0_scaled(KR*t);
}



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
					E12_r12 += (prefac * -s * rhoTerm * (2*M_PI)) * iG
						+ (prefac * c * rhoTermPrime * (rho12 ? 1./rho12 : 0.)) * (gInfo.RTR * rho12vec);
				}
				E += E12;
				a1.force -= E12_r12;
				a2.force += E12_r12;
			}
		}
		return E;
	}
};



CoulombWire::CoulombWire(const GridInfo& gInfo, const CoulombParams& params)
: Coulomb(gInfo, params), ws(gInfo.R), Vc(gInfo)
{	//Check orthogonality
	string dirName = checkOrthogonality(gInfo, params.iDir);
	
	//Select gauss-smoothing parameter:
	double maxBorderWidth = sqrt(0.5) * ws.inRadius(params.iDir);
	if(params.borderWidth > maxBorderWidth)
		die("Border width %lg bohrs must be less than %lg bohrs (Wigner-Seitz cell in-radius/sqrt(2)).\n",
			params.borderWidth, maxBorderWidth);
	double sigmaBorder = params.borderWidth / CoulombKernelDesc::nSigmasPerWidth;
	logPrintf("Selecting gaussian width %lg bohrs (for border width %lg bohrs).\n", sigmaBorder, params.borderWidth);

	//Create kernel description:
	vector3<bool> isTruncated(true, true, true); isTruncated[params.iDir] = false;
	vector3<> sigmaBorders(sigmaBorder, sigmaBorder, sigmaBorder);
	CoulombKernelDesc kernelDesc(gInfo.R, gInfo.S, isTruncated, sigmaBorders);
	//Load or compute the kernel:
	kernelDesc.computeKernel(Vc.data, ws, params.filename);
	Vc.set();
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

void setVcylindrical(size_t iStart, size_t iStop, vector3<int> S, const matrix3<> GGT, int iDir, double Rc, double* Vc)
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
