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

#include <core/CoulombPeriodic.h>
#include <core/CoulombSlab.h>
#include <core/CoulombWire.h>
#include <core/CoulombIsolated.h>
#include <core/Coulomb_internal.h>
#include <core/LatticeUtils.h>
#include <core/LoopMacros.h>
#include <core/BlasExtra.h>
#include <core/Thread.h>

std::shared_ptr<Coulomb> CoulombParams::createCoulomb(const GridInfo& gInfo) const
{	if(geometry != Periodic)
		logPrintf("\n---------- Setting up coulomb interaction ----------\n");
	switch(geometry)
	{	case Periodic:    return std::make_shared<CoulombPeriodic>(gInfo, *this);
		case Slab:        return std::make_shared<CoulombSlab>(gInfo, *this);
		case Wire:        return std::make_shared<CoulombWire>(gInfo, *this);
		case Cylindrical: return std::make_shared<CoulombCylindrical>(gInfo, *this);
		case Isolated:    return std::make_shared<CoulombIsolated>(gInfo, *this);
		case Spherical:   return std::make_shared<CoulombSpherical>(gInfo, *this);
		default: return 0; //never encountered (to suppress warning)
	}
}



//! Helper class for evaluating regularized Coulomb kernel for exchange
struct ExchangeEval
{
	ExchangeEval(const GridInfo& gInfo, const CoulombParams& params, const Coulomb& coulomb, double omega);
	~ExchangeEval();
	complexDataGptr operator()(complexDataGptr&& in, vector3<> kDiff) const;

private:
	const GridInfo& gInfo;
	const CoulombParams& params;
	const Coulomb& coulomb;
	double omega;
	
	//Shorthand for combinations of regularization method and geometry
	enum KernelMode
	{	PeriodicKernel, //regularization = None or AuxiliaryFunction, with geometry = Periodic
		SphericalKernel, //regularization = SphericalTruncated, or regularization = None and geometry = Spherical
		WignerSeitzKernel //regularization = WignerSeitzTruncated, or regularization = None and geometry = Isolated
	} kernelMode;
	
	//For all analytic modes (None, AuxiliaryFunction, SphericalTruncated):
	double Vzero; //G=0 term (0, the auxiliary correction, and the G=0 limit respectively)
	//For spherical mode
	double Rc;
	double* coeffErfcTilde; //quintic spline coefficients for radial fourier transform of truncated erfc/r
	double dGinvErfcTilde; //inverse of sample spacing in coeffErfcTilde
	size_t nSamplesErfcTilde; //number of samples in coeffErfcTilde spline
};



//--------------- class Coulomb ----------------

DataGptr Coulomb::operator()(const DataGptr& in) const
{	DataGptr out(in->clone()); //create destructible copy
	return (*this)((DataGptr&&)out);
}

complexDataGptr Coulomb::operator()(complexDataGptr&& in, vector3<> kDiff, double omega) const
{	auto exEvalOmega = exchangeEval.find(omega);
	assert(exEvalOmega != exchangeEval.end());
	return (*exEvalOmega->second)((complexDataGptr&&)in, kDiff);
}

complexDataGptr Coulomb::operator()(const complexDataGptr& in, vector3<> kDiff, double omega) const
{	complexDataGptr out(in->clone()); //create destructible copy
	return (*this)((complexDataGptr&&)out, kDiff, omega);
}

Coulomb::Coulomb(const GridInfo& gInfo, const CoulombParams& params)
: gInfo(gInfo), params(params)
{
}

void Coulomb::initExchangeEval()
{	//Initialize Exchange evaluators if required
	for(double omega: params.omegaSet)
		exchangeEval[omega]  = std::make_shared<ExchangeEval>(gInfo, params, *this, omega);
}


//-------- CPU implementation of Coulomb_internal.h --------

template<typename Coulomb_calc>
void coulombAnalytic_thread(size_t iStart, size_t iStop, vector3<int> S, const matrix3<>& GGT, const Coulomb_calc& calc, complex* data)
{	THREAD_halfGspaceLoop
	(	data[i] *= calc(iG, GGT);
	)
}
#define DECLARE_coulombAnalytic(Type) \
	void coulombAnalytic(vector3<int> S, const matrix3<>& GGT, const Coulomb##Type##_calc& calc, complex* data) \
	{	threadLaunch(coulombAnalytic_thread<Coulomb##Type##_calc>, S[0]*S[1]*(1+S[2]/2), S, GGT, calc, data); \
	}
DECLARE_coulombAnalytic(Periodic)
DECLARE_coulombAnalytic(Slab)
DECLARE_coulombAnalytic(Spherical)
#undef DECLARE_coulombAnalytic


template<typename Exchange_calc>
void exchangeAnalytic_thread(size_t iStart, size_t iStop, vector3<int> S, const matrix3<>& GGT, const Exchange_calc& calc,
	complex* data, const vector3<>& kDiff, double Vzero, double thresholdSq)
{	THREAD_fullGspaceLoop
	(	double kplusGsq = GGT.metric_length_squared(iG + kDiff);
		data[i] *= kplusGsq<thresholdSq ? Vzero : calc(kplusGsq);
	)
}
#define DECLARE_exchangeAnalytic(Type) \
	void exchangeAnalytic(vector3<int> S, const matrix3<>& GGT, const Exchange##Type##_calc& calc, \
		complex* data, const vector3<>& kDiff, double Vzero, double thresholdSq) \
	{	\
		threadLaunch(exchangeAnalytic_thread<Exchange##Type##_calc>, \
			S[0]*S[1]*S[2], S, GGT, calc, data, kDiff, Vzero, thresholdSq); \
	}
DECLARE_exchangeAnalytic(Periodic)
DECLARE_exchangeAnalytic(PeriodicScreened)
DECLARE_exchangeAnalytic(Spherical)
DECLARE_exchangeAnalytic(SphericalScreened)
#undef DECLARE_exchangeAnalytic


//--------------- class ExchangeEval ----------------


//-------------- Auxiliary Function method --------------------------
// [ P. Carrier, S. Rohra and A. Gorling, Phys. Rev. B 75, 205126 (2007) ]
// Reciprocal space function (with correct periodicity) singular at G=0 as 4 pi/G^2
// g is a dimensionless coordinate (coefficients to reciprocal lattice vectors)
//   When omega!=0, the G=0 behavior is (1-exp(-omega^2 G^2/4))/G^2, which is not singular,
//   but still benefits from this correction for brillouin zones larger than omega
inline double fSingular(const vector3<>& g, const matrix3<>& GGT, double omegaSq)
{	vector3<> sinPi; for(int i=0; i<3; i++) sinPi[i] = sin(M_PI*g[i]);
	vector3<> sin2Pi; for(int i=0; i<3; i++) sin2Pi[i] = sin(2*M_PI*g[i]);
	double effGsq =  (1./(M_PI*M_PI)) * //Periodic function from ref. which is >0 and ~ G^2 near G=0
		( sinPi[0]*sinPi[0]*GGT(0,0) + sinPi[1]*sinPi[1]*GGT(1,1) + sinPi[2]*sinPi[2]*GGT(2,2)
		+ 0.5*(sin2Pi[0]*sin2Pi[1]*GGT(0,1) + sin2Pi[1]*sin2Pi[2]*GGT(1,2) + sin2Pi[2]*sin2Pi[0]*GGT(2,0)));
	return (omegaSq ? (1.-exp(-0.25*effGsq/omegaSq)) : 1.) * (4*M_PI) / effGsq;
}
//Integrate fSingular on a box of size scale and center gCenter (in reciprocal lattice coordinates)
double fSingularIntegralBox(double scale, const vector3<> gCenter, const matrix3<>& GGT, double omegaSq)
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
	double ret = 0.0;
	double h = 0.5 * scale; 
	vector3<> g;
	for(int i0=-N; i0<=N; i0++)
	{	g[0] = gCenter[0] + h*x[N-abs(i0)]*(i0>0?1:-1);
		double w0 = w[N-abs(i0)];
		for(int i1=-N; i1<=N; i1++)
		{	g[1] = gCenter[1] + h*x[N-abs(i1)]*(i1>0?1:-1);
			double w01 = w0 * w[N-abs(i1)];
			for(int i2=-N; i2<=N; i2++)
			{	g[2] = gCenter[2] + h*x[N-abs(i2)]*(i2>0?1:-1);
				double w012 = w01 * w[N-abs(i2)];
				ret += w012 * fSingular(g, GGT, omegaSq);
			}
		}
	}
	return h*h*h * ret;
}
//Integrate fSingular between the box of size scale and scale/3 centered at the origin (in reciprocal lattice coordinates)
double fSingularIntegralBoxDiff(double scale, const matrix3<>& GGT, double omegaSq)
{	double scaleBy3 = scale/3.;
	double ret = 0.0;
	vector3<int> ig;
	for(ig[0]=-1; ig[0]<=1; ig[0]++)
	for(ig[1]=-1; ig[1]<=1; ig[1]++)
	for(ig[2]=-1; ig[2]<=1; ig[2]++)
		if(ig.length_squared()) //except the center box
			ret += fSingularIntegralBox(scaleBy3, scaleBy3*ig, GGT, omegaSq);
	return ret;
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



ExchangeEval::ExchangeEval(const GridInfo& gInfo, const CoulombParams& params, const Coulomb& coulomb, double omega)
: gInfo(gInfo), params(params), coulomb(coulomb), omega(omega)
{
	if(!omega) logPrintf("\n-------- Setting up exchange kernel --------\n");
	else logPrintf("\n--- Setting up screened exchange kernel (omega = %lg) ---\n", omega);

	const ExchangeRegularization& exReg = params.exchangeRegularization;
	
	if(exReg.method==ExchangeRegularization::None)
	{	assert(params.geometry==CoulombParams::Periodic
			|| params.geometry==CoulombParams::Isolated
			|| params.geometry==CoulombParams::Spherical);
	}
	
	if(exReg.method==ExchangeRegularization::None && params.geometry==CoulombParams::Periodic)
	{	Vzero = omega ? M_PI/(omega*omega) : 0.;
		kernelMode = PeriodicKernel;
		logPrintf("No singularity correction.\n");
	}
	
	if(exReg.method==ExchangeRegularization::AuxiliaryFunction)
	{	assert(params.geometry==CoulombParams::Periodic);
		double omegaSq = omega*omega;
		//Compute the integral of fSingular over the brillouin zone
		double fSingularIntegral = 0.;
		for(double scale=1.0; scale>1e-16; scale/=3.)
			fSingularIntegral += fSingularIntegralBoxDiff(scale, gInfo.GGT, omegaSq);
		//Compute the integral of fSingular as sampled by the k-point mesh:
		double fSingularSum = 0.;
		const std::vector<vector3<>>& kmesh = params.supercell->kmesh;
		for(unsigned i=1; i<kmesh.size(); i++) //ignore dk=0
			fSingularSum += fSingular(kmesh[i]-kmesh[0], gInfo.GGT, omegaSq);
		//Set difference to be the G=0 contribution:
		Vzero = fSingularIntegral * kmesh.size() - fSingularSum;
		kernelMode = PeriodicKernel;
		logPrintf("Singular Vxx(G=0) correction = %le\n", Vzero/(kmesh.size() * gInfo.detR));
	}
	
	if( (exReg.method==ExchangeRegularization::None && params.geometry==CoulombParams::Spherical)
		|| exReg.method==ExchangeRegularization::SphericalTruncated)
	{
		if(exReg.method==ExchangeRegularization::None)
			Rc = ((CoulombSpherical&)coulomb).Rc;
		else
		{	double Vsuper = fabs(det(params.supercell->Rsuper)); //supercell volume
			Rc = pow((3./(4*M_PI)) * Vsuper, 1./3);
		}
		if(omega) //Initialize look-up table for erfc fourier transform
		{	//Determine end of range for look-up table:
			double Gmax = 0.0;
			vector3<int> c;
			for(c[0]=-1; c[0]<=1; c[0]+=2) for(c[1]=-1; c[1]<=1; c[1]+=2) for(c[2]=-1; c[2]<=1; c[2]+=2)
			{	vector3<> f; for(int k=0; k<3; k++) f[k] = c[k]*(gInfo.S[k]/2 + 1);
				double G = sqrt(gInfo.GGT.metric_length_squared(f));
				if(G>Gmax) Gmax=G;
			}
			//Compute samples and spline coefficients:
			const double dG = 0.01 * (2*M_PI)/Rc;
			int nSamples = int(ceil(Gmax/dG)) + 10;
			std::vector<double> samples(nSamples);
			for(int i=0; i<nSamples; i++)
				samples[i] = truncatedErfcTilde(dG*i, omega, Rc);
			Vzero = samples[0];
			std::vector<double> coeff = QuinticSpline::getCoeff(samples);
			#ifdef GPU_ENABLED
			cudaMalloc(&coeffErfcTilde, sizeof(double)*coeff.size());
			cudaMemcpy(coeffErfcTilde, coeff.data(), sizeof(double)*coeff.size(), cudaMemcpyHostToDevice);
			gpuErrorCheck();
			#else
			coeffErfcTilde = new double[coeff.size()];
			memcpy(coeffErfcTilde, coeff.data(), sizeof(double)*coeff.size());
			#endif
			dGinvErfcTilde = 1.0/dG;
			nSamplesErfcTilde = nSamples;
		}
		else
			Vzero = (2*M_PI) * Rc*Rc;
		kernelMode = SphericalKernel;
		logPrintf("Truncated on a sphere of radius %lg bohrs.\n", Rc);
	}
	
	if( (exReg.method==ExchangeRegularization::None && params.geometry==CoulombParams::Isolated)
		|| exReg.method==ExchangeRegularization::WignerSeitzTruncated)
	{
		die("Not yet implemented.\n");
		
		kernelMode = SphericalKernel;
		logPrintf("Truncated on a sphere of radius %lg bohrs.\n", Rc);
	}
}

ExchangeEval::~ExchangeEval()
{	
	if(coeffErfcTilde)
	{
		#ifdef GPU_ENABLED
		cudaFree(coeffErfcTilde);
		#else
		delete[] coeffErfcTilde;
		#endif
	}
}


complexDataGptr ExchangeEval::operator()(complexDataGptr&& in, vector3<> kDiff) const
{
	const ExchangeRegularization& exReg = params.exchangeRegularization;
	
	#define CALL_exchangeAnalytic(calc) callPref(exchangeAnalytic)(gInfo.S, gInfo.GGT, calc, in->dataPref(), kDiff, Vzero, symmThresholdSq)
	switch(kernelMode)
	{	case PeriodicKernel:
		{	if(omega) CALL_exchangeAnalytic(ExchangePeriodicScreened_calc(omega));
			else CALL_exchangeAnalytic(ExchangePeriodic_calc());
			break;
		}
		case SphericalKernel:
		{	if(omega) CALL_exchangeAnalytic(ExchangeSphericalScreened_calc(coeffErfcTilde, dGinvErfcTilde, nSamplesErfcTilde));
			else CALL_exchangeAnalytic(ExchangeSpherical_calc(Rc));
			break;
		}
		case WignerSeitzKernel:
		{	die("Not yet implemented.\n");
			break;
		}
	}
	#undef CALL_exchangeAnalytic
	return in;
}
