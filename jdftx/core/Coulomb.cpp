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
#include <core/Coulomb_ExchangeEval.h>
#include <core/LoopMacros.h>
#include <core/BlasExtra.h>
#include <core/Thread.h>
#include <core/Operators.h>
#include "LatticeUtils.h"

CoulombParams::CoulombParams() : ionMargin(5.), embed(false), embedFluidMode(false)
{
}

std::shared_ptr<Coulomb> CoulombParams::createCoulomb(const GridInfo& gInfo) const
{	if(geometry != Periodic)
	{	logPrintf("\n---------- Setting up coulomb interaction ----------\n");
		Citations::add("Truncated Coulomb potentials", wsTruncationPaper);
	}
	
	if(embedFluidMode) //Use embedding box, but periodic Coulomb kernel (fluid response does the approx image separation)
	{	logPrintf("Fluid mode embedding: using embedded box, but periodic Coulomb kernel.\n");
		logPrintf("(Fluid response is responsible for (approximate) separation between periodic images.)\n");
		if(!embed)
			die("Fluids with coulomb truncation requires the use of command coulomb-truncation-embed.\n");
		return std::make_shared<CoulombPeriodic>(gInfo, *this);
	}
	
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

vector3<bool> CoulombParams::isTruncated() const
{	switch(geometry)
	{	case Isolated:
		case Spherical:
			return vector3<bool>(true, true, true);
		case Cylindrical:
		case Wire:
		{	vector3<bool> result(true, true, true);
			result[iDir] = false;
			return result;
		}
		case Slab:
		{	vector3<bool> result(false, false, false);
			result[iDir] = true;
			return result;
		}
		case Periodic:
		default:
			return vector3<bool>(false, false, false);
	}
}

void CoulombParams::splitEfield(const matrix3<>& R, vector3<>& RT_Efield_ramp, vector3< double >& RT_Efield_wave) const
{	vector3<bool> isTrunc = isTruncated();
	vector3<> RT_Efield = (~R) * Efield; //put in contravariant lattice coordinates
	for(int k=0; k<3; k++)
	{	//zero out component if it is within symmetry threshold:
		if(fabs(RT_Efield[k]/R.column(k).length()) < Efield.length() * symmThreshold)
			RT_Efield[k] = 0.;
		if(isTrunc[k])
		{	RT_Efield_ramp[k] = RT_Efield[k];
			RT_Efield_wave[k] = 0.;
		}
		else
		{	RT_Efield_ramp[k] = 0.;
			RT_Efield_wave[k] = RT_Efield[k];
		}
	}
}


//--------------- class Coulomb ----------------

template<typename scalar> void boundarySymmetrize(const std::vector< std::pair<int,IndexArray> >& symmIndex, scalar* data)
{	for(unsigned n=2; n<symmIndex.size(); n++)
		if(symmIndex[n].first)
			callPref(eblas_symmetrize)(symmIndex[n].first, n, symmIndex[n].second.dataPref(), data);
}

ScalarFieldTilde Coulomb::embedExpand(const ScalarFieldTilde& in) const
{	assert(params.embed);
	assert(&(in->gInfo) == &gInfoOrig);
	ScalarField out; nullToZero(out, gInfo);
	callPref(eblas_scatter_daxpy)(gInfoOrig.nr, 1., embedIndex.dataPref(), I(centerToO * in)->dataPref(), out->dataPref());
	boundarySymmetrize(symmIndex, out->dataPref());
	return J(out);
}

complexScalarFieldTilde Coulomb::embedExpand(complexScalarFieldTilde&& in) const
{	assert(params.embed);
	assert(&(in->gInfo) == &gInfoOrig);
	complexScalarField out; nullToZero(out, gInfo);
	callPref(eblas_scatter_zdaxpy)(gInfoOrig.nr, 1., embedIndex.dataPref(), I(Complex(centerToO) * ((complexScalarFieldTilde&&)in))->dataPref(), out->dataPref());
	boundarySymmetrize(symmIndex, out->dataPref());
	return J((complexScalarField&&)out);
}

ScalarFieldTilde Coulomb::embedShrink(const ScalarFieldTilde& in) const
{	assert(params.embed);
	assert(&(in->gInfo) == &gInfo);
	ScalarField Iin = I(in);
	boundarySymmetrize(symmIndex, Iin->dataPref());
	ScalarField out; nullToZero(out, gInfoOrig);
	callPref(eblas_gather_daxpy)(gInfoOrig.nr, 1., embedIndex.dataPref(), Iin->dataPref(), out->dataPref());
	return centerFromO * J(out);
}

complexScalarFieldTilde Coulomb::embedShrink(complexScalarFieldTilde&& in) const
{	assert(params.embed);
	assert(&(in->gInfo) == &gInfo);
	complexScalarField Iin = I((complexScalarFieldTilde&&)in);
	boundarySymmetrize(symmIndex, Iin->dataPref());
	complexScalarField out; nullToZero(out, gInfoOrig);
	callPref(eblas_gather_zdaxpy)(gInfoOrig.nr, 1., embedIndex.dataPref(), Iin->dataPref(), out->dataPref());
	return Complex(centerFromO) * J((complexScalarField&&)out);
}

ScalarFieldTilde Coulomb::operator()(ScalarFieldTilde&& in, PointChargeMode pointChargeMode) const
{	if(params.embed)
	{	ScalarFieldTilde outSR;
		if(pointChargeMode!=PointChargeNone) outSR = (*ionKernel) * in; //Special handling (range separation) to avoid Nyquist frequency issues
		if(pointChargeMode==PointChargeRight) in = gaussConvolve(in, ionWidth); //bandwidth-limit point charge to the right and compute long-ranged part below
		ScalarFieldTilde outLR = embedShrink(apply(embedExpand(in))); //Apply truncated Coulomb in expanded grid and shrink back
		if(pointChargeMode==PointChargeLeft) outLR = gaussConvolve(outLR, ionWidth); //since point-charge to left, bandwidth-limit long-range part computed above
		return outLR + outSR;
	}
	else return apply((ScalarFieldTilde&&)in);
}

ScalarFieldTilde Coulomb::operator()(const ScalarFieldTilde& in, PointChargeMode pointChargeMode) const
{	ScalarFieldTilde out(in->clone()); //create destructible copy
	return (*this)((ScalarFieldTilde&&)out, pointChargeMode);
}

complexScalarFieldTilde Coulomb::operator()(complexScalarFieldTilde&& in, vector3<> kDiff, double omega) const
{	auto exEvalOmega = exchangeEval.find(omega);
	assert(exEvalOmega != exchangeEval.end());
	if(params.embed) return embedShrink((*exEvalOmega->second)(embedExpand((complexScalarFieldTilde&&)in), kDiff));
	else return (*exEvalOmega->second)((complexScalarFieldTilde&&)in, kDiff);
}

complexScalarFieldTilde Coulomb::operator()(const complexScalarFieldTilde& in, vector3<> kDiff, double omega) const
{	complexScalarFieldTilde out(in->clone()); //create destructible copy
	return (*this)((complexScalarFieldTilde&&)out, kDiff, omega);
}

double Coulomb::energyAndGrad(std::vector<Atom>& atoms) const
{	if(!ewald) ((Coulomb*)this)->ewald = createEwald(gInfo.R, atoms.size());
	double Eewald = 0.;
	if(params.embed)
	{	matrix3<> embedScaleMat = Diag(embedScale);
		matrix3<> invEmbedScaleMat = inv(embedScaleMat);
		//Convert atom positions to embedding grid's lattice coordinates:
		for(unsigned i=0; i<atoms.size(); i++)
		{	Atom& a = atoms[i];
			//Center, restrict to WS-cell and check margins:
			a.pos = wsOrig->restrict(a.pos - xCenter);
			bool posValid = true;
			switch(params.geometry)
			{	case CoulombParams::Periodic:
					assert(!"Embedding not allowed in periodic geometry");
					break;
				case CoulombParams::Slab:
				{	double L = gInfoOrig.R.column(params.iDir).length();
					posValid = (fabs(a.pos[params.iDir]*L) < 0.5*L-params.ionMargin);
					break;
				}
				case CoulombParams::Cylindrical:
				case CoulombParams::Wire:
				{	posValid = (wsOrig->boundaryDistance(a.pos, params.iDir) > params.ionMargin);
					break;
				}
				case CoulombParams::Spherical:
				case CoulombParams::Isolated:
				{	posValid = (wsOrig->boundaryDistance(a.pos) > params.ionMargin);
					break;
				}
			}
			if(!posValid) die("Atom %d lies within the margin of %lg bohrs from the truncation boundary.\n" ionMarginMessage, i+1, params.ionMargin);
			//Scale to embedding-mesh lattice coords:
			a.pos = embedScaleMat * a.pos;
			a.force = invEmbedScaleMat * a.force;
		}
		//Compute on the unit cell of the embedding mesh:
		Eewald += ewald->energyAndGrad(atoms);
		//Convert positions and forces back to original mesh:
		for(Atom& a: atoms)
		{	a.pos = xCenter + invEmbedScaleMat * a.pos;
			a.force = embedScaleMat * a.force;
		}
	}
	else Eewald = ewald->energyAndGrad(atoms);
	//Electric field contributions if any:
	if(params.Efield.length_squared())
	{	vector3<> RT_Efield_ramp, RT_Efield_wave;
		params.splitEfield(gInfoOrig.R, RT_Efield_ramp, RT_Efield_wave);
		for(unsigned i=0; i<atoms.size(); i++)
		{	Atom& a = atoms[i];
			vector3<> x = a.pos - xCenter;
			for(int k=0; k<3; k++)
			{	//note that Z is negative of nuclear charge in present sign convention
				Eewald += a.Z * (RT_Efield_ramp[k]*x[k] + RT_Efield_wave[k]*sin(2*M_PI*x[k])/(2*M_PI));
				a.force[k] -= a.Z * (RT_Efield_ramp[k] + RT_Efield_wave[k]*cos(2*M_PI*x[k]));
			}
		}
	}
	return Eewald;
}

void getEfieldPotential_sub(size_t iStart, size_t iStop, const vector3<int>& S, const WignerSeitz* ws,
	const vector3<>& xCenter, const vector3<>& RT_Efield_ramp, const vector3<>& RT_Efield_wave, double* V)
{	matrix3<> invS = inv(Diag(vector3<>(S)));
	THREAD_rLoop
	(	vector3<> x = ws->restrict(invS * iv - xCenter);
		double Vcur = 0.;
		for(int k=0; k<3; k++)
			Vcur -= (RT_Efield_ramp[k]*x[k] + RT_Efield_wave[k]*sin(2*M_PI*x[k])/(2*M_PI));
		V[i] = Vcur;
	)
}
ScalarField Coulomb::getEfieldPotential() const
{	if(params.Efield.length_squared())
	{	vector3<> RT_Efield_ramp, RT_Efield_wave;
		params.splitEfield(gInfoOrig.R, RT_Efield_ramp, RT_Efield_wave);
		ScalarField V(ScalarFieldData::alloc(gInfoOrig));
		threadLaunch(getEfieldPotential_sub, gInfoOrig.nr, gInfoOrig.S, wsOrig, xCenter, RT_Efield_ramp, RT_Efield_wave, V->data());
		return V;
	}
	else return ScalarField();
}

void setEmbedIndex_sub(size_t iStart, size_t iStop, const vector3<int>& S, const vector3<int>& Sembed, const WignerSeitz* ws, int* embedIndex)
{	vector3<> invS; for(int k=0; k<3; k++) invS[k] = 1./S[k];
	THREAD_rLoop
	(	vector3<int> ivEmbed = ws->restrict(iv, S, invS); //wrapped coordinates within WS cell of original mesh
		for(int k=0; k<3; k++) //wrapped embedding mesh cooridnates within first fundamental domain:
		{	ivEmbed[k] = ivEmbed[k] % Sembed[k];
			if(ivEmbed[k]<0) ivEmbed[k] += Sembed[k];
		}
		embedIndex[i] = ivEmbed[2] + Sembed[2]*(ivEmbed[1] + Sembed[1]*ivEmbed[0]);
	)
}

//Initialize index maps for symmetric points on the boundary
//Note: in this case S is the double size mesh and Sorig is the smaller original mesh
void setEmbedBoundarySymm_sub(size_t iStart, size_t iStop, const vector3<int>& Sorig, const vector3<int>& S,
	const WignerSeitz* wsOrig, std::mutex* m, std::multimap<int,int>* boundaryMap)
{	//Compute partial index map per thread:
	std::multimap<int,int> bMap;
	vector3<> invS; for(int k=0; k<3; k++) invS[k] = 1./S[k];
	vector3<> invSorig; for(int k=0; k<3; k++) invSorig[k] = 1./Sorig[k];
	THREAD_rLoop
	(	vector3<int> ivWS = wsOrig->restrict(iv, S, invS); //this is a double-sized Wigner-Seitz cell restriction
		vector3<> xOrig; for(int k=0; k<3; k++) xOrig[k] = ivWS[k] * invSorig[k]; //original lattice coordinates
		if(wsOrig->onBoundary(xOrig))
		{	for(int k=0; k<3; k++)
				ivWS[k] = positiveRemainder(ivWS[k], Sorig[k]);
			bMap.insert(std::pair<int,int>(ivWS[2]+Sorig[2]*(ivWS[1]+Sorig[1]*ivWS[0]), i));
		}
	)
	//Accumulate index map over threads:
	m->lock();
	boundaryMap->insert(bMap.begin(), bMap.end());
	m->unlock();
}

inline void setIonKernel(int i, double Gsq, double expFac, double GzeroVal, double* kernel)
{	kernel[i] = (4*M_PI) * (Gsq ? (1.-exp(-expFac*Gsq))/Gsq : GzeroVal);
}

Coulomb::Coulomb(const GridInfo& gInfoOrig, const CoulombParams& params)
: gInfoOrig(gInfoOrig), params(params), gInfo(params.embed ? gInfoEmbed : gInfoOrig), xCenter(params.embedCenter), wsOrig(0)
{
	if(params.embed)
	{	//Initialize embedding grid:
		logPrintf("Setting up double-sized grid for truncated Coulomb potentials:\n");
		vector3<bool> isTruncated = params.isTruncated();
		for(int k=0; k<3; k++)
		{	int dimScale = isTruncated[k] ? 2 : 1;
			gInfoEmbed.S[k] = dimScale * gInfoOrig.S[k];
			gInfoEmbed.R.set_col(k, dimScale * gInfoOrig.R.column(k));
			embedScale[k] = 1./dimScale;
		}
		gInfoEmbed.initialize(true);
		//Report embedding center in various coordinate systems:
		vector3<int> ivCenter;
		for(int k=0; k<3; k++) //Convert to mesh coordinates (only used by createXSF to align scalar fields):
			ivCenter[k] = int(round(xCenter[k] * gInfoOrig.S[k]));
		vector3<> rCenter = gInfoOrig.R * xCenter; //Cartesian coordinates
		logPrintf("Integer grid location selected as the embedding center:\n");
		logPrintf("   Grid: "); ivCenter.print(globalLog, " %d ");
		logPrintf("   Lattice: "); xCenter.print(globalLog, " %lg ");
		logPrintf("   Cartesian: "); rCenter.print(globalLog, " %lg ");
		//Setup translation operators:
		nullToZero(centerToO, gInfoOrig);   initTranslation(centerToO, -rCenter);
		nullToZero(centerFromO, gInfoOrig); initTranslation(centerFromO, rCenter);
		//Setup Wigner-Seitz cell of original mesh and initialize index map:
		wsOrig = new WignerSeitz(gInfoOrig.R);
		embedIndex.init(gInfoOrig.nr);
		threadLaunch(setEmbedIndex_sub, gInfoOrig.nr, gInfoOrig.S, gInfo.S, wsOrig, embedIndex.data());
		//Setup boundary symmetrization:
		std::multimap<int,int> boundaryMap; std::mutex m;
		threadLaunch(setEmbedBoundarySymm_sub, gInfo.nr, gInfoOrig.S, gInfo.S, wsOrig, &m, &boundaryMap);
		std::map<int, std::vector<int> > symmEquiv; //for each n, n-fold symmetry equivalence classes of the boundary
		for(auto iter=boundaryMap.begin(); iter!=boundaryMap.end();)
		{	int count = 0;
			auto iterNext = iter;
			while(iterNext!=boundaryMap.end() && iterNext->first==iter->first)
			{	count++;
				iterNext++;
			}
			for(; iter!=iterNext; iter++) symmEquiv[count].push_back(iter->second);
		}
		symmIndex.resize(symmEquiv.rbegin()->first+1, std::make_pair(0,IndexArray()));
		for(const auto& entry: symmEquiv)
		{	int n = entry.first; if(n<=1) continue; //Ignore singletons
			int N = entry.second.size();
			const int* srcData = entry.second.data();
			assert(N > 0);
			assert(N % n == 0);
			symmIndex[n].first = N/n;
			symmIndex[n].second.init(N);
			memcpy(symmIndex[n].second.data(), srcData, sizeof(int)*N);
		}
		//Range-separation for ions:
		double hMax = std::max(gInfoOrig.h[0].length(), std::max(gInfoOrig.h[1].length(), gInfoOrig.h[2].length()));
		double Gnyq = M_PI / hMax; //Minimum Nyquist frequency on original grid
		ionWidth = sqrt(params.ionMargin / Gnyq); //Minimize real and reciprocal space errors
		logPrintf("Range-separation parameter for embedded mesh potentials due to point charges: %lg bohrs.\n", ionWidth);
		ionKernel = std::make_shared<RealKernel>(gInfoOrig);
		double expFac = 0.5*ionWidth*ionWidth;
		double GzeroVal = expFac;
		if(params.embedFluidMode)
			GzeroVal -= expFac * det(Diag(embedScale));
		applyFuncGsq(gInfoOrig, setIonKernel, expFac, GzeroVal, ionKernel->data());
	}
	
	//Check electric field:
	if(params.Efield.length_squared())
	{	vector3<> RT_Efield_ramp, RT_Efield_wave;
		params.splitEfield(gInfoOrig.R, RT_Efield_ramp, RT_Efield_wave);
		if(RT_Efield_ramp.length_squared() && !params.embed)
			die("Electric field with component along truncated direction requires coulomb-truncation-embed.");
		if(!wsOrig) wsOrig = new WignerSeitz(gInfoOrig.R);
	}
}

Coulomb::~Coulomb()
{	if(wsOrig) delete wsOrig;
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
//Specialization for slab mode:
template<> void exchangeAnalytic_thread<ExchangeSlab_calc>(size_t iStart, size_t iStop,
	vector3<int> S, const matrix3<>& GGT, const ExchangeSlab_calc& calc,
	complex* data, const vector3<>& kDiff, double Vzero, double thresholdSq)
{	THREAD_fullGspaceLoop
	(	data[i] *= calc(iG, GGT, kDiff, Vzero, thresholdSq);
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
DECLARE_exchangeAnalytic(Slab)
#undef DECLARE_exchangeAnalytic


void multRealKernel_thread(size_t iStart, size_t iStop,
	vector3<int> S, const double* kernel, complex* data)
{	THREAD_fullGspaceLoop( multRealKernel_calc(i, iG, S, kernel, data); )
}
void multRealKernel(vector3<int> S, const double* kernel, complex* data)
{	threadLaunch(multRealKernel_thread, S[0]*S[1]*S[2], S, kernel, data);
}

//Multiply a complexScalarFieldTilde by a kernel sampled with offset and rotation by rot
void multTransformedKernel_thread(size_t iStart, size_t iStop,
	vector3<int> S, const double* kernel, complex* data, const vector3<int>& offset)
{	THREAD_fullGspaceLoop( multTransformedKernel_calc(i, iG, S, kernel, data, offset); )
}
void multTransformedKernel(vector3<int> S, const double* kernel, complex* data, const vector3<int>& offset)
{	threadLaunch(multTransformedKernel_thread, S[0]*S[1]*S[2], S, kernel, data, offset);
}
