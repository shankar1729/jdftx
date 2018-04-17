/*-------------------------------------------------------------------
Copyright 2017 Ravishankar Sundararaman

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

#include <electronic/TetrahedralDOS.h>
#include <core/LatticeUtils.h>
#include <core/Util.h>
#include <algorithm>
#include <cfloat>
#include <map>
#include <set>

TetrahedralDOS::TetrahedralDOS(std::vector<vector3<>> kmesh, std::vector<int> iReduced,
	const matrix3<>& R, const matrix3<int>& super, int nSpins, int nBands, int nWeights, double weightSum)
: nSpins(nSpins), nBands(nBands), nWeights(nWeights),
nReduced(iReduced.size()
	? *std::max_element(iReduced.begin(),iReduced.end())+1 //reduced k-pt mapping provided
	: kmesh.size()), //no mapping to reduced; use full mesh
nStates(nReduced*nSpins), eigs(nStates*nBands, 0.), weights(nStates*nBands*nWeights, 1.)
{
	//Pick optimum tetrahedral tesselation of each parallelopiped cell in BZ:
	vector3<> dk[8]; //k-point offsets of basis parallelopiped relative to its first vertex
	matrix3<> GT = ~((2*M_PI) * inv(R));
	matrix3<> GTsup = ~((2*M_PI) * inv(R*super));
	matrix3<> GGT = (~GT) * GT; //G-space metric for kmesh
	{	//There are 4 distinct tesselations listed by signs of each basis vector (upto an overall sign)
		//Pick the one with the shortest common body-diagonal to yield the most suitable aspect ratios
		vector3<int> sBest;
		double dMin = DBL_MAX;
		vector3<int> s;
		s[0] = 1;
		for(s[1]=-1; s[1]<=1; s[1]+=2)
		for(s[2]=-1; s[2]<=1; s[2]+=2)
		{	double d = (GTsup * s).length();
			if(d < dMin)
			{	dMin = d;
				sBest = s;
			}
		}
		//Compute the offsets:
		matrix3<> invGT_GTsup = inv(GT) * GTsup;
		vector3<int> i;
		int ik = 0;
		for(i[0]=0; abs(i[0])<2; i[0]+=sBest[0])
		for(i[1]=0; abs(i[1])<2; i[1]+=sBest[1])
		for(i[2]=0; abs(i[2])<2; i[2]+=sBest[2])
			dk[ik++] = invGT_GTsup * i;
	}
	//--- parallelopiped vertex indices of the 6 tetrahedra in its tesselation
	int cellVerts[6][4] = {
		{0, 2, 3, 7},
		{0, 3, 1, 7},
		{0, 1, 5, 7},
		{0, 5, 4, 7},
		{0, 4, 6, 7},
		{0, 6, 2, 7} };
	
	//Construct Brillouin zone triangulation in an TetrahedralDOS object:
	PeriodicLookup<vector3<>> plook(kmesh, GGT);
	tetrahedra.resize(6*kmesh.size()); //6 tetrahedra per parallelopiped cell
	double Vtot = 0.;
	for(unsigned i=0; i<kmesh.size(); i++)
	{	const vector3<>& v0 = kmesh[i];
		//initialize the parallelopiped vertices and look up their indices:
		vector3<> v[8]; size_t iv[8];
		for(int j=0; j<8; j++)
		{	v[j] = v0 + dk[j];
			iv[j] = plook.find(v[j]);
			assert(iv[j] != string::npos);
		}
		//loop over tetrahedra:
		for(unsigned t=0; t<6; t++)
		{	Tetrahedron& tet = tetrahedra[6*i+t];
			for(int p=0; p<4; p++)
				tet.q[p] = iReduced.size() ? iReduced[iv[cellVerts[t][p]]] : iv[cellVerts[t][p]];
			tet.V = box(
				GT * (v[cellVerts[t][1]] - v[cellVerts[t][0]]),
				GT * (v[cellVerts[t][2]] - v[cellVerts[t][0]]),
				GT * (v[cellVerts[t][3]] - v[cellVerts[t][0]]) );
			Vtot += tet.V;
		}
	}
	double VnormFac = weightSum/Vtot;
	for(Tetrahedron& t: tetrahedra)
		t.V *= VnormFac; //normalize volume of tetrahedra to add up to qWeightSum/nSpins
}

void TetrahedralDOS::setEigs(const std::vector<diagMatrix>& E)
{	assert(int(E.size())==nStates);
	for(int q=0; q<nStates; q++)
		for(int b=0; b<nBands; b++)
			e(q,b) = E[q][b];
}

void TetrahedralDOS::setWeights(int iWeight, const std::vector<diagMatrix>& weights)
{	assert(int(weights.size())==nStates);
	for(int q=0; q<nStates; q++)
		for(int b=0; b<nBands; b++)
			w(iWeight,q,b) = weights[q][b];
}


//Replace clusters of eigenvalues that differ by less than Etol, to a single value
void TetrahedralDOS::weldEigenvalues(double Etol)
{	std::multimap<double,size_t> eigMap;
	for(size_t i=0; i<eigs.size(); i++)
		eigMap.insert(std::make_pair(eigs[i],i));
	for(auto i=eigMap.begin(); i!=eigMap.end();)
	{	auto j = i;
		double ePrev = i->first;
		j++;
		while(j!=eigMap.end() && (j->first <= ePrev+Etol)) //Stop when j'th eig differs by more than Etol from previous one
		{	ePrev = j->first;
			j++;
		}
		//Replace cluster [i,j) by its mean:
		double eMean=0.; size_t eCount=0;
		for(auto k=i; k!=j; k++) { eMean += k->first; eCount++; }
		eMean /= eCount;
		for(auto k=i; k!=j; k++) eigs[k->second] = eMean;
		//Move iterator past end of cluster:
		i = j;
	}
}

//Apply gaussian smoothing of width Esigma
TetrahedralDOS::Lspline TetrahedralDOS::gaussSmooth(const Lspline& in, double Esigma) const
{	assert(Esigma > 0.);
	//Initialize uniform energy grid:
	double Emin = in.front().first - 10*Esigma;
	double Emax = in.back().first + 10*Esigma;
	double dE = 0.2*Esigma;
	size_t nE = ceil((Emax-Emin)/dE);
	if(nE > 1000000) logPrintf(
		"WARNING: very fine energy grid for DOS. If this takes too long /\n"
		"         results in too large a file, either increase Esigma or\n"
		"         set it\n to zero (raw output of tetrahedron method).\n" );
	Lspline out(nE, std::make_pair(0., std::vector<double>(nWeights, 0.)));
	for(size_t iE=0; iE<nE; iE++) out[iE].first = Emin + iE*dE;
	//Apply the gaussian smoothing to each channel:
	const double EsigmaDen = 1./(Esigma*sqrt(2.));
	for(size_t iIn=0; iIn+1<in.size(); iIn++)
	{	const double& E0 = in[ iIn ].first; const std::vector<double>& w0 = in[ iIn ].second;
		const double& E1 = in[iIn+1].first; const std::vector<double>& w1 = in[iIn+1].second;
		if(E1==E0) continue;
		size_t iEstart = std::max(floor((E0-10*Esigma-Emin)/dE), 0.);
		size_t iEstop = std::min(ceil((E1+10*Esigma-Emin)/dE), double(nE));
		for(size_t iE=iEstart; iE<iEstop; iE++)
		{	double E = out[iE].first;
			double e0 = (E-E0)*EsigmaDen;
			double e1 = (E1-E)*EsigmaDen;
			double gaussTerm = (exp(-e0*e0) - exp(-e1*e1)) / (2*sqrt(M_PI) * (e0 + e1));
			double erfTerm = (erf(e0) + erf(e1)) / (2 * (e0 + e1));
			for(int iW=0; iW<nWeights; iW++)
				out[iE].second[iW] += gaussTerm*(w1[iW] - w0[iW]) + erfTerm*(e0*w1[iW] + e1*w0[iW]);
		}
	}
	return out;
}

//Write the density of states to a file, for a given state offset:
void TetrahedralDOS::printDOS(const Lspline& dos, string filename, string header)
{	logPrintf("Dumping '%s' ... ", filename.c_str()); logFlush();
	FILE* fp = fopen(filename.c_str(), "w");
	if(!fp) die("Could not open '%s' for writing.\n", filename.c_str());
	//Header:
	if(header.length())
		fprintf(fp, "%s\n", header.c_str());
	//Data:
	for(Lspline::const_iterator iter=dos.begin(); iter!=dos.end(); iter++)
	{	fprintf(fp, "%.15le", iter->first); 
		for(int i=0; i<nWeights; i++)
			fprintf(fp, "\t%.15le", iter->second[i]);
		fprintf(fp, "\n");
	}
	fclose(fp);
	logPrintf("done.\n"); logFlush();
}


//---- Internal implementation and private functions -----

//An interval in energy
struct Interval
{	double eStart, eStop;
	Interval(double eStart=0.0, double eStop=0.0):eStart(eStart),eStop(eStop) {}
	bool operator<(const Interval& i) const { if(eStart!=i.eStart) return eStart<i.eStart; else return eStop<i.eStop; }
};

//Coefficient set in an interval for an array of cubic splines, corresponding to
//b[0] (1-t)^3 + b[1] 3t(1-t)^2 + b[2] 3(1-t)t^2 + b[3] t^3 for each b in bArr
struct CsplineElem
{	typedef std::array<double,4> double4;
	std::vector<double4> bArr;
	
	//Allocate and initializie to zero if unallocated, otherwise leave unchanged
	void nullToZero(int nWeights)
	{	double4 zero4 = {{0.,0.,0.,0.}};
		if(!bArr.size()) bArr.assign(nWeights, zero4);
	}
	
	double value(double t, int i) const
	{	double4 b = bArr[i];
		double c[3], d[2];
		//deCasteljau's algorithm for cubic bezier:
		for(int k=0; k<3; k++) c[k] = b[k] + t*(b[k+1]-b[k]);
		for(int k=0; k<2; k++) d[k] = c[k] + t*(c[k+1]-c[k]);
		return d[0]+t*(d[1]-d[0]);
	}
	
	double deriv(double t, int i) const
	{	double4 b = bArr[i];
		double c[3], d[2];
		//Derivative is a bezier spline of degree 2 with coefficients:
		for(int k=0; k<3; k++) c[k] = 3.0*(b[k+1]-b[k]);
		//deCasteljau's algorithm for the quadtratic bezier:
		for(int k=0; k<2; k++) d[k] = c[k] + t*(c[k+1]-c[k]);
		return d[0]+t*(d[1]-d[0]);
	}
};

struct Cspline : public std::map<Interval,CsplineElem> //array of piecewise cubic splines (one for each weight function)
{	std::map<double, std::vector<double> > deltas; //additional delta functions (from tetrahedra with same energy for all vertices)
};

//Accumulate contribution from one tetrahedron (exactly a cubic spline for linear interpolation)
//to the weighted DOS for all weight functions (from a single band)
void TetrahedralDOS::accumTetrahedron(const Tetrahedron& t, int iBand, int iSpin, Cspline& wdos) const
{	//sort vertices in ascending order of energy:
	std::array<int,4> q = t.q;
	struct EnergyCmp
	{	const TetrahedralDOS& td; int iBand, stateOffset;
		EnergyCmp(const TetrahedralDOS& td, int iBand, int iSpin) : td(td), iBand(iBand), stateOffset(iSpin*td.nReduced) {}
		const double& e(int q) { return td.e(stateOffset+q, iBand); }
		const double* w(int q) { return &td.w(0, stateOffset+q, iBand); }
		bool operator()(int q1, int q2) { return e(q1)<e(q2); }
	};
	EnergyCmp eCmp(*this, iBand, iSpin); //accessor object and functor to sort vertices by energy
	std::sort(q.begin(), q.end(), eCmp);
	//Load energies and weights from memory:
	double e0=eCmp.e(q[0]), e1=eCmp.e(q[1]), e2=eCmp.e(q[2]), e3=eCmp.e(q[3]);
	const double *w0=eCmp.w(q[0]), *w1=eCmp.w(q[1]), *w2=eCmp.w(q[2]), *w3=eCmp.w(q[3]);
	//Area coefficient
	if(e3==e0)
	{	//Implies e0=e1=e2=e3, and the corresponding density of states is a delta function
		std::vector<double>& wDelta = wdos.deltas[e0];
		if(!wDelta.size()) wDelta.resize(nWeights, 0.);
		for(int i=0; i<nWeights; i++)
			wDelta[i] += t.V * (1./4) * (w0[i] + w1[i] + w2[i] + w3[i]);
		return;
	}
	double inv_e30 = 1.0/(e3-e0);
	double A = inv_e30 * t.V;
	//Intersection points and weights at constant planes:
	double E13_0 = (e1-e0)*inv_e30;
	double E20_3 = (e3-e2)*inv_e30;
	double E12_0 = 0., E21_3 = 0.;
	if(e2>e0) E12_0 = (e1-e0)/(e2-e0);
	if(e3>e1) E21_3 = (e3-e2)/(e3-e1);
	//Create the coefficients:
	CsplineElem *c01=0, *c12=0, *c23=0;
	if(e1>e0) { c01 = &wdos[Interval(e0,e1)]; c01->nullToZero(nWeights); }
	if(e2>e1) { c12 = &wdos[Interval(e1,e2)]; c12->nullToZero(nWeights); }
	if(e3>e2) { c23 = &wdos[Interval(e2,e3)]; c23->nullToZero(nWeights); }
	for(int i=0; i<nWeights; i++)
	{	double w0i=w0[i], w1i=w1[i], w2i=w2[i], w3i=w3[i];
		double wai = w0i + (w3i-w0i)*E13_0;
		double wbi = w0i + (w2i-w0i)*E12_0;
		double wci = w3i + (w0i-w3i)*E20_3;
		double wdi = w3i + (w1i-w3i)*E21_3;
		if(c01)
		{	CsplineElem::double4& b = c01->bArr[i];
			b[2] += A*E12_0*w0i;
			b[3] += A*E12_0*(w1i+wbi+wai);
		}
		if(c12)
		{	CsplineElem::double4& b = c12->bArr[i];
			b[0] += A*E12_0*(w1i+wbi+wai);
			b[1] += A*(w1i + (1.0/3)*(2*wai+wbi + E12_0*(2*w2i+wci)));
			b[2] += A*(w2i + (1.0/3)*(2*wci+wdi + E21_3*(2*w1i+wai)));
			b[3] += A*E21_3*(w2i+wci+wdi);
		}
		if(c23)
		{	CsplineElem::double4& b = c23->bArr[i];
			b[0] += A*E21_3*(w2i+wci+wdi);
			b[1] += A*E21_3*w3i;
		}
	}
}

//Coalesce overlapping splines: convert an arbitrary set of spline pieces into a regular ordered piecewise spline
void TetrahedralDOS::coalesceIntervals(Cspline& cspline) const
{	//Create a list of nodes for the splines:
	std::set<double> e;
	for(Cspline::const_iterator i=cspline.begin(); i!=cspline.end(); i++)
	{	e.insert(i->first.eStart);
		e.insert(i->first.eStop);
	}
	for(auto iter: cspline.deltas)
		e.insert(iter.first);
	//Now break each original spline piece down to the intervals in e
	std::vector<double> v(nWeights), d(nWeights); //temporary storage for values and derivatives
	std::map<Interval, CsplineElem> combined;
	for(Cspline::const_iterator cIter=cspline.begin(); cIter!=cspline.end(); cIter++)
	{	const double& eStart = cIter->first.eStart;
		const double& eStop = cIter->first.eStop;
		double inv_de = 1.0/(eStop-eStart);
		const CsplineElem& celem = cIter->second;
		std::set<double>::iterator j=e.find(eStart), jStop=e.find(eStop);
		double e = *j;
		double t = (e-eStart)*inv_de;
		for(int i=0; i<nWeights; i++)
		{	v[i] = celem.value(t, i);
			d[i] = celem.deriv(t, i)*inv_de;
		}
		while(j!=jStop)
		{	j++;
			double eNext = *j;
			double tNext = (eNext-eStart)*inv_de;
			CsplineElem& c = combined[Interval(e,eNext)];
			c.nullToZero(nWeights);
			for(int i=0; i<nWeights; i++)
			{	double viNext = celem.value(tNext, i);
				double diNext = celem.deriv(tNext, i)*inv_de;
				CsplineElem::double4& bi = c.bArr[i];
				bi[0] += v[i];
				bi[1] += v[i] + (1.0/3)*(eNext-e)*d[i];
				bi[2] += viNext - (1.0/3)*(eNext-e)*diNext;
				bi[3] += viNext;
				v[i] = viNext;
				d[i] = diNext;
			}
			e=eNext;
			t=tNext;
		}
	}
	//Fix interior discontinuities and incorporate interior deltas:
	auto cIter = combined.begin();
	double h = cIter->first.eStop - cIter->first.eStart;
	auto cIterNext = cIter; cIterNext++;
	while(cIterNext != combined.end())
	{	double hNext = cIterNext->first.eStop - cIterNext->first.eStart;
		assert(cIterNext->first.eStart == cIter->first.eStop);
		//Check for a delta at the intersection of the two intervals:
		std::vector<double>* deltaPtr = 0;
		auto dIter = cspline.deltas.find(cIter->first.eStop);
		if(dIter != cspline.deltas.end())
			deltaPtr = &(dIter->second);
		//Update weights at the interval intersection:
		for(int i=0; i<nWeights; i++)
		{	double integral = h*cIter->second.bArr[i][3] + hNext*cIterNext->second.bArr[i][0];
			if(deltaPtr) integral += 4*deltaPtr->at(i);
			double wMean = integral/(h + hNext);
			cIter->second.bArr[i][3] = wMean;
			cIterNext->second.bArr[i][0] = wMean;
		}
		//Advance to next interval
		h = hNext;
		cIter = cIterNext;
		cIterNext++;
	}
	//Incorporate deltas at the extrema:
	if(cspline.deltas.size())
	{	//Delta at bottom of band:
		auto cIter = combined.begin();
		auto dIter = cspline.deltas.begin();
		if(dIter->first == cIter->first.eStart)
		{	double h = cIter->first.eStop - cIter->first.eStart;
			for(int i=0; i<nWeights; i++)
				cIter->second.bArr[i][0] += dIter->second[i]*(4./h);
		}
		//Delta at top of band:
		cIter = combined.end(); cIter--; //iterator to last entry
		dIter = cspline.deltas.end(); dIter--; //iterator to last entry
		if(dIter->first == cIter->first.eStop)
		{	double h = cIter->first.eStop - cIter->first.eStart;
			for(int i=0; i<nWeights; i++)
				cIter->second.bArr[i][3] += dIter->second[i]*(4./h);
		}
	}
	//Set results:
	cspline.clear();
	cspline.insert(combined.begin(), combined.end());
	cspline.deltas.clear(); //deltas incorporated in above
}

static bool LsplineCmp(const TetrahedralDOS::LsplineElem& l1, const TetrahedralDOS::LsplineElem& l2) { return l1.first < l2.first; }

//Convert cubic splines to integrated linear splines which handle discontinuities and singularities better:
//The Cspline object must be coalesced before passing to this function
TetrahedralDOS::Lspline TetrahedralDOS::convertLspline(const Cspline& cspline) const
{	//Convert each cubic-spline interval into three linear spline ones:
	int nIntervals = cspline.size();
	Lspline lspline(1+3*nIntervals);
	auto lIter = lspline.begin();
	for(auto cIter = cspline.begin(); cIter != cspline.end(); cIter++) //loop over c-spline intervals
	{	const double& e0 = cIter->first.eStart;
		const double& e3 = cIter->first.eStop;
		const CsplineElem& ce = cIter->second;
		//Macro to add interval:
		#define ADD_interval(e, wiCode) \
		{	lIter->first = e; \
			lIter->second.resize(nWeights); \
			for(int i=0; i<nWeights; i++) \
				lIter->second[i] = (wiCode); \
			lIter++; \
		}
		//Add intervals:
		const double e1 = (1./3)*(2*e0 + e3);
		const double e2 = (1./3)*(e0 + 2*e3);
		const double a10 = 13./60;
		const double a11 = 0.60;
		const double a12 = 0.15;
		const double a13 = 1./30;
		if(cIter==cspline.begin()) ADD_interval(e0, ce.bArr[i][0]) //only need at start of spline
		ADD_interval(e1, a10*ce.bArr[i][0] + a11*ce.bArr[i][1] + a12*ce.bArr[i][2] + a13*ce.bArr[i][3])
		ADD_interval(e2, a10*ce.bArr[i][3] + a11*ce.bArr[i][2] + a12*ce.bArr[i][1] + a13*ce.bArr[i][0])
		ADD_interval(e3, ce.bArr[i][3])
	}
	assert(lIter == lspline.end());
	return lspline;
}

//Collect contributions from multiple linear splines (one for each band)
TetrahedralDOS::Lspline TetrahedralDOS::mergeLsplines(const std::vector<Lspline>& lsplines) const
{	//Collect list of energy nodes and interval boundaries:
	std::set<double> eSet, leftBoundarySet, rightBoundarySet;
	for(const Lspline& lspline: lsplines)
	{	leftBoundarySet.insert(lspline.front().first);
		rightBoundarySet.insert(lspline.back().first);
		for(unsigned j=0; j<lspline.size(); j++)
			eSet.insert(lspline[j].first);
	}
	std::set<double> boundarySet; //points that are in left or right boundaries, but not both
	std::set_symmetric_difference( //this eliminates edge-discontinuities in DOS from bands that meet at one energy
		leftBoundarySet.begin(), leftBoundarySet.end(),
		rightBoundarySet.begin(), rightBoundarySet.end(),
		std::inserter(boundarySet, boundarySet.end()) );
	//Compute merging weights for the left and right boundary weights:
	std::map<double,double> startScale, stopScale;
	for(double e: leftBoundarySet)
	{	if(boundarySet.count(e)) startScale[e] = 1.; //unmerged boundary
		else //merge weighted by interval length (preserves integral under trapezoidal rule):
		{	auto eIter = eSet.find(e);
			auto eIterMinus = eIter; eIterMinus--; double eMinus = *eIterMinus;
			auto eIterPlus = eIter; eIterPlus++; double ePlus = *eIterPlus;
			startScale[e] = (ePlus-e) / (ePlus-eMinus);
		}
	}
	for(double e: rightBoundarySet)
	{	if(boundarySet.count(e)) stopScale[e] = 1.; //unmerged boundary
		else //merge weighted by interval length (preserves integral under trapezoidal rule):
		{	auto eIter = eSet.find(e);
			auto eIterMinus = eIter; eIterMinus--; double eMinus = *eIterMinus;
			auto eIterPlus = eIter; eIterPlus++; double ePlus = *eIterPlus;
			stopScale[e] = (e-eMinus) / (ePlus-eMinus);
		}
	}
	//Create output linear spline with duplicate nodes for each boundary node (to represent discontinuities)
	Lspline combined(eSet.size() + boundarySet.size(), std::make_pair(0., std::vector<double>(nWeights,0.)));
	Lspline::iterator cIter = combined.begin();
	for(double e: eSet)
	{	(cIter++)->first = e;
		if(boundarySet.count(e))
			(cIter++)->first = e;
	}
	assert(cIter == combined.end());
	//Collect contributions from each linear spline:
	for(const Lspline& lspline: lsplines)
	{	Lspline::iterator cIter = std::upper_bound(combined.begin(), combined.end(), lspline.front(), LsplineCmp);
		cIter--; //Now this points to the last of the possibly duplicate entries for the starting energy
		double eStart = lspline.front().first; //start energy of current input spline
		double eStop = lspline.back().first; //end energy of current input spline
		Lspline::const_iterator leftIter = lspline.begin();
		Lspline::const_iterator rightIter = leftIter + 1;
		double e = eStart;
		while(true)
		{	double t = (e - leftIter->first) / (rightIter->first - leftIter->first);
			double endpointScale = 1.;
			if(e == eStart) endpointScale = startScale[eStart];
			if(e == eStop) endpointScale = stopScale[eStop];
			for(int i=0; i<nWeights; i++)
				cIter->second[i] += endpointScale * ((1.-t) * leftIter->second[i] + t * rightIter->second[i]);
			//Check for end of input spline:
			if(e == eStop) break; //note this ends at the first of the possibly duplicate entries for the ending energy
			//Advance and update input interval if required:
			cIter++;
			e = cIter->first;
			if(e > rightIter->first)
			{	leftIter++;
				rightIter++;
			}
		}
	}
	return combined;
}

//Generate the density of states for a given state offset:
TetrahedralDOS::Lspline TetrahedralDOS::getDOS(int iSpin, double Etol) const
{	std::vector<Lspline> lsplines(nBands);
	for(int iBand=0; iBand<nBands; iBand++)
	{	Cspline wdos;
		for(const Tetrahedron& t: tetrahedra)
			accumTetrahedron(t, iBand, iSpin, wdos);
		if(wdos.size()==0 && wdos.deltas.size()==1) // band is a single delta function
		{	double eDelta = wdos.deltas.begin()->first;
			const std::vector<double>& wDelta = wdos.deltas.begin()->second;
			lsplines[iBand].resize(3, std::make_pair(eDelta, std::vector<double>(nWeights, 0.)));
			lsplines[iBand][0].first = eDelta-0.5*Etol;
			lsplines[iBand][2].first = eDelta+0.5*Etol;
			for(int i=0; i<nWeights; i++)
				lsplines[iBand][1].second[i] = wDelta[i] * (2./Etol);
		}
		else
		{	coalesceIntervals(wdos);
			lsplines[iBand] = convertLspline(wdos);
		}
	}
	return mergeLsplines(lsplines);
}
