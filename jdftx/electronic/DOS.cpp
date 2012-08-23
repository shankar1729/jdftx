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

#include <electronic/DOS.h>
#include <electronic/Everything.h>
#include <electronic/SphericalHarmonics.h>
#include <electronic/ColumnBundle.h>
#include <electronic/operators.h>
#include <core/DataIO.h>
#include <array>
#include <CGAL/Plane_3.h>

#ifdef TRIANGULATE_ENABLED
#include <core/opt/Triangulate.h>
#endif

DOS::DOS() : occupied(false), Etol(1e-6)
{
}

void DOS::setup(const Everything& everything)
{	e = &everything;
	if(!weights.size()) //Add the default of just total density of states:
	{	Weight weight; weight.type = Weight::Total;
		weights.push_back(weight);
	}
}


#ifdef TRIANGULATE_ENABLED

//Density of states evaluation helper class
struct EvalDOS
{
	//Tetrahedron in Brillouin zone triangulation:
	struct Tetrahedron
	{	std::array<int,4> q; //state indices for the vertices
		double V; //volume of tetrahedron
	};
	std::vector<Tetrahedron> tetrahedra;
	
	int nWeights, nStates, nBands; double Etol;
	std::vector<double> eigs; //flat array of eigenvalues (inner index state, outer index bands)
	std::vector<double> weights; //flat array of DOS weights (inner index weight function, middle index state, and outer index bands)
	double& e(int iState, int iBand) { return eigs[iState + nStates*iBand]; } //access eigenvalue
	const double& e(int iState, int iBand) const { return eigs[iState + nStates*iBand]; } //access eigenvalue (const version)
	double& w(int iWeight, int iState, int iBand) { return weights[iWeight + nWeights*(iState + nStates*iBand)]; } //access weight
	const double& w(int iWeight, int iState, int iBand) const { return weights[iWeight + nWeights*(iState + nStates*iBand)]; } //access weight (const version)
	
	EvalDOS(int nCells, int nWeights, int nStates, int nBands, double Etol)
	: tetrahedra(nCells),
	nWeights(nWeights), nStates(nStates), nBands(nBands), Etol(Etol),
	eigs(nStates*nBands), weights(nWeights*nStates*nBands)
	{
	}
	
	//Replace clusters of eigenvalues that differ by less than Etol, to a single value
	void weldEigenvalues()
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
		{	if(!bArr.size()) bArr.assign(nWeights, double4({{0.,0.,0.,0.}}));
		}
		
		double value(double t, int i) const
		{	register double4 b = bArr[i];
			register double c[3], d[2];
			//deCasteljau's algorithm for cubic bezier:
			for(int k=0; k<3; k++) c[k] = b[k] + t*(b[k+1]-b[k]);
			for(int k=0; k<2; k++) d[k] = c[k] + t*(c[k+1]-c[k]);
			return d[0]+t*(d[1]-d[0]);
		}
		
		double deriv(double t, int i) const
		{	register double4 b = bArr[i];
			register double c[3], d[2];
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
	inline void accumTetrahedron(const Tetrahedron& t, int iBand, int stateOffset, Cspline& wdos) const
	{	//sort vertices in ascending order of energy:
		std::array<int,4> q = t.q;
		struct EnergyCmp
		{	const EvalDOS& eval; int iBand, stateOffset;
			EnergyCmp(const EvalDOS& eval, int iBand, int stateOffset) : eval(eval), iBand(iBand), stateOffset(stateOffset) {}
			const double& e(int q) { return eval.e(stateOffset+q, iBand); }
			const double* w(int q) { return &eval.w(0, stateOffset+q, iBand); }
			bool operator()(int q1, int q2) { return e(q1)<e(q2); }
		};
		EnergyCmp eCmp(*this, iBand, stateOffset); //accessor object and functor to sort vertices by energy
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
	void coalesceIntervals(Cspline& cspline) const
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
		for(Cspline::const_iterator i=cspline.begin(); i!=cspline.end(); i++)
		{	const double& eStart = i->first.eStart;
			const double& eStop = i->first.eStop;
			double inv_de = 1.0/(eStop-eStart);
			const CsplineElem& celem = i->second;
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
		//Replace continuous part with combined cubics:
		cspline.clear();
		cspline.insert(combined.begin(), combined.end());
	}
	
	typedef std::vector< std::pair<double, std::vector<double> > > Lspline; //set of linear splines (energies, and vector of DOS values for each weight function)

	//Convert cubic splines to integrated linear splines which handle discontinuities and singularities better:
	//The Cspline object must be coalesced before passing to this function
	Lspline convertLspline(const Cspline& cspline) const
	{	//Collect all the energy points:
		std::set<double> eSet;
		for(auto iter: cspline)
		{	eSet.insert(iter.first.eStart);
			eSet.insert(iter.first.eStop);
		}
		for(auto iter: cspline.deltas)
			eSet.insert(iter.first);
		//Allocate linear spline and set its nodes:
		Lspline lspline(eSet.size());
		auto vecIter = lspline.begin();
		for(double e: eSet)
			(vecIter++)->first = e;
		//Find left- and right-linear-weighted integrals for each interval:
		lspline[0].second.assign(nWeights, 0.);
		for(unsigned j=0; j<lspline.size()-1; j++)
		{	double eStart = lspline[j].first;
			double eStop = lspline[j+1].first;
			double h = eStop - eStart;
			Cspline::const_iterator cIter = cspline.find(Interval(eStart, eStop));
			if(cIter == cspline.end()) assert(!"Brillouin zone triangulation must have holes!\n");
			lspline[j+1].second.assign(nWeights, 0.);
			for(int i=0; i<nWeights; i++)
			{	const CsplineElem::double4& b = cIter->second.bArr[i];
				lspline[ j ].second[i] += h * (0.20*b[0] + 0.15*b[1] + 0.10*b[2] + 0.05*b[3]);
				lspline[j+1].second[i] += h * (0.05*b[0] + 0.10*b[1] + 0.15*b[2] + 0.20*b[3]);
			}
		}
		//Include contributions from deltas:
		//TODO
		//Convert interval integrals to linear spline coefficients:
		for(unsigned j=0; j<lspline.size(); j++)
		{	double eStart = (j==0) ? lspline[j].first : lspline[j-1].first;
			double eStop = (j+1==lspline.size()) ? lspline[j].first : lspline[j+1].first;
			double normFac = 2./(eStop-eStart);
			for(int i=0; i<nWeights; i++)
				lspline[j].second[i] *= normFac;
		}
		return lspline;
	}
	
	//Collect contributions from multiple linear splines (one for each band)
	Lspline mergeLsplines(const std::vector<Lspline>& lsplines) const
	{	//Collect list of energy nodes and interval boundaries:
		std::set<double> eSet, boundarySet;
		for(const Lspline& lspline: lsplines)
		{	boundarySet.insert(lspline.front().first);
			boundarySet.insert(lspline.back().first);
			for(unsigned j=0; j<lspline.size(); j++)
				eSet.insert(lspline[j].first);
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
		struct LsplineCmp
		{	typedef std::pair<double, std::vector<double> > LsplineElem;
			bool operator()(const LsplineElem& l1, const LsplineElem& l2) { return l1.first < l2.first; }
		};
		for(const Lspline& lspline: lsplines)
		{	Lspline::iterator cIter = std::upper_bound(combined.begin(), combined.end(), lspline.front(), LsplineCmp());
			cIter--; //Now this points to the second of the duplicate entries for the starting energy
			double eStop = lspline.back().first; //end energy of current input spline
			Lspline::const_iterator leftIter = lspline.begin();
			Lspline::const_iterator rightIter = leftIter + 1;
			while(true)
			{	double t = (cIter->first - leftIter->first) / (rightIter->first - leftIter->first);
				for(int i=0; i<nWeights; i++)
					cIter->second[i] += (1.-t) * leftIter->second[i] + t * rightIter->second[i];
				//Check for end of input spline:
				if(cIter->first == eStop) break; //note this ends at the earlier of the duplicate entries for the ending energy
				//Advance and update input interval if required:
				cIter++;
				if(cIter->first > rightIter->first)
				{	leftIter++;
					rightIter++;
				}
			}
		}
		return combined;
	}
	
	//Generate the density of states for a given state offset:
	Lspline getDOS(int stateOffset) const
	{	std::vector<Lspline> lsplines(nBands);
		for(int iBand=0; iBand<nBands; iBand++)
		{	Cspline wdos;
			for(const Tetrahedron& t: tetrahedra)
				accumTetrahedron(t, iBand, stateOffset, wdos);
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
	
	//Write the density of states to a file, for a given state offset:
	void printDOS(int stateOffset, string filename)
	{	logPrintf("Dumping '%s' ... ", filename.c_str()); logFlush();
		//Compute DOS:
		Lspline wdos = getDOS(stateOffset);
		//Output DOS:
		FILE* fp = fopen(filename.c_str(), "w");
		if(!fp) die("Could not open '%s' for writing.\n", filename.c_str());
		for(Lspline::const_iterator iter=wdos.begin(); iter!=wdos.end(); iter++)
		{	fprintf(fp, "%.15le", iter->first); 
			for(int i=0; i<nWeights; i++)
				fprintf(fp, "\t%.15le", iter->second[i]);
			fprintf(fp, "\n");
		}
		fclose(fp);
		logPrintf("Done.\n"); logFlush();
	}
	
	//Reduce an integer vector by its gcd:
	inline static int gcd(int x, int y)
	{	while(y != 0)
		{	int yPrev = y;
			y = x % y;
			x = yPrev;
		}
		return x;
	}
	inline static vector3<int> gcdReduce(const vector3<int>& d)
	{	int g = gcd(gcd(d[0], d[1]), d[2]);
		return vector3<int>(d[0]/g, d[1]/g, d[2]/g);
	}
	
	//Thread function for setting fourier transform of slice of half-width R centered at r0 parallel to lattice-plane d:
	inline static void sliceWeight_thread(int iStart, int iStop, const vector3<int>& S, const matrix3<>& GGT,
		complex* w, const vector3<>& r0, double R, const vector3<int>& d)
	{
		double planeSpacing = (2*M_PI) / sqrt(GGT.metric_length_squared(gcdReduce(d))); //lattice-plane spacing
		double prefac = (2*R) / planeSpacing;
		int dSq = d.length_squared();
		THREAD_halfGspaceLoop
		(	int iGsq = iG.length_squared();
			int iGdotd = dot(iG, d);
			if(iGsq * dSq == iGdotd * iGdotd) // => iG || d (Cauchy-Schwarz)
			{	double GR = R * sqrt(GGT.metric_length_squared(iG));
				w[i] = prefac * cis(-dot(iG,r0)) * bessel_jl(0,GR);
			}
			else w[i] = 0.;
		)
	}
	
	//Thread function for setting fourier transform of sphere of radius R centered at r0:
	inline static void sphereWeight_thread(int iStart, int iStop, const vector3<int>& S, const matrix3<>& GGT,
		complex* w, const vector3<>& r0, double R, double cellVolume)
	{
		double prefac = (4./3)*M_PI * pow(R,3) / cellVolume;
		THREAD_halfGspaceLoop
		(	double GR = R * sqrt(GGT.metric_length_squared(iG));
			w[i] = prefac * cis(-dot(iG,r0)) * (bessel_jl(0,GR) + bessel_jl(2,GR));
		)
	}
};


void DOS::dump()
{
	const ElecInfo& eInfo = e->eInfo;
	int nSpins = (eInfo.spinType==SpinNone) ? 1 : 2;
	int qCount = eInfo.qnums.size()/nSpins;
	
	//Create full-list of k-points (from those in symmetry-irreducible wedge):
	std::map<vector3<>, int> kpointMap; //map from unreduced kpoints to state index
	const std::vector<matrix3<int> >& sym = e->symm.getMatrices();
	const double kpointTol = 1e-4, kpointTolSq = kpointTol*kpointTol; //threshold for identifying kpoints
	for(int q=0; q<qCount; q++)
	{	vector3<> kState = eInfo.qnums[q].k;
		for(const matrix3<int>& mat: sym)
		{	vector3<> k = (~mat) * kState;
			//Reduce to centered zone (each reciprocal lattice coord in [-0.5,0.5))
			for(int i=0; i<3; i++)
			{	k[i] -= floor(k[i]);
				if(k[i]>=0.5)
					k[i] -= 1.;
			}
			//Check if this k-vector has already been encountered:
			bool found = false;
			for(const auto& kIndexPair: kpointMap)
				if(circDistanceSquared(k, kIndexPair.first) < kpointTolSq)
				{	found = true;
					break;
				}
			//Add to map if not yet encountered:
			if(!found) kpointMap[k] = q;
		}
	}
	
	//Triangulate the Brillouin zone
	std::vector< vector3<> > kpointVec; //flat list of kpoints
	std::vector<int> stateIndexVec; //flat list of state indices (same order as kpointVec)
	for(const auto& kIndexPair: kpointMap)
	{	kpointVec.push_back(kIndexPair.first);
		stateIndexVec.push_back(kIndexPair.second);
	}
	std::vector<Triangulation::Cell> cells;
	Triangulation::triangulate(kpointVec, cells); //wraps CGAL's Periodic 3D Delaunay triangulation
	
	//Convert the triangulation to an EvalDOS object:
	EvalDOS eval(cells.size(), weights.size(), eInfo.nStates, eInfo.nBands, Etol);
	double Vtot = 0.;
	for(unsigned c=0; c<cells.size(); c++)
	{	const Triangulation::Cell& cell = cells[c];
		EvalDOS::Tetrahedron& t = eval.tetrahedra[c];
		for(int i=0; i<4; i++) t.q[i] = stateIndexVec[cell.vertex[i].index];
		t.V = box(
			e->gInfo.GT * (cell.vertex[1].pos - cell.vertex[0].pos),
			e->gInfo.GT * (cell.vertex[2].pos - cell.vertex[0].pos),
			e->gInfo.GT * (cell.vertex[3].pos - cell.vertex[0].pos) );
		Vtot += t.V;
	}
	double VnormFac = 2./(nSpins*Vtot);
	for(EvalDOS::Tetrahedron& t: eval.tetrahedra)
		t.V *= VnormFac; //normalize volume of tetrahedra to add up to 2/nSpins
	
	//Compute the weight functions on the real-space grid:
	std::vector<DataRptr> weightFuncs(weights.size());
	for(unsigned iWeight=0; iWeight<weights.size(); iWeight++)
	{	Weight& weight = weights[iWeight];
		DataRptr& weightFunc = weightFuncs[iWeight];
		switch(weight.type)
		{	case Weight::Total:
				weightFunc = 0; //uniform weights, no field needed
				break;
			case Weight::Slice:
			case Weight::Sphere:
			case Weight::AtomSlice:
			case Weight::AtomSphere:
			{	//Get the center from atom position for Atom* modes:
				if(weight.type == Weight::AtomSlice || weight.type == Weight::AtomSphere)
					weight.center = e->iInfo.species[weight.specieIndex]->atpos[weight.atomIndex];
				//Compute weight function in Fourier space:
				const GridInfo& gInfo = e->gInfo;
				DataGptr weightFuncTilde(DataG::alloc(gInfo));
				complex* wData = weightFuncTilde->data();
				if(weight.type == Weight::Slice || weight.type == Weight::AtomSlice)
					threadLaunch(EvalDOS::sliceWeight_thread, gInfo.nG, gInfo.S, gInfo.GGT,
						wData, weight.center, weight.radius, weight.direction);
				else // (weight.type == Weight::Sphere || weight.type == Weight::Atomphere)
					threadLaunch(EvalDOS::sphereWeight_thread, gInfo.nG, gInfo.S, gInfo.GGT,
						wData, weight.center, weight.radius, fabs(gInfo.detR));
				//Store weight function in real space:
				weightFunc = I(weightFuncTilde, true);
				break;
			}
			case Weight::File:
			{	nullToZero(weightFunc, e->gInfo);
				loadRawBinary(weightFunc, weight.filename.c_str());
				break;
			}
			case Weight::Delim: break; //never encountered, just to suppress compiler warning
		}
	}
	
	//Compute the overlap of the density for each states and band, with each of the weight functions:
	for(int iState=0; iState<eInfo.nStates; iState++)
	{	for(int iBand=0; iBand<eInfo.nBands; iBand++)
		{	ColumnBundle C = e->eVars.C[iState].getSub(iBand, iBand+1);
			diagMatrix F(1, occupied ? e->eVars.F[iState][iBand] : 1.); //weight by filling if occupied, by 1 if total
			//Store the eigenvalue for this state and band:
			eval.e(iState, iBand) = e->eVars.Hsub_eigs[iState][iBand];
			//Compute the density for this state and band:
			DataRptr n = diagouterI(F, C);
			e->iInfo.augmentDensity(F, C, n); //pseudopotential contribution
			//Compute the weights:
			for(unsigned iWeight=0; iWeight<weights.size(); iWeight++)
				eval.w(iWeight, iState, iBand) = weightFuncs[iWeight]
					? e->gInfo.dV * dot(weightFuncs[iWeight], n) //mode with explicit wieght function
					: 1.; //Total DOS
		}
	}
	eval.weldEigenvalues();
	
	//Compute and print density of states:
	for(int iSpin=0; iSpin<nSpins; iSpin++)
		eval.printDOS(iSpin*qCount, e->dump.getFilename(nSpins==1 ? "dos" : (iSpin==0 ? "dosUp" : "dosDn")));
}

#else // TRIANGULATE_ENABLED

void DOS::dump()
{	string filename = e->dump.getFilename(e->eInfo.spinType==SpinNone ? "dos" : "dos*");
	logPrintf("Dumping '%s' ... Unavailable: please recompile JDFTx with CGAL enabled.\n", filename.c_str());
}

#endif // TRIANGULATE_ENABLED
