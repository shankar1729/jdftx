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
#include <core/SphericalHarmonics.h>
#include <electronic/ColumnBundle.h>
#include <core/ScalarFieldIO.h>
#include <core/LatticeUtils.h>
#include <array>

DOS::DOS() : Etol(1e-6), Esigma(0)
{
}

void DOS::setup(const Everything& everything)
{	e = &everything;
	if(!weights.size()) //Add the default of just total density of states:
	{	Weight weight;
		weight.type = Weight::Total;
		weight.fillingMode = Weight::Complete;
		weights.push_back(weight);
	}
	//Check compatibility of orbital modes with pseudopotentials:
	for(const Weight& weight: weights)
		if(weight.type == Weight::Orbital)
		{	const SpeciesInfo& sp = *(e->iInfo.species[weight.specieIndex]);
			const Weight::OrbitalDesc& oDesc = weight.orbitalDesc;
			int lMax = sp.lMaxAtomicOrbitals();
			if(lMax < 0)
				die("Species %s has no atomic orbitals and cannot be used for orbital-projected density of states.\n", sp.name.c_str());
			if(oDesc.l > lMax)
				die("Angular momentum of density-of-states projection orbital %s exceeds lMax=%d of species %s.\n",
					string(oDesc).c_str(), lMax, sp.name.c_str());
			int nAOl = sp.nAtomicOrbitals(oDesc.l);
			if(oDesc.n >= unsigned(nAOl))
				die("Principal quantum number (%d) of density-of-states projection orbital %s exceeds maximum value (%d) for l=%d of species %s.\n",
					oDesc.n+1, string(oDesc).c_str(), nAOl, oDesc.l, sp.name.c_str());
		}
	Citations::add("Linear-tetrahedron sampling for density of states",
		"G. Lehmann and M. Taut, Phys. status solidi (b) 54, 469 (1972)");
}

string DOS::Weight::getDescription(const Everything& e) const
{	ostringstream oss;
	if(type==Total) oss << "Total";
	if(type==Slice || type==AtomSlice)
		oss << "(" << direction[0] << "," << direction[1] << "," << direction[2] << ") slice of half-width ";
	if(type==Sphere || type==AtomSphere) oss << "Sphere of radius ";
	if(type==Slice || type==AtomSlice || type==Sphere || type==AtomSphere)
		oss << radius << " bohr at ";
	if(type==Orbital)
		oss << orbitalDesc << " orbital at ";
	if(type==OrthoOrbital)
		oss << orbitalDesc << " orthonormalized-orbital at ";
	if(type==Sphere || type==Slice)
	{	vector3<> c;
		if(e.iInfo.coordsType == CoordsLattice)
		{	c = center;
			oss << "lattice";
		}
		else
		{	c = e.gInfo.invR * center;
			oss << "cartesian";
		}
		oss << " (" << c[0] << "," << c[1] << "," << c[2] << ")";
	}
	if(type==AtomSphere || type==AtomSlice || type==Orbital || type==OrthoOrbital)
		oss << e.iInfo.species[specieIndex]->name << " #" << (atomIndex+1);
	if(type==File)
		oss << "Weighted by '" << filename << "'";
	if(fillingMode == Occupied)
		oss << " (Occupied)";
	return oss.str();
}

//Connection between string representations of orbitals and (l,m) pairs:
const EnumStringMap<DOS::Weight::OrbitalDesc>& getOrbitalDescMap()
{	static EnumStringMap<DOS::Weight::OrbitalDesc> orbitalDescMap
	(	//Spin-less specification (allowed for all spin modes for nonrelativistic pseudopotentials)
		DOS::Weight::OrbitalDesc(0,0, 0, SpinNone), "s",
		DOS::Weight::OrbitalDesc(1,-1, 0, SpinNone), "py",
		DOS::Weight::OrbitalDesc(1, 0, 0, SpinNone), "pz",
		DOS::Weight::OrbitalDesc(1, 1, 0, SpinNone), "px",
		DOS::Weight::OrbitalDesc(1,2, 0, SpinNone), "p",
		DOS::Weight::OrbitalDesc(2,-2, 0, SpinNone), "dxy",
		DOS::Weight::OrbitalDesc(2,-1, 0, SpinNone), "dyz",
		DOS::Weight::OrbitalDesc(2, 0, 0, SpinNone), "dz2",
		DOS::Weight::OrbitalDesc(2, 1, 0, SpinNone), "dxz",
		DOS::Weight::OrbitalDesc(2, 2, 0, SpinNone), "dx2-y2",
		DOS::Weight::OrbitalDesc(2,3, 0, SpinNone), "d",
		DOS::Weight::OrbitalDesc(3,-3, 0, SpinNone), "fy(3x2-y2)",
		DOS::Weight::OrbitalDesc(3,-2, 0, SpinNone), "fxyz",
		DOS::Weight::OrbitalDesc(3,-1, 0, SpinNone), "fyz2",
		DOS::Weight::OrbitalDesc(3, 0, 0, SpinNone), "fz3",
		DOS::Weight::OrbitalDesc(3, 1, 0, SpinNone), "fxz2",
		DOS::Weight::OrbitalDesc(3, 2, 0, SpinNone), "fz(x2-y2)",
		DOS::Weight::OrbitalDesc(3, 3, 0, SpinNone), "fx(x2-3y2)",
		DOS::Weight::OrbitalDesc(3,4, 0, SpinNone), "f",
		//Up/dn specification: noncollinear modes with nonrelativistic pseudopotentials only
		DOS::Weight::OrbitalDesc(0,0, 0, SpinZ), "sUp",
		DOS::Weight::OrbitalDesc(0,0, 1, SpinZ), "sDn",
		DOS::Weight::OrbitalDesc(1,-1, 0, SpinZ), "pyUp",
		DOS::Weight::OrbitalDesc(1,-1, 1, SpinZ), "pyDn",
		DOS::Weight::OrbitalDesc(1, 0, 0, SpinZ), "pzUp",
		DOS::Weight::OrbitalDesc(1, 0, 1, SpinZ), "pzDn",
		DOS::Weight::OrbitalDesc(1, 1, 0, SpinZ), "pxUp",
		DOS::Weight::OrbitalDesc(1, 1, 1, SpinZ), "pxDn",
		DOS::Weight::OrbitalDesc(1,2, 0, SpinZ), "pUp",
		DOS::Weight::OrbitalDesc(1,2, 1, SpinZ), "pDn",
		DOS::Weight::OrbitalDesc(2,-2, 0, SpinZ), "dxyUp",
		DOS::Weight::OrbitalDesc(2,-2, 1, SpinZ), "dxyDn",
		DOS::Weight::OrbitalDesc(2,-1, 0, SpinZ), "dyzUp",
		DOS::Weight::OrbitalDesc(2,-1, 1, SpinZ), "dyzDn",
		DOS::Weight::OrbitalDesc(2, 0, 0, SpinZ), "dz2Up",
		DOS::Weight::OrbitalDesc(2, 0, 1, SpinZ), "dz2Dn",
		DOS::Weight::OrbitalDesc(2, 1, 0, SpinZ), "dxzUp",
		DOS::Weight::OrbitalDesc(2, 1, 1, SpinZ), "dxzDn",
		DOS::Weight::OrbitalDesc(2, 2, 0, SpinZ), "dx2-y2Up",
		DOS::Weight::OrbitalDesc(2, 2, 1, SpinZ), "dx2-y2Dn",
		DOS::Weight::OrbitalDesc(2,3, 0, SpinZ), "dUp",
		DOS::Weight::OrbitalDesc(2,3, 1, SpinZ), "dDn",
		DOS::Weight::OrbitalDesc(3,-3, 0, SpinZ), "fy(3x2-y2)Up",
		DOS::Weight::OrbitalDesc(3,-3, 1, SpinZ), "fy(3x2-y2)Dn",
		DOS::Weight::OrbitalDesc(3,-2, 0, SpinZ), "fxyzUp",
		DOS::Weight::OrbitalDesc(3,-2, 1, SpinZ), "fxyzDn",
		DOS::Weight::OrbitalDesc(3,-1, 0, SpinZ), "fyz2Up",
		DOS::Weight::OrbitalDesc(3,-1, 1, SpinZ), "fyz2Dn",
		DOS::Weight::OrbitalDesc(3, 0, 0, SpinZ), "fz3Up",
		DOS::Weight::OrbitalDesc(3, 0, 1, SpinZ), "fz3Dn",
		DOS::Weight::OrbitalDesc(3, 1, 0, SpinZ), "fxz2Up",
		DOS::Weight::OrbitalDesc(3, 1, 1, SpinZ), "fxz2Dn",
		DOS::Weight::OrbitalDesc(3, 2, 0, SpinZ), "fz(x2-y2)Up",
		DOS::Weight::OrbitalDesc(3, 2, 1, SpinZ), "fz(x2-y2)Dn",
		DOS::Weight::OrbitalDesc(3, 3, 0, SpinZ), "fx(x2-3y2)Up",
		DOS::Weight::OrbitalDesc(3, 3, 1, SpinZ), "fx(x2-3y2)Dn",
		DOS::Weight::OrbitalDesc(3,4, 0, SpinZ), "fUp",
		DOS::Weight::OrbitalDesc(3,4, 1, SpinZ), "fDn",
		//j,mj specification (only option for spin-split channels of relativsitic pseudopotentials)
		DOS::Weight::OrbitalDesc(0,-1, 0, SpinOrbit), "s(-1/2)",
		DOS::Weight::OrbitalDesc(0,+0, 0, SpinOrbit), "s(+1/2)",
		DOS::Weight::OrbitalDesc(1,-2, 0, SpinOrbit), "p+(-3/2)",
		DOS::Weight::OrbitalDesc(1,-1, 0, SpinOrbit), "p+(-1/2)",
		DOS::Weight::OrbitalDesc(1,+0, 0, SpinOrbit), "p+(+1/2)",
		DOS::Weight::OrbitalDesc(1,+1, 0, SpinOrbit), "p+(+3/2)",
		DOS::Weight::OrbitalDesc(1,2, 0, SpinOrbit), "p+",
		DOS::Weight::OrbitalDesc(1,+0, 1, SpinOrbit), "p-(-1/2)",
		DOS::Weight::OrbitalDesc(1,+1, 1, SpinOrbit), "p-(+1/2)",
		DOS::Weight::OrbitalDesc(1,2, 1, SpinOrbit), "p-",
		DOS::Weight::OrbitalDesc(2,-3, 0, SpinOrbit), "d+(-5/2)",
		DOS::Weight::OrbitalDesc(2,-2, 0, SpinOrbit), "d+(-3/2)",
		DOS::Weight::OrbitalDesc(2,-1, 0, SpinOrbit), "d+(-1/2)",
		DOS::Weight::OrbitalDesc(2,+0, 0, SpinOrbit), "d+(+1/2)",
		DOS::Weight::OrbitalDesc(2,+1, 0, SpinOrbit), "d+(+3/2)",
		DOS::Weight::OrbitalDesc(2,+2, 0, SpinOrbit), "d+(+5/2)",
		DOS::Weight::OrbitalDesc(2,3, 0, SpinOrbit), "d+",
		DOS::Weight::OrbitalDesc(2,-1, 1, SpinOrbit), "d-(-3/2)",
		DOS::Weight::OrbitalDesc(2,+0, 1, SpinOrbit), "d-(-1/2)",
		DOS::Weight::OrbitalDesc(2,+1, 1, SpinOrbit), "d-(+1/2)",
		DOS::Weight::OrbitalDesc(2,+2, 1, SpinOrbit), "d-(+3/2)",
		DOS::Weight::OrbitalDesc(2,3, 1, SpinOrbit), "d-",
		DOS::Weight::OrbitalDesc(3,-4, 0, SpinOrbit), "f+(-7/2)",
		DOS::Weight::OrbitalDesc(3,-3, 0, SpinOrbit), "f+(-5/2)",
		DOS::Weight::OrbitalDesc(3,-2, 0, SpinOrbit), "f+(-3/2)",
		DOS::Weight::OrbitalDesc(3,-1, 0, SpinOrbit), "f+(-1/2)",
		DOS::Weight::OrbitalDesc(3,+0, 0, SpinOrbit), "f+(+1/2)",
		DOS::Weight::OrbitalDesc(3,+1, 0, SpinOrbit), "f+(+3/2)",
		DOS::Weight::OrbitalDesc(3,+2, 0, SpinOrbit), "f+(+5/2)",
		DOS::Weight::OrbitalDesc(3,+3, 0, SpinOrbit), "f+(+7/2)",
		DOS::Weight::OrbitalDesc(3,4, 0, SpinOrbit), "f+",
		DOS::Weight::OrbitalDesc(3,-2, 1, SpinOrbit), "f-(-5/2)",
		DOS::Weight::OrbitalDesc(3,-1, 1, SpinOrbit), "f-(-3/2)",
		DOS::Weight::OrbitalDesc(3,+0, 1, SpinOrbit), "f-(-1/2)",
		DOS::Weight::OrbitalDesc(3,+1, 1, SpinOrbit), "f-(+1/2)",
		DOS::Weight::OrbitalDesc(3,+2, 1, SpinOrbit), "f-(+3/2)",
		DOS::Weight::OrbitalDesc(3,+3, 1, SpinOrbit), "f-(+5/2)",
		DOS::Weight::OrbitalDesc(3,4, 1, SpinOrbit), "f-"
	);
	return orbitalDescMap;
}

void DOS::Weight::OrbitalDesc::parse(string desc)
{	//SPlit into principal quantum number and orbital code:
	size_t orbStart = desc.find_first_not_of("0123456789");
	if(orbStart == string::npos)
		throw "Orbital description '" + desc + "' does not have an orbital code.";
	string orbCode = desc.substr(orbStart);
	//Convert orbital code:
	if(!getOrbitalDescMap().getEnum(orbCode.c_str(), *this))
		throw "'" + orbCode + "' is not a valid orbital code for an orbital description.";
	//Extract the principal quantum number:
	if(orbStart)
	{	int nPlus1 = atoi(desc.substr(0, orbStart).c_str());
		if(nPlus1<=0) throw string("Principal quantum number in orbital description must be a positive integer");
		n = nPlus1 - 1;
	}
}

DOS::Weight::OrbitalDesc::operator string() const
{	ostringstream oss;
	if(n>0) oss << (n+1); //optional pseudo-principal quantum number 
	oss << getOrbitalDescMap().getString(*this);
	return oss.str();
}

bool DOS::Weight::OrbitalDesc::operator<(const DOS::Weight::OrbitalDesc& other) const
{	if(l != other.l) return l < other.l;
	if(m != other.m) return m < other.m;
	if(s != other.s) return s < other.s;
	if(n != other.n) return n < other.n;
	return spinType < other.spinType;
}

//Density of states evaluation helper class
struct EvalDOS
{
	//Tetrahedron in Brillouin zone triangulation:
	struct Tetrahedron
	{	std::array<int,4> q; //state indices for the vertices
		double V; //volume of tetrahedron
	};
	std::vector<Tetrahedron> tetrahedra;
	
	int nWeights, nStates, nBands; double Etol, Esigma;
	std::vector<double> eigs; //flat array of eigenvalues (inner index state, outer index bands)
	std::vector<double> weights; //flat array of DOS weights (inner index weight function, middle index state, and outer index bands)
	double& e(int iState, int iBand) { return eigs[iState + nStates*iBand]; } //access eigenvalue
	const double& e(int iState, int iBand) const { return eigs[iState + nStates*iBand]; } //access eigenvalue (const version)
	double& w(int iWeight, int iState, int iBand) { return weights[iWeight + nWeights*(iState + nStates*iBand)]; } //access weight
	const double& w(int iWeight, int iState, int iBand) const { return weights[iWeight + nWeights*(iState + nStates*iBand)]; } //access weight (const version)
	
	EvalDOS(int nCells, int nWeights, int nStates, int nBands, double Etol, double Esigma)
	: tetrahedra(nCells),
	nWeights(nWeights), nStates(nStates), nBands(nBands), Etol(Etol), Esigma(Esigma),
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
	
	typedef std::pair<double, std::vector<double> > LsplineElem;
	typedef std::vector<LsplineElem> Lspline; //set of linear splines (energies, and vector of DOS values for each weight function)
	static bool LsplineCmp(const LsplineElem& l1, const LsplineElem& l2) { return l1.first < l2.first; }
	
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
		for(const LsplineElem& delta: cspline.deltas)
		{	Lspline::iterator lIter = std::lower_bound(lspline.begin(), lspline.end(), delta, LsplineCmp);
			for(int i=0; i<nWeights; i++)
				lIter->second[i] += delta.second[i];
		}
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
	
	//Apply gaussian smoothing (if Esigma is non-zero)
	Lspline gaussSmooth(Lspline in) const
	{	if(!Esigma) return in;
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
		return gaussSmooth(mergeLsplines(lsplines));
	}
	
	//Write the density of states to a file, for a given state offset:
	void printDOS(int stateOffset, string filename, string header)
	{	logPrintf("Dumping '%s' ... ", filename.c_str()); logFlush();
		//Compute DOS:
		Lspline wdos = getDOS(stateOffset);
		//Output DOS:
		FILE* fp = fopen(filename.c_str(), "w");
		if(!fp) die("Could not open '%s' for writing.\n", filename.c_str());
		//Header:
		fprintf(fp, "%s\n", header.c_str());
		//Data:
		for(Lspline::const_iterator iter=wdos.begin(); iter!=wdos.end(); iter++)
		{	fprintf(fp, "%.15le", iter->first); 
			for(int i=0; i<nWeights; i++)
				fprintf(fp, "\t%.15le", iter->second[i]);
			fprintf(fp, "\n");
		}
		fclose(fp);
		logPrintf("done.\n"); logFlush();
	}
	
	//Thread function for setting fourier transform of slice of half-width R centered at r0 parallel to lattice-plane d:
	inline static void sliceWeight_thread(size_t iStart, size_t iStop, const vector3<int>& S, const matrix3<>& GGT,
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
				w[i] = prefac * cis(-(2.*M_PI)*dot(iG,r0)) * bessel_jl(0,GR);
			}
			else w[i] = 0.;
		)
	}
	
	//Thread function for setting fourier transform of sphere of radius R centered at r0:
	inline static void sphereWeight_thread(size_t iStart, size_t iStop, const vector3<int>& S, const matrix3<>& GGT,
		complex* w, const vector3<>& r0, double R, double cellVolume)
	{
		double prefac = (4./3)*M_PI * pow(R,3) / cellVolume;
		THREAD_halfGspaceLoop
		(	double GR = R * sqrt(GGT.metric_length_squared(iG));
			w[i] = prefac * cis(-(2.*M_PI)*dot(iG,r0)) * (bessel_jl(0,GR) + bessel_jl(2,GR));
		)
	}
};


//Check that a 3-vector i sintegral within tolerance and return integer version
inline vector3<int> round(const vector3<> v, const double tol)
{	double err;
	vector3<int> vInt = round(v, &err);
	assert(err < tol);
	return vInt;
}

void DOS::dump()
{
	const ElecInfo& eInfo = e->eInfo;
	int nSpins = eInfo.nSpins();
	int qCount = eInfo.qnums.size()/nSpins;
	const Supercell& supercell = *(e->coulombParams.supercell);
	
	//Pick optimum tetrahedral tesselation of each parallelopiped cell in BZ:
	vector3<> dk[8]; //k-point offsets of basis parallelopiped relative to its first vertex
	{	//There are 4 distinct tesselations listed by signs of each basis vector (upto an overall sign)
		//Pick the one with the shortest common body-diagonal to yield the most suitable aspect ratios
		vector3<int> sBest;
		double dMin = DBL_MAX;
		matrix3<> GTsup = ~((2*M_PI) * inv(supercell.Rsuper));
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
		matrix3<> invGT_GTsup = inv(e->gInfo.GT) * GTsup;
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
	
	//Construct Brillouin zone triangulation in an EvalDOS object:
	//--- lookup table of k-points using (integral) reciprocal superlattice coordinates
	std::map<vector3<int>,int> kpointMap;
	const std::vector<vector3<>>& kmesh = supercell.kmesh;
	for(unsigned i=0; i<kmesh.size(); i++)
	{	//Add kRange and its periodic images if it's too close to the edge (to catch roundoff issues in wrapping):
		vector3<> kRange[2];
		kRange[0] = kRange[1] = kmesh[i];
		for(int j=0; j<3; j++)
		{	if(fabs(+0.5-kRange[0][j])<symmThreshold) kRange[1][j] -= 1.;
			if(fabs(-0.5-kRange[0][j])<symmThreshold) kRange[1][j] += 1.;
		}
		for(int s0=0; s0<2; s0++)
		for(int s1=0; s1<2; s1++)
		for(int s2=0; s2<2; s2++)
			kpointMap[round((vector3<>(kRange[s0][0],kRange[s1][1],kRange[s2][2])-kmesh[0])*supercell.super, symmThreshold)] = i;
	}
	//--- add 6 tetrahedra per parallelopiped cell
	EvalDOS eval(6*kmesh.size(), weights.size(), eInfo.nStates, eInfo.nBands, Etol, Esigma);
	double Vtot = 0.;
	for(unsigned i=0; i<kmesh.size(); i++)
	{	const vector3<>& v0 = kmesh[i];
		//initialize the parallelopiped vertices and look up their indices:
		vector3<> v[8]; int iv[8];
		for(int j=0; j<8; j++)
		{	v[j] = v0 + dk[j];
			vector3<> vWrapped; for(int k=0; k<3; k++) vWrapped[k] = v[j][k] - floor(0.5+v[j][k]);
			auto ivIter = kpointMap.find(round((vWrapped-kmesh[0])*supercell.super, symmThreshold));
			assert(ivIter != kpointMap.end());
			iv[j] = ivIter->second;
		}
		//loop over tetrahedra:
		for(unsigned t=0; t<6; t++)
		{	EvalDOS::Tetrahedron& tet = eval.tetrahedra[6*i+t];
			for(int p=0; p<4; p++)
			{	const Supercell::KmeshTransform& kTransform = supercell.kmeshTransform[iv[cellVerts[t][p]]];
				tet.q[p] = kTransform.iReduced;
			}
			tet.V = box(
				e->gInfo.GT * (v[cellVerts[t][1]] - v[cellVerts[t][0]]),
				e->gInfo.GT * (v[cellVerts[t][2]] - v[cellVerts[t][0]]),
				e->gInfo.GT * (v[cellVerts[t][3]] - v[cellVerts[t][0]]) );
			Vtot += tet.V;
		}
	}
	double VnormFac = eInfo.qWeightSum/(nSpins*Vtot);
	for(EvalDOS::Tetrahedron& t: eval.tetrahedra)
		t.V *= VnormFac; //normalize volume of tetrahedra to add up to qWeightSum/nSpins
	
	//Set uniform weights for Total mode:
	for(int iState=eInfo.qStart; iState<eInfo.qStop; iState++)
		for(int iBand=0; iBand<eInfo.nBands; iBand++)
			for(unsigned iWeight=0; iWeight<weights.size(); iWeight++)
				if(weights[iWeight].type == Weight::Total)
					eval.w(iWeight, iState, iBand) = 1.;
	
	//Compute the weight functions on the real-space grid for weighted-density modes:
	std::vector<ScalarFieldArray> weightFuncs(weights.size());
	bool needDensity = false; //check if any of the modes need the density
	for(unsigned iWeight=0; iWeight<weights.size(); iWeight++)
	{	Weight& weight = weights[iWeight];
		ScalarFieldArray& weightFunc = weightFuncs[iWeight];
		switch(weight.type)
		{	case Weight::Slice:
			case Weight::Sphere:
			case Weight::AtomSlice:
			case Weight::AtomSphere:
			{	//Get the center from atom position for Atom* modes:
				if(weight.type == Weight::AtomSlice || weight.type == Weight::AtomSphere)
					weight.center = e->iInfo.species[weight.specieIndex]->atpos[weight.atomIndex];
				//Compute weight function in Fourier space:
				const GridInfo& gInfo = e->gInfo;
				ScalarFieldTilde weightFuncTilde(ScalarFieldTildeData::alloc(gInfo));
				complex* wData = weightFuncTilde->data();
				if(weight.type == Weight::Slice || weight.type == Weight::AtomSlice)
					threadLaunch(EvalDOS::sliceWeight_thread, gInfo.nG, gInfo.S, gInfo.GGT,
						wData, weight.center, weight.radius, weight.direction);
				else // (weight.type == Weight::Sphere || weight.type == Weight::Atomphere)
					threadLaunch(EvalDOS::sphereWeight_thread, gInfo.nG, gInfo.S, gInfo.GGT,
						wData, weight.center, weight.radius, fabs(gInfo.detR));
				//Store weight function in real space:
				weightFunc.resize(1);
				weightFunc[0] = I(weightFuncTilde);
				needDensity = true;
				break;
			}
			case Weight::File:
			{	nullToZero(weightFunc, e->gInfo, 1);
				loadRawBinary(weightFunc[0], weight.filename.c_str());
				needDensity = true;
				break;
			}
			default: //Total and orbital modes do not weight the density:
			{	break;
			}
		}
		if(eInfo.spinType == SpinZ && weightFunc.size())
		{	weightFunc.resize(eInfo.nDensities);
			weightFunc[1] = weightFunc[0];
		}
		if(eInfo.spinType == SpinVector)
		{	if(weight.Mhat.length() && weight.type==Weight::Total)
			{	//treat spin-projected total DOS as a density-weighted DOS with spatially uniform weight
				nullToZero(weightFunc, e->gInfo, 1);
				weightFunc[0] += 1.;
				needDensity = true;
			}
			if(weightFunc.size()) //convert to the appropriate 4-array depending on magnetization
			{	ScalarField w; std::swap(w, weightFunc[0]);
				weightFunc.resize(eInfo.nDensities);
				const vector3<>& M = weight.Mhat;
				weightFunc[0] = w * (1. + M[2]); //UpUp
				weightFunc[1] = w * (1. - M[2]); //DnDn
				weightFunc[2] = w * (+2.*M[0]); //Re(UpDn)
				weightFunc[3] = w * (-2.*M[1]); //Im(UpDn)
				if(M.length_squared()) weightFunc *= 0.5; //count only one spin
			}
		}
	}
	//Compute the overlap of the density for each states and band, with each of the density weight functions:
	if(needDensity)
	{	//Direct contributions from bands:
		for(int iState=eInfo.qStart; iState<eInfo.qStop; iState++)
		{	for(int iBand=0; iBand<eInfo.nBands; iBand++)
			{	//Compute the density for this state and band:
				diagMatrix F(1, 1.); //compute density with filling=1; incorporate fillings later per weight function if required
				ColumnBundle C = e->eVars.C[iState].getSub(iBand, iBand+1);
				ScalarFieldArray n = diagouterI(F, C, eInfo.nDensities, &e->gInfo);
				//Compute the weights:
				for(unsigned iWeight=0; iWeight<weights.size(); iWeight++)
					if(weightFuncs[iWeight].size())
						eval.w(iWeight, iState, iBand) = e->gInfo.dV * dot(weightFuncs[iWeight], n);
			}
		}
		//Ultrasoft augmentation:
		for(unsigned iWeight=0; iWeight<weights.size(); iWeight++)
			if(weightFuncs[iWeight].size())
			{	//Project weight functions from grid to atom-centered spherical functions:
				e->iInfo.augmentDensityGridGrad(weightFuncs[iWeight]);
				//Propagate from spherical functions to bands:
				diagMatrix Fq(eInfo.nBands, 1.); //actual fillings incorporated later if required
				for(int iState=eInfo.qStart; iState<eInfo.qStop; iState++)
				{	const std::vector<matrix>& VdagCq = e->eVars.VdagC[iState];
					std::vector<matrix> wVdagCq(e->iInfo.species.size());
					const QuantumNumber& qnum = e->eInfo.qnums[iState];
					e->iInfo.augmentDensitySphericalGrad(qnum, Fq, VdagCq, wVdagCq);
					//Sum diagonal contributions per band over species
					diagMatrix wAug(e->eInfo.nBands, 0.);
					for(size_t sp=0; sp<VdagCq.size(); sp++)
						if(wVdagCq[sp]) wAug += diag(dagger(VdagCq[sp]) * wVdagCq[sp]);
					//Augment with appropriate prefactors:
					for(int iBand=0; iBand<eInfo.nBands; iBand++)
						eval.w(iWeight, iState, iBand) += e->gInfo.dV * wAug[iBand];
				}
			}
	}
	
	//Check which orbital modes are needed:
	bool needOrbitals = false;
	bool needOrthoOrbitals = false;
	for(const Weight& weight: weights)
	{	if(weight.type == Weight::Orbital) needOrbitals = true;
		if(weight.type == Weight::OrthoOrbital) needOrthoOrbitals = true;
	}
	
	//Compute projections for orbital mode:
	if(needOrbitals || needOrthoOrbitals)
	{	for(int iState=eInfo.qStart; iState<eInfo.qStop; iState++)
		{	const ColumnBundle& C = e->eVars.C[iState];
			//Count atomic orbitals:
			int nOrbitals = 0;
			std::vector<int> spOffset(e->iInfo.species.size()); //species offset into atomic orbitals list
			for(unsigned sp=0; sp<e->iInfo.species.size(); sp++)
			{	spOffset[sp] = nOrbitals;
				nOrbitals += e->iInfo.species[sp]->nAtomicOrbitals();
			}
			//Atomic-orbital projections if needed:
			matrix CdagOpsi;
			if(needOrbitals)
				CdagOpsi = C ^ e->iInfo.getAtomicOrbitals(iState, true); //directly calculate overlap using Opsi
			//Ortho-orbital projections if needed:
			matrix CdagOpsiOrtho;
			if(needOrthoOrbitals)
			{	ColumnBundle psi = e->iInfo.getAtomicOrbitals(iState, false);
				//Orthogonalize:
				ColumnBundle Opsi = O(psi);
				matrix orthoMat = invsqrt(psi ^ Opsi); //orthonormalizing matrix
				CdagOpsiOrtho = (C ^ Opsi) * orthoMat;
			}
			for(unsigned iWeight=0; iWeight<weights.size(); iWeight++)
			{	const Weight& weight = weights[iWeight];
				if(weight.type == Weight::Orbital || weight.type == Weight::OrthoOrbital)
				{	const Weight::OrbitalDesc& oDesc = weight.orbitalDesc;
					const matrix& CdagOpsiRel = weight.type==Weight::Orbital ? CdagOpsi : CdagOpsiOrtho;
					int l = oDesc.l;
					//Check relativistic psp compatibility:
					if(e->iInfo.species[weight.specieIndex]->isRelativistic() && oDesc.l)
					{	if(!(oDesc.spinType==SpinOrbit || (oDesc.spinType==SpinNone && oDesc.m==l+1)))
							die("Individual orbital projections for l>0 in relativsitic pseudopotentials must use the j,mj specification (eg. 'p+(+1/2)').\n");
					}
					else
					{	if(oDesc.spinType==SpinOrbit)
							die("Only l>0 projections in relativistic pseudopotentials may use the j,mj specification.\n");
					}
					//Accumulate weights:
					int sStart, sStop;
					if(oDesc.spinType==SpinNone)
					{	sStart=0;
						sStop=eInfo.spinorLength();
					}
					else
					{	sStart=oDesc.s;
						sStop=oDesc.s+1;
					}
					for(int s=sStart; s<sStop; s++)
					{	int mMin, mMax;
						if(oDesc.m==l+1) //multi-orbital mode
						{	mMin = -l;
							mMax = +l;
							if(e->iInfo.species[weight.specieIndex]->isRelativistic())
								mMin -= (s ? -1 : +1);
						}
						else //single orbital mode
						{	mMin = oDesc.m;
							mMax = oDesc.m;
						}
						for(int m=mMin; m<=mMax; m++)
						{	int iCol = spOffset[weight.specieIndex] + e->iInfo.species[weight.specieIndex]->atomicOrbitalOffset(weight.atomIndex, oDesc.n, l, m, s);
							for(int iBand=0; iBand<eInfo.nBands; iBand++)
								eval.w(iWeight, iState, iBand) += (CdagOpsiRel.data()[CdagOpsiRel.index(iBand, iCol)]).norm();
						}
					}
				}
			}
		}
	}
	
	//Apply the filling weights and set the eigenvalues:
	for(int iState=eInfo.qStart; iState<eInfo.qStop; iState++)
		for(int iBand=0; iBand<eInfo.nBands; iBand++)
		{	for(unsigned iWeight=0; iWeight<weights.size(); iWeight++)
				if(weights[iWeight].fillingMode == Weight::Occupied)
					eval.w(iWeight, iState, iBand) *= e->eVars.F[iState][iBand];
			eval.e(iState, iBand) = e->eVars.Hsub_eigs[iState][iBand];
		}
		
	//Synchronize eigenvalues and weights between processes:
	if(mpiUtil->nProcesses()>1)
	{	if(mpiUtil->isHead())
		{	for(int iSrc=1; iSrc<mpiUtil->nProcesses(); iSrc++)
			{	int qStart = eInfo.qStartOther(iSrc);
				int qStop = eInfo.qStopOther(iSrc);
				std::vector<double> message((qStop-qStart)*eInfo.nBands*(weights.size()+1));
				mpiUtil->recv(message.data(), message.size(), iSrc, 0);
				const double* messagePtr = message.data();
				for(int iState=qStart; iState<qStop; iState++)
					for(int iBand=0; iBand<eInfo.nBands; iBand++)
					{	for(unsigned iWeight=0; iWeight<weights.size(); iWeight++)
							eval.w(iWeight, iState, iBand) = *(messagePtr++);
						eval.e(iState, iBand) = *(messagePtr++);
					}
			}
		}
		else
		{	//Pack data into a single message:
			std::vector<double> message((eInfo.qStop-eInfo.qStart)*eInfo.nBands*(weights.size()+1));
			double* messagePtr = message.data();
			for(int iState=eInfo.qStart; iState<eInfo.qStop; iState++)
				for(int iBand=0; iBand<eInfo.nBands; iBand++)
				{	for(unsigned iWeight=0; iWeight<weights.size(); iWeight++)
						*(messagePtr++) = eval.w(iWeight, iState, iBand);
					*(messagePtr++) = eval.e(iState, iBand);
				}
			mpiUtil->send(message.data(), message.size(), 0, 0);
		}
	}
	if(!mpiUtil->isHead()) return;
	
	//Compute and print density of states (head only):
	string header = "\"Energy\"";
	for(const Weight& weight: weights)
		header += ("\t\"" + weight.getDescription(*e) + "\"");
	eval.weldEigenvalues();
	for(int iSpin=0; iSpin<nSpins; iSpin++)
		eval.printDOS(iSpin*qCount, e->dump.getFilename(nSpins==1 ? "dos" : (iSpin==0 ? "dosUp" : "dosDn")), header);
}
