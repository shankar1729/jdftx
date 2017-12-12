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

#include <electronic/RadialSchrodinger.h>
#include <core/Util.h>
#include <stack>
#include <cfloat>
#include <algorithm>

RadialSchrodinger::RadialSchrodinger(const std::vector<double>& rArr, const std::vector<double>& drArr,
	const std::vector<double>& V, double Z, size_t iMatch)
: rArr(rArr), drArr(drArr), Z(Z), iMatch(iMatch ? iMatch : rArr.size()/2)
{
	assert(rArr.size()==drArr.size());
	assert(drArr.size()==V.size());
	assert(iMatch<rArr.size());
	
	//Set up cubic spline  triadiagonal system for V with natural B.C.
	std::vector<double> a(rArr.size(), 0.); //sub-diagonal (first entry ignored)
	std::vector<double> b(rArr.size(), 0.); //diagonal
	std::vector<double> c(rArr.size(), 0.); //super-diagonal (last entry ignored)
	std::vector<double> d(rArr.size(), 0.); //rArr.h.s
	for(size_t i=0; i<rArr.size()-1; i++)
	{	//Contribution from (rArr[i],rArr[i+1]) interval:
		double hInv = 1./(rArr[i+1]-rArr[i]);
		double rhs = 3*hInv*hInv*(V[i+1]-V[i]);
		b[i] += 2*hInv;   c[i] += hInv;     d[i] += rhs; //to i'th equation
		a[i+1] += hInv; b[i+1] += 2*hInv; d[i+1] += rhs; //to i+1'th equation
	}
	//Solve tridiagonal matrix equation to fully determine spline V(rArr):
	for(size_t i=1; i<rArr.size(); i++)
	{	double frac = a[i]/b[i-1];
		b[i] -= frac*c[i-1];
		d[i] -= frac*d[i-1];
	}
	std::vector<double> Vp(rArr.size()); //values of dV/dr at each node
	Vp[rArr.size()-1] = d[rArr.size()-1]/b[rArr.size()-1];
	for(size_t i=rArr.size()-1; i>0; i--)
		Vp[i-1] = (d[i-1] - c[i-1]*Vp[i])/b[i-1];
	//Sub-sample the smooth part of the potential using cubic splines (required for RK4):
	Vsub.resize(2*V.size()-1);
	for(size_t i=0; i<V.size()-1; i++)
	{	Vsub[2*i] = V[i];
		Vsub[2*i+1] = 0.5*(V[i]+V[i+1]) - 0.125*(Vp[i]-Vp[i-1])*(rArr[i]-rArr[i-1]);
	}
	Vsub.back() = V.back();
}

// Determine the optimum eigenfunction fillings (indexed by l and then nNodes)
std::vector< std::vector<double> > RadialSchrodinger::getFillings(double nElectrons,
	std::vector< std::vector<double> >* Eptr)
{
	std::vector< std::vector<double> > F;
	std::vector< std::vector<double> > E;
	std::vector<size_t> nNodesNext; //lowest number of nodes of as-yet unfilled orbitals, per l
	
	double nLeft = nElectrons; //number of electrons yet to be filled
	while(nLeft)
	{	//Search for lowest unoccupied eigenvalue in the l's explored so far
		double Eopt = DBL_MAX; int lOpt = -1;
		for(size_t l=0; l<nNodesNext.size(); l++)
		{	if(E[l].size()==nNodesNext[l]) //don't have unoccupied energy info
				E[l].push_back(getEig(l,nNodesNext[l]));
			if(E[l].back()<Eopt)
			{	Eopt = E[l].back();
				lOpt = l;
			}
		}
		//Look at the next higher l if lOpt is at the highest encountered l so far:
		if(size_t(lOpt+1)==nNodesNext.size())
		{	size_t l = lOpt+1;
			nNodesNext.push_back(0);
			E.resize(l+1);
			E[l].push_back(getEig(l,0));
			if(E[l].back()<Eopt)
			{	Eopt = E[l].back();
				lOpt = l;
			}
		}
		//Occupy the lowest unoccupied orbital:
		double f = std::min(nLeft, 2.*(2*lOpt+1));
		nLeft -= f;
		if(size_t(lOpt+1)>F.size()) F.resize(lOpt+1);
		F[lOpt].push_back(f);
		nNodesNext[lOpt]++;
	}
	
	if(Eptr)
	{	//Drop the unfilled orbital eigenvalues:
		E.resize(F.size());
		for(unsigned l=0; l<F.size(); l++)
			E[l].resize(F[l].size());
		*Eptr = E;
	}
	return F;
}

// Compute total non-interacting energy of atom, given its fillings, and collect any optional outputs
double RadialSchrodinger::compute(const std::vector< std::vector<double> >& F, Outputs outputs)
{	//Prepare outputs:
	double *n=0, *Dn=0, *tau=0;
	if(outputs.n) { outputs.n->assign(rArr.size(), 0.); n = &outputs.n->front(); }
	if(outputs.Dn) { outputs.Dn->assign(rArr.size(), 0.); Dn = &outputs.Dn->front(); }
	if(outputs.tau) { outputs.tau->assign(rArr.size(), 0.); tau = &outputs.tau->front(); }
	if(outputs.z) outputs.z->assign(F.size(), std::vector< std::vector<complex> >());
	if(outputs.E) outputs.E->assign(F.size(), std::vector<double>());
	//Loop over all occupied states:
	std::vector<complex> z(rArr.size()); //space for single eigenfunction
	double Etot = 0.;
	for(size_t l=0; l<F.size(); l++)
		for(size_t nNodes=0; nNodes<F[l].size(); nNodes++)
		{	double f = F[l][nNodes];
			double E = getEig(l,nNodes);
			solveSchEqn(l, E, &z.front());
			//Accumulate orbital-summed quantities:
			Etot += f*E;
			for(size_t i=0; i<rArr.size(); i++)
			{	//Radial wavefunction and its derivative:
				double u = z[i].real(), uPrime = z[i].imag();
				double r = rArr[i], rl = pow(r,l), rlm1 = l ? pow(r,l-1) : 0.; //r^(l-1) terms are non-zero only for l>0
				double R = rl*u, Rbyr = rlm1*u;
				double Rprime = rl*uPrime + l*Rbyr;
				//Compute required quantities:
				if(n) n[i] += f*(R*R);
				if(Dn) Dn[i] += f*(2*R*Rprime);
				if(tau) tau[i] += 0.5*f*(Rprime*Rprime + l*(l+1)*Rbyr*Rbyr);
			}
			//Store individual orbital quantities:
			if(outputs.E) outputs.E->at(l).push_back(E);
			if(outputs.z) outputs.z->at(l).push_back(z);
		}
	return Etot;
}


// Inward and outward integration of the Schrodinger equation
int RadialSchrodinger::solveSchEqn(int l, double E, complex* zSave)
{	//Cached integration weights
	static std::vector< std::vector<double> > wArr;
	if(size_t(l+1)>wArr.size()) wArr.resize(l+1);
	if(!wArr[l].size())
	{	wArr[l].resize(rArr.size());
		for(size_t i=0; i<rArr.size(); i++)
			wArr[l][i] = (4*M_PI) * drArr[i] * pow(rArr[i],2*l+2);
	}
	
	//Rescaling to avoid overflows: 
	const double rangeMin = 1e-50;
	const double rangeMax = 1e+50;
	const double rangeScale = rangeMin/rangeMax;
	std::stack<size_t> iScaleLeft, iScaleRight; //keep track of the locations where rescales happened
	iScaleLeft.push(0); iScaleRight.push(rArr.size());
	
	//Outward integration ('left' of the match point i.e. r<rArr[iMatch])
	double leftIntegral = 0.;
	double logDeriv0 = -Z*1./(l+1); //logarithmic derivative at r=0 (assuming V is finite)
	double r = rArr[0];
	complex z (1.+r*logDeriv0, logDeriv0);
	int prevSign = std::signbit(z.real()), nNodes = 0;
	for(size_t i=0; ; i++)
	{	if(i==iMatch) break;
		leftIntegral += z.real()*z.real() * wArr[l][i];
		if(zSave) zSave[i] = z;
		//RK4 step:
		double rNext = rArr[i+1];
		double h = rNext - r;
		complex k1 = schEqn(r,       2*i,   z,            l, E);
		complex k2 = schEqn(r+0.5*h, 2*i+1, z+(0.5*h)*k1, l, E);
		complex k3 = schEqn(r+0.5*h, 2*i+1, z+(0.5*h)*k2, l, E);
		complex k4 = schEqn(r+h,     2*i+2, z+h*k3,       l, E);
		z += ((1./6)*h) * (k1+k4+2.*(k2+k3));
		r = rNext;
		//Check for nodes:
		int curSign = std::signbit(z.real());
		if(prevSign ^ curSign) nNodes++;
		prevSign = curSign;
		//Rescale:
		if(fabs(z.real()) > rangeMax)
		{	z *= rangeScale;
			leftIntegral *= (rangeScale*rangeScale);
			iScaleLeft.push(i);
		}
	}
	complex zLeft = z;
	
	//Inward integration ('right' of the match point i.e. r>rArr[iMatch])
	double rightIntegral = 0.;
	r = rArr.back();
	z = complex(+0., -rangeMin);
	prevSign = std::signbit(z.real());
	for(size_t i=rArr.size()-1; ; i--)
	{	//Process outputs:
		rightIntegral += z.real()*z.real() * wArr[l][i];
		if(zSave) zSave[i] = z;
		if(i==iMatch) break;
		//RK4 step:
		double rNext = rArr[i-1];
		double h = rNext - r;
		complex k1 = schEqn(r,       2*i,   z,            l, E);
		complex k2 = schEqn(r+0.5*h, 2*i-1, z+(0.5*h)*k1, l, E);
		complex k3 = schEqn(r+0.5*h, 2*i-1, z+(0.5*h)*k2, l, E);
		complex k4 = schEqn(r+h,     2*i-2, z+h*k3,       l, E);
		z += ((1./6)*h) * (k1+k4+2.*(k2+k3));
		r = rNext;
		//Check for nodes:
		int curSign = std::signbit(z.real());
		if(prevSign ^ curSign) nNodes++;
		prevSign = curSign;
		//Rescale:
		if(fabs(z.real()) > rangeMax)
		{	z *= rangeScale;
			rightIntegral *= (rangeScale*rangeScale);
			iScaleRight.push(i);
		}
	}
	complex zRight = z;
	
	if(zSave)
	{	//Pick scale factors that match leftAbs and righAbs, and yield a normalized eigenfunction:
		double normFac = sqrt(1./(zLeft.norm()*rightIntegral + zRight.norm()*leftIntegral));
		double rightScale = normFac * zLeft.abs();
		double leftScale = normFac * zRight.abs();
		//Apply the scale factors to the saved wavefunction:
		size_t i=iMatch;
		while(iScaleRight.size())
		{	for(; i<iScaleRight.top(); i++)
				zSave[i] *= rightScale;
			rightScale *= rangeScale;
			iScaleRight.pop();
		}
		i=iMatch;
		while(iScaleLeft.size())
		{	for(; i>iScaleLeft.top(); i--)
				zSave[i-1] *= leftScale;
			leftScale *= rangeScale;
			iScaleLeft.pop();
		}
	}
	
	//Matching error:
	double err = (zRight.imag()*zLeft.real() - zLeft.imag()*zRight.real()) / sqrt(zRight.norm()*zLeft.norm());
	if(std::signbit(err) == std::signbit(zLeft.real()*zRight.real())) nNodes++;
	
	//Update error vs eigenvalue maps: (enables efficient eigenvalue root-finding)
	if(size_t(l+1)>cachedEerr.size()) cachedEerr.resize(l+1);
	cachedEerr[l][E] = err;
	//Update node count vs eigenvalue maps: (enables efficient eigenvalue bracketing)
	if(size_t(l+1)>nodesEmin.size()) nodesEmin.resize(l+1);
	if(size_t(l+1)>nodesEmax.size()) nodesEmax.resize(l+1);
	auto i = nodesEmin[l].find(nNodes);
	if(i!=nodesEmin[l].end())
	{	//Previously encountered nNodes, update range
		if(E < i->second) i->second = E;
		i=nodesEmax[l].find(nNodes);
		if(E > i->second) i->second = E;
	}
	else
	{	//First time encountering nNodes, set min and max to current E
		nodesEmin[l][nNodes] = E;
		nodesEmax[l][nNodes] = E;
	}
	return nNodes;
}


// Binary-search (interval halving or doubling as appropriate!)
// till given nNodes is encountered at particular l
// and update nodesEmin[l] and nodesEmax[l]
// NOTE: The maps must be seeded with two calls to solveSchEqn(l, E)
// with distinct E's at this l, before calling this routine.
void RadialSchrodinger::locateNodeCount(int l, int nNodes)
{
	auto lbound = nodesEmin[l].lower_bound(nNodes);

	if(lbound==nodesEmin[l].begin())
	{	//the required nNodes is below the explored E range:
		double Elo = nodesEmin[l].begin()->second;
		double Ehi = nodesEmax[l].rbegin()->second;
		do { Elo = 2*Elo - Ehi; } //double the interval by moving the lower end
		while(solveSchEqn(l, Elo) > nNodes);
		lbound = nodesEmin[l].lower_bound(nNodes); //update the iterator (since map changed)
	}
	if(lbound==nodesEmin[l].end())
	{	//the required nNodes is above the explored E range:
		double Elo = nodesEmin[l].begin()->second;
		double Ehi = nodesEmax[l].rbegin()->second;
		do { Ehi = 2*Ehi - Elo; } //double the interval by moving the upper end
		while(solveSchEqn(l, Ehi) < nNodes);
		lbound = nodesEmin[l].lower_bound(nNodes); //update the iterator (since map changed)
	}
	//At this point we've ensured that the explored E range includes nNodes

	if(lbound->first==nNodes) return; //exact match found for nNodes, so job done
	else
	{	//bisect until the required nNodes is encountered:
		auto ubound = nodesEmax[l].lower_bound(nNodes); ubound--;
		double Elo = ubound->second;
		double Ehi = lbound->second;
		while(true)
		{	double Emid = 0.5*(Elo+Ehi);
			int nNodesMid = solveSchEqn(l, Emid);
			if(nNodesMid<nNodes) Elo = Emid;
			else if (nNodesMid>nNodes) Ehi = Emid;
			else return; //exact match found
		}
	}
}

// Compute the error function whose roots are the eigenvalues
double RadialSchrodinger::getEerr(int l, double E)
{	if(size_t(l+1)>cachedEerr.size()) cachedEerr.resize(l+1);
	if(cachedEerr[l].find(E) == cachedEerr[l].end()) //E not yet encountered
		solveSchEqn(l, E); //this solve adds it to the cache
	return cachedEerr[l][E]; //retrieve result from cache
}

// Zero-in on a bracketed eigenvalue using Ridder's method
double RadialSchrodinger::findE(int l, double Elo, double Ehi, double tol)
{
	double errLo = getEerr(l, Elo); if(errLo==0.0) return Elo;
	double errHi = getEerr(l, Ehi); if(errHi==0.0) return Ehi;
	assert(std::signbit(errLo) != std::signbit(errHi));
	while(true)
	{	double Emid = 0.5*(Ehi+Elo);
		double errMid = getEerr(l, Emid); if(errMid==0.0 || fabs(Ehi-Elo)<tol) return Emid;
		double Enew = Emid + (Emid-Elo)*errMid*copysign(1./sqrt(errMid*errMid-errHi*errLo), errLo-errHi);
		double errNew = getEerr(l, Enew); if(errNew==0.0) return Enew;
		if(copysign(errMid,errNew) != errMid)
		{	Elo = Emid; errLo = errMid;
			Ehi = Enew; errHi = errNew;
		}
		else if(copysign(errLo,errNew) != errLo)
		{	Ehi = Enew; errHi = errNew;
		}
		else
		{	Elo = Enew; errLo = errNew;
		}
	}
}

// Find eigenvalue (and optionally eigenfunction) for a specific node count and angular momentum
double RadialSchrodinger::getEig(int l, int nNodes, complex* zEvec)
{	//Seed the node count maps:
	if(size_t(l+1)>=nodesEmin.size() || !nodesEmin[l].size())
	{	double Vmin = *std::min_element(Vsub.begin(),Vsub.end());
		solveSchEqn(l, Vmin-0.5*Z*Z/((l+1)*(l+1))); //lower bound on eigenvalue
		solveSchEqn(l, Vmin+1.); //typical energy range for eigenvalue spacing
		//These don't need to bracket all the eigenvalues, the explored range will be expanded if required
	}
	//First check if already computed:
	if(size_t(l+1)>cachedE.size()) cachedE.resize(l+1);
	double E;
	if(cachedE[l].find(nNodes) != cachedE[l].end())
	{	E = cachedE[l][nNodes];  //already computed
	}
	else
	{	//Locate energies with nNodes and nNodes+1 nodes: that will bound an eigenstate of nNodes nodes
		locateNodeCount(l, nNodes);
		locateNodeCount(l, nNodes+1);
		//Zero in on the bracketed eigenvalue (pick the tightest bracket possible):
		double Emin = nodesEmax[l][nNodes];
		double Emax = nodesEmin[l][nNodes+1];
		E = findE(l, Emin, Emax, 1e-15*hypot(Emin,Emax));  //zero-in till machine precision
		cachedE[l][nNodes] = E; //cache for future queries
	}
	if(zEvec) solveSchEqn(l, E, zEvec); //get the eigenfunction if requested
	return E;
}


//------------------------------ class InvertKS --------------------------------------

InvertKS::InvertKS(const RadialFunctionR& nTarget) : r(nTarget.r), dr(nTarget.dr), n0(nTarget.f)
{	V.assign(r.size(), 0.);
	nElectrons = 0.;
	for(size_t i=0; i<r.size(); i++)
		nElectrons += 4*M_PI*r[i]*r[i]*dr[i] * n0[i];
}

void InvertKS::step(const diagMatrix& dV, double alpha)
{	axpy(alpha, dV, V);
}

double InvertKS::compute(diagMatrix* E_V, diagMatrix* KE_V)
{	RadialSchrodinger atom(r, dr, V, 0);
	F = atom.getFillings(nElectrons);
	double Etot = atom.compute(F, RadialSchrodinger::Outputs(&n));
	if(E_V) E_V->resize(r.size());
	double Vdotn0 = 0.;
	for(size_t i=0; i<r.size(); i++)
	{	double wr = 4*M_PI*r[i]*r[i]*dr[i];
		Vdotn0 += wr * V[i] * n0[i];
		if(E_V) E_V->at(i) = wr*(n0[i] - n[i]);
	}
	if(KE_V) *KE_V = precondition(*E_V);
	return Vdotn0 - Etot;
}

diagMatrix InvertKS::precondition(const diagMatrix& grad)
{	const double kappa = 1./r.back();
	//Solve radial Helmholtz equation with constant kappa: (equivalent to 1/(G^2+kappa^2) in PW)
	diagMatrix Kgrad(r.size());
	//Outward contributions:
	double Qin = 0.;
	for(size_t i=0; i<r.size(); i++)
	{	double expc = r[i] ? exp(-kappa*r[i])/r[i] : 0.;
		double sinch = (kappa && r[i]) ? sinh(kappa*r[i])/(kappa*r[i]) : 1.;
		double dQin = grad[i] * sinch;
		Qin += 0.5*dQin;
		Kgrad[i] = Qin * expc;
		Qin += 0.5*dQin;
	}
	//Inward contributions:
	double Qout = 0.;
	for(size_t i=r.size()-1; ; i--)
	{	double expc = r[i] ? exp(-kappa*r[i])/r[i] : 0.;
		double sinch = (kappa && r[i]) ? sinh(kappa*r[i])/(kappa*r[i]) : 1.;
		double dQout = grad[i] * expc;
		Qout += 0.5*dQout;
		Kgrad[i] += Qout * sinch;
		Qout += 0.5*dQout;
		if(!i) break;
	}
	//Hold end of grid potential at 0
	double KgradOffset = -Kgrad.back();
	for(size_t i=0; i<r.size(); i++)
		Kgrad[i] += KgradOffset;
	return Kgrad;
}

//Get the orbital KE density
RadialFunctionR InvertKS::getTau()
{	std::vector<double> tau;
	RadialSchrodinger atom(r, dr, V, 0);
	atom.compute(F, RadialSchrodinger::Outputs(0, 0, &tau));
	
	RadialFunctionR tauRadial(r.size());
	tauRadial.set(r, dr);
	bool tailCorrect = false;
	for(size_t i=0; i<r.size(); i++)
	{	double iPlus = i+1<r.size() ? i+1 : i;
		double iMinus = i ? i-1 : i;
		double Dn0byn0 = log(n0[iPlus]/n0[iMinus]) / (r[i] * log(r[iPlus]/r[iMinus]));
		double tauVW0 = 0.125 * Dn0byn0 * Dn0byn0 * n0[i];
		if(tauVW0 > tau[i]) tailCorrect=true;
		tauRadial.f[i] = (tailCorrect ? tauVW0 : tau[i]);
		//tauRadial.f[i] = (tauVW0 + 0.3*n0[i]*pow(3*M_PI*M_PI*n0[i], 2./3));
	}
	return tauRadial;
}
