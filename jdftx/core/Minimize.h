/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Kendra Letchworth Weaver

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

#ifndef JDFTX_CORE_MINIMIZE_H
#define JDFTX_CORE_MINIMIZE_H

#include <cmath>
#include <algorithm>
#include <deque>
#include <core/Util.h>
#include <core/MinimizeParams.h>

//! @addtogroup optimization
//! @{

/**
@file Minimize.h
@brief Nonlinear minimization templates
*/

/** Interface (abstract base class) for the minimization algorithm template
	@tparam Vector A data type that represents a direction in the tangent space of the parameter manifold,
	which must have the following functions/operators defined: \n
		- Vector& operator*=(double); (scalar multiply) \n
		- void axpy(double alpha, const Vector& X, Vector& Y); (accumulate operation: Y += alpha*X) \n
		- double dot(const Vector& X, const Vector& Y); (inner product) \n
		- Vector clone(const Vector& X); (create a copy) \n
		- void randomize(Vector&); (initialize with random numbers, used for auto-fdtests)
*/
template<typename Vector> struct Minimizable
{
	//! Move the state in parameter space along direction dir with scale alpha
	virtual void step(const Vector& dir, double alpha)=0;
	
	//! Returns the objective function at the current state and store the gradient in grad, if non-null.
	virtual double compute(Vector* grad)=0; 
	
	//! Override to enable preconditioning: return the preconditioned gradient, given the gradient
	//! The last call to compute() is guaranteed to be at the same position, so a cached result may be returned
	virtual Vector precondition(const Vector& grad) { return clone(grad); }
	
	//! Override for optional processing/reporting after each/every few iterations
	//! It should return whether the state was modified
	virtual bool report(int iter) { return false; }
	
	//! Constrain an arbitrary vector to the space of free directions for minimize.
	//! Used only to generate a random direction for fdTest within the valid minimization subspace.
	virtual void constrain(Vector&) {}
	
	//! Minimize this objective function with algorithm controlled by params and return the minimized value
	double minimize(const MinimizeParams& params);
	
	//! Checks the consistency of the value and gradient returned by compute.
	//! params is used primarily to control output
	void fdTest(const MinimizeParams& params);
};

/** Interface (abstract base class) for linear conjugate gradients template which
	minimizes (1/2) dot(state, hessian*state) - dot(state, rhs) which is equivalent
	to the linear solve hessian * state == rhs
	@tparam Vector Same requirements as the Vector for #Minimizable
*/
template<typename Vector> struct LinearSolvable
{
	Vector state; //!< the location of the minimum, obtained by solving hessian * state == rhs
	
	//! Return vector multiplied by the hessian of the objective function
	virtual Vector hessian(const Vector&)=0;
	
	//! Override to enable preconditioning: return the preconditioned vector, given a vector
	virtual Vector precondition(const Vector& v) { return clone(v); }
	
	//! Solve the linear system hessian * state == rhs using conjugate gradients:
	//! @return the number of iterations taken to achieve target tolerance
	int solve(const Vector& rhs, const MinimizeParams& params);
};

//! @}
//---------------------- Implementation ----------------------------
//!@cond

namespace MinimizePrivate
{
	//Each of these linmin methods advance the parameters in obj along direction d
	//updating the energy E, gradient g, and the step-size alpha
	//The return value specifies if the step succeeded at reducing E
	//If the step fails, alpha MUST contain the total progress along dir
	//made by this step, so that minimize may reset it back to the original value


	//Equation-of-motion / Relaxation method stepping
	//NOTE: Criterion for success of this method is different from the others
	// It only ensures that the energy is not NaN/Inf.
	template<typename Vector>
	bool linminRelax(Minimizable<Vector>& obj, const MinimizeParams& p,
		const Vector& d, double alphaT, double& alpha, double& E, Vector& g)
	{
		alpha = alphaT; //constant step-size equal to the starting value
		obj.step(d, alpha);
		E = obj.compute(&g);
		if(!std::isfinite(E))
		{	fprintf(p.fpLog, "%s\tRelax step failed with %s = %le\n.", p.linePrefix, p.energyLabel, E);
			return false;
		}
		else return true;
	}


	//Quadratic line minimization, a more robust version of the one from the old electronic code
	template<typename Vector>
	bool linminQuad(Minimizable<Vector>& obj, const MinimizeParams& p,
		const Vector& d, double alphaT, double& alpha, double& E, Vector& g)
	{
		double alphaPrev = 0.0; //the progress made so far along d
		double Eorig = E;
		double gdotd = dot(g,d); //directional derivative at starting point
		if(gdotd >= 0.0)
		{	fprintf(p.fpLog, "%s\tBad step direction: g.d > 0.\n", p.linePrefix);
			alpha = alphaPrev;
			return false;
		}

		//Test step and step size prediction:
		double ET = 0.0; //test step energy
		for(int s=0; s<p.nAlphaAdjustMax; s++)
		{	if(alphaT < p.alphaTmin)
			{	fprintf(p.fpLog, "%s\talphaT below threshold %le. Quitting step.\n", p.linePrefix, p.alphaTmin);
				alpha = alphaPrev;
				return false;
			}
			//Try the test step:
			obj.step(d, alphaT-alphaPrev); alphaPrev = alphaT;
			ET = obj.compute(0);
			//Check if step crossed domain of validity of parameter space:
			if(!std::isfinite(ET))
			{	alphaT *= p.alphaTreduceFactor;
				fprintf(p.fpLog, "%s\tTest step failed with %s = %le, reducing alphaT to %le.\n",
					p.linePrefix, p.energyLabel, ET, alphaT);
				continue;
			}
			//Predict step size:
			alpha = 0.5*pow(alphaT,2)*gdotd/(alphaT*gdotd + E - ET);
			//Check reasonableness of predicted step size:
			if(alpha<0)
			{	//Curvature has the wrong sign
				//That implies ET < E, so accept test step as the real step
				alpha = alphaT;
				fprintf(p.fpLog, "%s\tWrong curvature in test step, using alpha=alphaT.\n", p.linePrefix);
				E = obj.compute(&g);
				return true;
			}
			if(alpha/alphaT > p.alphaTincreaseFactor)
			{	alphaT *= p.alphaTincreaseFactor;
				fprintf(p.fpLog, "%s\tPredicted alpha/alphaT>%lf, increasing alphaT to %le.\n",
					p.linePrefix, p.alphaTincreaseFactor, alphaT);
				continue;
			}
			if(alphaT/alpha < p.alphaTreduceFactor)
			{	alphaT *= p.alphaTreduceFactor;
				fprintf(p.fpLog, "%s\tPredicted alpha/alphaT<%lf, reducing alphaT to %le.\n",
					p.linePrefix, p.alphaTreduceFactor, alphaT);
				continue;
			}
			//Successful test step:
			break;
		}
		if(!std::isfinite(E))
		{	fprintf(p.fpLog, "%s\tTest step failed %d times. Quitting step.\n", p.linePrefix, p.nAlphaAdjustMax);
			alpha = alphaPrev;
			return false;
		}

		//Actual step:
		for(int s=0; s<p.nAlphaAdjustMax; s++)
		{	//Try the step:
			obj.step(d, alpha-alphaPrev); alphaPrev=alpha;
			E = obj.compute(&g);
			if(!std::isfinite(E))
			{	alpha *= p.alphaTreduceFactor;
				fprintf(p.fpLog, "%s\tStep failed with %s = %le, reducing alpha to %le.\n",
					p.linePrefix, p.energyLabel, E, alpha);
				continue;
			}
			if(E > Eorig)
			{	alpha *= p.alphaTreduceFactor;
				fprintf(p.fpLog, "%s\tStep increased %s by %le, reducing alpha to %le.\n",
					p.linePrefix, p.energyLabel, E-Eorig, alpha);
				continue;
			}
			//Step successful:
			break;
		}
		if(!std::isfinite(E) || E>Eorig)
		{	fprintf(p.fpLog, "%s\tStep failed to reduce %s after %d attempts. Quitting step.\n",
				p.linePrefix, p.energyLabel, p.nAlphaAdjustMax);
			return false;
		}
		return true;
	}


	//Cubic line minimization, designed to handle fluids which can be highly non-quadratic
	template<typename Vector>
	bool linminCubic(Minimizable<Vector>& obj, const MinimizeParams& p,
		const Vector& d, double alphaT, double& alpha, double& E, Vector& g)
	{
		double alphaPrev = 0;
		double Eprev = E;
		double gdotdPrev = dot(g,d); //directional derivative at starting point
		if(gdotdPrev >= 0.0)
		{	fprintf(p.fpLog, "%s\tBad step direction: g.d > 0.\n", p.linePrefix);
			alpha = 0;
			return false;
		}

		//Don't need to distinguish test and actual steps for cubic linmin!
		alpha = alphaT; //this is the initial tentative step
		double alphaState = 0.0; //location of obj's state along search direction since the start of this linmin
		bool alphaOnRoot = false; //whether current alpha has been set to a root from a previous step
		for(int s=0; s<p.nAlphaAdjustMax; s++)
		{	//Move by alpha:
			obj.step(d, alpha-alphaState); alphaState=alpha;
			E = obj.compute(&g);
			//Check for domain error:
			if(!std::isfinite(E))
			{	alpha *= p.alphaTreduceFactor;
				fprintf(p.fpLog, "%s\tStep failed with %s = %le, reducing alpha to %le.\n",
					p.linePrefix, p.energyLabel, E, alpha);
				continue;
			}
			//Check if we're done:
			if(alphaOnRoot)
			{	if(E<=Eprev) return true;
				else fprintf(p.fpLog, "%s\tStep increased %s by %le, trying another cubic\n",
					p.linePrefix, p.energyLabel, E-Eprev);
			}
			double gdotd = dot(g,d);
			//Have a valid cubic spline interval [alphaPrev,alpha]
			//with values Eprev, E and derivatives gdotdPrev, gdotd
			//Change coordinates to make that the unit interval in t:
			double Ep0 = gdotdPrev*(alpha-alphaPrev);
			double Ep1 = gdotd*(alpha-alphaPrev);
			double deltaE = E-Eprev;
			//dE/dt is a quadratic At^2-2Bt+C, with coefficients:
			double A = 3*(Ep0 + Ep1 - 2*deltaE);
			double B = 2*Ep0 + Ep1 - 3*deltaE;
			double C = Ep0;
			assert(Ep0<=0);
			//No maxima/minima or 0-slope inflection point:
			if(B*B-A*C<=0)
			{	assert(E<=Eprev);
				assert(gdotdPrev<=0);
				//This implies energy is monotonically decreasing (we know Ep0<0)
				double alphaNew = alphaPrev + p.alphaTincreaseFactor*(alpha-alphaPrev);
				alphaPrev = alpha; Eprev=E; gdotdPrev=gdotd;
				alpha = alphaNew; alphaOnRoot = false;
				//Since E<Eprev, scoot the interval over to start at alpha (future failures won't affect this progress)
				alphaPrev -= alphaState;
				alpha -= alphaState;
				alphaState = 0;
				continue;
			}
			//Locate minimum of E(t):
			double disc = sqrt(B*B-A*C);
			double tOpt = B>0 ? (B+disc)/A : C/(B-disc);
			double Eopt = Eprev + tOpt*(C + tOpt*(-B + tOpt*A/3));
			if(tOpt>=0 && tOpt<p.alphaTincreaseFactor && Eopt<=E && Eopt<=Eprev)
			{	//This is a valid energy reducing, not too far off solution
				double alphaNew = alphaPrev + tOpt*(alpha-alphaPrev);
				if(tOpt<=1)
				{	alpha = alphaNew;
				}
				else
				{	assert(E<=Eprev);
					assert(gdotdPrev<=0);
					//E<Eprev, so never have to revisit left of alpha
					alphaPrev = alpha; Eprev=E; gdotdPrev=gdotd;
					alpha = alphaNew;
					//Commit the progress thus far, protecting it from future failures in step:
					alphaPrev -= alphaState;
					alpha -= alphaState;
					alphaState = 0;
				}
				alphaOnRoot = true;
				continue;
			}
			else
			{	//If E>Eprev, then must have had tOpt in (0,1) with Eopt<Eprev<E.
				//Therefore definitely have E<Eprev here:
				assert(E<=Eprev);
				assert(gdotdPrev<=0);
				//Scoot interval over to start at E:
				double alphaNew = alphaPrev + p.alphaTincreaseFactor*(alpha-alphaPrev);
				alphaPrev = alpha; Eprev=E; gdotdPrev=gdotd;
				alpha = alphaNew; alphaOnRoot = false;
				//Commit the progress thus far, protecting it from future failures in step:
				alphaPrev -= alphaState;
				alpha -= alphaState;
				alphaState = 0;
				continue;
			}
		}
		//Check if current state is invalid or worse than prev:
		if(!std::isfinite(E) || E>Eprev)
		{	alpha = alphaState; //minimize will roll back to the last known good state
			return false;
		}
		else return true;
	}
}

//Energy difference convergence check
class EdiffCheck : std::deque<double>
{	unsigned nDiff;
	double threshold;
public:
	EdiffCheck(unsigned nDiff, double threshold) : nDiff(nDiff), threshold(fabs(threshold))
	{
	}
	bool checkConvergence(double E)
	{	if(!size()) { push_back(E); return false; } //first element
		if(E >= back()) { clear(); push_back(E); return false; } //energy increased, reset converge list
		//Have atleast one energy difference in list:
		push_back(E);
		if(size()==nDiff+2) pop_front(); //discard old unneeded elements
		if(size()==nDiff+1)
		{	for(unsigned i=0; i<nDiff; i++)
				if(at(i+1) < at(i)-threshold)
					return false;
			return true;
		}
		else return false;
	}
};

//Search direction reset check
class DirResetCheck : std::deque<bool>
{	unsigned num, den;
public:
	DirResetCheck(unsigned num, unsigned den) : num(num), den(den)
	{
	}
	bool checkReset(bool sdReset) //returns true if recent reset count exceeds threshold
	{	push_back(sdReset);
		if(size()==den+1) pop_front(); //discard old unneeded elements
		if(size()==den) return std::count(begin(), end(), true) > num;
		return false;
	}
};

template<typename Vector> double Minimizable<Vector>::minimize(const MinimizeParams& p)
{	using namespace MinimizePrivate;

	if(p.fdTest) fdTest(p); // finite difference test

	Vector g, gPrev; //current and previous gradient
	double E = compute(&g); //get initial energy and gradient
	EdiffCheck ediffCheck(p.nEnergyDiff, p.energyDiffThreshold); //list of past energies
	DirResetCheck dirResetCheck(p.nDirResetNum, p.nDirResetDen); //reste history
	
	Vector d = clone(g); //step direction (will be reset in first iteration)
	bool forceGradDirection = true; //whether current direction is along the gradient
	MinimizeParams::DirectionUpdateScheme currentDirUpdateScheme = p.dirUpdateScheme; //initially use the specified scheme, may switch to SD on trouble
	
	double alphaT = p.alphaTstart; //test step size
	double alpha = alphaT; //actual step size
	double beta = 0.0; //CG prev search direction mix factor
	double gKNorm = 0.0, gKNormPrev = 0.0; //current and previous norms of the preconditioned gradient

	//Select the linmin method:
	bool (*linmin)(Minimizable<Vector>&, const MinimizeParams&, const Vector&, double, double&, double&, Vector&)=0;
	switch(p.linminMethod)
	{	case MinimizeParams::Relax: linmin = linminRelax<Vector>; break;
		case MinimizeParams::Quad:  linmin = linminQuad<Vector>; break;
		case MinimizeParams::Cubic: linmin = linminCubic<Vector>; break;
	}

	//Iterate until convergence, max iteration count or kill signal
	int iter=0;
	for(iter=0; !killFlag; iter++)
	{	
		if(report(iter)) //optional reporting/processing
		{	E = compute(&g); //update energy and gradient if state was modified
			fprintf(p.fpLog, "%s\tState modified externally: resetting search direction.\n", p.linePrefix);
			fflush(p.fpLog);
			forceGradDirection = true; //reset search direction
		}
		
		Vector Kg = precondition(g);
		gKNorm = dot(g,Kg);
		fprintf(p.fpLog, "%sIter: %3d  %s: %22.15le  |grad|_K: %10.3le  alpha: %10.3le",
			p.linePrefix, iter, p.energyLabel, E, sqrt(gKNorm/p.nDim), alpha);

		//Print prev step stats and set CG direction parameter if necessary
		beta = 0.0;
		if(!forceGradDirection)
		{	double dotgd = dot(g,d);
			double dotgPrevKg = dot(gPrev, Kg);

			double linmin = dotgd/sqrt(dot(g,g)*dot(d,d));
			double cgtest = dotgPrevKg/sqrt(gKNorm*gKNormPrev);
			fprintf(p.fpLog, ", linmin = %10.3le", linmin);
			fprintf(p.fpLog, ", cgtest = %10.3le", cgtest);

			//Update beta:
			switch(currentDirUpdateScheme)
			{	case MinimizeParams::FletcherReeves:  beta = gKNorm/gKNormPrev; break;
				case MinimizeParams::PolakRibiere:    beta = (gKNorm-dotgPrevKg)/gKNormPrev; break;
				case MinimizeParams::HestenesStiefel: beta = (gKNorm-dotgPrevKg)/(dotgd-dot(d,gPrev)); break;
				case MinimizeParams::SteepestDescent: beta = 0.0; break;
			}
			if(beta<0.0)
			{	fprintf(p.fpLog, "\n%sEncountered beta<0, resetting CG.", p.linePrefix);
				beta=0.0;
			}
		}
		forceGradDirection = false;
		fprintf(p.fpLog, "\n"); fflush(p.fpLog);
		if(sqrt(gKNorm/p.nDim) < p.knormThreshold)
		{	fprintf(p.fpLog, "%sConverged (|grad|_K<%le).\n", p.linePrefix, p.knormThreshold);
			fflush(p.fpLog); return E;
		}
		if(ediffCheck.checkConvergence(E))
		{	fprintf(p.fpLog, "%sConverged (|Delta %s|<%le for %d iters).\n",
				p.linePrefix, p.energyLabel, p.energyDiffThreshold, p.nEnergyDiff);
			fflush(p.fpLog); return E;
		}
		if(!std::isfinite(gKNorm))
		{	fprintf(p.fpLog, "%s|grad|_K=%le. Stopping ...\n", p.linePrefix, gKNorm);
			fflush(p.fpLog); return E;
		}
		if(!std::isfinite(E))
		{	fprintf(p.fpLog, "%sE=%le. Stopping ...\n", p.linePrefix, E);
			fflush(p.fpLog); return E;
		}
		if(iter>=p.nIterations) break;
		gPrev = g;
		gKNormPrev = gKNorm;

		//Update search direction
		d *= beta; axpy(-1.0, Kg, d);  // d = beta*d - Kg

		//Line minimization
		if(linmin(*this, p, d, alphaT, alpha, E, g))
		{	//linmin succeeded:
			if(p.updateTestStepSize)
			{	alphaT = alpha;
				if(alphaT<p.alphaTmin) //bad step size
					alphaT = p.alphaTstart; //make sure next test step size is not too bad
			}
			fflush(p.fpLog);
			dirResetCheck.checkReset(false);
		}
		else
		{	//linmin failed:
			fprintf(p.fpLog, "%s\tUndoing step.\n", p.linePrefix);
			step(d, -alpha);
			E = compute(&g);
			if(beta)
			{	//Failed, but not along the gradient direction:
				fprintf(p.fpLog, "%s\tStep failed: resetting search direction.\n", p.linePrefix);
				fflush(p.fpLog);
				forceGradDirection = true; //reset search direction
				if(dirResetCheck.checkReset(true))
				{	fprintf(p.fpLog,
						"%s\tSearch direction was reset %d of the last %d iterations, switching to Steepest Descents\n",
						p.linePrefix, p.nDirResetNum, p.nDirResetDen);
					currentDirUpdateScheme = MinimizeParams::SteepestDescent;
				}
			}
			else
			{	//Failed along the gradient direction
				fprintf(p.fpLog, "%s\tStep failed along negative gradient direction.\n", p.linePrefix);
				fprintf(p.fpLog, "%sProbably at roundoff error limit. (Stopping)\n", p.linePrefix);
				fflush(p.fpLog);
				return E;
			}
		}
	}
	fprintf(p.fpLog, "%sNone of the convergence criteria satisfied after %d iterations.\n", p.linePrefix, iter);
	return E;
}


template<typename Vector> void Minimizable<Vector>::fdTest(const MinimizeParams& p)
{	
	const double deltaMin = 1e-9;
	const double deltaMax = 1e+1;
	const double deltaScale = 1e+1;
	string fdPrefixString = p.linePrefix + string("fdTest: ");
	const char* fdPrefix = fdPrefixString.c_str();
	fprintf(p.fpLog, "%s--------------------------------------\n", fdPrefix);
	Vector g;
	double E0 = compute(&g);
	
	Vector dx;
	{	// Set the direction to be a random vector of the same norm
		// as the preconditioned gradient times the initial test step size
		Vector Kg = precondition(g);
		dx = clone(Kg);
		randomize(dx);
		constrain(dx);
		dx *= p.alphaTstart * sqrt(dot(Kg,Kg)/dot(dx,dx));
	}
	double dE_ddelta = dot(dx, g); //directional derivative at delta=0

	double deltaPrev=0;
	for(double delta=deltaMin; delta<=deltaMax; delta*=deltaScale)
	{	double dE = dE_ddelta*delta;
		step(dx, delta-deltaPrev); deltaPrev=delta;
		double deltaE = compute(0) - E0;
		fprintf(p.fpLog, "%s   delta=%le:\n", fdPrefix, delta);
		fprintf(p.fpLog, "%s      d%s Ratio: %19.16lf\n", fdPrefix, p.energyLabel, deltaE/dE);
		fprintf(p.fpLog, "%s      d%s Error: %19.16lf\n", fdPrefix, p.energyLabel, sqrt(p.nDim)*1.1e-16/fabs(dE));
	}
	fprintf(p.fpLog, "%s--------------------------------------\n", fdPrefix);
	step(dx, -deltaPrev); //restore state to original value
}


template<typename Vector> int LinearSolvable<Vector>::solve(const Vector& rhs, const MinimizeParams& p)
{	//Initialize:
	Vector r = clone(rhs); axpy(-1.0, hessian(state), r); //residual r = rhs - A.state;
	Vector z = precondition(r), d = r; //the preconditioned residual and search direction
	double beta=0.0, rdotzPrev=0.0, rdotz = dot(r, z);

	//Check initial residual
	double rzNorm = sqrt(fabs(rdotz)/p.nDim);
	fprintf(p.fpLog, "%sInitial:  sqrt(|r.z|): %12.6le\n", p.linePrefix, rzNorm);
	if(rzNorm<p.knormThreshold) { fprintf(p.fpLog, "%sConverged sqrt(r.z)<%le\n", p.linePrefix, p.knormThreshold); return 0; }

	//Main loop:
	int iter;
	for(iter=0; iter<p.nIterations && !killFlag; iter++)
	{	//Update search direction:
		if(rdotzPrev)
		{	beta = rdotz/rdotzPrev;
			d *= beta; axpy(1.0, z, d); // d = z + beta*d
		}
		else d = clone(z); //fresh search direction (along gradient)
		//Step:
		Vector w = hessian(d);
		double alpha = rdotz/dot(w,d);
		axpy(alpha, d, state);
		axpy(-alpha, w, r);
		z = precondition(r);
		rdotzPrev = rdotz;
		rdotz = dot(r, z);
		//Print info:
		double rzNorm = sqrt(fabs(rdotz)/p.nDim);
		fprintf(p.fpLog, "%sIter: %3d  sqrt(|r.z|): %12.6le  alpha: %12.6le  beta: %13.6le\n",
			p.linePrefix, iter, rzNorm, alpha, beta);
		//Check convergence:
		if(rzNorm<p.knormThreshold) { fprintf(p.fpLog, "%sConverged sqrt(r.z)<%le\n", p.linePrefix, p.knormThreshold); return iter; }
	}
	fprintf(p.fpLog, "%sGradient did not converge within threshold in %d iterations\n", p.linePrefix, iter);
	return iter;
}

//!@endcond
#endif // JDFTX_CORE_MINIMIZE_H
