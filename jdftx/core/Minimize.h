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
	
	//! Override to synchronize scalars over MPI processes (if the same minimization is happening in sync over many processes)
	virtual double sync(double x) const { return x; }
	
	//! Minimize this objective function with algorithm controlled by params and return the minimized value
	double minimize(const MinimizeParams& params);
	
	//! Checks the consistency of the value and gradient returned by compute.
	//! params is used primarily to control output
	void fdTest(const MinimizeParams& params);
	
private:
	typedef bool (*Linmin)(Minimizable<Vector>&, const MinimizeParams&, const Vector&, double, double&, double&, Vector&);
	Linmin getLinmin(const MinimizeParams& params) const;
	double lBFGS(const MinimizeParams& params); //limited memory BFGS implementation (differs sufficiently from CG to be justify a separate implementation)
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
	virtual Vector hessian(const Vector&) const=0;
	
	//! Override to enable preconditioning: return the preconditioned vector, given a vector
	virtual Vector precondition(const Vector& v) const { return clone(v); }
	
	//! Override to synchronize scalars over MPI processes (if the same minimization is happening in sync over many processes)
	virtual double sync(double x) const { return x; }
	
	//! Solve the linear system hessian * state == rhs using conjugate gradients:
	//! @return the number of iterations taken to achieve target tolerance
	int solve(const Vector& rhs, const MinimizeParams& params);
};

//! @}
//---------------------- Implementation ----------------------------
//!@cond

#include <core/Minimize_linmin.h>
#include <core/Minimize_lBFGS.h>

template<typename Vector> double Minimizable<Vector>::minimize(const MinimizeParams& p)
{	if(p.fdTest) fdTest(p); // finite difference test
	if(p.dirUpdateScheme == MinimizeParams::LBFGS) return lBFGS(p);
	
	Vector g, gPrev; //current and previous gradient
	double E = sync(compute(&g)); //get initial energy and gradient
	EdiffCheck ediffCheck(p.nEnergyDiff, p.energyDiffThreshold); //list of past energies
	
	Vector d = clone(g); //step direction (will be reset in first iteration)
	bool forceGradDirection = true; //whether current direction is along the gradient
	MinimizeParams::DirectionUpdateScheme currentDirUpdateScheme = p.dirUpdateScheme; //initially use the specified scheme, may switch to SD on trouble
	
	double alphaT = p.alphaTstart; //test step size
	double alpha = alphaT; //actual step size
	double beta = 0.0; //CG prev search direction mix factor
	double gKNorm = 0.0, gKNormPrev = 0.0; //current and previous norms of the preconditioned gradient

	//Select the linmin method:
	Linmin linmin = getLinmin(p);
	
	//Iterate until convergence, max iteration count or kill signal
	int iter=0;
	for(iter=0; !killFlag; iter++)
	{
		if(report(iter)) //optional reporting/processing
		{	E = sync(compute(&g)); //update energy and gradient if state was modified
			fprintf(p.fpLog, "%s\tState modified externally: resetting search direction.\n", p.linePrefix);
			fflush(p.fpLog);
			forceGradDirection = true; //reset search direction
		}
		
		Vector Kg = precondition(g);
		gKNorm = sync(dot(g,Kg));
		fprintf(p.fpLog, "%sIter: %3d  %s: ", p.linePrefix, iter, p.energyLabel);
		fprintf(p.fpLog, p.energyFormat, E);
		fprintf(p.fpLog, "  |grad|_K: %10.3le  alpha: %10.3le", sqrt(gKNorm/p.nDim), alpha);

		//Print prev step stats and set CG direction parameter if necessary
		beta = 0.0;
		if(!forceGradDirection)
		{	double dotgd = sync(dot(g,d));
			double dotgPrevKg = sync(dot(gPrev, Kg));

			double linmin = dotgd/sqrt(sync(dot(g,g))*sync(dot(d,d)));
			double cgtest = dotgPrevKg/sqrt(gKNorm*gKNormPrev);
			fprintf(p.fpLog, "  linmin: %10.3le", linmin);
			fprintf(p.fpLog, "  cgtest: %10.3le", cgtest);

			//Update beta:
			switch(currentDirUpdateScheme)
			{	case MinimizeParams::FletcherReeves:  beta = gKNorm/gKNormPrev; break;
				case MinimizeParams::PolakRibiere:    beta = (gKNorm-dotgPrevKg)/gKNormPrev; break;
				case MinimizeParams::HestenesStiefel: beta = (gKNorm-dotgPrevKg)/(dotgd-sync(dot(d,gPrev))); break;
				case MinimizeParams::SteepestDescent: beta = 0.0; break;
				case MinimizeParams::LBFGS: break; //Should never encounter since LBFGS handled separately; just to eliminate compiler warnings
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
		}
		else
		{	//linmin failed:
			fprintf(p.fpLog, "%s\tUndoing step.\n", p.linePrefix);
			step(d, -alpha);
			E = sync(compute(&g));
			if(beta)
			{	//Failed, but not along the gradient direction:
				fprintf(p.fpLog, "%s\tStep failed: resetting search direction.\n", p.linePrefix);
				fflush(p.fpLog);
				forceGradDirection = true; //reset search direction
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
	double E0 = sync(compute(&g));
	
	Vector dx;
	{	// Set the direction to be a random vector of the same norm
		// as the preconditioned gradient times the initial test step size
		Vector Kg = precondition(g);
		dx = clone(Kg);
		randomize(dx);
		constrain(dx);
		dx *= p.alphaTstart * sqrt(sync(dot(Kg,Kg))/sync(dot(dx,dx)));
	}
	double dE_ddelta = sync(dot(dx, g)); //directional derivative at delta=0

	double deltaPrev=0;
	for(double delta=deltaMin; delta<=deltaMax; delta*=deltaScale)
	{	double dE = dE_ddelta*delta;
		step(dx, delta-deltaPrev); deltaPrev=delta;
		double deltaE = sync(compute(0)) - E0;
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
	double beta=0.0, rdotzPrev=0.0, rdotz = sync(dot(r, z));

	//Check initial residual
	double rzNorm = sqrt(fabs(rdotz)/p.nDim);
	fprintf(p.fpLog, "%sInitial:  sqrt(|r.z|): %12.6le\n", p.linePrefix, rzNorm); fflush(p.fpLog);
	if(rzNorm<p.knormThreshold) { fprintf(p.fpLog, "%sConverged sqrt(r.z)<%le\n", p.linePrefix, p.knormThreshold); fflush(p.fpLog); return 0; }

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
		double alpha = rdotz/sync(dot(w,d));
		axpy(alpha, d, state);
		axpy(-alpha, w, r);
		z = precondition(r);
		rdotzPrev = rdotz;
		rdotz = sync(dot(r, z));
		//Print info:
		double rzNorm = sqrt(fabs(rdotz)/p.nDim);
		fprintf(p.fpLog, "%sIter: %3d  sqrt(|r.z|): %12.6le  alpha: %12.6le  beta: %13.6le\n",
			p.linePrefix, iter, rzNorm, alpha, beta); fflush(p.fpLog);
		//Check convergence:
		if(rzNorm<p.knormThreshold) { fprintf(p.fpLog, "%sConverged sqrt(r.z)<%le\n", p.linePrefix, p.knormThreshold); fflush(p.fpLog); return iter; }
	}
	fprintf(p.fpLog, "%sGradient did not converge within threshold in %d iterations\n", p.linePrefix, iter); fflush(p.fpLog);
	return iter;
}

//!@endcond
#endif // JDFTX_CORE_MINIMIZE_H
