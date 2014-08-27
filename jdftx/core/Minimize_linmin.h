/*-------------------------------------------------------------------
Copyright 2014 Ravishankar Sundararaman

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

#ifndef JDFTX_CORE_MINIMIZE_LINMIN_H
#define JDFTX_CORE_MINIMIZE_LINMIN_H

#include <core/Minimize.h> //this include is present only to aid IDE's autocompletion (does nothing since this file is always used via Minimize.h)
#include <deque>

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
		E = obj.sync(obj.compute(&g));
		if(!std::isfinite(E))
		{	fprintf(p.fpLog, "%s\tRelax step failed with %s = %le\n.", p.linePrefix, p.energyLabel, E); fflush(p.fpLog);
			return false;
		}
		else return true;
	}


	//Quadratic line minimization
	template<typename Vector>
	bool linminQuad(Minimizable<Vector>& obj, const MinimizeParams& p,
		const Vector& d, double alphaT, double& alpha, double& E, Vector& g)
	{
		double alphaPrev = 0.0; //the progress made so far along d
		double Eorig = E;
		double gdotd = obj.sync(dot(g,d)); //directional derivative at starting point
		if(gdotd >= 0.0)
		{	fprintf(p.fpLog, "%s\tBad step direction: g.d > 0.\n", p.linePrefix); fflush(p.fpLog);
			alpha = alphaPrev;
			return false;
		}

		//Test step and step size prediction:
		double ET = 0.0; //test step energy
		for(int s=0; s<p.nAlphaAdjustMax; s++)
		{	if(alphaT < p.alphaTmin)
			{	fprintf(p.fpLog, "%s\talphaT below threshold %le. Quitting step.\n", p.linePrefix, p.alphaTmin); fflush(p.fpLog);
				alpha = alphaPrev;
				return false;
			}
			//Try the test step:
			obj.step(d, alphaT-alphaPrev); alphaPrev = alphaT;
			ET = obj.sync(obj.compute(0));
			//Check if step crossed domain of validity of parameter space:
			if(!std::isfinite(ET))
			{	alphaT *= p.alphaTreduceFactor;
				fprintf(p.fpLog, "%s\tTest step failed with %s = %le, reducing alphaT to %le.\n",
					p.linePrefix, p.energyLabel, ET, alphaT); fflush(p.fpLog);
				continue;
			}
			//Predict step size:
			alpha = 0.5*pow(alphaT,2)*gdotd/(alphaT*gdotd + E - ET);
			//Check reasonableness of predicted step size:
			if(alpha<0)
			{	//Curvature has the wrong sign
				//That implies ET < E, so accept step for now, and try descending further next time
				alphaT *= p.alphaTincreaseFactor;
				fprintf(p.fpLog, "%s\tWrong curvature in test step, increasing alphaT to %le.\n", p.linePrefix, alphaT); fflush(p.fpLog);
				E = obj.sync(obj.compute(&g));
				return true;
			}
			if(alpha/alphaT > p.alphaTincreaseFactor)
			{	alphaT *= p.alphaTincreaseFactor;
				fprintf(p.fpLog, "%s\tPredicted alpha/alphaT>%lf, increasing alphaT to %le.\n",
					p.linePrefix, p.alphaTincreaseFactor, alphaT); fflush(p.fpLog);
				continue;
			}
			if(alphaT/alpha < p.alphaTreduceFactor)
			{	alphaT *= p.alphaTreduceFactor;
				fprintf(p.fpLog, "%s\tPredicted alpha/alphaT<%lf, reducing alphaT to %le.\n",
					p.linePrefix, p.alphaTreduceFactor, alphaT); fflush(p.fpLog);
				continue;
			}
			//Successful test step:
			break;
		}
		if(!std::isfinite(E))
		{	fprintf(p.fpLog, "%s\tTest step failed %d times. Quitting step.\n", p.linePrefix, p.nAlphaAdjustMax); fflush(p.fpLog);
			alpha = alphaPrev;
			return false;
		}

		//Actual step:
		for(int s=0; s<p.nAlphaAdjustMax; s++)
		{	//Try the step:
			obj.step(d, alpha-alphaPrev); alphaPrev=alpha;
			E = obj.sync(obj.compute(&g));
			if(!std::isfinite(E))
			{	alpha *= p.alphaTreduceFactor;
				fprintf(p.fpLog, "%s\tStep failed with %s = %le, reducing alpha to %le.\n",
					p.linePrefix, p.energyLabel, E, alpha); fflush(p.fpLog);
				continue;
			}
			if(E > Eorig)
			{	alpha *= p.alphaTreduceFactor;
				fprintf(p.fpLog, "%s\tStep increased %s by %le, reducing alpha to %le.\n",
					p.linePrefix, p.energyLabel, E-Eorig, alpha); fflush(p.fpLog);
				continue;
			}
			//Step successful:
			break;
		}
		if(!std::isfinite(E) || E>Eorig)
		{	fprintf(p.fpLog, "%s\tStep failed to reduce %s after %d attempts. Quitting step.\n",
				p.linePrefix, p.energyLabel, p.nAlphaAdjustMax); fflush(p.fpLog);
			return false;
		}
		return true;
	}


	//Cubic line minimization, designed to handle fluids which can be highly non-quadratic
	template<typename Vector>
	bool linminCubicWolfe(Minimizable<Vector>& obj, const MinimizeParams& p,
		const Vector& d, double alphaT, double& alpha, double& E, Vector& g)
	{
		double Eprev = E;
		double gdotdPrev = obj.sync(dot(g,d)); //directional derivative at starting point
		if(gdotdPrev >= 0.0)
		{	fprintf(p.fpLog, "%s\tBad step direction: g.d > 0.\n", p.linePrefix); fflush(p.fpLog);
			alpha = 0;
			return false;
		}
		double E0 = Eprev, gdotd0 =gdotdPrev; //Always use initial energy and gradient for Wolfe test (even if part of line search has been committed)
		
		alpha = alphaT; //this is the initial tentative step
		double alphaPrev = 0; //alpha of the other point in the cubic interval
		double alphaState = 0.; //alpha that correspodns to state of the objective function
		for(int s=0;; s++)
		{	//Move by alpha:
			obj.step(d, alpha-alphaState); alphaState=alpha;
			E = obj.sync(obj.compute(&g));
			double gdotd = obj.sync(dot(g,d));
			if(s > p.nAlphaAdjustMax) break; //step limit
			//Check for domain error:
			if(!std::isfinite(E))
			{	alpha = alphaPrev + p.alphaTreduceFactor*(alpha-alphaPrev);
				fprintf(p.fpLog, "%s\tStep failed with %s = %le, reducing alpha to %le.\n",
					p.linePrefix, p.energyLabel, E, alpha); fflush(p.fpLog);
				continue;
			}
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
			double tMin = NAN; //location of minimum (NAN if none)
			if(B*B-A*C >= 0) //E'(t) has at least one root
			{	double disc = sqrt(B*B-A*C);
				double tOpt = B>0 ? (B+disc)/A : C/(B-disc); //only root of E'(t) for which E''(t) = 2*(A*t-B) >= 0
				double Eopt = Eprev + tOpt*(C + tOpt*(-B + tOpt*A/3)); //interpolated E(tOpt)
				if(std::isfinite(tOpt) && Eopt<=E && Eopt<=Eprev) //well-defined local minimum lower than current endpoints
					tMin = tOpt;
			}
			if(std::isfinite(tMin))
			{	tMin = std::min(tMin, p.alphaTincreaseFactor); //cap on forward-in-t step
				tMin = std::max(tMin, 1.-p.alphaTincreaseFactor); //cap on backward-in-t step
			}
			else
			{	//no local minimum lower than current endpoints within interval
				//therefore at least one of the endpoints must have function value decreasing in the outwards direction
				if(Ep1<=0) tMin = p.alphaTincreaseFactor; //E(t) decreases above t=1, take max forward-in-t step
				else tMin = 1.-p.alphaTincreaseFactor; //E(t) decreases below t=0, take max backward-in-t step
			}
			double alphaNew = alphaPrev + tMin*(alpha-alphaPrev);
			//Check if we're done (Wolfe conditions):
			if(E <= E0 + p.wolfeEnergy*alpha*gdotd0 && gdotd >= p.wolfeGradient*gdotd0)
			{	if(p.updateTestStepSize) alphaT = alphaNew; //save computed optimum step size for next time
				return true;
			}
			fprintf(p.fpLog, "%s\tWolfe criterion not satisfied: alpha: %lg  (E-E0)/|gdotd0|: %lg  gdotd/gdotd0: %lg (taking cubic step)\n",
				p.linePrefix, alpha, (E-E0)/fabs(gdotd0), gdotd/gdotd0); fflush(p.fpLog);
			//Prepare for next step:
			if(E<Eprev) { alphaPrev = alpha; Eprev=E; gdotdPrev=gdotd; } //pick the lower energy points amongst E and Eprev as the next 'prev' point
			alpha = alphaNew; //set the alpha chosen above as the other endpoint for the next cubic
		}
		//Check if current state is invalid or worse than starting point:
		if(!std::isfinite(E) || E>E0) return false; //minimize will roll back to the last known good state
		else return true;
	}
}

template<typename Vector> typename Minimizable<Vector>::Linmin Minimizable<Vector>::getLinmin(const MinimizeParams& p) const
{	using namespace MinimizePrivate;
	switch(p.linminMethod)
	{	case MinimizeParams::DirUpdateRecommended:
		{	switch(p.dirUpdateScheme)
			{	case MinimizeParams::SteepestDescent: return linminRelax<Vector>;
				case MinimizeParams::LBFGS: return linminCubicWolfe<Vector>;
				default: return linminQuad<Vector>; //Default for all nonlinear CG methods
			}
		}
		case MinimizeParams::Relax: return linminRelax<Vector>;
		case MinimizeParams::Quad: return linminQuad<Vector>;
		case MinimizeParams::CubicWolfe: return linminCubicWolfe<Vector>;
	}
	return 0;
}

#endif //JDFTX_CORE_MINIMIZE_LINMIN_H
