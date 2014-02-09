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
		double gdotd = obj.sync(dot(g,d)); //directional derivative at starting point
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
			ET = obj.sync(obj.compute(0));
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
				E = obj.sync(obj.compute(&g));
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
			E = obj.sync(obj.compute(&g));
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
		double gdotdPrev = obj.sync(dot(g,d)); //directional derivative at starting point
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
			E = obj.sync(obj.compute(&g));
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
			double gdotd = obj.sync(dot(g,d));
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

template<typename Vector> typename Minimizable<Vector>::Linmin Minimizable<Vector>::getLinmin(const MinimizeParams& p) const
{	using namespace MinimizePrivate;
	switch(p.linminMethod)
	{	case MinimizeParams::Relax: return linminRelax<Vector>; break;
		case MinimizeParams::Quad:  return linminQuad<Vector>; break;
		case MinimizeParams::Cubic: return linminCubic<Vector>; break;
	}
	return 0;
}

#endif //JDFTX_CORE_MINIMIZE_LINMIN_H
