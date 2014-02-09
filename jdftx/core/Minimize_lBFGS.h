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

#ifndef JDFTX_CORE_MINIMIZE_LBFGS_H
#define JDFTX_CORE_MINIMIZE_LBFGS_H

#include <core/Minimize_linmin.h>
#include <core/Minimize.h> //this include is present only to aid IDE's autocompletion (does nothing since this file is always used via Minimize.h)
#include <memory>
#include <list>
#include <stack>

//Following Nocedal and Liu, Math Prog 45, 503 (1989)
template<typename Vector> double Minimizable<Vector>::lBFGS(const MinimizeParams& p)
{	
	Vector g; //gradient
	double E = sync(compute(&g)); //get initial energy and gradient
	EdiffCheck ediffCheck(p.nEnergyDiff, p.energyDiffThreshold); //list of past energies

	double alphaT = p.alphaTstart; //test step size
	double alpha = alphaT; //actual step size
	double linminTest = 0.;
	
	//History of variable and residual changes:
	struct History
	{	Vector s; //change in variable (= alpha d)
		Vector y; //change in residual (= g - gPrev)
		double rho; //= 1/dot(s,y)
		double gamma; //= dot(s,y)/dot(y,y)
	};
	std::list< std::shared_ptr<History> > history;
	
	//Select the linmin method:
	Linmin linmin = getLinmin(p);
	
	//Iterate until convergence, max iteration count or kill signal
	int iter=0;
	for(iter=0; !killFlag; iter++)
	{	
		if(report(iter)) //optional reporting/processing
		{	E = sync(compute(&g)); //update energy and gradient if state was modified
			fprintf(p.fpLog, "%s\tState modified externally: resetting history.\n", p.linePrefix);
			fflush(p.fpLog);
			history.clear();
		}
		
		double gNorm = sync(dot(g,g));
		fprintf(p.fpLog, "%sIter: %3d  %s: %22.15le  |grad|: %10.3le  alpha: %10.3le",
			p.linePrefix, iter, p.energyLabel, E, sqrt(gNorm/p.nDim), alpha);
		if(linminTest) fprintf(p.fpLog, ", linmin = %10.3le", linminTest);
		
		//Check stopping conditions:
		fprintf(p.fpLog, "\n"); fflush(p.fpLog);
		if(sqrt(gNorm/p.nDim) < p.knormThreshold)
		{	fprintf(p.fpLog, "%sConverged (|grad|<%le).\n", p.linePrefix, p.knormThreshold);
			fflush(p.fpLog); return E;
		}
		if(ediffCheck.checkConvergence(E))
		{	fprintf(p.fpLog, "%sConverged (|Delta %s|<%le for %d iters).\n",
				p.linePrefix, p.energyLabel, p.energyDiffThreshold, p.nEnergyDiff);
			fflush(p.fpLog); return E;
		}
		if(!std::isfinite(gNorm))
		{	fprintf(p.fpLog, "%s|grad|=%le. Stopping ...\n", p.linePrefix, gNorm);
			fflush(p.fpLog); return E;
		}
		if(!std::isfinite(E))
		{	fprintf(p.fpLog, "%sE=%le. Stopping ...\n", p.linePrefix, E);
			fflush(p.fpLog); return E;
		}
		if(iter>=p.nIterations) break;
		
		//Prepare container that will be committed to history below: (and avoid copies)
		auto h = std::make_shared<History>();
		Vector& gPrev = h->y;
		Vector& d = h->s;
		gPrev = clone(g);
		
		//Compute search direction:
		d = clone(g);
		std::stack<double> a; //alpha in the reference renamed to 'a' here to not clash with step size
		for(auto h=history.rbegin(); h!=history.rend(); h++)
		{	a.push( (*h)->rho * sync(dot((*h)->s, d)) );
			axpy(-a.top(), (*h)->y, d);
		}
		d = precondition(d);
		if(history.size()) //H0 = K at first iteration, but incorporate scaling on subsequent ones (scheme M3 from the paper)
			d *= history.back()->gamma;
		for(auto h=history.begin(); h!=history.end(); h++)
		{	double b = (*h)->rho * sync(dot((*h)->y, d));
			axpy(a.top()-b, (*h)->s, d);
			a.pop();
		}
		d *= -1;
		
		//Line minimization
		if(linmin(*this, p, d, alphaT, alpha, E, g))
		{	//linmin succeeded:
			if(p.updateTestStepSize)
			{	alphaT = alpha;
				if(alphaT<p.alphaTmin) //bad step size
					alphaT = p.alphaTstart; //make sure next test step size is not too bad
			}
			fflush(p.fpLog);
		}
		else
		{	//linmin failed:
			fprintf(p.fpLog, "%s\tUndoing step.\n", p.linePrefix);
			step(d, -alpha);
			E = sync(compute(&g));
			if(history.size())
			{	//Failed, but not along the gradient direction:
				fprintf(p.fpLog, "%s\tStep failed: resetting history.\n", p.linePrefix);
				fflush(p.fpLog);
				history.clear();
				linminTest = 0.;
				continue;
			}
			else
			{	//Failed along the gradient direction
				fprintf(p.fpLog, "%s\tStep failed along negative gradient direction.\n", p.linePrefix);
				fprintf(p.fpLog, "%sProbably at roundoff error limit. (Stopping)\n", p.linePrefix);
				fflush(p.fpLog);
				return E;
			}
		}
		
		//Update history:
		linminTest = sync(dot(g,d))/sqrt(sync(dot(g,g))*sync(dot(d,d)));
		d *= alpha; //d -> alpha * d, which is change of state, and it resides in h->s
		gPrev *= -1.; axpy(1., g, gPrev); //gPrev -> g - gPrev, and it resides in h->y
		double ydots = sync(dot(h->y, h->s));
		h->rho = 1./ydots;
		h->gamma = ydots / sync(dot(h->y, precondition(h->y)));
		history.push_back(h);
		while(history.size() > p.history) history.pop_front(); //dicard oldest entries
	}
	fprintf(p.fpLog, "%sNone of the convergence criteria satisfied after %d iterations.\n", p.linePrefix, iter);
	return E;
}

#endif //JDFTX_CORE_MINIMIZE_LBFGS_H
