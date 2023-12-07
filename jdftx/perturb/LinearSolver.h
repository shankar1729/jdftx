/*-------------------------------------------------------------------
Copyright 2023 Brandon Li

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

#ifndef PERTURB_LINEARSOLVER_H_
#define PERTURB_LINEARSOLVER_H_

#include <cstdio>

//! @addtogroup Algorithms
//! @{

//! Parameters to control the minimization algorithm
struct SolverParams
{
	int nIterations;
	FILE* fpLog; //!< Stream to log iterations to
	const char* linePrefix; //!< prefix for each output line of minimizer, useful for nested minimizations (default "CG\t")
	
	enum Algorithms
	{	CG,
		MINRES
	} algorithm;
	
	double residualTol;
	double residualDiffThreshold;
	bool CGBypass;
	bool recomputeResidual;
	
	//! Set the default values
	SolverParams() 
	: nIterations(0), fpLog(stdout), linePrefix("Linear solver:\t"), algorithm(MINRES), residualTol(1e-4), residualDiffThreshold(1e-4), CGBypass(false), recomputeResidual(false) {}
};

/** Interface for symmetric (possibly indefinite) linear systems of the form Ax+b = 0.
	Choose between conjugate gradient or MINRES algorithms. Used in variational perturbation solver.
	CG implementation more memory-efficient compared to LinearSolvable. MINRES uses more memory but converges smoothly.
	@tparam Vector Same requirements as the Vector for #Minimizable
*/
template<typename Vector> struct LinearSolvableIndefinite 
{
	Vector state; //!< Solution to the system A.state+b = 0
	
	//! Linear operator A applied to the vector v.
	virtual void hessian(Vector& Av, const Vector& v) = 0;
	
	//! Apply precondioner to vector v. Required for CG algorithm
	virtual void precondition(Vector& v) = 0;
	
	//! Apply square root of preconditioner to v. Required for MINRES algorithm
	virtual void Ksqrt(Vector& v) = 0;
	
	//! Override to synchronize scalars over MPI processes (if the same minimization is happening in sync over many processes)
	virtual double sync(double x) const { return x; }
	
	//! Solve the linear system Ax+b=0 using algorithm specified by params
	//! @return the number of iterations taken to achieve target tolerance
	int solve(const Vector& b, const SolverParams& params);
	
	//! Solve Ax+b=0 using linear conjugate gradients.
	//! @return the number of iterations taken to achieve target tolerance
	int CG(const Vector& b, const SolverParams& params);
	
	//! Solve Ax+b=0 using MINRES algorithm.
	//! @return the number of iterations taken to achieve target tolerance
	int MINRES(const Vector& b, const SolverParams& params);
};

//! @}
//---------------------- Implementation ----------------------------
//!@cond


template<typename Vector> int LinearSolvableIndefinite<Vector>::solve(const Vector& b, const SolverParams& p)
{
	switch(p.algorithm) {
		case SolverParams::CG: return CG(b, p);
		case SolverParams::MINRES: return MINRES(b,p);
		default: return 0;
	}
}

template<typename Vector> int LinearSolvableIndefinite<Vector>::CG(const Vector& b, const SolverParams& p)
{	fprintf(p.fpLog, "\nBeginning conjugate gradient minimization.\n");
	
	//Initialize:
	Vector r; hessian(r, state); axpy(1.0, b, r); r *= -1.0; //residual r = -(A.state + b);
	Vector z = clone(r); precondition(z);
	Vector d, w;
	double beta=0.0, rdotzPrev=0.0, rdotz = sync(dot(r, z));

	Vector optstate; double optresidual = rdotz; int optiter = 0;
	int consecErrIncreases = 0;
	
	double rKr;
	{
		Vector Kb = clone(b);
		precondition(Kb);
		rKr = sync(dot(b, Kb));
	}
	
	//Check initial residual
	double relErr = sqrt(fabs(rdotz)/fabs(rKr));
	double relErrPrev = relErr;
	fprintf(p.fpLog, "%sInitial:  ||r||_k/||b||_k: %12.6le\n", p.linePrefix, relErr); fflush(p.fpLog);
	if(relErr<p.residualTol) { fprintf(p.fpLog, "%sConverged ||r||_k/||b||_k<%le\n", p.linePrefix, p.residualTol); fflush(p.fpLog); return 0; }

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
		hessian(w, d);
		double alpha = rdotz/sync(dot(w,d));
		axpy(alpha, d, state);
		axpy(-alpha, w, r);
		
		if (iter > 1 && iter % 10 == 0 && p.recomputeResidual) {
			hessian(r, state); axpy(1.0, b, r); r *= -1.0;
			fprintf(p.fpLog, "Recomputing residual.\n");
		}
		
		z = clone(r);
		precondition(z);
		rdotzPrev = rdotz;
		rdotz = sync(dot(r, z));
		
		//Print info:
		double relErr = sqrt(fabs(rdotz)/fabs(rKr));
		//double relErr = sqrt(fabs(sync(dot(z,z)))/fabs(rKr));
		fprintf(p.fpLog, "%sIter: %3d  ||r||_k/||b||_k: %12.6le  alpha: %12.6le  beta: %13.6le  t[s]: %9.2lf\n",
			p.linePrefix, iter, relErr, alpha, beta, clock_sec()); fflush(p.fpLog);

		//Check convergence:
		if(relErr<p.residualTol) { fprintf(p.fpLog, "%sConverged ||r||_k/||b||_k<%le\n", p.linePrefix, p.residualTol); fflush(p.fpLog); return iter; }
		
		if(relErr > relErrPrev) {
			if (consecErrIncreases == 0 && rdotzPrev < optresidual) {
				axpy(-alpha, d, state);
				optstate = state;
				optiter = iter-1;
				optresidual = rdotzPrev;
				axpy(alpha, d, state);
			}
			consecErrIncreases++;
		}
		else
			consecErrIncreases = 0;
		
		if (consecErrIncreases >= 3 && !p.CGBypass) { fprintf(p.fpLog, "Terminating solver as relative error has increased 3 times in a row.\n"); ; fflush(p.fpLog); break; }
		relErrPrev = relErr;
	}
	fprintf(p.fpLog, "%sGradient did not converge within threshold in %d iterations\n", p.linePrefix, iter); fflush(p.fpLog);
	if (optiter > 0) {
		state = optstate;
		fprintf(p.fpLog, "Setting final state to output of iteration %d with |r||_k/||b||_k = %g.\n", optiter, sqrt(optresidual)/fabs(rKr));
	}
	return iter;
}



template<typename Vector> int LinearSolvableIndefinite<Vector>::MINRES(const Vector& b, const SolverParams& p)
{	fprintf(p.fpLog, "\nBeginning MINRES solve.\n");
	Vector r = clone(state);
	
	Ksqrt(r); hessian(r, r); r *= -1.0; axpy(-1.0, b, r); Ksqrt(r); //r = -Ksqrt.b - Ksqrt.A.Ksqrt.x
	double rho=sqrt(dot(r,r));
	Vector v = clone(r); v *= (1/rho);
	double beta=0; double betatilde=0; double c=-1; double s=0;
	Vector vold=clone(state); Vector w=clone(state); Vector wtt=clone(v);
	Vector vhat, wt;
	
	int consecErrIncreases = 0;
	
	double rKr;
	{
		Vector Kb = clone(b);
		precondition(Kb);
		rKr = sync(dot(b, Kb));
	}
	double relErr = rho/sqrt(fabs(rKr));
	double relErrPrev = relErr;
	
	int iter;
	for(iter=0; iter<p.nIterations && !killFlag; iter++)
	{
		vhat = clone(v);
		Ksqrt(vhat); hessian(vhat, vhat); Ksqrt(vhat);
		axpy(-beta, vold, vhat);
		
		double alpha=dot(v,vhat);
		axpy(-alpha, v, vhat);
		beta=sqrt(dot(vhat,vhat));
		vold=clone(v);
		v = clone(vhat); v *= (1/beta);
		double l1=s*alpha-c*betatilde;
		double l2=s*beta;
		double alphatilde=-s*betatilde-c*alpha;
		betatilde=c*beta;
		double l0=sqrt(alphatilde*alphatilde+beta*beta);
		c=alphatilde/l0;
		s=beta/l0;
		wt=clone(wtt); axpy(-l1, w, wt);
		wtt=clone(v); axpy(-l2,w, wtt);
		w = clone(wt); w *= (1/l0);
		axpy(rho*c,w,state);
		rho=s*rho;
		
		double relErr = rho/sqrt(fabs(rKr));
		fprintf(p.fpLog, "%sIter: %3d  ||r||_k/||b||_k: %12.6le t[s]: %9.2lf\n",
			p.linePrefix, iter, relErr, clock_sec()); fflush(p.fpLog);
		if (relErr < p.residualTol) { fprintf(p.fpLog, "%sConverged ||r||_k/||b||_k<%le\n", p.linePrefix, p.residualTol); fflush(p.fpLog); Ksqrt(state); return iter; }
		
		if((relErrPrev - relErr)/relErr < p.residualDiffThreshold)
			consecErrIncreases++;
		else
			consecErrIncreases = 0;
		
		if (consecErrIncreases >= 2) { fprintf(p.fpLog, "Terminating solver since (||r_prev||-||r||)/||r||<%g for two iterations in a row.\n", p.residualTol); fflush(p.fpLog); break; }
		relErrPrev = relErr;
	}
	
	Ksqrt(state);
	fprintf(p.fpLog, "%sGradient did not converge within threshold in %d iterations\n", p.linePrefix, iter); fflush(p.fpLog);
	return iter;
}

//!@endcond
#endif /* PERTURB_LINEARSOLVER_H_ */
