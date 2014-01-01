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

#include <electronic/Everything.h>
#include <electronic/matrix.h>
#include <electronic/SpeciesInfo.h>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <limits>
#include <list>

int ElecInfo::whose(int q) const
{	if(mpiUtil->nProcesses()>1) return std::upper_bound(qStopArr.begin(),qStopArr.end(), q) - qStopArr.begin();
	else return 0;
}

int ElecInfo::findHOMO(int q) const
{	int HOMO = 0;
	for(int n=(e->eVars.F[q].size()-1); n>=0; n--)
		if(e->eVars.F[q][n] > 1e-3){ HOMO = n; break;	}
	return HOMO;
}

ElecInfo::ElecInfo()
: nBands(0), nStates(0), qStart(0), qStop(0), spinType(SpinNone), spinRestricted(false), nElectrons(0), 
fillingsUpdate(ConstantFillings), kT(1e-3), mu(std::numeric_limits<double>::quiet_NaN()),
mixInterval(0), subspaceRotation(false), hasU(false), nBandsOld(0),
fillingMixFraction(0.5), dnPrev(mu), muMeasuredPrev(mu), //set dnPrev and muMeasuredPrev to NaN
Cmeasured(1.), Cweight(0.), dnMixFraction(0.7)
{
}

void ElecInfo::setup(const Everything &everything, std::vector<diagMatrix>& F, Energies& ener)
{	e = &everything;
	betaBy2 = 0.5/kT;
	
	logPrintf("\n---------- Setting up k-points, bands, fillings ----------\n");
	if(spinRestricted) assert(spinType==SpinZ);
	
	//k-points are folded before symmetry setup, now reduce them under symmetries:
	kpointsReduce();
	nStates = qnums.size();
	
	//Determine distribution amongst processes: (TODO: weight division by node capability in heterogeneous contexts)
	qStart = (nStates * mpiUtil->iProcess()) / mpiUtil->nProcesses();
	qStop = (nStates * (mpiUtil->iProcess()+1)) / mpiUtil->nProcesses();
	qStopArr.resize(mpiUtil->nProcesses());
	for(int iProc=0; iProc<mpiUtil->nProcesses(); iProc++)
		qStopArr[iProc] = (nStates * (iProc+1)) / mpiUtil->nProcesses();
	
	// allocate the fillings matrices.
	F.resize(nStates);

	// Convert initial charge(s) qNet to the correct size for spin / no-spin:
	if(spinType==SpinNone && qNet.size()==2) { qNet[0]+=qNet[1]; qNet.pop_back(); }
	if(spinType==SpinZ && qNet.size()==1) { qNet.push_back(qNet[0]*=0.5); }
	if(!qNet.size()) qNet.assign(spinType==SpinNone ? 1 : 2, 0.);
	double wInv = spinType==SpinNone ? 0.5 : 1.0; //normalization factor from external to internal fillings
	
	// Figure out the number of bands and number of electrons
	logPrintf("Computing the number of bands and number of electrons\n");
	
	//--- Compute nElectrons and nBands based on atomic valencies:
	std::vector<double> nsElectrons = qNet; //manually specified excess charge(s)
	nElectrons = 0;
	for(unsigned s=0; s<nsElectrons.size(); s++) //each spin channel
	{	for(auto sp: e->iInfo.species)
			nsElectrons[s] += sp->Z * sp->atpos.size() * (1./nsElectrons.size());
		if(nsElectrons[s]<0)
			die("Number of electrons in spin channel %d is negative (= %lg)\n", s, nsElectrons[s]);
		nElectrons += nsElectrons[s];
	}
	
	//--- Calculate nBands if elec-n-bands not given:
	if(!nBands)
	{	if(fillingsUpdate == ConstantFillings) //pick nBands to just accommodate spin channel with most electrons
		{	double nsElectronsMax = *std::max_element(nsElectrons.begin(), nsElectrons.end());
			nBands = std::max(1, (int)ceil(nsElectronsMax*wInv));
		}
		else //set to number of atomic orbitals (which will ensure that complete bands are included)
			nBands = e->iInfo.nAtomicOrbitals(); //this estimate is usually on the high side, but it leads to better convergence than a stingier value
	}
	
	//--- No initial fillings, fill the lowest orbitals in each spin channel:
	if(!initialFillingsFilename.length())
	{	logPrintf("Calculating initial fillings.\n");
		for(int q=qStart; q<qStop; q++)
		{	F[q].assign(nBands, 0.);
			double ftot = nsElectrons[qnums[q].index()] * wInv; //total of all fillings in this state
			for(int b = 0; b<nBands; b++)
			{	if(ftot>1) F[q][b]=1;
				else if(ftot>0) { F[q][b]=ftot; break; }
				ftot -= F[q][b];
			}
		}
	}
	else
	{	logPrintf("Reading initial fillings from file %s.\n", initialFillingsFilename.c_str());
		if(nBandsOld <= 0) nBandsOld=nBands;
		read(F, initialFillingsFilename.c_str(), nBandsOld);
		
		for(int q=qStart; q<qStop; q++)
		{	F[q] *= wInv; //NOTE: fillings are always 0 to 1 internally, but read/write 0 to 2 for SpinNone
			
			//Check fillings:
			for(int b=0; b<nBandsOld; b++)
			{	if(F[q][b]<0) die("Filling of state %d band %d is negative.\n", q, b);
				if(F[q][b]>1) die("Filling of state %d band %d exceeds maximum.\n", q, b);
			}
			
			//Determine the total change in sum(F[q][0:nBands-1]):
			double sumFchange = qNet[qnums[q].index()] * wInv; //manual excess fillings
			for(int b=nBands;  b<nBandsOld; b++)
				sumFchange += F[q][b]; //discarded bands
			
			F[q].resize(nBands, 0.); //drop extra entries, or pad with 0s
			
			//If depleting electrons, do it from the highest (least occupied) bands
			if(sumFchange<0.)
			{	double qRemaining = -sumFchange;
				for(int b=nBands-1; b>=0; b--)
				{	if(F[q][b]<qRemaining) { qRemaining-=F[q][b]; F[q][b]=0.; }
					else { F[q][b]-=qRemaining; qRemaining=0.; break; }
				}
				if(qRemaining>0.)
					die("elec-initial-charge wants to deplete more electrons than available in state %d.\n"
						"   Note: automatic fillings charging does not redistribute fillings between\n"
						"   kpoints. You may do that manually, but maybe that's not a good idea.\n", q);
			}
			
			//If adding electrons, saturate the lowest (most occupied) bands
			if(sumFchange>0.)
			{	double qRemaining = sumFchange;
				for(int b=0; b<nBands; b++)
				{	if(F[q][b]+qRemaining>1) { qRemaining-=(1.-F[q][b]); F[q][b]=1.; }
					else { F[q][b]+=qRemaining; qRemaining=0.; break; }
				}
				if(qRemaining>0.)
					die("elec-initial-charge / nBands change over-saturates state %d.\n"
						"   Note: automatic fillings charging / band change does not redistribute fillings\n"
						"   between kpoints. You may do that manually, but maybe that's not a good idea.\n", q);
			}
		}
		
		//Compute electron count (override the previous atom-valence + qNet based value):
		nElectrons=0.;
		for(int q=qStart; q<qStop; q++)
			nElectrons += qnums[q].weight * trace(F[q]);
		mpiUtil->allReduce(nElectrons, MPIUtil::ReduceSum, true);
	}
	
	//subspace_rotation is false by default
	if(fillingsUpdate != ConstantFillings)
	{	//make sure there are extra bands:
		int nBandsMin = (int)ceil(nElectrons/2.0);
		if(nBands <= nBandsMin)
			die("%d bands insufficient for fermi fillings with %lg electrons (need at least %d, recommend > %d)\n",
				nBands, nElectrons, nBandsMin+1, nBandsMin+5);
		
		subspaceRotation = true;
		logPrintf("Turning on subspace rotations for fermi fillings.\n");
	}
	else
	{	//make sure there are sufficient bands:
		int nBandsMin = (int)ceil(nElectrons/2.0);
		if(nBands < nBandsMin)
			die("%d bands insufficient for %lg electrons (need at least %d)\n",
				nBands, nElectrons, nBandsMin);
		
		if(!e->cntrl.fixed_n) //Band structure never requires subspace-rotations
		{	//Check if fillings are scalar (then subspace rotations are not needed as [Hsub,F]=0)
			bool scalarFillings = true;
			for(int q=qStart; q<qStop; q++)
				scalarFillings &= F[q].isScalar();
			mpiUtil->allReduce(scalarFillings, MPIUtil::ReduceLAnd);
			if(!scalarFillings)
			{	subspaceRotation = true;
				logPrintf("Turning on subspace rotations due to non-diagonal fillings.\n");
			}
		}
	}
	
	//Set the Legendre multipliers corresponding to the initial fillings
	updateFillingsEnergies(F, ener);
	
	// Print out the current status of the electronic info before leaving
	logPrintf("nElectrons: %10.6f   nBands: %d   nStates: %d", nElectrons, nBands, nStates);
	if(e->cntrl.shouldPrintEigsFillings)
	{	logPrintf("   initial fillings:\n");
		printFillings(globalLog);
	}
	else logPrintf("\n");
}

void ElecInfo::printFillings(FILE* fp) const
{	//NOTE: fillings are always 0 to 1 internally, but read/write 0 to 2 for SpinNone
	if(mpiUtil->isHead())
		for(int q=0; q<nStates; q++)
		{	diagMatrix const* Fq = &e->eVars.F[q];
			diagMatrix FqTemp;
			if(!isMine(q))
			{	FqTemp.resize(nBands);
				FqTemp.recv(whose(q));
				Fq = &FqTemp;
			}
			((*Fq) * (spinType==SpinNone ? 2 : 1)).print(fp, "%20.14le ");
		}
	else
		for(int q=qStart; q<qStop; q++)
			e->eVars.F[q].send(0);
}

void ElecInfo::mixFillings(std::vector<diagMatrix>& F, Energies& ener)
{	const ElecVars &eVars = e->eVars;
	double muMix; //A chemical potential whose corresponding fermi distribution is mixed in
	
	if(!std::isnan(mu)) //Fixed mu:
	{	//Fit a fermi-dirac distribution to the current fillings:
		double g0, &C = Cmeasured; //density of states at current mu, and measured electrostatic capacitance
		double muMeasured = fitMu(F, eVars.Hsub_eigs, &g0); //this is the "measured" mu that fits the fillings

		//If not the first iteration, update estimate of Ceff
		if(!std::isnan(dnPrev))  //improve the capacitance estimate
		{	double curC = dnPrev/(muMeasured - muMeasuredPrev);  //observed total effective capacitance
			double curWeight = fabs(dnPrev); //weight of current measurement
			C = (Cweight + curWeight)/(Cweight/C + curWeight/curC);
			Cweight += curWeight;
			Cweight *= 0.9; //gradually forget the older measurements
		}
		muMeasuredPrev = muMeasured;
	
		//Determine the amount to change the charge by, and adjust chemical potential to get that charge
		double dn = dnMixFraction * C * (mu-muMeasured);
		muMix = findMu(eVars.Hsub_eigs, nElectrons + dn/fillingMixFraction); //this is the mu with which the fillings are mixed in
		nElectrons += dn;
		dnPrev = dn;
		logPrintf("FillingsMix:  MeasuredMu: %.15le  g0: %le  C: %le  dn: %.15le  n: %.15le  SetMu: %.15le\n",
			muMeasured, g0, C, dn, nElectrons, muMix);
	}
	else
	{	//Fixed charge, just bisect to find mu:
		muMix = findMu(eVars.Hsub_eigs, nElectrons);
		logPrintf("FillingsMix:  mu: %.15le  nElectrons: %.15le", muMix, nElectrons);
		if(spinType == SpinZ)
		{	double spinPol = integral(e->eVars.n[0] - e->eVars.n[1]);
			logPrintf("  magneticMoment: %.5f", spinPol);
		}
		logPrintf("\n");
	}
	
	//Mix in fillings at new mu:
	for(int q=qStart; q<qStop; q++)
		F[q] = (1-fillingMixFraction)*F[q] + fillingMixFraction*fermi(muMix, eVars.Hsub_eigs[q]);
	
	if(e->cntrl.shouldPrintEigsFillings)
	{	logPrintf("\nCorresponding fillings:\n");
		printFillings(globalLog);
	}
	
	// Update fillings contribution to free energy:
	updateFillingsEnergies(F, ener);
}


// Fermi legendre multipliers (TS and optionally muN)
void ElecInfo::updateFillingsEnergies(const std::vector<diagMatrix>& F, Energies& ener) const
{
	ener.TS = 0.0;
	for(int q=qStart; q<qStop; q++)
	{	double Sq = 0.0;
		for(double Fqi: F[q])
		{	if(Fqi>1e-300) Sq += Fqi*log(Fqi);
			if(1-Fqi>1e-300) Sq += (1-Fqi)*log(1-Fqi);
		}
		ener.TS -= kT * qnums[q].weight * Sq;
	}
	mpiUtil->allReduce(ener.TS, MPIUtil::ReduceSum);

	//Grand canonical multiplier if fixed mu:
	if(!std::isnan(mu)) ener.muN = mu * nElectrons;
}

//-------------------- Fermi function utilities --------------------------

diagMatrix ElecInfo::fermi(double mu, const diagMatrix& eps) const
{	diagMatrix ret(eps);
	for(unsigned i=0; i<eps.size(); i++) ret[i] = fermi(mu, eps[i]);
	return ret;
}

diagMatrix ElecInfo::fermiPrime(double mu, const diagMatrix& eps) const
{	diagMatrix ret(eps);
	for(unsigned i=0; i<eps.size(); i++) ret[i] = fermiPrime(mu, eps[i]);
	return ret;
}

matrix ElecInfo::fermiGrad(double mu, const diagMatrix& eps, const matrix& gradF) const
{	matrix gradEps(gradF); //copy input
	complex* gradEpsData = gradEps.data();
	//calculate result in place:
	for(int i=0; i<gradF.nRows(); i++)
		for(int j=0; j<gradF.nCols(); j++)
		{	double deps = eps[i] - eps[j];
			gradEpsData[gradEps.index(i,j)]
				*= (fabs(deps)<1e-6) ? fermiPrime(mu,eps[i]) : (fermi(mu,eps[i])-fermi(mu,eps[j]))/deps;
		}
	return gradEps;
}


//Number of electrons in a fermi distribution of given mu and eigenvalues:
double ElecInfo::nElectronsFermi(double mu, const std::vector<diagMatrix>& eps) const
{	double N = 0.0;
	for(int q=qStart; q<qStop; q++)
		for(double epsCur: eps[q])
			N += qnums[q].weight*fermi(mu, epsCur);
	mpiUtil->allReduce(N, MPIUtil::ReduceSum, true);
	return N;
}

//Return the mu that would match Ntarget at the current eigenvalues (bisection method)
double ElecInfo::findMu(const std::vector<diagMatrix>& eps, double nElectrons) const
{	const bool& verbose = e->cntrl.shouldPrintMuSearch;
	if(verbose) logPrintf("\nBisecting to find mu(nElectrons=%.15le)\n", nElectrons);
	//Find a range which is known to bracket the result:
	double muMin=-0.1; while(nElectronsFermi(muMin,eps)>=nElectrons) muMin-=0.1;
	double muMax=+0.0; while(nElectronsFermi(muMax,eps)<=nElectrons) muMax+=0.1;
	//Bisect:
	while(muMax-muMin>=1e-10*kT)
	{	double mu = 0.5*(muMin + muMax);
		double N = nElectronsFermi(mu, eps);
		if(verbose) logPrintf("MUBISECT: mu = [ %.15le %.15le %.15le ]  N = %le\n", muMin, mu, muMax, N);
		if(N>nElectrons) muMax = mu;
		else muMin = mu;
	}
	return 0.5*(muMin + muMax);
}


//-------------- Fermi function fit -------------------

#include <gsl/gsl_multifit_nlin.h>

struct FitParams
{	const ElecInfo& eInfo;
	const std::vector<diagMatrix>& F;
	const std::vector<diagMatrix>& eps;
};

//Compute residuals (error in each filling for the guessed mu)
int fitMu_fdf(const gsl_vector* x, void* params, gsl_vector* residual, gsl_matrix* jacobian)
{	const FitParams& p = *((const FitParams*)params);
	double mu = gsl_vector_get(x, 0);
	int resIndex = p.eInfo.nBands*p.eInfo.qStart;
	for(int q=p.eInfo.qStart; q<p.eInfo.qStop; q++)
	{	for(int b=0; b<p.eInfo.nBands; b++)
		{	//Compute the fermi function and derivative
			double f = p.eInfo.fermi(mu, p.eps[q][b]);
			double fPrime = -p.eInfo.fermiPrime(mu, p.eps[q][b]); //note: fprime is df/deps (we need df/dmu)
			double weight = p.eInfo.qnums[q].weight;
			if(residual) gsl_vector_set(residual, resIndex, weight*(f - p.F[q][b]));
			if(jacobian) gsl_matrix_set(jacobian, resIndex, 0, weight*fPrime);
			resIndex++;
		}
	}
	if(mpiUtil->nProcesses()>1)
	{	for(int iSrc=0; iSrc<mpiUtil->nProcesses(); iSrc++)
		{	int qStart = p.eInfo.qStartOther(iSrc);
			int qCount = p.eInfo.qStopOther(iSrc) - qStart;
			mpiUtil->bcast(gsl_vector_ptr(residual, qStart*p.eInfo.nBands),    qCount*p.eInfo.nBands, iSrc);
			mpiUtil->bcast(gsl_matrix_ptr(jacobian, qStart*p.eInfo.nBands, 0), qCount*p.eInfo.nBands, iSrc);
		}
	}
	return 0;
}
int fitMu_f(const gsl_vector* x, void* params, gsl_vector* residual) { return fitMu_fdf(x,params,residual,0); }
int fitMu_df(const gsl_vector* x, void* params, gsl_matrix* jacobian) { return fitMu_fdf(x,params,0,jacobian); }

double nrm2(gsl_vector* f) { return eblas_dnrm2(f->size, f->data, f->stride); }

double ElecInfo::fitMu(const std::vector<diagMatrix>& F, const std::vector<diagMatrix>& eps, double* dndmu) const
{	const bool& verbose = e->cntrl.shouldPrintMuSearch;
	//Set up GSL fit utility:
	unsigned nResiduals = nStates*nBands;
	gsl_multifit_fdfsolver *solver = gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, nResiduals, 1);
	FitParams params = {*this, F, eps};
	gsl_multifit_function_fdf fitFunc = {&fitMu_f, &fitMu_df, &fitMu_fdf, nResiduals, 1, (void*)&params};
	gsl_vector_const_view xInit = gsl_vector_const_view_array(&mu, 1);
	gsl_multifit_fdfsolver_set(solver, &fitFunc, &xInit.vector);
	//Fit:
	if(verbose) logPrintf("\nFitting Fermi-dirac distribution to current fillings and eigenvalues\n");
	bool converged=false;
	for(int iter=0; iter<500; iter++)
	{
		gsl_multifit_fdfsolver_iterate(solver);
		if(verbose)
			logPrintf("MUFIT %3d: mu = %.15le Residual = %le\n", iter, gsl_vector_get(solver->x, 0), nrm2(solver->f));
		if(gsl_multifit_test_delta(solver->dx, solver->x, 1e-14, 1e-14) != GSL_CONTINUE)
		{	converged=true;
			break;
		}
	}
	if(!converged) die("Convergence failure in Fermi-dirac distribution fit to the fillings.\n");
	double mu = gsl_vector_get(solver->x, 0);
	if(dndmu) //dn/dmu is just the sum of w df/dmu, i.e. the sum of all elements of the Jacobian
	{	*dndmu = 0.0;
		for(unsigned i=0; i<nResiduals; i++) *dndmu += gsl_matrix_get(solver->J, i, 0);
	}
	gsl_multifit_fdfsolver_free(solver);
	return mu;
}

// Positive remainder
double fmodPositive(double x, double y)
{	double ret = fmod(x, y);
	if(ret<0) ret+=fabs(y);
	return ret;
}

void ElecInfo::kpointsFold()
{	//Check and regularize the supplied k-points:
	double wSum = 0.0;
	for(unsigned j=0; j<qnums.size(); j++)
	{	wSum += qnums[j].weight;
		for(int k=0; k<3; k++)
			qnums[j].k[k] = fmodPositive(qnums[j].k[k], 1.0);
	}
	if(fabs(wSum-1.0)>1e-12) die("Specified kpoint weights do not add to 1\n");
	//Create the folded array:
	logPrintf("\nFolded %lu k-points by %dx%dx%d to ", qnums.size(), kfold[0], kfold[1], kfold[2]);
	int nFold = kfold[0]*kfold[1]*kfold[2];
	std::vector<QuantumNumber> qnumsFolded;
	qnumsFolded.reserve(nFold*qnums.size());
	vector3<int> i;
	for(unsigned j=0; j<qnums.size(); j++)
	{	QuantumNumber q = qnums[j];
		q.weight /= nFold;
		for(i[0]=0; i[0]<kfold[0]; i[0]++)
			for(i[1]=0; i[1]<kfold[1]; i[1]++)
				for(i[2]=0; i[2]<kfold[2]; i[2]++)
				{	for(int k=0; k<3; k++)
					{	q.k[k] = (qnums[j].k[k] + i[k]) / kfold[k];
						if(q.k[k]>0.5) q.k[k]-=1.; //reduce to centered unit interval
					}
					qnumsFolded.push_back(q);
				}
	}
	//Exchange the folded and unfolded k-points:
	std::swap(qnums, qnumsFolded);
	//Output the folded kpoint coordinates
	logPrintf("%lu k-points.\n", qnums.size());
}


void ElecInfo::kpointsReduce()
{	//Reduce under symmetries:
	std::list<QuantumNumber> qlist = e->symm.reduceKmesh(qnums);
	if(qnums.size()==qlist.size())
		logPrintf("No reducable k-points. ");
	else
		logPrintf("Reduced to %lu k-points under symmetry. ", qlist.size());
	qnums.assign(qlist.begin(), qlist.end()); //Convert back to a vector for efficient addressing
	//Handle spin states / spin dgeneracy:
	if(spinType==SpinNone)
	{	for(unsigned i=0; i<qnums.size(); i++)
			qnums[i].weight *= 2; //double the weight for each state to account for spin degeneracy
	}
	else //SpinZ
	{	unsigned nkPoints = qnums.size();
		qnums.insert(qnums.end(), qlist.begin(), qlist.end()); //add second copy for other spin
		for(unsigned ik=0; ik<nkPoints; ik++)
		{	//Set spin values:
			qnums[ik].spin = +1;
			qnums[ik+nkPoints].spin = -1;
		}
	}
	//Output the reduced kpoints:
	logPrintf("States including spin/spin-weights:\n");
	kpointsPrint(true);
}

void ElecInfo::kpointsPrint(bool printSpin) const
{	for(unsigned q=0; q<qnums.size(); q++)
	{	kpointPrint(q); logPrintf("\n");
	}
}

void ElecInfo::kpointPrint(int q, bool printSpin) const
{
	logPrintf("%5d\t[ %10.6f %10.6f %10.6f ] %8.6f", q,
			qnums[q].k[0], qnums[q].k[1], qnums[q].k[2], qnums[q].weight);
		if(printSpin) logPrintf("  spin %2d", qnums[q].spin);
}


//-------------- diagMatrix/matrix array parallel I/O ------------------

void ElecInfo::read(std::vector<diagMatrix>& M, const char *fname, int nRowsOverride) const
{	int nRows = nRowsOverride ? nRowsOverride : nBands;
	M.resize(nStates);
	MPIUtil::File fp; mpiUtil->fopenRead(fp, fname, nStates*nRows*sizeof(double));
	mpiUtil->fseek(fp, qStart*nRows*sizeof(double), SEEK_SET);
	for(int q=qStart; q<qStop; q++)
	{	M[q].resize(nBands);
		mpiUtil->fread(M[q].data(), sizeof(double), nRows, fp);
	}
	mpiUtil->fclose(fp);
}

void ElecInfo::read(std::vector<matrix>& M, const char *fname, int nRowsOverride, int nColsOverride) const
{	int nRows = nRowsOverride ? nRowsOverride : nBands;
	int nCols = nColsOverride ? nColsOverride : nBands;
	M.resize(nStates);
	MPIUtil::File fp; mpiUtil->fopenRead(fp, fname, nStates*nRows*nCols*sizeof(complex));
	mpiUtil->fseek(fp, qStart*nRows*nCols*sizeof(complex), SEEK_SET);
	for(int q=qStart; q<qStop; q++)
	{	M[q].init(nRows, nCols);
		mpiUtil->fread(M[q].data(), sizeof(complex), M[q].nData(), fp);
	}
	mpiUtil->fclose(fp);
}

void ElecInfo::write(const std::vector<diagMatrix>& M, const char *fname, int nRowsOverride) const
{	int nRows = nRowsOverride ? nRowsOverride : nBands;
	assert(int(M.size())==nStates);
	MPIUtil::File fp; mpiUtil->fopenWrite(fp, fname);
	mpiUtil->fseek(fp, qStart*nRows*sizeof(double), SEEK_SET);
	for(int q=qStart; q<qStop; q++)
	{	assert(M[q].nRows()==nRows);
		mpiUtil->fwrite(M[q].data(), sizeof(double), nRows, fp);
	}
	mpiUtil->fclose(fp);
}

void ElecInfo::write(const std::vector<matrix>& M, const char *fname, int nRowsOverride, int nColsOverride) const
{	int nRows = nRowsOverride ? nRowsOverride : nBands;
	int nCols = nColsOverride ? nColsOverride : nBands;
	assert(int(M.size())==nStates);
	MPIUtil::File fp; mpiUtil->fopenWrite(fp, fname);
	mpiUtil->fseek(fp, qStart*nRows*nCols*sizeof(complex), SEEK_SET);
	for(int q=qStart; q<qStop; q++)
	{	assert(M[q].nRows()==nRows);
		assert(M[q].nCols()==nCols);
		mpiUtil->fwrite(M[q].data(), sizeof(complex), M[q].nData(), fp);
	}
	mpiUtil->fclose(fp);
}
