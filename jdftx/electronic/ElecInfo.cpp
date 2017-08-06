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

#include <electronic/ElecInfo.h>
#include <electronic/Everything.h>
#include <electronic/SpeciesInfo.h>
#include <core/matrix.h>
#include <fluid/Euler.h>
#include <algorithm>
#include <limits>
#include <list>

int ElecInfo::findHOMO(int q) const
{	int HOMO = 0;
	for(int n=(e->eVars.F[q].size()-1); n>=0; n--)
		if(e->eVars.F[q][n] > 0.5){ HOMO = n; break; }
	return HOMO;
}

ElecInfo::ElecInfo()
: nBands(0), nStates(0), qStart(0), qStop(0), spinType(SpinNone), nElectrons(0), 
fillingsUpdate(FillingsConst), scalarFillings(true), smearingType(SmearingFermi), smearingWidth(1e-3), mu(NAN), muLoop(false),
hasU(false), nBandsOld(0),
Qinitial(0.), Minitial(0.), Mconstrain(false)
{
}

void ElecInfo::setup(const Everything &everything, std::vector<diagMatrix>& F, Energies& ener)
{	e = &everything;
	
	switch(spinType)
	{	case SpinNone:   nDensities = 1; spinWeight = 2; break;
		case SpinZ:      nDensities = 2; spinWeight = 1; break;
		case SpinVector: nDensities = 4; spinWeight = 1; break;
		case SpinOrbit:  nDensities = 1; spinWeight = 1; break;
	}
	
	logPrintf("\n---------- Setting up k-points, bands, fillings ----------\n");
	
	//k-points are folded before symmetry setup, now reduce them under symmetries:
	kpointsReduce();
	nStates = qnums.size();
	
	//Determine distribution amongst processes:
	qDivision.init(nStates, mpiUtil);
	qDivision.myRange(qStart, qStop);
	
	//Allocate the fillings matrices.
	F.resize(nStates);

	//Calculate initial charges per spin channel:
	std::vector<double> qNet(nSpins());
	for(int s=0; s<int(qNet.size()); s++)
		qNet[s] = (Qinitial + (s ? -1 : +1)*Minitial) / qNet.size();
	double wInv = 1./spinWeight; //normalization factor from external to internal fillings
	qWeightSum = qNet.size() * spinWeight;
	
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
	{	double nsElectronsMax = *std::max_element(nsElectrons.begin(), nsElectrons.end());
		int nBandsMin = std::max(1, (int)ceil(nsElectronsMax*wInv));
		if(fillingsUpdate == FillingsConst) //pick nBands to just accommodate spin channel with most electrons
		{	nBands = nBandsMin;
		}
		else //set to number of atomic orbitals (which will ensure that complete bands are included)
			nBands = std::max(nBandsMin+1, e->iInfo.nAtomicOrbitals()); //this estimate is usually on the high side, but it leads to better convergence than a stingier value
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
			const double fMin = 0.;
			const double fMax = 1.0833154705877; //maximum allowed for Cold smearing
			for(int b=0; b<nBandsOld; b++)
			{	if(F[q][b]<fMin) die("Filling of state %d band %d is negative.\n", q, b);
				if(F[q][b]>fMax) die("Filling of state %d band %d exceeds maximum.\n", q, b);
			}
			
			//Determine the total change in sum(F[q][0:nBands-1]):
			double sumFchange = 0.; 
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
	
	//Check that number of bands is sufficient:
	int nBandsMin = (int)ceil(nElectrons/qWeightSum);
	if(fillingsUpdate == FillingsHsub)
	{	//make sure there are extra bands:
		if(nBands*qWeightSum <= nElectrons)
			die("%d bands insufficient for variable fillings with %lg electrons (need at least %d, recommend > %d)\n",
				nBands, nElectrons, nBandsMin+1, nBandsMin+5);
	}
	else //FillingsConst:
	{	if(nBands < nBandsMin)
			die("%d bands insufficient for %lg electrons (need at least %d)\n",
				nBands, nElectrons, nBandsMin);
		
		//Check for non-scalar fillings in a variational minimize calculation:
		if(!e->cntrl.fixed_H && !e->cntrl.scf)
		{	scalarFillings = true;
			for(int q=qStart; q<qStop; q++)
				scalarFillings &= F[q].isScalar();
			mpiUtil->allReduce(scalarFillings, MPIUtil::ReduceLAnd);
			if(!scalarFillings)
				logPrintf("Turning on subspace rotations due to non-scalar fillings.\n");
		}
	}
	
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
			((*Fq) * spinWeight).print(fp, "%.15lf ");
		}
	else
		for(int q=qStart; q<qStop; q++)
			e->eVars.F[q].send(0);
}

void ElecInfo::smearReport(const double* muOverride) const
{	double Bz = 0.;
	logPrintf("\tFillingsUpdate:  mu: %+.9lf  nElectrons: %.6lf",
		muOverride ? *muOverride //use override value if provided
			: ( !std::isnan(mu) ? mu //use target mu if in fixed-mu mode
				: findMu(e->eVars.Hsub_eigs, nElectrons, Bz) ), //determine from eigenvalues otherwise
		nElectrons);
	if(spinType == SpinZ)
	{	ScalarField Mfield = e->eVars.n[0] - e->eVars.n[1];
		double Mabs = integral(sqrt(Mfield*Mfield));
		double Mtot = integral(Mfield);
		logPrintf("  magneticMoment: [ Abs: %7.5f  Tot: %+8.5f ]", Mabs, Mtot);
	}
	if(spinType == SpinVector)
	{	VectorField Mfield;
		Mfield[0] = (+2.)*e->eVars.n[2];
		Mfield[1] = (-2.)*e->eVars.n[3];
		Mfield[2] = e->eVars.n[0] - e->eVars.n[1];
		//Absolute magnetization:
		double Mabs = integral(sqrt(lengthSquared(Mfield)));
		//Total magentziation and direction
		vector3<> Mtot = e->gInfo.dV * sumComponents(Mfield);
		vector3<> euler;
		if(Mtot.length()) getEulerAxis(Mtot, euler);
		euler *= 180./M_PI; //convert to degrees
		logPrintf("  magneticMoment: [ Abs: %7.5f  Tot: %7.5f  theta: %6.2f  phi: %+7.2f ]", Mabs, Mtot.length(), euler[1], euler[0]);
	}
	logPrintf("\n"); logFlush();
}


// Fermi legendre multipliers (TS and optionally muN)
void ElecInfo::updateFillingsEnergies(const std::vector<diagMatrix>& eps, Energies& ener) const
{
	//Determine relevant Legenedre multipliers:
	double Bz = 0.;
	double mu = !std::isnan(this->mu) ? this->mu : findMu(eps, nElectrons, Bz);
	
	//Calculate entropy contributions:
	ener.TS = 0.0;
	for(int q=qStart; q<qStop; q++)
		ener.TS += smearingWidth * qnums[q].weight * trace(smearEntropy(muEff(mu,Bz,q), eps[q]));
	mpiUtil->allReduce(ener.TS, MPIUtil::ReduceSum);

	//Grand canonical multiplier if fixed mu:
	if(!std::isnan(this->mu)) ener.muN = this->mu * nElectrons;
}

//-------------------- Fermi function utilities --------------------------

double ElecInfo::smear(double mu, double eps) const
{	double x = (eps-mu)/(2.*smearingWidth);
	switch(smearingType)
	{	case SmearingFermi: return 0.5*(1.-tanh(x));
		case SmearingGauss: return 0.5*erfc(x);
		case SmearingCold: return 0.5*erfc(x+sqrt(0.5)) + exp(-std::pow(x+sqrt(0.5),2))/sqrt(2*M_PI);
		default: return NAN; //to suppress warning
	}
}

double ElecInfo::smearPrime(double mu, double eps) const
{	double x = (eps-mu)/(2.*smearingWidth);
	switch(smearingType)
	{	case SmearingFermi: return -0.25/(smearingWidth * std::pow(cosh(x), 2));
		case SmearingGauss: return -exp(-x*x) / (2.*sqrt(M_PI)*smearingWidth);
		case SmearingCold: return -exp(-std::pow(x+sqrt(0.5),2)) * (2.+x*sqrt(2.)) / (2.*sqrt(M_PI)*smearingWidth);
		default: return NAN; //to suppress warning
	}
}

double ElecInfo::smearEntropy(double mu, double eps) const
{	double x = (eps-mu)/(2.*smearingWidth);
	switch(smearingType)
	{	case SmearingFermi: 
		{	double f = 0.5*(1.-tanh(x));
			double S = 0.;
			if(f>1e-300) S -= f*log(f);
			if(1-f>1e-300) S -= (1-f)*log(1-f);
			return S;
		}
		case SmearingGauss: return exp(-x*x) / sqrt(M_PI);
		case SmearingCold: return exp(-std::pow(x+sqrt(0.5),2)) * (1.+x*sqrt(2.)) / sqrt(M_PI);
		default: return NAN; //to suppress warning
	}
}

diagMatrix ElecInfo::smear(double mu, const diagMatrix& eps) const
{	diagMatrix ret(eps);
	for(unsigned i=0; i<eps.size(); i++) ret[i] = smear(mu, eps[i]);
	return ret;
}

diagMatrix ElecInfo::smearPrime(double mu, const diagMatrix& eps) const
{	diagMatrix ret(eps);
	for(unsigned i=0; i<eps.size(); i++) ret[i] = smearPrime(mu, eps[i]);
	return ret;
}

diagMatrix ElecInfo::smearEntropy(double mu, const diagMatrix& eps) const
{	diagMatrix ret(eps);
	for(unsigned i=0; i<eps.size(); i++) ret[i] = smearEntropy(mu, eps[i]);
	return ret;
}


matrix ElecInfo::smearGrad(double mu, const diagMatrix& eps, const matrix& gradF) const
{	matrix gradEps(gradF); //copy input
	complex* gradEpsData = gradEps.data();
	//calculate result in place:
	for(int i=0; i<gradF.nRows(); i++)
		for(int j=0; j<gradF.nCols(); j++)
		{	double deps = eps[i] - eps[j];
			gradEpsData[gradEps.index(i,j)]
				*= (fabs(deps)<1e-6) ? smearPrime(mu,eps[i]) : (smear(mu,eps[i])-smear(mu,eps[j]))/deps;
		}
	return gradEps;
}


//Number of electrons and magnetization in a fermi distribution of given mu, Bz and eigenvalues:
double ElecInfo::magnetizationCalc(double mu, double Bz, const std::vector<diagMatrix>& eps, double& N) const
{	N = 0.;
	double M = 0.;
	for(int q=qStart; q<qStop; q++)
	{	double s = qnums[q].spin;
		double muEff = this->muEff(mu, Bz, q);
		for(double epsCur: eps[q])
		{	double wf = qnums[q].weight * smear(muEff, epsCur);
			N += wf;
			M += s * wf;
		}
	}
	mpiUtil->allReduce(N, MPIUtil::ReduceSum, true);
	mpiUtil->allReduce(M, MPIUtil::ReduceSum, true);
	return M;
}

//Calculate nElectrons at given mu, bisecting on Bz if M is constrained
double ElecInfo::nElectronsCalc(double mu, const std::vector< diagMatrix >& eps, double& Bz) const
{	Bz = 0.;
	double N = 0.;
	if(Mconstrain)
	{	const bool& verbose = e->cntrl.shouldPrintMuSearch;
		if(verbose) logPrintf("\nBisecting to find Bz(M=%5lf)\n", Minitial);
		//Find a range which is known to bracket the result:
		const double absTol = 1e-10, relTol = 1e-14;
		double Mtol = std::max(absTol, relTol*fabs(Minitial));
		double BzMin=-0.1, BzMax=+0.1;
		while(magnetizationCalc(mu,BzMin,eps,N)>=Minitial+Mtol) BzMin-=(BzMax-BzMin);
		while(magnetizationCalc(mu,BzMax,eps,N)<=Minitial-Mtol) BzMax+=(BzMax-BzMin);
		//Bisect:
		double BzTol = std::max(absTol*smearingWidth, relTol*std::max(fabs(BzMin),fabs(BzMax)));
		while(BzMax-BzMin>=BzTol)
		{	double Bz = 0.5*(BzMin + BzMax);
			double M = magnetizationCalc(mu, Bz, eps, N);
			if(verbose) logPrintf("BzBISECT: Bz = [ %.15le %.15le %.15le ]  M = %le\n", BzMin, Bz, BzMax, M);
			if(M>Minitial) BzMax = Bz;
			else BzMin = Bz;
		}
		Bz = 0.5*(BzMin + BzMax);
	}
	magnetizationCalc(mu, Bz, eps, N);
	return N;
}


//Return the mu that would match Ntarget at the current eigenvalues (bisection method)
double ElecInfo::findMu(const std::vector<diagMatrix>& eps, double nElectrons, double& Bz) const
{	const bool& verbose = e->cntrl.shouldPrintMuSearch;
	if(verbose) logPrintf("\nBisecting to find mu(nElectrons=%.15le)\n", nElectrons);
	//Find a range which is known to bracket the result:
	const double absTol = 1e-10, relTol = 1e-14;
	double nTol = std::max(absTol, relTol*fabs(nElectrons));
	double muMin=-0.1, muMax=+0.0;
	while(nElectronsCalc(muMin,eps,Bz)>=nElectrons+nTol) muMin-=(muMax-muMin);
	while(nElectronsCalc(muMax,eps,Bz)<=nElectrons-nTol) muMax+=(muMax-muMin);
	//Bisect:
	double muTol = std::max(absTol*smearingWidth, relTol*std::max(fabs(muMin),fabs(muMax)));
	while(muMax-muMin>=muTol)
	{	double mu = 0.5*(muMin + muMax);
		double N = nElectronsCalc(mu, eps, Bz);
		if(verbose) logPrintf("MUBISECT: mu = [ %.15le %.15le %.15le ]  N = %le\n", muMin, mu, muMax, N);
		if(N>nElectrons) muMax = mu;
		else muMin = mu;
	}
	return 0.5*(muMin + muMax);
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
	std::vector<QuantumNumber> qRed = e->symm.reduceKmesh(qnums);
	if(qnums.size()==qRed.size())
		logPrintf("No reducable k-points. ");
	else
		logPrintf("Reduced to %lu k-points under symmetry. ", qRed.size());
	qnums.assign(qRed.begin(), qRed.end());
	//Include spin weight factor in k-point weights:
	for(unsigned i=0; i<qnums.size(); i++)
		qnums[i].weight *= spinWeight;
	//Add extra states due to spin-degrees of freedom:
	if(spinType==SpinZ) //(Note: noncollinear modes use spinor ColumnBundles and only one state per k-point)
	{	unsigned nkPoints = qnums.size();
		qnums.insert(qnums.end(), qRed.begin(), qRed.end()); //add second copy for other spin
		for(unsigned ik=0; ik<nkPoints; ik++)
		{	//Set spin values:
			qnums[ik].spin = +1;
			qnums[ik+nkPoints].spin = -1;
		}
	}
	//Output the reduced kpoints:
	if(e->cntrl.shouldPrintKpointsBasis)
	{	logPrintf("States including spin/spin-weights:\n");
		kpointsPrint(globalLog, true);
	}
	else logPrintf("\n");
}

void ElecInfo::kpointsPrint(FILE* fp, bool printSpin) const
{	for(unsigned q=0; q<qnums.size(); q++)
	{	kpointPrint(fp, q, printSpin); fprintf(fp, "\n");
	}
}

void ElecInfo::kpointPrint(FILE* fp, int q, bool printSpin) const
{	fprintf(fp, "%5d  [ %+.7f %+.7f %+.7f ]  %.9f", q, qnums[q].k[0], qnums[q].k[1], qnums[q].k[2], qnums[q].weight);
	if(printSpin && qnums[q].spin) fprintf(fp, "  spin %+d", qnums[q].spin); //only for spin-polarized calculations
}


//-------------- diagMatrix/matrix array parallel I/O ------------------

void ElecInfo::read(std::vector<diagMatrix>& M, const char *fname, int nRowsOverride) const
{	int nRows = nRowsOverride ? nRowsOverride : nBands;
	M.resize(nStates);
	MPIUtil::File fp; mpiUtil->fopenRead(fp, fname, nStates*nRows*sizeof(double));
	mpiUtil->fseek(fp, qStart*nRows*sizeof(double), SEEK_SET);
	for(int q=qStart; q<qStop; q++)
	{	M[q].resize(nRows);
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

void ElecInfo::appendWrite(const std::vector<diagMatrix>& M, const char *fname, int nRowsOverride) const
{	int nRows = nRowsOverride ? nRowsOverride : nBands;
	assert(int(M.size())==nStates);
	MPIUtil::File fp; mpiUtil->fopenAppend(fp, fname);
	mpiUtil->fseek(fp, qStart*nRows*sizeof(double), SEEK_CUR);
	for(int q=qStart; q<qStop; q++)
	{	assert(M[q].nRows()==nRows);
		mpiUtil->fwrite(M[q].data(), sizeof(double), nRows, fp);
	}
	mpiUtil->fclose(fp);
}