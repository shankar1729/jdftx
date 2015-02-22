/*-------------------------------------------------------------------
Copyright 2014 Deniz Gunceler, Ravishankar Sundararaman

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

#include <electronic/Dump.h>
#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/operators.h>
#include <electronic/Dump_internal.h>
#include <core/WignerSeitz.h>
#include <core/Operators.h>
#include <core/Units.h>

//---------------------- Excitations -----------------------------------------------

void dumpExcitations(const Everything& e, const char* filename)
{
	const GridInfo& g = e.gInfo;

	struct excitation
	{	int q,o,u;
		double dE;
		double dreal, dimag, dnorm;
		
		excitation(int q, int o, int u, double dE, double dreal, double dimag, double dnorm): q(q), o(o), u(u), dE(dE), dreal(dreal), dimag(dimag), dnorm(dnorm){};
		
		inline bool operator<(const excitation& other) const {return dE<other.dE;}
		void print(FILE* fp) const { fprintf(fp, "%5i %3i %3i %12.5e %12.5e %12.5e %12.5e\n", q, o, u, dE, dreal, dimag, dnorm); }
	};
	std::vector<excitation> excitations;

	double maxHOMO=-DBL_MAX, minLUMO=DBL_MAX; // maximum (minimum) of all HOMOs (LUMOs) in all qnums
	int maxHOMOq=0, minLUMOq=0, maxHOMOn=0, minLUMOn=0; //Indices and energies for the indirect gap
	
	//Select relevant eigenvals:
	std::vector<diagMatrix> eigsQP;
	if(e.exCorr.orbitalDep && e.dump.count(std::make_pair(DumpFreq_End, DumpOrbitalDep)))
	{	//Search for an eigenvalsQP file:
		string fname = e.dump.getFilename("eigenvalsQP");
		FILE* fp = fopen(fname.c_str(), "r");
		if(fp)
		{	fclose(fp);
			eigsQP.resize(e.eInfo.nStates);
			e.eInfo.read(eigsQP, fname.c_str());
		}
	}
	const std::vector<diagMatrix>& eigs = eigsQP.size() ? eigsQP : e.eVars.Hsub_eigs;
	
	// Integral kernel's for Fermi's golden rule
	ScalarField r0, r1, r2;
	nullToZero(r0, g); 	nullToZero(r1, g); 	nullToZero(r2, g);
	applyFunc_r(g, Moments::rn_pow_x, 0, g.R, 1, vector3<>(0.,0.,0.), r0->data());
	applyFunc_r(g, Moments::rn_pow_x, 1, g.R, 1, vector3<>(0.,0.,0.), r1->data());
	applyFunc_r(g, Moments::rn_pow_x, 2, g.R, 1, vector3<>(0.,0.,0.), r2->data());
	
	//Find and cache all excitations in system (between same qnums)
	bool insufficientBands = false;
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{	//Find local HOMO and check band sufficiency:
		int HOMO = e.eInfo.findHOMO(q);
		if(HOMO+1>=e.eInfo.nBands) { insufficientBands=true; break; }
		
		//Update global HOMO and LUMO of current process:
		if(eigs[q][HOMO]   > maxHOMO) { maxHOMOq = q; maxHOMOn = HOMO;   maxHOMO = eigs[q][HOMO];   }
		if(eigs[q][HOMO+1] < minLUMO) { minLUMOq = q; minLUMOn = HOMO+1; minLUMO = eigs[q][HOMO+1]; }
		
		for(int o=HOMO; o>=0; o--)
		{	for(int u=(HOMO+1); u<e.eInfo.nBands; u++)
			{	complex x = integral(I(e.eVars.C[q].getColumn(u,0))*r0*I(e.eVars.C[q].getColumn(o,0)));
				complex y = integral(I(e.eVars.C[q].getColumn(u,0))*r1*I(e.eVars.C[q].getColumn(o,0)));
				complex z = integral(I(e.eVars.C[q].getColumn(u,0))*r2*I(e.eVars.C[q].getColumn(o,0)));
				vector3<> dreal(x.real(), y.real(),z.real());
				vector3<> dimag(x.imag(), y.imag(),z.imag());
				vector3<> dnorm(sqrt(x.norm()), sqrt(y.norm()),sqrt(z.norm()));
				double dE = eigs[q][u]-eigs[q][o]; //Excitation energy
				excitations.push_back(excitation(q, o, u, dE, dreal.length_squared(), dimag.length_squared(), dnorm.length_squared()));
			}
		}
	}
	mpiUtil->allReduce(insufficientBands, MPIUtil::ReduceLOr);
	if(insufficientBands)
	{	logPrintf("Insufficient bands to calculate excited states!\n");
		logPrintf("Increase the number of bands (elec-n-bands) and try again!\n");
		return;
	}
	
	//Transmit results to head process:
	if(mpiUtil->isHead())
	{	excitations.reserve(excitations.size() * mpiUtil->nProcesses());
		for(int jProcess=1; jProcess<mpiUtil->nProcesses(); jProcess++)
		{	//Receive data:
			size_t nExcitations; mpiUtil->recv(nExcitations, jProcess, 0);
			std::vector<int> msgInt(4 + nExcitations*3); 
			std::vector<double> msgDbl(2 + nExcitations*4);
			mpiUtil->recv(msgInt.data(), msgInt.size(), jProcess, 1);
			mpiUtil->recv(msgDbl.data(), msgDbl.size(), jProcess, 2);
			//Unpack:
			std::vector<int>::const_iterator intPtr = msgInt.begin();
			std::vector<double>::const_iterator dblPtr = msgDbl.begin();
			//--- globals:
			int j_maxHOMOq = *(intPtr++); int j_maxHOMOn = *(intPtr++); double j_maxHOMO = *(dblPtr++);
			int j_minLUMOq = *(intPtr++); int j_minLUMOn = *(intPtr++); double j_minLUMO = *(dblPtr++);
			if(j_maxHOMO > maxHOMO) { maxHOMOq=j_maxHOMOq; maxHOMOn=j_maxHOMOn; maxHOMO=j_maxHOMO; }
			if(j_minLUMO < minLUMO) { minLUMOq=j_minLUMOq; minLUMOn=j_minLUMOn; minLUMO=j_minLUMO; }
			//--- excitation array:
			for(size_t iExcitation=0; iExcitation<nExcitations; iExcitation++)
			{	int q = *(intPtr++); int o = *(intPtr++); int u = *(intPtr++);
				double dE = *(dblPtr++);
				double dreal = *(dblPtr++); double dimag = *(dblPtr++); double dnorm = *(dblPtr++);
				excitations.push_back(excitation(q, o, u, dE, dreal, dimag, dnorm));
			}
		}
	}
	else
	{	//Pack data:
		std::vector<int> msgInt; std::vector<double> msgDbl;
		size_t nExcitations = excitations.size();
		msgInt.reserve(4 + nExcitations*3);
		msgDbl.reserve(2 + nExcitations*4);
		msgInt.push_back(maxHOMOq); msgInt.push_back(maxHOMOn); msgDbl.push_back(maxHOMO);
		msgInt.push_back(minLUMOq); msgInt.push_back(minLUMOn); msgDbl.push_back(minLUMO);
		for(const excitation& e: excitations)
		{	msgInt.push_back(e.q); msgInt.push_back(e.o); msgInt.push_back(e.u);
			msgDbl.push_back(e.dE);
			msgDbl.push_back(e.dreal); msgDbl.push_back(e.dimag); msgDbl.push_back(e.dnorm);
		}
		//Send data:
		mpiUtil->send(nExcitations, 0, 0);
		mpiUtil->send(msgInt.data(), msgInt.size(), 0, 1);
		mpiUtil->send(msgDbl.data(), msgDbl.size(), 0, 2);
	}

	//Process and print excitations:
	if(!mpiUtil->isHead()) return;
	
	FILE* fp = fopen(filename, "w");
	if(!fp) die("Error opening %s for writing.\n", filename);
	
	std::sort(excitations.begin(), excitations.end());
	const excitation& opt = excitations.front();
	fprintf(fp, "Using %s eigenvalues.      HOMO: %.5f   LUMO: %.5f  \n", eigsQP.size() ? "discontinuity-corrected QP" : "KS", maxHOMO, minLUMO);
	fprintf(fp, "Optical (direct) gap: %.5e (from n = %i to %i in qnum = %i)\n", opt.dE, opt.o, opt.u, opt.q);
	fprintf(fp, "Indirect gap: %.5e (from (%i, %i) to (%i, %i))\n\n", minLUMO-maxHOMO, maxHOMOq, maxHOMOn, minLUMOq, minLUMOn);
	
	fprintf(fp, "Optical excitation energies and corresponding electric dipole transition strengths\n");
	fprintf(fp, "qnum   i   f      dE        |<psi1|r|psi2>|^2 (real, imag, norm)\n");
	for(const excitation& e: excitations) e.print(fp);
	fclose(fp);
}


//---------------------------- Moments ------------------------------------

namespace Moments
{
	inline double map_to_interval(double x)
	{	double temp =  modf(x, new double);
		if(temp>0.5)
			return temp-1.;
		else if(temp<-0.5)
			return temp+1.;
		else
			return temp;
	}
	inline vector3<> map_to_interval(vector3<> x)
	{	for(int j=0;j<3;j++)
			x[j] = map_to_interval(x[j]);
		return x;
	}

	void rn_pow_x(int i, vector3<> r, int dir, matrix3<> R, double moment, vector3<> r0, double* rx)
	{	vector3<> lx = map_to_interval(inv(R)*(r-r0));
		rx[i] = pow((R*lx)[dir], moment);
	}
	
	void dumpMoment(const Everything& e, const char* filename, int moment, vector3<> origin)
	{	const GridInfo& g = e.gInfo;
		
		FILE* fp = fopen(filename, "w");
		if(!fp) die("Error opening %s for writing.\n", filename);	
		
		// Calculates the electronic moment about origin
		ScalarField r0, r1, r2;
		nullToZero(r0, g); 	nullToZero(r1, g); 	nullToZero(r2, g);
		applyFunc_r(g, rn_pow_x, 0, g.R, moment, origin, r0->data());
		applyFunc_r(g, rn_pow_x, 1, g.R, moment, origin, r1->data());
		applyFunc_r(g, rn_pow_x, 2, g.R, moment, origin, r2->data());
		vector3<> elecMoment;
		elecMoment[0] = integral(e.eVars.n[0]*r0);
		elecMoment[1] = integral(e.eVars.n[0]*r1);
		elecMoment[2] = integral(e.eVars.n[0]*r2);
		fprintf(fp, "Electron moment of order %i: %f\t%f\t%f", moment, elecMoment[0], elecMoment[1], elecMoment[2]);
		
		// Calculates the ionic moment about the origin
		vector3<> ionMoment(0., 0., 0.);
		for(auto sp: e.iInfo.species)
			for(unsigned n=0; n < sp->atpos.size(); n++)
			{	vector3<> cartesianCoord = g.R * map_to_interval(sp->atpos[n]);
				for(int j=0; j<3; j++)
					ionMoment[j] += -sp->Z * pow(cartesianCoord[j], moment); 
			}
		fprintf(fp, "\nIon moment of order %i: %f\t%f\t%f", moment, ionMoment[0], ionMoment[1], ionMoment[2]);

		
		// Calculates the total (elec+ion) dipole moment
		fprintf(fp, "\nTotal moment of order %i: %f\t%f\t%f", moment, ionMoment[0]+elecMoment[0], ionMoment[1]+elecMoment[1], ionMoment[2]+elecMoment[2]);		
	}
}


//---------------------------- XC analysis --------------------------------

namespace XC_Analysis
{
	static const double cutoff = 1e-8;
	
	void regularize(int i, vector3<> r, double* tau)
	{	tau[i] = (tau[i]<cutoff ? 0. : tau[i]);
	}
	void spness_kernel(int i, vector3<> r, double* tauW, double* tau, double* spness)
	{	spness[i] = (tau[i]>cutoff ? tauW[i]/tau[i] : 1.);
	}
	
	ScalarFieldArray tauWeizsacker(const Everything& e)
	{	ScalarFieldArray tauW(e.eVars.n.size());
		for(size_t j=0; j<e.eVars.n.size(); j++)
			tauW[j] = (1./8.)*lengthSquared(gradient(e.eVars.n[j]))*pow(e.eVars.n[j], -1);
		return tauW;
	}
	
	ScalarFieldArray spness(const Everything& e)
	{
			const auto& tau = (e.exCorr.needsKEdensity() ? e.eVars.tau : e.eVars.KEdensity());
		
			ScalarFieldArray tauW = tauWeizsacker(e);
			ScalarFieldArray spness(e.eVars.n.size()); nullToZero(spness, e.gInfo);
			for(size_t j=0; j<e.eVars.n.size(); j++)
			{	applyFunc_r(e.gInfo, regularize, tau[j]->data());
				applyFunc_r(e.gInfo, regularize, tauW[j]->data());
				//spness[j] = tauW[j]*pow(tau[j], -1);
				applyFunc_r(e.gInfo, spness_kernel, tauW[j]->data(), tau[j]->data(), spness[j]->data());
			}

			return spness;
	}
	
	ScalarFieldArray sHartree(const Everything& e)
	{
		ScalarFieldArray sVh(e.eVars.n.size());
		ScalarFieldTildeArray nTilde = J(e.eVars.n);
		for(size_t s=0; s<e.eVars.n.size(); s++)
			sVh[s] = I((*(e.coulomb))(nTilde[s]));
		
		return sVh;
	}
}


//-------------------------- Solvation radii ----------------------------------

inline void set_rInv(size_t iStart, size_t iStop, const vector3<int>& S, const matrix3<>& RTR, const WignerSeitz* ws, double* rInv)
{	vector3<> invS; for(int k=0; k<3; k++) invS[k] = 1./S[k];
	matrix3<> meshMetric = Diag(invS) * RTR * Diag(invS);
	const double rWidth = 1.;
	THREAD_rLoop( 
		double r = sqrt(meshMetric.metric_length_squared(ws->restrict(iv, S, invS)));
		if(r < rWidth) //Polynomial with f',f" zero at origin and f,f',f" matched at rWidth
		{	double t = r/rWidth;
			rInv[i] = (1./rWidth) * (1. + 0.5*(1.+t*t*t*(-2.+t)) + (1./12)*(1.+t*t*t*(-4.+t*3.))*2. );
		}
		else rInv[i] = 1./r;
	)
}

void Dump::dumpRsol(ScalarField nbound, string fname)
{	
	//Compute normalization factor for the partition:
	int nAtomsTot = 0; for(const auto& sp: e->iInfo.species) nAtomsTot += sp->atpos.size();
	const double nFloor = 1e-5/nAtomsTot; //lower cap on densities to prevent Nyquist noise in low density regions
	ScalarField nAtomicTot;
	for(const auto& sp: e->iInfo.species)
	{	RadialFunctionG nRadial;
		logSuspend(); sp->getAtom_nRadial(0,0, nRadial, true); logResume();
		for(unsigned atom=0; atom<sp->atpos.size(); atom++)
		{	ScalarField nAtomic = radialFunction(e->gInfo, nRadial, sp->atpos[atom]);
			double nMin, nMax; callPref(eblas_capMinMax)(e->gInfo.nr, nAtomic->dataPref(), nMin, nMax, nFloor);
			nAtomicTot += nAtomic;
		}
	}
	ScalarField nboundByAtomic = (nbound*nbound) * inv(nAtomicTot);

	ScalarField rInv0(ScalarFieldData::alloc(e->gInfo));
	{	logSuspend(); WignerSeitz ws(e->gInfo.R); logResume();
		threadLaunch(set_rInv, e->gInfo.nr, e->gInfo.S, e->gInfo.RTR, &ws, rInv0->data());
	}
	
	//Compute bound charge 1/r and 1/r^2 expectation values weighted by atom-density partition:
	FILE* fp = fopen(fname.c_str(), "w");
	fprintf(fp, "#Species   rMean +/- rSigma [bohrs]   (rMean +/- rSigma [Angstrom])   sqrt(Int|nbound^2|) in partition\n");
	for(const auto& sp: e->iInfo.species)
	{	RadialFunctionG nRadial;
		logSuspend(); sp->getAtom_nRadial(0,0, nRadial, true); logResume();
		for(unsigned atom=0; atom<sp->atpos.size(); atom++)
		{	ScalarField w = radialFunction(e->gInfo, nRadial, sp->atpos[atom]) * nboundByAtomic;
			//Get r centered at current atom:
			ScalarFieldTilde trans; nullToZero(trans, e->gInfo); initTranslation(trans, e->gInfo.R*sp->atpos[atom]);
			ScalarField rInv = I(trans * J(rInv0), true);
			//Compute moments:
			double wNorm = integral(w);
			double rInvMean = integral(w * rInv) / wNorm;
			double rInvSqMean = integral(w * rInv * rInv) / wNorm;
			double rInvSigma = sqrt(rInvSqMean - rInvMean*rInvMean);
			double rMean = 1./rInvMean;
			double rSigma = rInvSigma / (rInvMean*rInvMean);
			//Print stats:
			fprintf(fp, "Rsol %s    %.2lf +/- %.2lf    ( %.2lf +/- %.2lf A )   Qrms: %.1le\n", sp->name.c_str(),
				rMean, rSigma, rMean/Angstrom, rSigma/Angstrom, sqrt(wNorm));
		}
	}
	fclose(fp);
}

