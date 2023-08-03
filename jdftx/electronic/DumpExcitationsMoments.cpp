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
#include <electronic/ColumnBundleTransform.h>
#include <electronic/Dump_internal.h>
#include <core/WignerSeitz.h>
#include <core/Operators.h>
#include <core/Units.h>
#include <core/LatticeUtils.h>

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
	VectorField r;
	{	nullToZero(r, g);
		logSuspend();
		WignerSeitz ws(e.gInfo.R);
		logResume();
		const vector3<>& x0 = e.coulombParams.embedCenter; //origin
		
		vector3<double*> rData(r.data());
		matrix3<> invS = inv(Diag(vector3<>(e.gInfo.S)));
		vector3<int> iv;
		size_t i = 0;
		for(iv[0]=0; iv[0]<e.gInfo.S[0]; iv[0]++)
		for(iv[1]=0; iv[1]<e.gInfo.S[1]; iv[1]++)
		for(iv[2]=0; iv[2]<e.gInfo.S[2]; iv[2]++)
		{	vector3<> x = x0 + ws.reduce(invS*iv - x0); //lattice coordinates wrapped to WS centered on x0
			storeVector(e.gInfo.R*x, rData, i++);
		}
	}
	
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
			{	vector3<> dreal, dimag, dnorm;
				for(int iDir=0; iDir<3; iDir++)
				{	complex xi = integral(I(e.eVars.C[q].getColumn(u,0)) * r[iDir] * I(e.eVars.C[q].getColumn(o,0)));
					dreal[iDir] = xi.real();
					dimag[iDir] = xi.imag();
					dnorm[iDir] = xi.abs();
				}
				double dE = eigs[q][u]-eigs[q][o]; //Excitation energy
				excitations.push_back(excitation(q, o, u, dE, dreal.length_squared(), dimag.length_squared(), dnorm.length_squared()));
			}
		}
	}
	mpiWorld->allReduce(insufficientBands, MPIUtil::ReduceLOr);
	if(insufficientBands)
	{	logPrintf("Insufficient bands to calculate excited states!\n");
		logPrintf("Increase the number of bands (elec-n-bands) and try again!\n");
		return;
	}
	
	//Transmit results to head process:
	if(mpiWorld->isHead())
	{	excitations.reserve(excitations.size() * mpiWorld->nProcesses());
		for(int jProcess=1; jProcess<mpiWorld->nProcesses(); jProcess++)
		{	//Receive data:
			size_t nExcitations; mpiWorld->recv(nExcitations, jProcess, 0);
			std::vector<int> msgInt(4 + nExcitations*3); 
			std::vector<double> msgDbl(2 + nExcitations*4);
			mpiWorld->recvData(msgInt, jProcess, 1);
			mpiWorld->recvData(msgDbl, jProcess, 2);
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
		mpiWorld->send(nExcitations, 0, 0);
		mpiWorld->sendData(msgInt, 0, 1);
		mpiWorld->sendData(msgDbl, 0, 2);
	}

	//Process and print excitations:
	if(!mpiWorld->isHead()) return;
	
	FILE* fp = fopen(filename, "w");
	if(!fp) die_alone("Error opening %s for writing.\n", filename);
	
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

//---------------------------- FCI ------------------------------------

//Split doublet index i12 that ends on N^2 or N(N+1)/2 depending on sym = false or true, into i1 and i2
inline void splitDoubletIndex(size_t i12, size_t N, bool sym, size_t& i1, size_t& i2)
{	if(sym) //N(N+1)/2 case
	{	i1 = 0; while((i1+1)*(i1+2)<=2*i12) i1++;
		i2 = i12 - (i1*(i1+1))/2;
	}
	else //N^2 case
	{	i1 = i12 / N;
		i2 = i12 - i1 * N;
	}
}

void dumpFCI(const Everything& e, const char* filename)
{	static StopWatch watch("dumpFCI"); watch.start();
	
	const size_t nBands = e.eInfo.nBands;
	const int nStates = e.eInfo.nStates;
	const int nOrbitals = nBands * nStates;
	const int nSpinor = e.eInfo.spinorLength();
	const double Kcut = 1e-8; //threshold for non-zero matrix element
	
	FILE* fp = NULL;
	if(mpiWorld->isHead())
	{	fp = fopen(filename, "w");
		if(!fp) die_alone("Error opening %s for writing.\n", filename);
		
		//Write header:
		fprintf(fp, " &FCI\n  NORB=%d,\n  NELEC=%d,\n  MS2=%d,\n  UHF=TRUE,\n  ORBSYM= ",
			nOrbitals, int(round(e.eInfo.nElectrons)), int(round(fabs(e.eInfo.Minitial))) );
		for(int i=0; i<nOrbitals; i++)
			fprintf(fp, "0, ");
		fprintf(fp, "\n  ISYM=1,\n  NPROP= 1 1 1,\n  PROPBITLEN= 15\n &END\n");
	}
	
	//Coulomb integrals:
	#define GETwfns(a) \
		ColumnBundle Ctmp##a; \
		if(not e.eInfo.isMine(q##a)) \
			Ctmp##a.init(nBands, e.basis[q##a].nbasis*nSpinor, &(e.basis[q##a]), &(e.eInfo.qnums[q##a]), isGpuEnabled()); \
		const ColumnBundle& C##a = e.eInfo.isMine(q##a) ? e.eVars.C[q##a] : Ctmp##a; \
		mpiWorld->bcastData((ColumnBundle&)C##a, e.eInfo.whose(q##a));
	#define UPDATE_conjIpsi(a) \
		for(int s=0; s<nSpinor; s++) \
			conjIpsi##a[s] = conj(I(C##a.getColumn(i##a,s)));
	#define UPDATE_nPair(a,b, dk) \
		{	complexScalarField nPair; \
			for(int s=0; s<nSpinor; s++) \
				nPair += conjIpsi##a[s] * I(C##b.getColumn(i##b,s)); \
			if((dk).length_squared() > symmThresholdSq) \
				multiplyBlochPhase(nPair, dk); /*correct for k-wrapping upto reciprocal lattice vector*/ \
			n##a##b = J(nPair); \
		}
	#define UPDATE_12 \
		{	complexScalarFieldTilde n12; \
			UPDATE_nPair(1,2, vector3<>()) \
			Kn12 = O((*e.coulombWfns)(n12, dk21, 0.)); \
		}
	#define UPDATE_34 \
		UPDATE_nPair(3,4, dk43-dk21)
	
	//--- Quadruple loop over k-spin indices such that 12 and 34 have same spin indices
	for(int q1=0; q1<nStates; q1++)
	{	GETwfns(1)
		for(int q2=0; q2<=q1; q2++) if(e.eInfo.qnums[q1].spin == e.eInfo.qnums[q2].spin)
		{	GETwfns(2)
			vector3<> dk21 = e.eInfo.qnums[q2].k - e.eInfo.qnums[q1].k;
			for(int q3=0; q3<nStates; q3++)
			{	GETwfns(3)
				for(int q4=0; q4<=q3; q4++) if(e.eInfo.qnums[q3].spin == e.eInfo.qnums[q4].spin)
				{	vector3<> dk43 = e.eInfo.qnums[q4].k - e.eInfo.qnums[q3].k;
					if(circDistanceSquared(dk21, dk43) > symmThresholdSq) continue; //momentum conservation
					GETwfns(4)
					
					//Split band quadruplets over MPI:
					size_t N12 = (q1==q2 ? (nBands*(nBands+1))/2 : nBands*nBands);
					size_t N34 = (q3==q4 ? (nBands*(nBands+1))/2 : nBands*nBands);
					size_t N1234 = N12 * N34;
					size_t i1234start=0, i1234stop=0;
					TaskDivision(N1234, mpiWorld).myRange(i1234start, i1234stop);
					
					//Initial indices on this process:
					size_t i1234 = i1234start; //quaruplet index
					size_t i12 = i1234 / N34; //12 doublet index
					size_t i34 = i1234 - i12*N34; //34 doublet index
					size_t i1, i2, i3, i4; //individual band indices
					splitDoubletIndex(i12, nBands, q1==q2, i1, i2);
					splitDoubletIndex(i34, nBands, q3==q4, i3, i4);
					
					//Initialize wavefunctions, pair densities for initial indices:
					std::vector<complexScalarField> conjIpsi1(nSpinor), conjIpsi3(nSpinor);
					complexScalarFieldTilde Kn12, n34;
					//--- 12 pair
					UPDATE_conjIpsi(1)
					UPDATE_12
					//--- 34 pair
					UPDATE_conjIpsi(3)
					UPDATE_34
					
					//Loop over band quadruplets for this process:
					ostringstream oss; //buffer containing output from this process
					while(i1234 < i1234stop)
					{	
						complex K1234 = dot(Kn12, n34);
						if(K1234.abs() > Kcut)
						{	if(fabs(K1234.imag()) > Kcut)
								oss << "( " << K1234.real() << ", " << K1234.imag() << ") "; //complex element
							else
								oss << K1234.real() << ' '; //real element
							oss << i1*nStates+q1+1 << ' '
								<< i2*nStates+q2+1 << ' '
								<< i3*nStates+q3+1 << ' '
								<< i4*nStates+q4+1 << '\n';
						}
						
						i1234++; if(i1234==i1234stop) break;
						i4++;
						if(i4==(q3==q4 ? i3+1 : nBands))
						{	i4=0;
							i3++;
							if(i3==nBands)
							{	i3=0;
								i2++;
								if(i2==(q1==q2 ? i1+1 : nBands))
								{	i2=0;
									i1++;
									UPDATE_conjIpsi(1)
								}
								UPDATE_12
							}
							UPDATE_conjIpsi(3)
						}
						UPDATE_34
					}
					
					//Write matrix elements from HEAD:
					if(mpiWorld->isHead())
					{	fputs(oss.str().c_str(), fp); //local data
						oss.clear();
						for(int jProcess=1; jProcess<mpiWorld->nProcesses(); jProcess++)
						{	string buf; mpiWorld->recv(buf, jProcess, 0);
							fputs(buf.c_str(), fp); //data from other processes
						}
						fflush(fp);
					}
					else mpiWorld->send(oss.str(), 0, 0);
				}
			}
		}
	}
	#undef GETwfns
	#undef UPDATE_conjIpsi
	#undef UPDATE_nPair
	#undef UPDATE_12
	#undef UPDATE_34
	
	//Compute and collect one-particle matrix elements / eigenvalues:
	std::vector<matrix> H0(nStates);
	std::vector<diagMatrix> E(nStates);
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{
		const ColumnBundle& Cq = e.eVars.C[q];
		ScalarFieldArray Vloc(e.eInfo.nDensities); nullToZero(Vloc, e.gInfo);
		for(unsigned s=0; s<Vloc.size(); s++)
			if(s<2)
			{	Vloc[s] = Jdag(O(e.iInfo.Vlocps));
				if(e.eVars.Vexternal.size())
					Vloc[s] += JdagOJ(e.eVars.Vexternal[s]);
			}
		ColumnBundle H0Cq = Idag_DiagV_I(Cq, Vloc) - 0.5*L(Cq); //local PS and KE
		{	//Add nonlocal PS contribution:
			std::vector<matrix> HVdagCq(e.iInfo.species.size());
			e.iInfo.EnlAndGrad(e.eInfo.qnums[q], e.eVars.F[q], e.eVars.VdagC[q], HVdagCq);
			e.iInfo.projectGrad(HVdagCq, Cq, H0Cq); //HCq += (nonlocalPS projectors) * Cq
		}
		H0[q] = Cq ^ H0Cq; //matrix elements of KE + net pseudopotential
		E[q] = e.eVars.Hsub_eigs[q];
		
		if(not mpiWorld->isHead())
		{	mpiWorld->sendData(E[q], 0, q);
			mpiWorld->sendData(H0[q], 0, q);
		}
	}
	if(mpiWorld->isHead())
	{	for(int q=0; q<nStates; q++)
			if(!e.eInfo.isMine(q))
			{	H0[q].init(nBands, nBands);
				E[q].resize(nBands);
				mpiWorld->recvData(E[q], e.eInfo.whose(q), q);
				mpiWorld->recvData(H0[q], e.eInfo.whose(q), q);
			}
	}
	
	//Output one-particle matrix elements / eigenvalues:
	if(mpiWorld->isHead())
	{	//One-particle
		for(int q=0; q<nStates; q++)
			for(size_t b1=0; b1<nBands; b1++)
				for(size_t b2=0; b2<=b1; b2++)
				{	complex M = H0[q](b1,b2);
					if(M.abs() > Kcut)
					{	if(fabs(M.imag()) > Kcut)
							fprintf(fp, "( %lg, %lg )", M.real(), M.imag()); //complex element
						else
							fprintf(fp, "%lg", M.real()); //real element
						fprintf(fp, " %lu %lu 0 0\n", b1*nStates+q+1, b2*nStates+q+1);
					}
				}
		//Eigenvalues
		for(int q=0; q<nStates; q++)
			for(size_t b=0; b<nBands; b++)
				fprintf(fp, "%lg %lu 0 0 0\n", E[q][b], b*nStates+q+1);
		//Ewald energy
		double Eewald = e.ener.E["Eewald"];
		if(fabs(Eewald) > Kcut)
			fprintf(fp, "%lg 0 0 0 0\n", Eewald);
		fclose(fp);
	}
	watch.stop();
}

//---------------------------- Moments ------------------------------------

void dumpMoment(const Everything& e, const char* filename)
{	logSuspend();
	WignerSeitz ws(e.gInfo.R);
	logResume();
	const vector3<>& x0 = e.coulombParams.embedCenter; //origin
	
	//Compute electronic moment in lattice coordinates:
	vector3<> elecMoment;
	ScalarField n = e.eVars.get_nTot();
	double* nData = n->data();
	matrix3<> invS = inv(Diag(vector3<>(e.gInfo.S)));
	vector3<int> iv;
	for(iv[0]=0; iv[0]<e.gInfo.S[0]; iv[0]++)
	for(iv[1]=0; iv[1]<e.gInfo.S[1]; iv[1]++)
	for(iv[2]=0; iv[2]<e.gInfo.S[2]; iv[2]++)
	{	vector3<> x = x0 + ws.reduce(invS*iv - x0); //lattice coordinates wrapped to WS centered on x0
		elecMoment += x * (*(nData)++); //collect in lattice coordinates
	}
	elecMoment *= e.gInfo.dV; //integration weight
	
	//Compute ionic moment in lattice coordinates:
	vector3<> ionMoment;
	for(auto sp: e.iInfo.species)
		for(vector3<> pos: sp->atpos)
		{	vector3<> x = x0 + ws.reduce(pos - x0); //wrap atom position to WS centered on x0
			ionMoment -= sp->Z * x;
		}
	
	//Zero contributions in periodic directions:
	for(int iDir=0; iDir<3; iDir++)
		if(not e.coulombParams.isTruncated()[iDir])
		{	elecMoment[iDir] = 0.;
			ionMoment[iDir] = 0.;
		}
	
	//Convert to Cartesian and compute total:
	elecMoment = e.gInfo.R * elecMoment;
	ionMoment = e.gInfo.R * ionMoment;
	vector3<> totalMoment = elecMoment + ionMoment;
	
	//Output:
	FILE* fp = fopen(filename, "w");
	if(!fp) die("Error opening %s for writing.\n", filename);	
	fprintf(fp, "#Contribution %4s %12s %12s\n", "x", "y", "z");
	fprintf(fp, "%10s %12.6lf %12.6lf %12.6lf\n", "Electronic", elecMoment[0], elecMoment[1], elecMoment[2]);
	fprintf(fp, "%10s %12.6lf %12.6lf %12.6lf\n", "Ionic", ionMoment[0], ionMoment[1], ionMoment[2]);
	fprintf(fp, "%10s %12.6lf %12.6lf %12.6lf\n", "Total", totalMoment[0], totalMoment[1], totalMoment[2]);
	fclose(fp);
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
		double r = sqrt(meshMetric.metric_length_squared(ws->reduce(iv, S, invS)));
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
			ScalarField rInv = I(trans * J(rInv0));
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

void Dump::dumpUnfold()
{	const GridInfo& gInfo = e->gInfo;
	const ElecInfo& eInfo = e->eInfo;
	//Header
	//--- unit cell grid:
	const matrix3<int>& M = Munfold;
	matrix3<> invM = inv(matrix3<>(M));
	int nUnits = abs(det(M));
	GridInfo gInfoUnit;
	gInfoUnit.R = gInfo.R * invM;
	gInfoUnit.Gmax = sqrt(2*e->cntrl.Ecut); //Ecut = 0.5 Gmax^2
	logSuspend();
	gInfoUnit.initialize();
	logResume();
	//--- list of unit cell k-points for each supercell k:
	std::vector<std::vector<vector3<>>> kUnit(eInfo.nStates);
	for(int q=0; q<eInfo.nStates; q++)
	{	PeriodicLookup<vector3<>> plook(kUnit[q], gInfoUnit.GGT, nUnits);
		for(const vector3<int>& iG: e->basis[q].iGarr)
		{	vector3<> k = (eInfo.qnums[q].k + iG) * invM; //corresponding unit cell k
			for(int iDir=0; iDir<3; iDir++) k[iDir] -= ceil(-0.5+k[iDir]); //wrap to (-0.5,0.5]
			if(plook.find(k) == string::npos) //not encountered yet
			{	plook.addPoint(kUnit[q].size(), k); //update periodic look-up table
				kUnit[q].push_back(k);
				if(int(kUnit[q].size())==nUnits) break; //all unit cell k for this supercell k found
			}
		}
		std::sort(kUnit[q].begin(), kUnit[q].end());
	}
	//--- write header:
	string fname = getFilename("bandUnfoldHeader");
	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
	if(mpiWorld->isHead())
	{	FILE* fp = fopen(fname.c_str(), "w");
		for(int q=0; q<eInfo.nStates; q++)
		{	fprintf(fp, "Supercell kpoint ");
			eInfo.kpointPrint(fp, q, true);
			fprintf(fp, "\n");
			for(vector3<>& k:  kUnit[q])
				fprintf(fp, "%+.7f %+.7f %+.7f\n", k[0], k[1], k[2]);
			fprintf(fp, "\n");
		}
		fclose(fp);
	}
	logPrintf("done\n"); logFlush();
	//Weights of each unit cell k:
	fname = getFilename("bandUnfold");
	logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
	std::vector<diagMatrix> unitWeights(eInfo.nStates);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
	{	unitWeights[q].reserve(nUnits * eInfo.nBands);
		for(const vector3<>& k: kUnit[q])
		{	//Create basis and qnum for unit cell k:
			Basis basisUnit;
			logSuspend();
			basisUnit.setup(gInfoUnit, e->iInfo, e->cntrl.Ecut, k);
			logResume();
			QuantumNumber qnumUnit; qnumUnit.k = k;
			//Setup wave function transformation:
			int nSpinor = eInfo.spinorLength();
			ColumnBundleTransform cbt(k, basisUnit, eInfo.qnums[q].k, e->basis[q], nSpinor, SpaceGroupOp(), +1, M);
			//Extract unit cell wavefunctions:
			const ColumnBundle& C = e->eVars.C[q];
			ColumnBundle OC = O(C);
			ColumnBundle Cunit(1, basisUnit.nbasis*nSpinor, &basisUnit, &qnumUnit, isGpuEnabled());
			ColumnBundle OCunit = Cunit.similar();
			for(int b=0; b<eInfo.nBands; b++)
			{	Cunit.zero(); cbt.gatherAxpy(1., C,b, Cunit,0);
				OCunit.zero(); cbt.gatherAxpy(1., OC,b, OCunit,0);
				unitWeights[q].push_back(trace(Cunit^OCunit).real());
			}
		}
	}
	eInfo.write(unitWeights, fname.c_str(), nUnits*eInfo.nBands);
	logPrintf("done\n"); logFlush();
}


void dumpProjections(const Everything& e, const char* filename, bool ortho, bool norm)
{	const ElecVars& eVars = e.eVars;
	const ElecInfo& eInfo = e.eInfo;
	const IonInfo& iInfo = e.iInfo;
	FILE* fp = 0;
	if(mpiWorld->isHead())
	{	fp = fopen(filename, "w");
		//Write header:
		fprintf(fp, "%d states, %d bands, %d %sorbital-projections, %lu species\n",
			eInfo.nStates, eInfo.nBands, iInfo.nAtomicOrbitals(),
			(ortho ? "ortho-" : ""), iInfo.species.size());
		fprintf(fp, "# Symbol nAtoms nOrbitalsPerAtom lMax nShells(l=0) ... nShells(l=lMax)\n");
		for(const auto& sp: iInfo.species)
		{	int nAtoms = sp->atpos.size();
			int nOrbitalsPerAtom = sp->nAtomicOrbitals() / nAtoms;
			int lMax = sp->lMaxAtomicOrbitals();
			fprintf(fp, "%s %d %d %d", sp->name.c_str(), nAtoms, nOrbitalsPerAtom, lMax);
			for(int l=0; l<=lMax; l++)
				fprintf(fp, " %d", sp->nAtomicOrbitals(l));
			fprintf(fp, "\n");
		}
	}
	for(int q=0; q<eInfo.nStates; q++)
	{	matrix proj; //orbitals: nOrbitals x nBands
		if(eInfo.isMine(q))
		{	if(ortho)
			{	ColumnBundle psi = iInfo.getAtomicOrbitals(q, false);
				ColumnBundle Opsi = O(psi);
				matrix orthoMat = invsqrt(psi ^ Opsi); //orthonormalizing matrix
				proj = orthoMat * (Opsi ^ eVars.C[q]);
			}
			else proj = iInfo.getAtomicOrbitals(q, true) ^ eVars.C[q]; //dagger(Opsi).Cq
			if(not mpiWorld->isHead()) mpiWorld->sendData(proj, 0, q); //send to head for writing
		}
		if(mpiWorld->isHead())
		{	if(not eInfo.isMine(q)) //recv from process that stored q
			{	proj.init(iInfo.nAtomicOrbitals(), eInfo.nBands);
				mpiWorld->recvData(proj, eInfo.whose(q), q);
			}
			//Write projections:
			fprintf(fp, "# ");
			eInfo.kpointPrint(fp, q, true);
			fprintf(fp, "; lines per band, with %s per orbital:\n",
				norm ? "|projection|^2" : "projectionRe, projectionIm");
			complex* projData = proj.data();
			for(int b=0; b<proj.nCols(); b++) //bands
			{	for(int a=0; a<proj.nRows(); a++) //orbitals
				{	if(norm)
						fprintf(fp, "%9.7lf ", projData->norm()); //write |projection|^2
					else
						fprintf(fp, "%9.7lf %9.7lf ", projData->real(), projData->imag()); //write complex projection
					projData++;
				}
				fprintf(fp, "\n");
			}
		}
	}
	if(mpiWorld->isHead()) fclose(fp);
}
