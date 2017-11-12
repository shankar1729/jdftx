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

#include <electronic/ExCorr_OrbitalDep_GLLBsc.h>
#include <electronic/ColumnBundle.h>
#include <electronic/Everything.h>
#include <core/ScalarFieldIO.h>
#include <core/Units.h>

ExCorr_OrbitalDep_GLLBsc::ExCorr_OrbitalDep_GLLBsc(const Everything& e) : ExCorr::OrbitalDep(e)
{	smearingWidth = (e.eInfo.fillingsUpdate==ElecInfo::FillingsConst) ? 0. : e.eInfo.smearingWidth;
}


std::vector<double> ExCorr_OrbitalDep_GLLBsc::getExtremalEnergy(bool HOMO) const
{	int nSpins = e.eVars.n.size();
	if(smearingWidth) //smeared version
	{	std::vector<double> eExtremal(nSpins, 0.), wExtremal(nSpins, 0.);
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		{	int s = e.eInfo.qnums[q].index();
			for(int b=0; b<e.eInfo.nBands; b++)
			{	double f = HOMO ? e.eVars.F[q][b] : 1.-e.eVars.F[q][b];
				double w = e.eInfo.qnums[q].weight * (f<=0.5 ? 0. : (1.-f)*(f-0.5)*(f-0.5)); //a smooth function w(f) such that w(1/2) = w'(1/2) = w(1) = 0
				eExtremal[s] += w * e.eVars.Hsub_eigs[q][b];
				wExtremal[s] += w;
			}
		}
		mpiWorld->allReduce(eExtremal.data(), eExtremal.size(), MPIUtil::ReduceSum);
		mpiWorld->allReduce(wExtremal.data(), wExtremal.size(), MPIUtil::ReduceSum);
		for(int s=0; s<nSpins; s++) eExtremal[s] /= wExtremal[s];
		return eExtremal;
	}
	else //direct max/min version
	{	std::vector<double> eExtremal(nSpins, HOMO ? -DBL_MAX : DBL_MAX);
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
		{	int s = e.eInfo.qnums[q].index();
			for(int b=0; b<e.eInfo.nBands; b++)
			{	if(HOMO)
				{	if(e.eVars.F[q][b]>=0.5) eExtremal[s] = std::max(eExtremal[s], e.eVars.Hsub_eigs[q][b]);
				}
				else
				{	if(e.eVars.F[q][b]<=0.5) eExtremal[s] = std::min(eExtremal[s], e.eVars.Hsub_eigs[q][b]);
				}
			}
		}
		mpiWorld->allReduce(eExtremal.data(), eExtremal.size(), HOMO ? MPIUtil::ReduceMax : MPIUtil::ReduceMin);
		return eExtremal;
	}
}

ScalarFieldArray ExCorr_OrbitalDep_GLLBsc::getPotential() const
{	int nSpins = e.eVars.n.size();
	if(!e.eVars.Hsub_eigs[e.eInfo.qStart].size()) return ScalarFieldArray(nSpins); //no eigenvalues yet
	std::vector<double> eHOMO = getExtremalEnergy(true);
	return getPotential(eHOMO);
}

void ExCorr_OrbitalDep_GLLBsc::dump() const
{	int nSpins = e.eVars.n.size();
	if(!e.eVars.Hsub_eigs[e.eInfo.qStart].size()) return; //no eigenvalues yet
	
	//Detect HOMO and LUMO, and compute discontinuity potential:
	std::vector<double> eHOMO = getExtremalEnergy(true);
	std::vector<double> eLUMO = getExtremalEnergy(false);
	for(double e: eLUMO)
		if(e==+DBL_MAX)
		{	logPrintf("OrbitalDep dump: (GLLBsc) No unoccupied states in the calculation; can't determine discontinuity potential.\n");
			return;
		}
	ScalarFieldArray VsclocDisc(nSpins);
	if(e.cntrl.fixed_H)
	{	assert(e.eVars.VFilenamePattern.length());
		for(int s=0; s<nSpins; s++)
		{	string fname = e.eVars.VFilenamePattern;
			size_t pos = fname.find("$VAR");
			assert(pos != string::npos);
			fname.replace(pos,4, string("VsclocDisc") + (nSpins==1 ? "" : (s==0 ? "_up" : "_dn")));
			nullToZero(VsclocDisc[s], e.gInfo);
			logPrintf("Reading discontinuity potential from '%s' ... ", fname.c_str()); logFlush();
			loadRawBinary(VsclocDisc[s], fname.c_str());
			logPrintf("done\n"); logFlush();
		}
		logPrintf("WARNING: eigenvalsQP from a fixed potential run needs the correct fillings to be set / read in.\n");
		logPrintf("  Automatic fillings would be wrong for metals, but the discontinuity correction is zero then anyway.\n");
	}
	else
	{	VsclocDisc = getPotential(eHOMO, &eLUMO);
		for(ScalarField& V: VsclocDisc) V = JdagOJ(V); //include the dV factor (like Vscloc)
		#define StartDump(prefix) \
			string fname = e.dump.getFilename(prefix); \
			logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
		#define EndDump \
			logPrintf("done\n"); logFlush();
		#define DUMP(object, prefix) \
			{	StartDump(prefix) \
				if(mpiWorld->isHead()) saveRawBinary(object, fname.c_str()); \
				EndDump \
			}
		if(e.eInfo.spinType == SpinZ)
		{	DUMP(VsclocDisc[0], "VsclocDisc_up")
			DUMP(VsclocDisc[1], "VsclocDisc_dn")
		}
		else DUMP(VsclocDisc[0], "VsclocDisc")
		#undef DUMP
	}
	
	//Output quasiparticle energies using first-order pertubation theory:
	std::vector<diagMatrix> eigsQP(e.eInfo.nStates);
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{	matrix Hpert = e.eVars.C[q] ^ Idag_DiagV_I(e.eVars.C[q], VsclocDisc);
		diagMatrix emptyMask = eye(e.eInfo.nBands) - e.eVars.F[q]; //mask to select unoccupied states
		matrix Hqp = e.eVars.Hsub_eigs[q] + emptyMask * Hpert * emptyMask;
		matrix Hqp_evars;
		Hqp.diagonalize(Hqp_evars, eigsQP[q]);
	}
	StartDump("eigenvalsQP")
	e.eInfo.write(eigsQP, fname.c_str());
	EndDump
	#undef StartDump
	#undef EndDump
	#undef DUMP
	
	//Print quasiparticle gap:
	std::vector<double> eHOMOqp(nSpins, -DBL_MAX), eLUMOqp(nSpins, +DBL_MAX);
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{	int s = e.eInfo.qnums[q].index();
		for(int b=0; b<e.eInfo.nBands; b++)
		{	if(e.eVars.F[q][b]>=0.5) eHOMOqp[s] = std::max(eHOMOqp[s], eigsQP[q][b]);
			if(e.eVars.F[q][b]<=0.5) eLUMOqp[s] = std::min(eLUMOqp[s], eigsQP[q][b]);
		}
	}
	mpiWorld->allReduce(eHOMOqp.data(), eHOMOqp.size(), MPIUtil::ReduceMax);
	mpiWorld->allReduce(eLUMOqp.data(), eLUMOqp.size(), MPIUtil::ReduceMin);
	#define PRINT_GAP(s, sName) \
		{	double Egap = eLUMOqp[s] - eHOMOqp[s]; \
			double Edisc = eLUMOqp[s] - eLUMO[s]; \
			logPrintf("\tQuasiparticle gap%s: %lg Eh (%lg eV)  [Discontinuity contribution: %lg Eh (%lg eV)]\n", \
				sName, Egap, Egap/eV, Edisc, Edisc/eV); \
		}
	if(e.eInfo.spinType == SpinZ)
	{	PRINT_GAP(0, " (up)")
		PRINT_GAP(1, " (dn)")
	}
	else PRINT_GAP(0, "")
}

double smoothedSqrt(double de, double width)
{	if(width)
	{	double X = de/width;
		if(X<-5.) return 0.;
		else if(X<+5.) return sqrt(width * 0.5*(exp(-X*X)/sqrt(M_PI) + X*(1. + erf(X))));
		else return sqrt(de);
	}
	else return sqrt(std::max(0., de));
}

ScalarFieldArray ExCorr_OrbitalDep_GLLBsc::getPotential(std::vector<double> eHOMO, std::vector<double>* eLUMO) const
{	int nSpins = eHOMO.size();
	const double Kx = 8*sqrt(2)/(3*M_PI*M_PI);
	ScalarFieldArray V(nSpins);
	e.iInfo.augmentDensityInit();
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{	const QuantumNumber& qnum = e.eInfo.qnums[q];
		int s = qnum.index();
		diagMatrix Feff(e.eInfo.nBands);
		for(int b=0; b<e.eInfo.nBands; b++)
		{	double deTerm = smoothedSqrt(eHOMO[s]-e.eVars.Hsub_eigs[q][b], smearingWidth); //orbital-dep potential
			if(eLUMO) deTerm = smoothedSqrt((*eLUMO)[s]-e.eVars.Hsub_eigs[q][b], smearingWidth) - deTerm; //convert to the discontinuity contribution
			Feff[b] = e.eVars.F[q][b] * Kx * deTerm;
		}
		V += qnum.weight * diagouterI(Feff, e.eVars.C[q], V.size(), &e.gInfo); //without the 1/n(r) denominator
		e.iInfo.augmentDensitySpherical(qnum, Feff, e.eVars.VdagC[q]); //ultrasoft contribution
	}
	e.iInfo.augmentDensityGrid(V);
	for(int s=0; s<nSpins; s++)
	{	nullToZero(V[s], e.gInfo);
		V[s]->allReduce(MPIUtil::ReduceSum);
		e.symm.symmetrize(V[s]);
		V[s] *= inv(e.eVars.n[s]); //denominator added here
	}
	return V;
}
