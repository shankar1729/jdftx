/*-------------------------------------------------------------------
Copyright 2026 Ravishankar Sundararaman, Ben Rich.

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

#include <electronic/VanDerWaalsD3.h>
#include <electronic/VanDerWaalsD3_data.h>
#include <electronic/Everything.h>

//! Fast integer powers, bypassing slower std::pow
template<int m> double intPow(double x)
{	if(m < 0) return intPow<-m>(1./x);
	double result = 1., xPow = x; //xPow = x^(2^k) for bit k of m
	for(int k=m; k>0; k/=2)
	{	if(k % 2) result *= xPow;
		xPow *= xPow;
	}
	return result;
}

namespace D3
{
	AtomParams getAtomParams(int Z)
	{	if((Z < 1) or (Z > D3::Zmax))
			die("\nAtomic number %d not supported in DFT-D3\n\n", Z);
		AtomParams ap;
		ap.Z = Z;
		ap.sqrtQ = sqrtQ[Z - 1]; //Note: per-atom arrays have 0-based index
		ap.k2Rcov = k2Rcov[Z - 1];
		for(int iCN=0; iCN<D3::C6::nCN; iCN++)
			if(Z == D3::C6::Z[iCN])
			{	ap.CN.push_back(D3::C6::CN[iCN]);
				ap.iCN.push_back(iCN);
			}
		return ap;
	}
	
	matrix AtomParams::getL(double observedCN, matrix& Lprime) const
	{	matrix L(CN.size(), 1);
		Lprime.init(CN.size(), 1);
		complex* Ldata = L.data();
		complex* LprimeData = Lprime.data();
		double Lsum = 0., LprimeSum = 0.;
		for(double CNi: CN)
		{	double Li = exp(-D3::k3 * intPow<2>(observedCN - CNi));
			double LiPrime = (-2. * D3::k3) * (observedCN - CNi) * Li;
			Lsum += Li;
			LprimeSum += LiPrime;
			*(Ldata++) = Li;
			*(LprimeData++) = LiPrime;
		}
		//Check for underflow issues:
		if(Lsum < 1E-40)
		{	//Only happens when observedCN > max(CN) + 5.6
			//If this is the case, contributions other than from max(CN) negligible at 1E-16 level
			int iCNmax = std::max_element(CN.begin(), CN.end()) - CN.begin();
			Lprime.zero();
			L.zero();
			L.set(iCNmax, 0, 1.);
			return L;
		}
		//Normalize:
		double invLsum = 1./Lsum;
		Lprime = invLsum*(Lprime - (invLsum*LprimeSum) * L);
		return L * invLsum;
	}
	
	PairParams getPairParams(const AtomParams& ap1, const AtomParams& ap2)
	{	PairParams pp;
		//Radius:
		int iZmin = std::min(ap1.Z, ap2.Z) - 1; //index into arrays by lower of two atomic numbers
		int iZmax = std::max(ap1.Z, ap2.Z) - 1;
		int iZ12 = (iZmax*(iZmax+1))/2 + iZmin; //index into unique Z pairs
		pp.R0 = D3::R0_Angstrom[iZ12] * Angstrom;
		//C6 coefficients:
		pp.C6.init(ap1.CN.size(), ap2.CN.size());
		complex* C6data = pp.C6.data();
		for(int iCN2: ap2.iCN)
			for(int iCN1: ap1.iCN)
			{	int iCNmin = std::min(iCN1, iCN2);
				int iCNmax = std::max(iCN1, iCN2);
				int iCN12 = (iCNmax*(iCNmax+1))/2 + iCNmin; //index into unique iCN pairs
				*(C6data++) = D3::C6::coeff[iCN12];
			}
		return pp;
	}
}

//----- Implementation of class VanDerWaalsD3 -----

VanDerWaalsD3::VanDerWaalsD3(const Everything& e)
: VanDerWaals(e), dampBJ(e.iInfo.vdWstyle == VDW_D3BJ)
{
	logPrintf("\nInitializing DFT-D3%s calculator:\n", dampBJ ? "(BJ)" : "");
	
	//Get parameters for exchange-correlation functional
	string xcName = e.exCorr.getName();
	D3::setXCscale(xcName, dampBJ, s6, sr6, s8, sr8);
	logPrintf("\tParameters set for %s functional\n", xcName.c_str());
	if(dampBJ)
	{	logPrintf("\ts6: %6.3lf  s8: %6.3lf\n", s6, s8);
		logPrintf("\ta1: %6.3lf  a2: %6.3lf\n", sr6, sr8);
	}
	else
	{	logPrintf("\ts6: %6.3lf  s_r6: %6.3lf\n", s6, sr6);
		logPrintf("\ts8: %6.3lf  s_r8: %6.3lf\n", s8, sr8);
	}

	//Get per-atom parameters:
	logPrintf("\tPer-atom parameters loaded for:\n");
	for(size_t spIndex=0; spIndex<e.iInfo.species.size(); spIndex++)
	{	const auto& sp = e.iInfo.species[spIndex];
		assert(sp->atomicNumber);
		D3::AtomParams ap = D3::getAtomParams(sp->atomicNumber);
		atomParams.push_back(ap);
		logPrintf("\t%2s:  sqrtQ[a0]: %6.3f  Rcov[a0]: %6.3f  CN: [",
			sp->name.c_str(), ap.sqrtQ, ap.k2Rcov / D3::k2);
		for(double CNi: ap.CN) logPrintf(" %.2f", CNi);
		logPrintf(" ]\n");
	}

	//Get per-pair paramaters:
	pairParams.resize(atomParams.size());
	for(size_t iSp1=0; iSp1<atomParams.size(); iSp1++)
	{	pairParams[iSp1].resize(atomParams.size());
		for(size_t iSp2=0; iSp2<atomParams.size(); iSp2++)
			pairParams[iSp1][iSp2] = D3::getPairParams(atomParams[iSp1], atomParams[iSp2]);
	}

	Citations::add(
		"DFT-D3 dispersion correction",
		dampBJ
		? "S. Grimme, S. Ehrlich and L. Goerigk, J. Comp. Chem. 32, 1456 (2011)"
		: "S. Grimme, J. Antony, S. Ehrlich and H. Krieg, J. Chem. Phys. 132, 154104 (2010)"
	);
}


double VanDerWaalsD3::getScaleFactor(string exCorrName, double scaleOverride) const
{	return 0.;  //global scale factor not used in D3
}


//! Return r^-n pair-potential energy and set its derivative E_r.
//! The potential is damped at length scale R0 and exponent alpha_n.
template<int n, int alpha_n> double vdWpotential(double invr, double R0, double& E_r)
{	//Main r^-n potential (and r derivative):
	double pot = intPow<n>(invr);
	double pot_r = (-n) * pot * invr;
	//Exponent in denominator of damping factor (and r derivative):
	double damp_term = 6 * intPow<alpha_n>(R0 * invr);
	double damp_term_r = (-alpha_n) * damp_term * invr;
	//Overall damping factor (and r derivative):
	double damp = 1./(1. + damp_term);
	double damp_r = -damp_term_r * damp * damp;
	//Combine to get final energy and derivative:
	E_r = pot_r * damp + pot * damp_r;
	return pot * damp;
}

//! Return r^-n pair-potential energy and set its derivative E_r for the BJ-damping form.
template<int n> double vdWpotentialBJ(double r, double R0, double a1, double a2, double& E_r)
{	//Main r^-n potential (and r derivative):
	double subcdamp = a1*R0 + a2; // eqn 6 of Grimme 2011
	double cdamp = intPow<n>(subcdamp); // raised to n (constant wrt r)
	double rn = intPow<n>(r);
	double E = 1 / (rn + cdamp);
	E_r = -n * intPow<n-1>(r) * E * E;
	return E;
}

//! Return trace(transpose(L1) * C6 * L2), with C6 the real parts of a matrix in its column-major storage order.
inline double contractC6(const std::vector<double>& C6, const std::vector<double>& L1, const std::vector<double>& L2)
{	const double* C6col = C6.data();
	double result = 0.;
	for(double L2j: L2)
	{	double L1C6 = 0.;
		for(double L1i: L1)
			L1C6 += L1i * *(C6col++);
		result += L1C6 * L2j;
	}
	return result;
}


//! Thread-parallelize energy calculations over first atom
void energyAndGrad_sub(int c1, int nAtoms, const Atom* atoms, const D3::AtomParams* atomParams,
	const std::vector<D3::PairParams>* pairParams, const std::vector<double>* Lall, const std::vector<double>* LprimeAll,
	const std::vector<std::vector<double>>* C6all, bool dampBJ, double s6, double s8, double sr6, double sr8,
	const GridInfo* gInfo, double rCutSq, vector3<int> S, size_t iStart, size_t iStop, double* diagC6,
	double* E6tot, double* E8tot, double* E_CNtot, vector3<>* forcesTot, matrix3<>* E_RRTtot, std::mutex* lock)
{	const D3::AtomParams& ap1 = atomParams[atoms[c1].sp];
	const std::vector<double>& L1 = Lall[c1], &L1prime = LprimeAll[c1]; //C6 interpolation weights for atom 1
	double E6 = 0., E8 = 0.; std::vector<double> E_CN(nAtoms); //contributions of the pairs of atom c1,
	std::vector<vector3<>> forces(nAtoms); matrix3<> E_RRT;      //added to the totals (E6tot etc.) at the end

	for(int c2=0; c2<nAtoms; c2++)
	{	const D3::AtomParams& ap2 = atomParams[atoms[c2].sp];
		const std::vector<double>& L2 = Lall[c2], &L2prime = LprimeAll[c2]; //C6 interpolation weights for atom 2

		//Compute C6 and C8 for this pair:
		const D3::PairParams& pp = pairParams[atoms[c1].sp][atoms[c2].sp];
		double ratio8by6 = 3. * ap1.sqrtQ * ap2.sqrtQ;
		double C6 = contractC6(C6all[atoms[c1].sp][atoms[c2].sp], L1, L2);
		double C8 = C6 * ratio8by6;
		double R0 = dampBJ ? sqrt(C8/C6) : pp.R0; //BJ damping uses C8/C6 to get R0
		if(c1 == c2) diagC6[c1] = C6; //only for reporting
		
		//Sum energy, force and stress contributions over periodic images:
		double E12_C6 = 0., E12_C8 = 0.; //energy contributions upto C6 and C8 prefactors
		THREAD_halfGspaceLoop(
			const vector3<int>& iR = iG;
			vector3<> x = iR + (atoms[c1].pos - atoms[c2].pos);
			double rSq = gInfo->RTR.metric_length_squared(x);
			if(rSq and rSq<=rCutSq)
			{	double r = sqrt(rSq);
				double cellWeight = (iR[2] ? 1. : 0.5); //account for double-counting in half-space cut plane
				double invr = 1./r;
				double term6_r; double term6; double term8_r; double term8;
				if(dampBJ) 
				{	term6 = (vdWpotentialBJ<6>(r, R0, sr6, sr8, term6_r));
					term8 = (vdWpotentialBJ<8>(r, R0, sr6, sr8, term8_r));
				} 
				else
				{	term6 = (vdWpotential<6, D3::alpha6>(invr, sr6 * R0, term6_r));
					term8 = (vdWpotential<8, D3::alpha8>(invr, sr8 * R0, term8_r));
				}
				E12_C6 -= cellWeight * s6 * term6;
				E12_C8 -= cellWeight * s8 * term8;
				//Colect forces and/or stresses:
				double E_r_by_r = (-invr * cellWeight) * (C6*s6*term6_r + C8*s8*term8_r);
				vector3<> E_x = E_r_by_r * (gInfo->RTR * x); 
				forces[c1] -= E_x;
				forces[c2] += E_x;
				if(E_RRTtot)
				{	const vector3<> rVec = gInfo->R * x;
					E_RRT += E_r_by_r * outer(rVec, rVec);
				}
			}
		)
		E6 += E12_C6 * C6;
		E8 += E12_C8 * C8;
		
		//Propagate gradients to CN:
		double E12_C6_tot = E12_C6 + E12_C8 * ratio8by6; //total derivative w.r.t C6
		E_CN[c1] += E12_C6_tot * contractC6(C6all[atoms[c1].sp][atoms[c2].sp], L1prime, L2);
		E_CN[c2] += E12_C6_tot * contractC6(C6all[atoms[c1].sp][atoms[c2].sp], L1, L2prime);
	}
	std::lock_guard<std::mutex> guard(*lock);
	*E6tot += E6; *E8tot += E8;
	if(E_RRTtot) *E_RRTtot += E_RRT;
	for(int c=0; c<nAtoms; c++)
	{	E_CNtot[c] += E_CN[c];
		forcesTot[c] += forces[c];
	}
}


double VanDerWaalsD3::energyAndGrad(std::vector<Atom>& atoms, const double scaleFac, matrix3<>* E_RRTptr) const
{	static StopWatch watch("VanDerWaalsD3::energyAndGrad"); watch.start();
	const double rCut = e.iInfo.ljOverride ? e.iInfo.ljOverride : 200.; //Truncate summation at 1/r^6 ~ 10^-16
	const double rCutSq = rCut*rCut;
	vector3<int> S; size_t iStart, iStop; //Number of unit cells in supercell to reach rCut and MPI division
	setNeighborSampling(rCut, S, iStart, iStop);
	logPrintf("\nComputing DFT-D3 correction:\n");

	//Get coordination numbers:
	std::vector<double> CN;
	computeCN(atoms, CN);
	
	//Compute energy and direct force/stress contributions :
	double E6 = 0., E8 = 0.; //!< r^-6 and r^-8 energies
	std::vector<double> E_CN(atoms.size()); //coordination number gradients
	std::vector<vector3<>> forces(atoms.size()); //VDW forces per atom
	matrix3<> E_RRT; //Stress * volume (updated only if E_RRTptr is non-null)
	std::vector<double> diagC6(atoms.size()); //diagonal C6 for reporting
	//C6 interpolation weights per atom and reference C6 per pair of species, as plain arrays for contractC6:
	std::vector<std::vector<double>> Lall(atoms.size()), LprimeAll(atoms.size());
	std::vector<std::vector<std::vector<double>>> C6all(pairParams.size());
	for(size_t c=0; c<atoms.size(); c++)
	{	matrix Lprime, L = atomParams[atoms[c].sp].getL(CN[c], Lprime);
		for(int i=0; i<L.nRows(); i++)
		{	Lall[c].push_back(L.data()[i].real());
			LprimeAll[c].push_back(Lprime.data()[i].real());
		}
	}
	for(size_t iSp1=0; iSp1<pairParams.size(); iSp1++)
		for(const D3::PairParams& pp: pairParams[iSp1])
		{	C6all[iSp1].push_back(std::vector<double>());
			for(size_t i=0; i<pp.C6.nData(); i++)
				C6all[iSp1].back().push_back(pp.C6.data()[i].real());
		}
	std::mutex lock;
	if(atoms.size()) threadedLoop(energyAndGrad_sub, atoms.size(), int(atoms.size()), atoms.data(), atomParams.data(),
		pairParams.data(), Lall.data(), LprimeAll.data(), C6all.data(), dampBJ, s6, s8, sr6, sr8, &e.gInfo, rCutSq, S,
		iStart, iStop, diagC6.data(), &E6, &E8, E_CN.data(), forces.data(), E_RRTptr ? &E_RRT : NULL, &lock);
	report(diagC6, "diagonal-C6", atoms, " %.2f");
	mpiWorld->allReduce(E6, MPIUtil::ReduceSum);
	mpiWorld->allReduce(E8, MPIUtil::ReduceSum);
	mpiWorld->allReduceData(E_CN, MPIUtil::ReduceSum);
	logPrintf("EvdW_6 = %11.6lf\n", E6);
	logPrintf("EvdW_8 = %11.6lf\n", E8);
	
	//Propagate gradients w.r.t CN to forces/stresses
	propagateCNgradient(atoms, E_CN, forces, E_RRTptr ? &E_RRT : NULL);

	//Collect forces and stresses:
	mpiWorld->allReduceData(forces, MPIUtil::ReduceSum, true);
	for(int c=0; c<int(atoms.size()); c++)
		atoms[c].force += forces[c];
	if(E_RRTptr)
	{	mpiWorld->allReduce(E_RRT, MPIUtil::ReduceSum, true);
		*E_RRTptr += E_RRT;
	}
	watch.stop();
	return E6 + E8;
}


//! Thread-parallelize coordination number calculations over first atom
void computeCN_sub(int c1, int nAtoms, const Atom* atoms, const D3::AtomParams* atomParams,
	const GridInfo* gInfo, double rCutSq, vector3<int> S, size_t iStart, size_t iStop, double* CNtot, std::mutex* lock)
{	const D3::AtomParams& ap1 = atomParams[atoms[c1].sp];
	std::vector<double> CN(nAtoms); //contributions of the pairs of atom c1, added to CNtot at the end
	for(int c2=0; c2<nAtoms; c2++)
	{	const D3::AtomParams& ap2 = atomParams[atoms[c2].sp];
		double k2RcovSum = ap1.k2Rcov + ap2.k2Rcov;
		THREAD_halfGspaceLoop(
			const vector3<int>& iR = iG;
			vector3<> x = iR + (atoms[c1].pos - atoms[c2].pos);
			double rSq = gInfo->RTR.metric_length_squared(x);
			if(rSq and rSq<=rCutSq)
			{	double r = sqrt(rSq);
				double cellWeight = (iR[2] ? 1. : 0.5); //account for double-counting in half-space cut plane
				double CNterm = cellWeight/(1. + exp(-D3::k1*(k2RcovSum/r - 1.)));
				CN[c1] += CNterm;
				CN[c2] += CNterm;
			}
		)
	}
	std::lock_guard<std::mutex> guard(*lock);
	for(int c=0; c<nAtoms; c++)
		CNtot[c] += CN[c];
}


//Compute local coordination number
void VanDerWaalsD3::computeCN(const std::vector<Atom>& atoms, std::vector<double>& CN) const
{	const double rCut = 50., rCutSq=rCut*rCut; //Damping factor in CN calculation drops off more quickly than dispersion term
	vector3<int> S; size_t iStart, iStop; //Number of unit cells in supercell to reach rCut and MPI division
	setNeighborSampling(rCut, S, iStart, iStop);
	CN.assign(atoms.size(), 0.);
	std::mutex lock;
	if(atoms.size()) threadedLoop(computeCN_sub, atoms.size(), int(atoms.size()), atoms.data(), atomParams.data(),
		&e.gInfo, rCutSq, S, iStart, iStop, CN.data(), &lock);
	mpiWorld->allReduceData(CN, MPIUtil::ReduceSum);
	report(CN, "coordination-number", atoms);
}


//! Thread-parallelize coordination number gradient propagation with respect to first atom
void propagateCNgradient_sub(int c1, int nAtoms, const Atom* atoms, const D3::AtomParams* atomParams, const double* E_CN,
	const GridInfo* gInfo, double rCutSq, vector3<int> S, size_t iStart, size_t iStop,
	vector3<>* forcesTot, matrix3<>* E_RRTtot, std::mutex* lock)
{	const D3::AtomParams& ap1 = atomParams[atoms[c1].sp];
	std::vector<vector3<>> forces(nAtoms); matrix3<> E_RRTsum; //contributions of the pairs of atom c1,
	matrix3<>* E_RRT = E_RRTtot ? &E_RRTsum : NULL;            //added to the totals (forcesTot, E_RRTtot) at the end
	for(int c2=0; c2<nAtoms; c2++)
	{	const D3::AtomParams& ap2 = atomParams[atoms[c2].sp];
		double k2RcovSum = ap1.k2Rcov + ap2.k2Rcov;
		THREAD_halfGspaceLoop(
			const vector3<int>& iR = iG;
			vector3<> x = iR + (atoms[c1].pos - atoms[c2].pos);
			double rSq = gInfo->RTR.metric_length_squared(x);
			if(rSq and rSq<=rCutSq)
			{	double r = sqrt(rSq);
				double invr = 1./r;
				double cellWeight = (iR[2] ? 1. : 0.5); //account for double-counting in half-space cut plane
				double E_CNterm = E_CN[c1] + E_CN[c2];
				double expTerm = exp(-D3::k1*(k2RcovSum*invr - 1.));
				double expTerm_r = expTerm * D3::k1 * (k2RcovSum * invr * invr);
				double E_r_by_r = (-invr * E_CNterm * cellWeight * expTerm_r) / intPow<2>(1+expTerm);
				//Collect forces and/or stresses:
				vector3<> E_x = E_r_by_r * (gInfo->RTR * x); 
				forces[c1] -= E_x;
				forces[c2] += E_x;
				if(E_RRT)
				{	const vector3<> rVec = gInfo->R * x;
					*E_RRT += E_r_by_r * outer(rVec, rVec);
				}
			}
		)
	}
	std::lock_guard<std::mutex> guard(*lock);
	for(int c=0; c<nAtoms; c++)
		forcesTot[c] += forces[c];
	if(E_RRTtot) *E_RRTtot += E_RRTsum;
}


//Propagate coordination-number gradient to forces, and optionally, stresses:
void VanDerWaalsD3::propagateCNgradient(const std::vector<Atom>& atoms, const std::vector<double>& E_CN,
	std::vector<vector3<>>& forces, matrix3<>* E_RRT) const
{	//Set up same neighbor cells as used in computeCN()
	const double rCut = 50., rCutSq=rCut*rCut;
	vector3<int> S; size_t iStart, iStop;
	setNeighborSampling(rCut, S, iStart, iStop);
	//Propagate gradients corresponding to computeCN()
	std::mutex lock;
	if(atoms.size()) threadedLoop(propagateCNgradient_sub, atoms.size(), int(atoms.size()), atoms.data(), atomParams.data(),
		E_CN.data(), &e.gInfo, rCutSq, S, iStart, iStop, forces.data(), E_RRT, &lock);
}


void VanDerWaalsD3::report(const std::vector<double>& result, string name, const std::vector<Atom>& atoms, const char* fmt) const
{	size_t c = 0;
	for(int sp=0; sp<int(atomParams.size()); sp++)
	{	logPrintf("# %s %s", name.c_str(), e.iInfo.species[sp]->name.c_str());
		while((c < atoms.size()) and (atoms[c].sp == sp))
		{	logPrintf(fmt, result[c]);
			c++;
		}
		logPrintf("\n");
	}
	logFlush();
}


void VanDerWaalsD3::setNeighborSampling(double rCut, vector3<int>& S, size_t& iStart, size_t& iStop) const
{	vector3<bool> isTruncated = e.coulombParams.isTruncated();
	for(int k=0; k<3; k++)
		S[k] = 1 + 2*(isTruncated[k] ? 0 : (int)ceil(rCut/e.gInfo.R.column(k).length()));
	size_t nCellsHlf = S[0] * S[1] * (S[2]/2+1);; //similar to the half-G-space used for FFTs
	TaskDivision(nCellsHlf, mpiWorld).myRange(iStart, iStop);
}
