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

#include <wannier/WannierMinimizerRS.h>
#include <core/WignerSeitz.h>

inline void setRealSpaceMeasures(size_t iStart, size_t iStop, const vector3<int>& S, const matrix3<>& R, const WignerSeitz* ws, vector3<double*> r, double* rSq)
{	vector3<> invS(1./S[0], 1./S[1], 1./S[2]);
	THREAD_rLoop
	(	vector3<> x; for(int k=0; k<3; k++) x[k] = iv[k] * invS[k]; //lattice coordinates
		x = ws->restrict(x); //wrap to Wigner-Seitz cell
		vector3<> rVec = ws->boundaryDistance(x) ? R * x : vector3<>(); //convert to Cartesian (and zero boundaries)
		storeVector(rVec, r, i);
		rSq[i] = rVec.length_squared();
	)
}

WannierMinimizerRS::WannierMinimizerRS(const Everything& e, const Wannier& wannier)
: WannierMinimizer(e, wannier, true)
{
	//Initialize the r and r^2 scalar fields:
	nullToZero(r, gInfoSuper);
	nullToZero(rSq, gInfoSuper);
	WignerSeitz ws(gInfoSuper.R);
	threadLaunch(setRealSpaceMeasures, gInfoSuper.nr, gInfoSuper.S, gInfoSuper.R, &ws, r.data(), rSq->data());
	
	//Split centers over MPI processes:
	TaskDivision(nCenters, mpiWorld).myRange(nStart, nStop);
	
	//Create MPI communicators for rotations:
	std::vector<std::set<int>> ranksNeeded(kMesh.size()); //rank of other processes that need each rotation
	for(size_t i=0; i<kMesh.size(); i++)
	{	int iProc = whose(i);
		for(int iSpin=0; iSpin<nSpins; iSpin++)
		{	int jProc = whose_q(i, iSpin);
			if(jProc != iProc) ranksNeeded[i].insert(jProc);
		}
	}
	for(size_t i=0; i<kMesh.size(); i++)
		if(ranksNeeded[i].size() //some communication is needed
			&& ( isMine(i) //this  process is the source
				|| ranksNeeded[i].count(mpiWorld->iProcess()) ) ) //this process is a consumer
		{
			std::vector<int> ranks(1, whose(i)); //source rank
			ranks.insert(ranks.end(), ranksNeeded[i].begin(), ranksNeeded[i].end());
			kMesh[i].mpi = std::make_shared<MPIUtil>(mpiWorld, ranks);
		}
}

void WannierMinimizerRS::initialize(int iSpin)
{	this->iSpin = iSpin;
}

static complexScalarField operator*(complex s, const complexScalarField& in)
{	complexScalarField out = clone(in);
	callPref(eblas_zscal)(out->nElem, s, out->dataPref(),1);
	return out;
}

double WannierMinimizerRS::getOmega(bool grad, bool invariant)
{	resumeOperatorThreading();
	//Compute supercell wave-functions:
	ColumnBundle Csuper(nCenters, basisSuper.nbasis*nSpinor, &basisSuper, &qnumSuper, isGpuEnabled());
	Csuper.zero();
	for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
	{	const KmeshEntry& ki = kMesh[i];
		axpyWfns(ki.point.weight, ki.U(0,nBands,0,nCenters), ki.point, iSpin, Csuper);
	}
	mpiWorld->allReduceData(Csuper, MPIUtil::ReduceSum);
	
	//Translation to handle high offsets
	std::vector<vector3<>> dr(nCenters), drReverse(nCenters);
	for(int n=0; n<nCenters; n++)
	{	drReverse[n] = gInfoSuper.invR * rPinned[n]; //translation from origin to rPinned
		dr[n] = -drReverse[n]; //and translation to rPinned as origin, both in supercell-lattice coords
	}
	translateColumns(Csuper, dr.data());
	
	//Compute spread (and optionally its gradient w.r.t Csuper):
	double Omega = 0.;
	rSqExpect.assign(nCenters, 0.);
	rExpect.assign(nCenters, vector3<>());
	ColumnBundle Omega_Csuper; if(grad) { Omega_Csuper = Csuper.similar(); Omega_Csuper.zero(); }
	for(int n=nStart; n<nStop; n++)
	{	std::vector<complexScalarField> psi(nSpinor), Omega_psi(nSpinor);
		//|r^2| contribution:
		for(int s=0; s<nSpinor; s++)
		{	psi[s] = I(Csuper.getColumn(n,s));
			complexScalarField rSq_psi = rSq * psi[s];
			rSqExpect[n] += gInfoSuper.dV * dot(psi[s], rSq_psi).real();
			if(grad) Omega_psi[s] += (2*gInfoSuper.dV) * rSq_psi;
		}
		//|r|^2 contribution:
		for(int dir=0; dir<3; dir++)
			for(int s=0; s<nSpinor; s++)
				rExpect[n][dir] += gInfoSuper.dV * dot(psi[s], r[dir]* psi[s]).real();
		//Accumulate contributions to Omega:
		bool shouldPin = pinned[n] && !invariant;
		const vector3<> r0_n = shouldPin ? vector3<>() : rExpect[n];
		Omega += (rSqExpect[n] - 2*dot(rExpect[n],r0_n) + r0_n.length_squared());
		//Shift back from rPinned to original origin for stored outputs:
		rSqExpect[n] -= rExpect[n].length_squared();
		rExpect[n] += rPinned[n]; 
		rSqExpect[n] += rExpect[n].length_squared();
		//Propagate gradient via rExpect:
		if(grad)
		{	for(int dir=0; dir<3; dir++)
				for(int s=0; s<nSpinor; s++)
					Omega_psi[s] += (-2.*r0_n[dir] * 2*gInfoSuper.dV) * r[dir]*psi[s];
		}
		//Off-diagonal corrections to get the invariant part (expensive):
		if(invariant)
		{	for(int m=0; m<nCenters; m++) if(m != n)
			{	std::vector<complexScalarField> psi_m(2), Omega_psi_m(2);
				for(int s=0; s<nSpinor; s++)
					psi_m[s] = I(Csuper.getColumn(m,s));
				for(int dir=0; dir<3; dir++)
				{	complex rExpect_mn;
					for(int s=0; s<nSpinor; s++)
						rExpect_mn += gInfoSuper.dV * dot(psi_m[s], r[dir] * psi[s]);
					Omega -= rExpect_mn.norm();
					if(grad)
						for(int s=0; s<nSpinor; s++)
						{	Omega_psi[s] -= 2*rExpect_mn * gInfoSuper.dV * (r[dir] * psi_m[s]);
							Omega_psi_m[s] -= 2*rExpect_mn.conj() * gInfoSuper.dV * (r[dir] * psi[s]);
						}
				}
				if(grad)
					for(int s=0; s<nSpinor; s++)
						Omega_Csuper.accumColumn(m,s, Idag(Omega_psi_m[s]));
			}
		}
		if(grad)
			for(int s=0; s<nSpinor; s++)
				Omega_Csuper.accumColumn(n,s, Idag(Omega_psi[s]));
	}
	mpiWorld->allReduce(Omega, MPIUtil::ReduceSum);
	mpiWorld->allReduceData(rSqExpect, MPIUtil::ReduceSum);
	mpiWorld->allReduceData(rExpect, MPIUtil::ReduceSum);
	
	if(grad)
	{	//Collect across processes:
		mpiWorld->allReduceData(Omega_Csuper, MPIUtil::ReduceSum);
		translateColumns(Omega_Csuper, drReverse.data());
		//Propagate to rotations:
		for(unsigned i=0; i<kMesh.size(); i++) if(isMine_q(i,iSpin))
		{	KmeshEntry& ki = kMesh[i];
			matrix Omega_U;
			axpyWfns_grad(ki.point.weight, Omega_U, ki.point, iSpin, Omega_Csuper);
			ki.Omega_UdotU += Omega_U(0,nCenters, 0,nBands) * ki.U;
		}
	}
	suspendOperatorThreading();
	return Omega;
}

double WannierMinimizerRS::getOmega(bool grad)
{	return getOmega(grad, false);
}

double WannierMinimizerRS::getOmegaI(bool grad)
{	return getOmega(grad, true);
}
