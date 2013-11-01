/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman

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

#include <electronic/SpeciesInfo.h>
#include <electronic/SpeciesInfo_internal.h>
#include <electronic/Everything.h>
#include <electronic/matrix.h>
#include <electronic/operators.h>
#include <electronic/ColumnBundle.h>

//------- primary SpeciesInfo functions involved in simple energy and gradient calculations (with norm-conserving pseudopotentials) -------


//Return non-local energy and optionally accumulate its electronic and/or ionic gradients for a given quantum number
double SpeciesInfo::EnlAndGrad(const diagMatrix& Fq, const ColumnBundle& Cq, ColumnBundle& HCq, std::vector< vector3<> >* forces) const
{	static StopWatch watch("EnlAndGrad"); watch.start();
	if(!atpos.size()) return 0.0; //unused species
	const GridInfo &gInfo = e->gInfo;
	const Basis& basis = *Cq.basis;
	
	double Enlq = 0.0;
	for(int l=0; l<int(VnlRadial.size()); l++)
	{	unsigned nProj = VnlRadial[l].size(); if(!nProj) continue; //skip l if no projectors
		//Copy Mnl into a block-diagonal form for all atoms:
		tiledBlockMatrix M(this->Mnl[l], atpos.size()); 
		//Allocate temporaries:
		ColumnBundle V = Cq.similar(nProj * atpos.size());
		ColumnBundle dV[3]; vector3<complex*> dVdata;
		if(forces)
			for(int k=0; k<3; k++)
			{	dV[k] = V.similar();
				dVdata[k] = dV[k].dataPref();
			}
		for(int m=-l; m<=l; m++)
		{	// Calculate the nonlocal projectors (and optionally their spatial derivatives):
			for(unsigned p=0; p<nProj; p++)
			{	size_t offs = p * basis.nbasis;
				size_t atomStride = nProj * basis.nbasis;
				callPref(Vnl)(basis.nbasis, atomStride, atpos.size(), l, m, V.qnum->k, basis.iGarrPref, gInfo.G,
					atposPref, VnlRadial[l][p], V.dataPref()+offs, forces, dVdata+offs);
			}
			// Compute contribution to Enl:
			matrix VdagC = V^Cq;
			matrix MVdagC = M*VdagC;
			Enlq += trace(Fq * dagger(VdagC) * MVdagC).real();
			//Update electronic gradient if requested
			if(HCq) HCq += V * MVdagC;
			//Update forces if requested:
			if(forces)
			{	for(int k=0; k<3; k++)
				{	diagMatrix fdiag = diag(MVdagC * Fq * dagger(dV[k]^Cq));
					for(unsigned atom=0; atom<atpos.size(); atom++)
						for(unsigned p=0; p<nProj; p++)
							(*forces)[atom][k] -= 2.0*Cq.qnum->weight * fdiag[p+atom*nProj];
				}
			}
		}
	}
	watch.stop();
	return Enlq;
}

//Compute DFT+U corrections:
double SpeciesInfo::computeU(const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C,
	std::vector<ColumnBundle>* HC, std::vector<vector3<> >* forces) const
{	if(!plusU.size()) return 0.; //no U for this species
	const ElecInfo& eInfo = e->eInfo;
	int nSpins = eInfo.spinType==SpinNone ? 1 : 2; //number of spins
	int qCount = eInfo.nStates/nSpins; //number of states of each spin
	double wSpinless = 0.5*nSpins; //factor multiplying state weights to get to spinless weights
	double Utot = 0.;
	for(int s=0; s<nSpins; s++)
		for(auto Uparams: plusU)
		{	int mCount = 2*Uparams.l+1; //number of m's at given l
			double prefac = 0.5 * Uparams.UminusJ / wSpinless;
			//Compute the density matrix:
			matrix rho;
			for(int q=s*qCount; q<(s+1)*qCount; q++)
			{	ColumnBundle Opsi(C[q].similar(atpos.size() * mCount));
				setOpsi(Opsi, Uparams.n, Uparams.l);
				matrix CdagOpsi = C[q] ^ Opsi;
				rho += (eInfo.qnums[q].weight*wSpinless) * dagger(CdagOpsi) * F[q] * CdagOpsi;
			}
			//Compute contributions to U and its derivative w.r.t density matrix rho:
			matrix U_rho = zeroes(rho.nRows(), rho.nCols());
			for(unsigned a=0; a<atpos.size(); a++)
			{	matrix rhoSub = rho(a,atpos.size(),rho.nRows(), a,atpos.size(),rho.nCols());
				Utot += prefac * trace(rhoSub - rhoSub*rhoSub).real();
				U_rho.set(a,atpos.size(),rho.nRows(), a,atpos.size(),rho.nCols(), prefac * (eye(mCount) - 2.*rhoSub));
			}
			//Propagate gradient from U_rho to wavefunctions or ionic positions if required:
			if(HC || forces)
			{	for(int q=s*qCount; q<(s+1)*qCount; q++)
				{	ColumnBundle Opsi(C[q].similar(atpos.size() * mCount));
					std::vector<ColumnBundle> dOpsi;
					setOpsi(Opsi, Uparams.n, Uparams.l, forces ? &dOpsi : 0);
					matrix CdagOpsi = C[q] ^ Opsi;
					if(HC) HC->at(q) += wSpinless * Opsi * (U_rho * dagger(CdagOpsi)); //gradient upto state weight and fillings
					if(forces)
					{	for(int k=0; k<3; k++)
						{	diagMatrix fMat = wSpinless * diag(U_rho * dagger(CdagOpsi) * F[q] * (C[q]^dOpsi[k]));
							for(unsigned a=0; a<atpos.size(); a++)
								(*forces)[a][k] -= 2.0*C[q].qnum->weight * trace(fMat(a,atpos.size(),fMat.nRows()));
						}
					}
				}
			}
		}
	return Utot;
}

void SpeciesInfo::updateLocal(DataGptr& Vlocps, DataGptr& rhoIon, DataGptr& nChargeball,
	DataGptr& nCore, DataGptr& tauCore) const
{	if(!atpos.size()) return; //unused species
	((SpeciesInfo*)this)->updateLatticeDependent(); //update lattice dependent quantities (if lattice vectors have changed)
	const GridInfo& gInfo = e->gInfo;

	//Prepare optional outputs:
	complex *nChargeballData=0, *nCoreData=0, *tauCoreData=0;
	if(Z_chargeball) { nullToZero(nChargeball, gInfo); nChargeballData = nChargeball->dataPref(); }
	if(nCoreRadial) { nullToZero(nCore, gInfo); nCoreData = nCore->dataPref(); }
	if(tauCoreRadial) { nullToZero(tauCore, gInfo); tauCoreData = tauCore->dataPref(); }
	
	//Calculate in half G-space:
	double invVol = 1.0/gInfo.detR;
	callPref(::updateLocal)(gInfo.S, gInfo.GGT,
		Vlocps->dataPref(), rhoIon->dataPref(), nChargeballData, nCoreData, tauCoreData,
		atpos.size(), atposPref, invVol, VlocRadial,
		Z, nCoreRadial, tauCoreRadial, Z_chargeball, width_chargeball);
}


std::vector< vector3<double> > SpeciesInfo::getLocalForces(const DataGptr& ccgrad_Vlocps,
	const DataGptr& ccgrad_rhoIon, const DataGptr& ccgrad_nChargeball,
	const DataGptr& ccgrad_nCore, const DataGptr& ccgrad_tauCore) const
{	
	if(!atpos.size()) return std::vector< vector3<double> >(); //unused species, return empty forces
	
	const GridInfo& gInfo = e->gInfo;
	complex* ccgrad_rhoIonData = ccgrad_rhoIon ? ccgrad_rhoIon->dataPref() : 0;
	complex* ccgrad_nChargeballData = (Z_chargeball && ccgrad_nChargeball) ? ccgrad_nChargeball->dataPref() : 0;
	complex* ccgrad_nCoreData = nCoreRadial ? ccgrad_nCore->dataPref() : 0;
	complex* ccgrad_tauCoreData = (tauCoreRadial && ccgrad_tauCore) ? ccgrad_tauCore->dataPref() : 0;
	
	//Propagate ccgrad* to gradient w.r.t structure factor:
	DataGptr ccgrad_SG(DataG::alloc(gInfo, isGpuEnabled())); //complex conjugate gradient w.r.t structure factor
	callPref(gradLocalToSG)(gInfo.S, gInfo.GGT,
		ccgrad_Vlocps->dataPref(), ccgrad_rhoIonData, ccgrad_nChargeballData,
		ccgrad_nCoreData, ccgrad_tauCoreData, ccgrad_SG->dataPref(), VlocRadial,
		Z, nCoreRadial, tauCoreRadial, Z_chargeball, width_chargeball);
	
	//Now propagate that gradient to each atom of this species:
	DataGptrVec gradAtpos; nullToZero(gradAtpos, gInfo);
	vector3<complex*> gradAtposData; for(int k=0; k<3; k++) gradAtposData[k] = gradAtpos[k]->dataPref();
	std::vector< vector3<> > forces(atpos.size());
	for(unsigned at=0; at<atpos.size(); at++)
	{	callPref(gradSGtoAtpos)(gInfo.S, atpos[at], ccgrad_SG->dataPref(), gradAtposData);
		for(int k=0; k<3; k++)
			forces[at][k] = -sum(gradAtpos[k]); //negative gradient
	}
	return forces;
}
