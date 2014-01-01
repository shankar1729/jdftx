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
double SpeciesInfo::EnlAndGrad(const QuantumNumber& qnum, const diagMatrix& Fq, const matrix& VdagCq, matrix& HVdagCq) const
{	static StopWatch watch("EnlAndGrad"); watch.start();
	if(!atpos.size()) return 0.; //unused species
	if(!MnlAll) return 0.; //purely local psp
	int nProj = MnlAll.nRows();
	
	matrix MVdagC = zeroes(VdagCq.nRows(), VdagCq.nCols());
	double Enlq = 0.0;
	for(unsigned atom=0; atom<atpos.size(); atom++)
	{	matrix atomVdagC = VdagCq(atom*nProj,(atom+1)*nProj, 0,VdagCq.nCols());
		matrix MatomVdagC = MnlAll * atomVdagC;
		MVdagC.set(atom*nProj,(atom+1)*nProj, 0,VdagCq.nCols(), MatomVdagC);
		Enlq += trace(Fq * dagger(atomVdagC) * MatomVdagC).real();
	}
	HVdagCq += MVdagC;
	watch.stop();
	return Enlq;
}

//Compute DFT+U corrections:
double SpeciesInfo::computeU(const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C,
	std::vector<ColumnBundle>* HC, std::vector<vector3<> >* forces, FILE* fpRhoAtom) const
{	static StopWatch watch("computeU"); watch.start();
	if(!plusU.size()) return 0.; //no U for this species
	const ElecInfo& eInfo = e->eInfo;
	int nSpins = eInfo.spinType==SpinNone ? 1 : 2; //number of spins
	double wSpinless = 0.5*nSpins; //factor multiplying state weights to get to spinless weights
	double Utot = 0.;
	for(auto Uparams: plusU)
	{	int mCount = 2*Uparams.l+1; //number of m's at given l
		int matSize = mCount * atpos.size(); //density matrix size
		double prefac = 0.5 * Uparams.UminusJ / wSpinless;
		//Compute the density matrix:
		std::vector<matrix> rho(nSpins), U_rho(nSpins);
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{	int s = eInfo.qnums[q].index();
			ColumnBundle Opsi(C[q].similar(matSize));
			setOpsi(Opsi, Uparams.n, Uparams.l);
			matrix CdagOpsi = C[q] ^ Opsi;
			rho[s] += (eInfo.qnums[q].weight*wSpinless) * dagger(CdagOpsi) * F[q] * CdagOpsi;
		}
		for(int s=0; s<nSpins; s++)
		{	//Collect contributions from all processes:
			if(!rho[s]) rho[s] = zeroes(matSize, matSize);
			rho[s].allReduce(MPIUtil::ReduceSum);
			//Symmetrize:
			for(unsigned sp=0; sp<e->iInfo.species.size(); sp++)
				if(e->iInfo.species[sp].get()==this)
					e->symm.symmetrizeSpherical(rho[s], sp);
			//Compute contributions to U and its derivative w.r.t density matrix rho:
			U_rho[s] = zeroes(matSize,matSize);
			for(unsigned a=0; a<atpos.size(); a++)
			{	matrix rhoSub = rho[s](a,atpos.size(),matSize, a,atpos.size(),matSize);
				Utot += prefac * trace(rhoSub - rhoSub*rhoSub).real();
				U_rho[s].set(a,atpos.size(),matSize, a,atpos.size(),matSize, prefac * (eye(mCount) - 2.*rhoSub));
				if(fpRhoAtom) //print atomic density matrices
				{	ostringstream oss; if(Uparams.n) oss << (Uparams.n + 1); oss << string("spdf")[Uparams.l];
					fprintf(fpRhoAtom, "%s  spin: %+d  orbital: %s  atom# %d\n",
						name.c_str(), nSpins-1-2*s, oss.str().c_str(), a+1);
					rhoSub.print_real(fpRhoAtom, "\t%+9.6lf");
					fprintf(fpRhoAtom, "\n");
				}
			}
		}
		//Propagate gradient from U_rho to wavefunctions or ionic positions if required:
		if(HC || forces)
		{	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
			{	int s = eInfo.qnums[q].index();
				ColumnBundle Opsi(C[q].similar(matSize));
				std::vector<ColumnBundle> dOpsi;
				setOpsi(Opsi, Uparams.n, Uparams.l);
				matrix CdagOpsi = C[q] ^ Opsi;
				if(HC) HC->at(q) += wSpinless * Opsi * (U_rho[s] * dagger(CdagOpsi)); //gradient upto state weight and fillings
				if(forces)
				{	diagMatrix fCartMat[3];
					for(int k=0; k<3; k++)
						fCartMat[k] = wSpinless * diag(U_rho[s] * dagger(CdagOpsi) * F[q] * (C[q]^D(Opsi,k)));
					for(unsigned a=0; a<atpos.size(); a++)
					{	vector3<> fCart; //proportional to Cartesian force
						for(int k=0; k<3; k++) fCart[k] = trace(fCartMat[k](a,atpos.size(),fCartMat[k].nRows()));
						(*forces)[a] += 2.*C[q].qnum->weight * (e->gInfo.RT * fCart);
					}
				}
			}
		}
	}
	watch.stop();
	return Utot;
}

void SpeciesInfo::updateLocal(DataGptr& Vlocps, DataGptr& rhoIon, DataGptr& nChargeball,
	DataGptr& nCore, DataGptr& tauCore) const
{	if(!atpos.size()) return; //unused species
	((SpeciesInfo*)this)->updateLatticeDependent(); //update lattice dependent quantities (if lattice vectors have changed)
	((SpeciesInfo*)this)->cachedV.clear(); //clear any cached projectors
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

void SpeciesInfo::accumNonlocalForces(const ColumnBundle& Cq, const matrix& VdagC, const matrix& E_VdagC, const matrix& grad_CdagOCq, std::vector<vector3<> >& forces) const
{	matrix DVdagC[3]; //cartesian gradient of VdagC
	{	auto V = getV(Cq);
		for(int k=0; k<3; k++)
			DVdagC[k] = D(*V,k)^Cq;
	}
	int nProj = MnlAll.nRows();
	//Loop over atoms:
	for(unsigned atom=0; atom<atpos.size(); atom++)
	{	matrix atomVdagC = VdagC(atom*nProj,(atom+1)*nProj, 0,E_VdagC.nCols());
		matrix E_atomVdagC = E_VdagC(atom*nProj,(atom+1)*nProj, 0,E_VdagC.nCols());
		if(QintAll) E_atomVdagC += QintAll * atomVdagC * grad_CdagOCq; //Contribution via overlap augmentation
		
		vector3<> fCart; //proportional to cartesian force
		for(int k=0; k<3; k++)
		{	matrix atomDVdagC = DVdagC[k](atom*nProj,(atom+1)*nProj, 0,E_VdagC.nCols());
			fCart[k] = trace(E_atomVdagC * dagger(atomDVdagC)).real();
		}
		forces[atom] += 2.*Cq.qnum->weight * (e->gInfo.RT * fCart);
	}
}

std::shared_ptr<ColumnBundle> SpeciesInfo::getV(const ColumnBundle& Cq) const
{	const QuantumNumber& qnum = *(Cq.qnum);
	const Basis& basis = *(Cq.basis);
	int nProj = MnlAll.nRows();
	if(!nProj) return 0; //purely local psp
	//First check cache
	if(e->cntrl.cacheProjectors)
	{	auto iter = cachedV.find(qnum.k);
		if(iter != cachedV.end()) //found
			return iter->second; //return cached value
	}
	//No cache / not found in cache; compute:
	std::shared_ptr<ColumnBundle> V = std::make_shared<ColumnBundle>(Cq.similar(nProj*atpos.size()));
	int iProj = 0;
	for(int l=0; l<int(VnlRadial.size()); l++)
		for(unsigned p=0; p<VnlRadial[l].size(); p++)
			for(int m=-l; m<=l; m++)
			{	size_t offs = iProj * basis.nbasis;
				size_t atomStride = nProj * basis.nbasis;
				callPref(Vnl)(basis.nbasis, atomStride, atpos.size(), l, m, qnum.k, basis.iGarrPref, basis.gInfo->G, atposPref, VnlRadial[l][p], V->dataPref()+offs);
				iProj++;
			}
	//Add to cache if necessary:
	if(e->cntrl.cacheProjectors)
		((SpeciesInfo*)this)->cachedV[qnum.k] = V;
	return V;
}
