/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman
Copyright 1996-2003 Sohrab Ismail-Beigi

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
#include <electronic/ColumnBundle.h>
#include <core/Util.h>
#include <fstream>
#include <sstream>


#ifdef GPU_ENABLED
void SpeciesInfo::sync_atposGpu()
{	if(!atpos.size()) return; //unused species
	//Transfer atomic positions to the GPU:
	cudaMemcpy(atposGpu, &atpos[0], sizeof(vector3<>)*atpos.size(), cudaMemcpyHostToDevice);
	gpuErrorCheck();
}
#endif

SpeciesInfo::SpeciesInfo()
{
	Z = 0.0;
	Z_chargeball = 0.0; width_chargeball = 0.0;
	dE_dnG = 0.0;
	mass = 0.0;

	lLocCpi = -1; //default: custom channel / highest l channel
	dGloc = 0.02; // default grid seperation for full G functions.
	dGnl  = 0.02; // default grid separation for reduced G operations
	pulayfilename ="none";
}

SpeciesInfo::~SpeciesInfo()
{	if(atpos.size())
	{
		#ifdef GPU_ENABLED
		cudaFree(atposGpu);
		#endif
		VlocRadial.free();
		nCoreRadial.free();
		for(auto& Vnl_l: VnlRadial) for(auto& Vnl_lp : Vnl_l) Vnl_lp.free();
		for(auto& Qijl: Qradial) Qijl.second.free();
		for(auto& proj_l: projRadial) for(auto& proj_lp: proj_l) proj_lp.free();
		for(auto& psi_l: psiRadial) for(auto& psi_lp: psi_l) psi_lp.free();
	}
}


void SpeciesInfo::setup(const Everything &everything)
{	e = &everything;
	if(!atpos.size()) return; //unused species
	
	//Read pseudopotential
	ifstream ifs(potfilename.c_str());
	if(!ifs.is_open()) die("Can't open pseudopotential file '%s' for reading.\n", potfilename.c_str());
	logPrintf("\nReading pseudopotential file '%s':\n",potfilename.c_str());
	switch(pspFormat)
	{	case Pot: readPot(ifs); break;
		case Cpi: readCpi(ifs); break;
		case Fhi: readFhi(ifs); break;
		case Uspp: readUspp(ifs); break;
	}
	
	setupPulay(); //Pulay info

	#ifdef GPU_ENABLED
	//Alloc and init GPU atomic positions:
	cudaMalloc(&atposGpu, sizeof(vector3<>)*atpos.size());
	sync_atposGpu();
	atposPref = atposGpu;
	#else
	atposPref = &atpos[0];
	#endif
}


void SpeciesInfo::print(FILE* fp) const
{	if(!atpos.size()) return; //unused species
	for(unsigned at=0; at<atpos.size(); at++)
	{	vector3<> pos = atpos[at]; //always in gInfo coordinates
		if(e->iInfo.coordsType == CoordsCartesian)
			pos = e->gInfo.R * pos; //convert to Cartesian coordinates
		fprintf(fp, "ion %s %19.15lf %19.15lf %19.15lf %lg\n",
			name.c_str(), pos[0], pos[1], pos[2], moveScale[at]);
	}
}

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

void SpeciesInfo::augmentOverlap(const ColumnBundle& Cq, ColumnBundle& OCq) const
{	static StopWatch watch("augmentOverlap"); watch.start();
	if(!atpos.size()) return; //unused species
	if(!Qint.size()) return; //no overlap augmentation
	const GridInfo &gInfo = e->gInfo;
	const Basis& basis = *Cq.basis;
	
	for(int l=0; l<int(VnlRadial.size()); l++)
	{	unsigned nProj = VnlRadial[l].size(); if(!nProj) continue; //skip l if no projectors
		//Copy Qint into a block-diagonal form for all atoms:
		tiledBlockMatrix Q(this->Qint[l], atpos.size());
		//Allocate temporaries:
		ColumnBundle V = Cq.similar(nProj * atpos.size());
		for(int m=-l; m<=l; m++)
		{	// Calculate the nonlocal projectors:
			for(unsigned p=0; p<nProj; p++)
			{	size_t offs = p * basis.nbasis;
				size_t atomStride = nProj * basis.nbasis;
				callPref(Vnl)(basis.nbasis, atomStride, atpos.size(), l, m, V.qnum->k, basis.iGarrPref, gInfo.G,
					atposPref, VnlRadial[l][p], V.dataPref()+offs, false, vector3<complex*>());
			}
			// Augment overlap:
			OCq += V * (Q * (V^Cq));
		}
	}
	watch.stop();
}


void SpeciesInfo::augmentDensity(const diagMatrix& Fq, const ColumnBundle& Cq, DataRptr& n) const
{	static StopWatch watch("augmentDensity"); watch.start();
	if(!atpos.size()) return; //unused species
	if(!Qint.size()) return; //no overlap augmentation
	const GridInfo &gInfo = e->gInfo;
	const Basis& basis = *Cq.basis;
	//Allocate temporaries:
	int nProj = 0;
	for(unsigned l=0; l<VnlRadial.size(); l++)
		nProj += (2*l+1)*VnlRadial[l].size();
	ColumnBundle V = Cq.similar(nProj);
	DataGptr Q(DataG::alloc(gInfo, isGpuEnabled()));
	DataGptr nAugTilde;
	//Loop over atoms:
	for(unsigned atom=0; atom<atpos.size(); atom++)
	{	//Initialize all projectors at this atom:
		int iProj = 0;
		for(int l=0; l<int(VnlRadial.size()); l++)
			for(unsigned p=0; p<VnlRadial[l].size(); p++)
				for(int m=-l; m<=l; m++)
				{	size_t offs = iProj * basis.nbasis;
					size_t atomStride = nProj * basis.nbasis;
					callPref(Vnl)(basis.nbasis, atomStride, 1, l, m, V.qnum->k, basis.iGarrPref, gInfo.G,
						atposPref+atom, VnlRadial[l][p], V.dataPref()+offs, false, vector3<complex*>());
					iProj++;
				}
		matrix VdagC = V ^ Cq;
		matrix Rho = VdagC * Fq * dagger(VdagC); //density matrix in projector basis
		
		int i1 = 0;
		//Triple loop over first projector:
		for(int l1=0; l1<int(VnlRadial.size()); l1++)
		for(int p1=0; p1<int(VnlRadial[l1].size()); p1++)
		for(int m1=-l1; m1<=l1; m1++)
		{	//Triple loop over second projector:
			int i2 = 0;
			for(int l2=0; l2<int(VnlRadial.size()); l2++)
			for(int p2=0; p2<int(VnlRadial[l2].size()); p2++)
			for(int m2=-l2; m2<=l2; m2++)
			{	if(i2<=i1) //rest handled by i1<->i2 symmetry 
				//NOTE: this => (l2,m2)<=(l1,m1) which is assumed in YlmProd implementation
					for(int l=abs(l1-l2); l<=l1+l2; l+=2) //total angular momentum
					{	QijIndex qIndex = { l1, p1, l2, p2, l };
						auto Qijl = Qradial.find(qIndex);
						if(Qijl==Qradial.end()) continue; //no entry at this l
						//Initialize Q function (returns false if channel is zero):
						if(callPref(Qr)(l1,m1, l2,m2, l, gInfo.S, gInfo.G, Qijl->second,
							atpos[atom], Q->dataPref(), 0, vector3<complex*>()))
						{	//Accumulate contribution in fourier space:
							nAugTilde += Q *
								( ((i1==i2 ? 1 : 2)/gInfo.detR)
								* (Rho.data()[Rho.index(i2,i1)] * cis(0.5*M_PI*(l2-l1))).real() );
						}
					}
				i2++;
			}
			i1++;
		}
	}
	n += Cq.qnum->weight * I(nAugTilde,true);
	watch.stop();
}

void SpeciesInfo::augmentDensityGrad(const diagMatrix& Fq, const ColumnBundle& Cq, const DataRptr& Vscloc,
	ColumnBundle& HCq, std::vector<vector3<> >* forces, const matrix& gradCdagOCq) const
{	static StopWatch watch("augmentDensityGrad"); watch.start();
	if(!atpos.size()) return; //unused species
	if(!Qint.size()) return; //no overlap augmentation
	const GridInfo &gInfo = e->gInfo;
	const Basis& basis = *Cq.basis;
	//Allocate temporaries:
	int nProj = 0;
	for(unsigned l=0; l<VnlRadial.size(); l++)
		nProj += (2*l+1)*VnlRadial[l].size();
	ColumnBundle V = Cq.similar(nProj);
	ColumnBundle dV[3];
	DataGptr gradAtpos[3];
	if(forces)
		for(int k=0; k<3; k++)
		{	dV[k] = V.similar();
			gradAtpos[k] = DataG::alloc(gInfo, isGpuEnabled());
		}
	DataGptr Q(DataG::alloc(gInfo, isGpuEnabled()));
	DataGptr IdagVscloc = Idag(Vscloc);
	//Loop over atoms:
	for(unsigned atom=0; atom<atpos.size(); atom++)
	{	vector3<complex*> dVdata, gradAtposData;
		if(forces)
		{	for(int k=0; k<3; k++)
			{	dVdata[k] = dV[k].dataPref();
				gradAtpos[k]->zero();
				gradAtposData[k] = gradAtpos[k]->dataPref();
			}
		}
		//Initialize all projectors at this atom:
		int iProj = 0;
		for(int l=0; l<int(VnlRadial.size()); l++)
			for(unsigned p=0; p<VnlRadial[l].size(); p++)
				for(int m=-l; m<=l; m++)
				{	size_t offs = iProj * basis.nbasis;
					size_t atomStride = nProj * basis.nbasis;
					callPref(Vnl)(basis.nbasis, atomStride, 1, l, m, V.qnum->k, basis.iGarrPref, gInfo.G,
						atposPref+atom, VnlRadial[l][p], V.dataPref()+offs, forces, dVdata+offs);
					iProj++;
				}
		matrix VdagC = V ^ Cq;
		matrix Rho = VdagC * Fq * dagger(VdagC); //density matrix in projector basis
		matrix Qint(nProj, nProj); Qint.zero(); //Full nProj x nProj version of this->Qint[l]
		matrix gradRho(nProj, nProj); gradRho.zero(); //gradient w.r.t density matrix in projector basis
		
		int i1 = 0;
		//Triple loop over first projector:
		for(int l1=0; l1<int(VnlRadial.size()); l1++)
		for(int p1=0; p1<int(VnlRadial[l1].size()); p1++)
		for(int m1=-l1; m1<=l1; m1++)
		{	//Triple loop over second projector:
			int i2 = 0;
			for(int l2=0; l2<int(VnlRadial.size()); l2++)
			for(int p2=0; p2<int(VnlRadial[l2].size()); p2++)
			for(int m2=-l2; m2<=l2; m2++)
			{	if(i2<=i1) //rest handled by i1<->i2 symmetry 
				{	//NOTE: this => (l2,m2)<=(l1,m1) which is assumed in YlmProd implementation
					for(int l=abs(l1-l2); l<=l1+l2; l+=2) //total angular momentum
					{	QijIndex qIndex = { l1, p1, l2, p2, l };
						auto Qijl = Qradial.find(qIndex);
						if(Qijl==Qradial.end()) continue; //no entry at this l
						//Initialize Q function (returns false if channel is zero):
						if(callPref(Qr)(l1,m1, l2,m2, l, gInfo.S, gInfo.G, Qijl->second,
							atpos[atom], Q->dataPref(), forces ? IdagVscloc->dataPref() : 0, gradAtposData,
							((i1==i2 ? 1 : 2)/gInfo.detR) * Cq.qnum->weight * 
								(Rho.data()[Rho.index(i2,i1)] * cis(0.5*M_PI*(l2-l1))).real() ))
						{	//Accumulate contribution to density matrix gradient:
							complex gradRho_ij = (dot(Q, IdagVscloc)/gInfo.detR) * cis(0.5*M_PI*(l2-l1));
							gradRho.data()[gradRho.index(i2,i1)] += gradRho_ij.conj();
							if(i1!=i2)
								gradRho.data()[gradRho.index(i1,i2)] += gradRho_ij;
						}
					}
				}
				if(l1==l2 && m1==m2)
				{	Qint.data()[Qint.index(i1,i2)] = this->Qint[l1].data()[this->Qint[l1].index(p1,p2)];
					Qint.data()[Qint.index(i2,i1)] = this->Qint[l1].data()[this->Qint[l1].index(p2,p1)];
				}
				i2++;
			}
			i1++;
		}
		matrix gradRhoVdagC = gradRho * VdagC;
		if(HCq) HCq += V * gradRhoVdagC;
		if(forces)
		{	for(int k=0; k<3; k++)
			{	matrix dVdagC = dV[k]^Cq;
				(*forces)[atom][k] -= sum(gradAtpos[k]) //Contribution via dQ
					+ 2.*Cq.qnum->weight *
						( trace(gradRhoVdagC * Fq * dagger(dVdagC)).real() //Contribution via dV
						+ trace(Qint * VdagC * gradCdagOCq * dagger(dVdagC)).real() ); //Contribution via overlap
			}
		}
	}
	watch.stop();
}

//Set atomic orbitals in column bundle from radial functions (almost same operation as setting Vnl)
void SpeciesInfo::setAtomicOrbitals(ColumnBundle& Y, int colOffset) const
{	if(!atpos.size()) return;
	const Basis& basis = *Y.basis;
	int iCol = colOffset; //current column
	for(int l=0; l<int(psiRadial.size()); l++)
		for(unsigned p=0; p<psiRadial[l].size(); p++)
		{	for(int m=-l; m<=l; m++)
			{	// Calculate the nonlocal projectors (and optionally their spatial derivatives):
				size_t offs = iCol * basis.nbasis;
				size_t atomStride = basis.nbasis;
				callPref(Vnl)(basis.nbasis, atomStride, atpos.size(), l, m, Y.qnum->k, basis.iGarrPref, e->gInfo.G,
					atposPref, psiRadial[l][p], Y.dataPref()+offs, false, vector3<complex*>());
				iCol += atpos.size();
			}
		}
}
int SpeciesInfo::nAtomicOrbitals() const
{	int nOrbitals = 0;
	for(int l=0; l<int(psiRadial.size()); l++)
		nOrbitals += (2*l+1)*psiRadial[l].size();
	return nOrbitals * atpos.size();
}

// Binary-write the PAW projector matrices (if any) for a particular state, looped over atoms, l and then m
void SpeciesInfo::writeProjectors(const ColumnBundle& Cq, FILE* fp) const
{	if(!atpos.size()) return; //unused species
	const GridInfo &gInfo = e->gInfo;
	const Basis& basis = *Cq.basis;
	
	//count up the number of projectors for each atom (multiple projectors for each l, m)
	int nProj = 0;
	for(int l=0; l<int(projRadial.size()); l++)
		nProj += (2*l+1) * projRadial[l].size();
	ColumnBundle proj = Cq.similar(nProj);
	
	for(unsigned at=0; at<atpos.size(); at++)
	{	// Calculate the projector matrix for each atom:
		int iProj = 0;
		for(int l=0; l<int(projRadial.size()); l++)
			for(int m=-l; m<=l; m++)
				for(unsigned p=0; p<projRadial[l].size(); p++)
				{	size_t offs = (iProj++) * basis.nbasis;
					size_t atomStride = basis.nbasis;
					callPref(Vnl)(basis.nbasis, atomStride, 1, l, m, proj.qnum->k, basis.iGarrPref, gInfo.G,
						atposPref+at, projRadial[l][p], proj.dataPref()+offs, false, vector3<complex*>());
				}
		// Save the projection:
		(proj^Cq).write(fp);
	}
}


void SpeciesInfo::updateLocal(DataGptr& Vlocps, DataGptr& rhoIon, DataGptr& nChargeball,
	DataGptr& nCore, DataGptr& tauCore) const
{	if(!atpos.size()) return; //unused species
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
		Z, e->iInfo.ionChargeWidth, nCoreRadial, tauCoreRadial, Z_chargeball, width_chargeball);
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
	DataGptr grad_SG(DataG::alloc(gInfo, isGpuEnabled())); //gradient w.r.t structure factor
	callPref(gradLocalToSG)(gInfo.S, gInfo.GGT,
		ccgrad_Vlocps->dataPref(), ccgrad_rhoIonData, ccgrad_nChargeballData,
		ccgrad_nCoreData, ccgrad_tauCoreData, grad_SG->dataPref(), VlocRadial,
		Z, e->iInfo.ionChargeWidth, nCoreRadial, tauCoreRadial, Z_chargeball, width_chargeball);
	
	//Now propagate that gradient to each atom of this species:
	DataGptr gradAtpos[3]; vector3<complex*> gradAtposData;
	for(int k=0; k<3; k++)
	{	gradAtpos[k] = DataG::alloc(gInfo, isGpuEnabled());
		gradAtposData[k] = gradAtpos[k]->dataPref();
	}
	std::vector< vector3<> > forces(atpos.size());
	for(unsigned at=0; at<atpos.size(); at++)
	{	callPref(gradSGtoAtpos)(gInfo.S, atpos[at], grad_SG->dataPref(), gradAtposData);
		for(int k=0; k<3; k++)
			forces[at][k] = -sum(gradAtpos[k]); //negative gradient
	}
	return forces;
}

bool SpeciesInfo::QijIndex::operator<(const SpeciesInfo::QijIndex& other) const
{	//Bring both indices to the upper triangular part:
	QijIndex q1 = *this; q1.sortIndices();
	QijIndex q2 = other; q2.sortIndices();
	//Compare:
	if(q1.l1 < q2.l1) return true;
	if(q1.l1 == q2.l1)
	{	if(q1.p1 < q2.p1) return true;
		if(q1.p1 == q2.p1)
		{	if(q1.l2 < q2.l2) return true;
			if(q1.l2 == q2.l2)
			{	if(q1.p2 < q2.p2) return true;
				if(q1.p2 == q2.p2)
					return q1.l < q2.l;
			}
		}
	}
	return false;
}
void SpeciesInfo::QijIndex::sortIndices()
{	if(l1>l2)
	{	std::swap(l1,l2);
		std::swap(p1,p2);
	}
	else if(l1==l2 && p1>p2)
	{	std::swap(p1,p2);
	}
}
