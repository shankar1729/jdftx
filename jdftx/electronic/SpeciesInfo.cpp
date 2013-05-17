/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman, Deniz Gunceler
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
#include <electronic/operators.h>
#include <electronic/symbols.h>
#include <electronic/ColumnBundle.h>
#include <core/LatticeUtils.h>
#include <core/DataMultiplet.h>
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

inline bool isParallel(vector3<> x, vector3<> y)
{	return fabs(1.-fabs(dot(x, y)/(x.length() * y.length()))) < symmThreshold;
}

bool SpeciesInfo::Constraint::isEquivalent(const Constraint& otherConstraint, const matrix3<>& transform) const
{ 	if(moveScale != otherConstraint.moveScale) return false; //Ensure same moveSCale
	if(type != otherConstraint.type) return false; //Ensure same constraint type
	return (type==None) or isParallel(transform*d, otherConstraint.d);
}

int SpeciesInfo::Constraint::getDimension() const
{	if(not moveScale) return 0;
	switch(type)
	{	case Linear: return 1;
		case Planar: return 2;
		default: return 3;
	}
}

void SpeciesInfo::Constraint::print(FILE* fp, const Everything& e) const
{	vector3<> d = this->d; //in cartesian coordinates
	if(e.iInfo.coordsType == CoordsLattice)
		d = ~(e.gInfo.R) * d; //print in lattice coordinates
	fprintf(fp, "  %s %.14lg %.14lg %.14lg", constraintTypeMap.getString(type), d[0], d[1], d[2]);
}

vector3<> SpeciesInfo::Constraint::operator()(const vector3<>& grad)
{	vector3<> scaledGrad = moveScale * grad;
	switch(type)
	{	case Linear: return dot(scaledGrad, d)*d/ d.length_squared();
		case Planar: return  scaledGrad - dot(scaledGrad, d)*d/d.length_squared();	
		default: return scaledGrad;
	}
}


SpeciesInfo::SpeciesInfo()
{
	Z = 0.0;
	atomicNumber = 0;
	Z_chargeball = 0.0; width_chargeball = 0.0;
	dE_dnG = 0.0;
	mass = 0.0;
	coreRadius = 0.;
	
	lLocCpi = -1; //default: custom channel / highest l channel
	dGloc = 0.02; // default grid seperation for full G functions.
	dGnl  = 0.02; // default grid separation for reduced G operations
	pulayfilename ="none";
	OpsiRadial = &psiRadial;
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
		if(OpsiRadial != &psiRadial)
		{	for(auto& Opsi_l: *OpsiRadial) for(auto& Opsi_lp: Opsi_l) Opsi_lp.free();
			delete OpsiRadial;
		}
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
	estimateAtomEigs();
	if(coreRadius) logPrintf("  Core radius for overlap checks: %.2lf bohrs.\n", coreRadius);
	else if(!VnlRadial.size()) logPrintf("  Disabling overlap check for local pseudopotential.\n");
	else logPrintf("  Warning: could not determine core radius; disabling overlap check for this species.\n");
	setupPulay(); //Pulay info

	//Check for augmentation:
	if(Qint.size())
	{	bool needKE = e->exCorr.needsKEdensity();
		bool needEXX = (e->exCorr.exxFactor()!=0.);
		for(auto exCorr: e->exCorrDiff)
		{	needKE |= e->exCorr.needsKEdensity();
			needEXX |= (e->exCorr.exxFactor()!=0.);
		}
		if(needKE || needEXX)
			die("\nUltrasoft pseudopotentials do not currently support meta-GGA or hybrid functionals.\n");
	}
	
	//Generate atomic number from symbol, if not stored in pseudopotential:
	if(!atomicNumber)
	{	AtomicSymbol atSym;
		if(!atomicSymbolMap.getEnum(name.c_str(), atSym))
			die("\nCould not determine atomic number for species '%s'.\n"
				"Either use a pseudopotential which contains this information,\n"
				"or set the species name to be the chemical symbol for that atom type.\n", name.c_str());
		atomicNumber = int(atSym);
	}
	
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
		fprintf(fp, "ion %s %19.15lf %19.15lf %19.15lf %lg",
			name.c_str(), pos[0], pos[1], pos[2], constraints[at].moveScale);
		if(constraints[at].type != Constraint::None)
			constraints[at].print(fp, *e);
		fprintf(fp, "\n");
	}
	
	//Output magnetic moments for spin-polarized calculations:
	//--- Only supported for pseudopotentials with atomic orbitals
	if(e->eInfo.spinType == SpinZ && OpsiRadial->size())
	{	diagMatrix M(atpos.size(), 0.); //magnetic moments
		for(int q=0; q<e->eInfo.nStates; q++) //states
		{	const ColumnBundle& Cq = e->eVars.C[q];
			const diagMatrix& Fq = e->eVars.F[q];
			ColumnBundle Opsi = Cq.similar(atpos.size()); //space for atomic orbitals
			const QuantumNumber& qnum = *(Cq.qnum);
			const Basis& basis = *(Cq.basis);
			for(int l=0; l<int(OpsiRadial->size()); l++) //angular momenta
				for(const RadialFunctionG& curOpsiRadial: OpsiRadial->at(l)) //principal quantum number (shells of pseudo-atom)
					for(int m=-l; m<=+l; m++) //angular momentum directions
					{	//Compute the atomic orbitals:
						callPref(Vnl)(basis.nbasis, basis.nbasis, atpos.size(), l, m, qnum.k, basis.iGarrPref, e->gInfo.G,
							atposPref, curOpsiRadial, Opsi.dataPref(), false, vector3<complex*>());
						//Accumulate electron counts:
						matrix CdagOpsi = Cq ^ Opsi;
						M += (qnum.spin * qnum.weight) * diag(dagger(CdagOpsi) * Fq * CdagOpsi);
					}
		}
		fprintf(fp, "# magnetic-moments %s", name.c_str());
		for(double m: M) fprintf(fp, " %+lg", m);
		fprintf(fp, "\n");
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
		
		//Triple loop over first projector:
		int i1 = 0;
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


void SpeciesInfo::setOpsi(ColumnBundle& Opsi, unsigned n, int l, std::vector<ColumnBundle>* dOpsi) const
{	if(!atpos.size()) return;
	assert(Opsi.basis); assert(Opsi.qnum);
	assert((2*l+1)*int(atpos.size()) <= Opsi.nCols());
	const Basis& basis = *Opsi.basis;
	vector3<complex*> dOpsiData;
	if(dOpsi)
	{	dOpsi->resize(3);
		for(int k=0; k<3; k++)
		{	if(!dOpsi->at(k)) dOpsi->at(k) = Opsi.similar();
			assert((2*l+1)*int(atpos.size()) <= dOpsi->at(k).nCols());
			dOpsiData[k] = dOpsi->at(k).dataPref();
		}
	}
	int iCol = 0; //current column
	for(int m=-l; m<=l; m++)
	{	//Set atomic orbitals for all atoms at specified (n,l,m):
		size_t offs = iCol * basis.nbasis;
		size_t atomStride = basis.nbasis;
		callPref(Vnl)(basis.nbasis, atomStride, atpos.size(), l, m, Opsi.qnum->k, basis.iGarrPref, e->gInfo.G,
			atposPref, OpsiRadial->at(l)[n], Opsi.dataPref()+offs, dOpsi, dOpsiData+offs);
		iCol += atpos.size();
	}
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

void SpeciesInfo::accumulateAtomicDensity(DataGptrCollection& nTilde) const
{	//Collect list of distinct magnetizations and corresponding atom positions:
	std::map<double, std::vector< vector3<> > > Matpos;
	if(initialMagneticMoments.size())
	{	for(unsigned a=0; a<atpos.size(); a++)
			Matpos[initialMagneticMoments[a]].push_back(atpos[a]);
	}
	else Matpos[0.] = atpos;
	//Scratch space for atom positions on GPU:
	#ifdef GPU_ENABLED
	vector3<>* atposCurPref;
	cudaMalloc(&atposCurPref, sizeof(vector3<>)*atpos.size());
	#endif
	//For each magnetization:
	for(auto mapEntry: Matpos)
	{	double M = mapEntry.first; if(e->eInfo.spinType == SpinNone) assert(!M);
		const std::vector< vector3<> >& atposCur = mapEntry.second;
		#ifdef GPU_ENABLED
		cudaMemcpy(atposCurPref, atposCur.data(), atposCur.size()*sizeof(vector3<>), cudaMemcpyHostToDevice);
		#else
		const vector3<>* atposCurPref = atposCur.data();
		#endif
		//Compute structure factor for atoms with current magnetization:
		DataGptr SG; nullToZero(SG, e->gInfo);
		callPref(getSG)(e->gInfo.S, atposCur.size(), atposCurPref, 1./e->gInfo.detR, SG->dataPref());
		//Collect contributions:
		for(int s=0; s<(M ? 2 : 1); s++)
		{	RadialFunctionG nRadial;
			getAtom_nRadial(s, M, nRadial);
			DataGptr nTilde_s = nRadial * SG; //contribution per spin 
			nRadial.free();
			nTilde[s] += nTilde_s;
			if(!M) nTilde.back() += nTilde_s; //same contribution to both spins (irrespective of spinType)
		}
	}
	//Cleanup:
	#ifdef GPU_ENABLED
	cudaFree(atposCurPref);
	#endif
}

//Set atomic orbitals in column bundle from radial functions (almost same operation as setting Vnl)
void SpeciesInfo::setAtomicOrbitals(ColumnBundle& Y, int colOffset, std::vector<AtomConfig>* atomConfig) const
{	if(!atpos.size()) return;
	assert(Y.basis); assert(Y.qnum);
	//Check sizes:
	int colMax = colOffset; for(int l=0; l<int(psiRadial.size()); l++) colMax += psiRadial[l].size() * (2*l+1) * atpos.size();
	assert(colMax <= Y.nCols());
	if(atomConfig && atomConfig->size()<size_t(colMax)) atomConfig->resize(colMax);
	//Set orbitals and associated info if requested:
	const Basis& basis = *Y.basis;
	int iCol = colOffset; //current column
	for(int l=0; l<int(psiRadial.size()); l++)
		for(unsigned p=0; p<psiRadial[l].size(); p++)
		{	for(int m=-l; m<=l; m++)
			{	//Set atomic orbitals for all atoms at specified (n,l,m):
				size_t offs = iCol * basis.nbasis;
				size_t atomStride = basis.nbasis;
				callPref(Vnl)(basis.nbasis, atomStride, atpos.size(), l, m, Y.qnum->k, basis.iGarrPref, e->gInfo.G,
					atposPref, psiRadial[l][p], Y.dataPref()+offs, false, vector3<complex*>());
				if(atomConfig)
					for(unsigned a=0; a<atpos.size(); a++)
					{	AtomConfig& ac = atomConfig->at(iCol + a);
						ac.n = p;
						ac.l = l;
						ac.iAtom = a;
						//ac.F = atomF[l][p];
					}
				iCol += atpos.size();
			}
		}
}
void SpeciesInfo::setAtomicOrbital(ColumnBundle& Y, int col, unsigned iAtom, unsigned n, int l, int m) const
{	//Check inputs:
	assert(Y.basis); assert(Y.qnum);
	assert(col >= 0); assert(col < Y.nCols());
	assert(iAtom < atpos.size());
	assert(l >= 0); assert(unsigned(l) < psiRadial.size());
	assert(n < psiRadial[l].size());
	assert(m >= -l); assert(m <= l);
	//Call Vnl() to set the column:
	const Basis& basis = *Y.basis;
	size_t offs = col * basis.nbasis;
	size_t atomStride = basis.nbasis;
	callPref(Vnl)(basis.nbasis, atomStride, 1, l, m, Y.qnum->k, basis.iGarrPref, e->gInfo.G,
		atposPref+iAtom, psiRadial[l][n], Y.dataPref()+offs, false, vector3<complex*>());
}
int SpeciesInfo::nAtomicOrbitals() const
{	int nOrbitals = 0;
	for(int l=0; l<int(psiRadial.size()); l++)
		nOrbitals += (2*l+1)*psiRadial[l].size();
	return nOrbitals * atpos.size();
}
int SpeciesInfo::lMaxAtomicOrbitals() const
{	return int(psiRadial.size()) - 1;
}
int SpeciesInfo::nAtomicOrbitals(int l) const
{	assert(l >= 0); assert(unsigned(l) < psiRadial.size());
	return psiRadial[l].size();
}
int SpeciesInfo::atomicOrbitalOffset(unsigned int iAtom, unsigned int n, int l, int m) const
{	assert(iAtom < atpos.size());
	assert(l >= 0); assert(unsigned(l) < psiRadial.size());
	assert(n < psiRadial[l].size());
	assert(m >= -l); assert(m <= l);
	int iProj = l + m; //#projectors before this one at current l,n
	for(int L=0; L<=l; L++) //#projectors from previous l,n:
		iProj += (L==l ? n : psiRadial[l].size()) * (2*L+1);
	return iProj * atpos.size() + iAtom;
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
