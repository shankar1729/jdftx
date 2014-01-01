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

#include <electronic/IonInfo.h>
#include <electronic/Everything.h>
#include <electronic/SpeciesInfo.h>
#include <electronic/ExCorr.h>
#include <electronic/ColumnBundle.h>
#include <electronic/VanDerWaals.h>
#include <electronic/operators.h>
#include <fluid/FluidSolver.h>
#include <cstdio>
#include <cmath>
#include <core/Units.h>

#define MIN_ION_DISTANCE 1e-10

IonInfo::IonInfo()
{	shouldPrintForceComponents = false;
	vdWenable = false;
	vdWscale = 0.;
}

void IonInfo::setup(const Everything &everything)
{	e = &everything;

	//Force output in same coordinate system as input forces
	if(forcesOutputCoords==ForcesCoordsPositions)
		forcesOutputCoords = coordsType==CoordsLattice ? ForcesCoordsLattice : ForcesCoordsCartesian;

	logPrintf("\n---------- Setting up pseudopotentials ----------\n");
		
	//Choose width of the nuclear gaussian:
	switch(ionWidthMethod)
	{	case IonWidthManual: break; //manually specified value
		case IonWidthFFTbox:
			//set to a multiple of the maximum grid spacing
			ionWidth = 0.0;
			for (int i=0; i<3; i++)
			{	double dRi = e->gInfo.h[i].length();
				ionWidth = std::max(ionWidth, 1.6*dRi);
			}
			break;
		case IonWidthEcut:
			ionWidth = 0.8*M_PI / sqrt(2*e->cntrl.Ecut);
			break;
	}
	logPrintf("Width of ionic core gaussian charges (only for fluid interactions / plotting) set to %lg\n", ionWidth);
	
	// Call the species setup routines
	int nAtomsTot=0;
	for(auto sp: species)
	{	nAtomsTot += sp->atpos.size();
		sp->setup(*e);
	}
	if(!nAtomsTot) logPrintf("Warning: no atoms in the calculation.\n");
	
	if(not checkPositions())
		die("\nAtoms are too close, have overlapping pseudopotential cores.\n\n");
	
	if(ionWidth && (e->eVars.fluidParams.fluidType != FluidNone) && e->eVars.fluidParams.ionicScreening())
		logPrintf("\nCorrection to mu due to finite nuclear width = %lg\n", ionWidthMuCorrection());
}

void IonInfo::printPositions(FILE* fp) const
{	fprintf(fp, "# Ionic positions in %s coordinates:\n", coordsType==CoordsLattice ? "lattice" : "cartesian");
	for(auto sp: species) sp->print(fp);
}

// Check for overlapping atoms, returns true if okay
bool IonInfo::checkPositions() const
{	bool okay = true;
	double sizetest = 0;
	vector3<> vtest[2];

	for(auto sp: species)
		for(unsigned n=0; n < sp->atpos.size(); n++)
		{	if(sp->coreRadius == 0.) continue;
			vtest[0] = sp->atpos[n];
			for(auto sp1: species)
			{	if(sp1->coreRadius == 0.) continue;
				for (unsigned n1 = ((sp1==sp) ? (n+1) : 0); n1 < sp1->atpos.size(); n1++)
				{	vtest[1] = vtest[0] - sp1->atpos[n1];
					for(int i=0; i<3; i++) // Get periodic distance
						vtest[1][i] = vtest[1][i] - floor(vtest[1][i] + 0.5);
					sizetest = sqrt(dot(e->gInfo.R*vtest[1], e->gInfo.R*vtest[1]));
					if (coreOverlapCondition==additive and (sizetest < (sp->coreRadius + sp1->coreRadius)))
					{	logPrintf("\nWARNING: %s #%d and %s #%d are closer than the sum of their core radii.",
							sp->name.c_str(), n, sp1->name.c_str(), n1);
						okay = false;
					}
					else if (coreOverlapCondition==vector and (sizetest < sqrt(pow(sp->coreRadius, 2) + pow(sp1->coreRadius, 2))))
					{	logPrintf("\nWARNING: %s #%d and %s #%d are closer than the vector-sum of their core radii.",
							sp->name.c_str(), n, sp1->name.c_str(), n1);
						okay = false;
					}
					else if(sizetest < MIN_ION_DISTANCE)
					{	die("\nERROR: Ions %s #%d and %s #%d are on top of eachother.\n\n", sp->name.c_str(), n, sp1->name.c_str(), n1);
					}
				}
			}
		}
		
	if(not okay) // Add another line after printing core overlap warnings
		logPrintf("\n");
		
	return okay;
}

double IonInfo::getZtot() const
{	double Ztot=0.;
	for(auto sp: species)
		Ztot += sp->Z * sp->atpos.size();
	return Ztot;
}

double IonInfo::ionWidthMuCorrection() const
{	return (4*M_PI/ e->gInfo.detR) * (-0.5*ionWidth*ionWidth) * getZtot();
}


void IonInfo::update(Energies& ener)
{	const GridInfo &gInfo = e->gInfo;

	if(not checkPositions())
		ener.E["core-overlap"] = NAN;
	else
		ener.E["core-overlap"] = 0.;

	//----------- update Vlocps, rhoIon, nCore and nChargeball --------------
	initZero(Vlocps, gInfo);
	initZero(rhoIon, gInfo);
	if(nChargeball) nChargeball->zero();
	DataGptr nCoreTilde, tauCoreTilde;
	for(auto sp: species) //collect contributions to the above from all species
		sp->updateLocal(Vlocps, rhoIon, nChargeball, nCoreTilde, tauCoreTilde);
	//Add long-range part to Vlocps and smoothen rhoIon:
	Vlocps += (*e->coulomb)(rhoIon, Coulomb::PointChargeRight);
	rhoIon = gaussConvolve(rhoIon, ionWidth);
	//Process partial core density:
	if(nCoreTilde) nCore = I(nCoreTilde, true); // put in real space
	if(tauCoreTilde)
	{	tauCore = I(tauCoreTilde, true); // put in real space
		double tauMin, tauMax; //cap negative KE densities at 0:
		callPref(eblas_capMinMax)(gInfo.nr, tauCore->dataPref(), tauMin, tauMax, 0.);
	}
	
	//---------- energies dependent on ionic positions alone ----------------
	
	//Energies due to partial electronic cores:
	ener.E["Exc_core"] = nCore ? -e->exCorr(nCore, 0, false, &tauCore) : 0.0;
	
	//Energies from pair-potential terms (Ewald etc.):
	pairPotentialsAndGrad(&ener);
	
	//Pulay corrections:
	double dEtot_dnG = 0.0; //derivative of Etot w.r.t nG  (G-vectors/unit volume)
	for(auto sp: species)
		dEtot_dnG += sp->atpos.size() * sp->dE_dnG;
	
	double nbasisAvg = 0.0;
	for(int q=e->eInfo.qStart; q<e->eInfo.qStop; q++)
		nbasisAvg += 0.5*e->eInfo.qnums[q].weight * e->basis[q].nbasis;
	mpiUtil->allReduce(nbasisAvg, MPIUtil::ReduceSum);
	
	ener.E["Epulay"] = dEtot_dnG * 
		( sqrt(2.0)*pow(e->cntrl.Ecut,1.5)/(3.0*M_PI*M_PI) //ideal nG
		-  nbasisAvg/e->gInfo.detR ); //actual nG
}


double IonInfo::ionicEnergyAndGrad(IonicGradient& forces) const
{	const ElecInfo &eInfo = e->eInfo;
	const ElecVars &eVars = e->eVars;
	
	//---------- Forces from pair potential terms (Ewald etc.) ---------
	IonicGradient forcesPairPot; forcesPairPot.init(*this);
	pairPotentialsAndGrad(0, &forcesPairPot);
	e->symm.symmetrize(forcesPairPot);
	forces = forcesPairPot;
	if(shouldPrintForceComponents)
		forcesPairPot.print(*e, globalLog, "forcePairPot");
	
	//---------- local part: Vlocps, chargeball, partial core etc.: --------------
	//compute the complex-conjugate gradient w.r.t the relevant densities/potentials:
	const DataGptr ccgrad_Vlocps = J(eVars.get_nTot()); //just the electron density for Vlocps
	const DataGptr ccgrad_nChargeball = eVars.V_cavity; //cavity potential for chargeballs
	DataGptr ccgrad_rhoIon = (*e->coulomb)(ccgrad_Vlocps, Coulomb::PointChargeLeft); //long-range portion of Vlocps for rhoIon
	if(eVars.d_fluid) //and electrostatic potential due to fluid (if any):
		ccgrad_rhoIon += gaussConvolve(eVars.d_fluid, ionWidth);
	DataGptr ccgrad_nCore, ccgrad_tauCore;
	if(nCore) //cavity potential and exchange-correlation coupling to electron density for partial cores:
	{	ccgrad_nCore = eVars.V_cavity;
		//Exchange-correlation:
		DataRptr VxcCore, VtauCore;
		DataRptrCollection Vxc(eVars.n.size()), Vtau;
		e->exCorr(nCore, &VxcCore, false, &tauCore, &VtauCore);
		e->exCorr(eVars.get_nXC(), &Vxc, false, &eVars.tau, &Vtau);
		DataRptr VxcAvg = (Vxc.size()==1) ? Vxc[0] : 0.5*(Vxc[0]+Vxc[1]); //spin-avgd potential
		ccgrad_nCore += J(VxcAvg - VxcCore);
		//Contribution through tauCore (metaGGAs only):
		if(e->exCorr.needsKEdensity() || e->eVars.fluidParams.useTau)
		{	DataRptr VtauAvg = (Vtau.size()==1) ? Vtau[0] : 0.5*(Vtau[0]+Vtau[1]);
			if(VtauAvg) ccgrad_tauCore += J(VtauAvg - VtauCore);
		}
	}
	//Propagate those gradients to forces:
	IonicGradient forcesLoc; forcesLoc.init(*this);
	for(unsigned sp=0; sp<species.size(); sp++)
		forcesLoc[sp] = species[sp]->getLocalForces(ccgrad_Vlocps, ccgrad_rhoIon,
			ccgrad_nChargeball, ccgrad_nCore, ccgrad_tauCore);
	forcesLoc += e->eVars.fluidForces; //include extra fluid forces (if any)
	e->symm.symmetrize(forcesLoc);
	forces += forcesLoc;
	if(shouldPrintForceComponents)
		forcesLoc.print(*e, globalLog, "forceLoc");
	
	//--------- Forces due to nonlocal pseudopotential contributions ---------
	IonicGradient forcesNL; forcesNL.init(*this);
	computeU(eVars.F, eVars.C, 0, &forcesNL); //Include DFT+U contribution if any
	augmentDensityGridGrad(eVars.Vscloc, &forcesNL);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++)
	{	const QuantumNumber& qnum = e->eInfo.qnums[q];
		//Collect gradients with respect to VdagCq (not including fillings and state weight):
		std::vector<matrix> HVdagCq(species.size()); 
		EnlAndGrad(qnum, eVars.F[q], eVars.VdagC[q], HVdagCq);
		augmentDensitySphericalGrad(qnum, eVars.F[q], eVars.VdagC[q], HVdagCq);
		//Propagate to atomic positions:
		for(unsigned sp=0; sp<species.size(); sp++) if(HVdagCq[sp])
			species[sp]->accumNonlocalForces(eVars.C[q], eVars.VdagC[q][sp], HVdagCq[sp]*eVars.F[q], eVars.grad_CdagOC[q], forcesNL[sp]);
	}
	for(auto& force: forcesNL) //Accumulate contributions over processes
		mpiUtil->allReduce((double*)force.data(), 3*force.size(), MPIUtil::ReduceSum);
	e->symm.symmetrize(forcesNL);
	forces += forcesNL;
	if(shouldPrintForceComponents)
		forcesNL.print(*e, globalLog, "forceNL");
	
	return relevantFreeEnergy(*e);
}

double IonInfo::EnlAndGrad(const QuantumNumber& qnum, const diagMatrix& Fq, const std::vector<matrix>& VdagCq, std::vector<matrix>& HVdagCq) const
{	double Enlq = 0.0;
	for(unsigned sp=0; sp<species.size(); sp++)
		Enlq += species[sp]->EnlAndGrad(qnum, Fq, VdagCq[sp], HVdagCq[sp]);
	return Enlq;
}

void IonInfo::augmentOverlap(const ColumnBundle& Cq, ColumnBundle& OCq, std::vector<matrix>* VdagCq) const
{	if(VdagCq) VdagCq->resize(species.size());
	for(unsigned sp=0; sp<species.size(); sp++)
		species[sp]->augmentOverlap(Cq, OCq, VdagCq ? &VdagCq->at(sp) : 0);
}

void IonInfo::augmentDensityInit() const
{	for(auto sp: species) ((SpeciesInfo&)(*sp)).augmentDensityInit();
}
void IonInfo::augmentDensitySpherical(const QuantumNumber& qnum, const diagMatrix& Fq, const std::vector<matrix>& VdagCq) const
{	for(unsigned sp=0; sp<species.size(); sp++)
		((SpeciesInfo&)(*species[sp])).augmentDensitySpherical(qnum, Fq, VdagCq[sp]);
}
void IonInfo::augmentDensityGrid(DataRptrCollection& n) const
{	for(auto sp: species) sp->augmentDensityGrid(n);
}
void IonInfo::augmentDensityGridGrad(const DataRptrCollection& E_n, IonicGradient* forces) const
{	for(unsigned sp=0; sp<species.size(); sp++)
		((SpeciesInfo&)(*species[sp])).augmentDensityGridGrad(E_n, forces ? &forces->at(sp) : 0);
}
void IonInfo::augmentDensitySphericalGrad(const QuantumNumber& qnum, const diagMatrix& Fq, const std::vector<matrix>& VdagCq, std::vector<matrix>& HVdagCq) const
{	for(unsigned sp=0; sp<species.size(); sp++)
		species[sp]->augmentDensitySphericalGrad(qnum, Fq, VdagCq[sp], HVdagCq[sp]);
}

void IonInfo::project(const ColumnBundle& Cq, std::vector<matrix>& VdagCq, matrix* rotExisting) const
{	VdagCq.resize(species.size());
	for(unsigned sp=0; sp<e->iInfo.species.size(); sp++)
	{	if(rotExisting && VdagCq[sp]) VdagCq[sp] = VdagCq[sp] * (*rotExisting); //rotate and keep the existing projections
		else
		{	auto V = e->iInfo.species[sp]->getV(Cq);
			if(V) VdagCq[sp] = (*V) ^ Cq;
		}
	}
}

void IonInfo::projectGrad(const std::vector<matrix>& HVdagCq, const ColumnBundle& Cq, ColumnBundle& HCq) const
{	for(unsigned sp=0; sp<species.size(); sp++)
		if(HVdagCq[sp]) HCq += *(species[sp]->getV(Cq)) * HVdagCq[sp];
}

double IonInfo::computeU(const std::vector<diagMatrix>& F, const std::vector<ColumnBundle>& C,
	std::vector<ColumnBundle>* HC, IonicGradient* forces, FILE* fpRhoAtom) const
{	double U = 0.;
	for(unsigned sp=0; sp<species.size(); sp++)
		U += species[sp]->computeU(F, C, HC, forces ? &forces->at(sp) : 0, fpRhoAtom);
	return U;
}

int IonInfo::nAtomicOrbitals() const
{	int nAtomic = 0;
	for(auto sp: species)
		nAtomic += sp->nAtomicOrbitals();
	return nAtomic;
}

ColumnBundle IonInfo::getAtomicOrbitals(int q, int extraCols) const
{	ColumnBundle psi(nAtomicOrbitals()+extraCols, e->basis[q].nbasis, &e->basis[q], &e->eInfo.qnums[q], isGpuEnabled());
	int iCol=0;
	for(auto sp: species)
	{	sp->setAtomicOrbitals(psi, iCol);
		iCol += sp->nAtomicOrbitals();
	}
	return psi;
}


void IonInfo::pairPotentialsAndGrad(Energies* ener, IonicGradient* forces) const
{
	//Obtain the list of atomic positions and charges:
	std::vector<Atom> atoms;
	for(auto sp: species)
		for(const vector3<>& pos: sp->atpos)
			atoms.push_back(Atom(sp->Z, pos, vector3<>(0.,0.,0.), sp->atomicNumber));
	//Compute Ewald sum and gradients (this also moves each Atom::pos into fundamental zone)
	double Eewald = e->coulomb->energyAndGrad(atoms);
	//Compute optional pair-potential terms:
	double EvdW = 0.;
	if(vdWenable)
	{	double scaleFac = e->vanDerWaals->getScaleFactor(e->exCorr.getName(), vdWscale);
		EvdW = e->vanDerWaals->energyAndGrad(atoms, scaleFac); //vanDerWaals energy+force
	}
	//Store energies and/or forces if requested:
	if(ener)
	{	ener->E["Eewald"] = Eewald;
		ener->E["EvdW"] = EvdW;
	}
	if(forces)
	{	auto atom = atoms.begin();
		for(unsigned sp=0; sp<species.size(); sp++)
			for(unsigned at=0; at<species[sp]->atpos.size(); at++)
				(*forces)[sp][at] = (atom++)->force;
	}
	
}
