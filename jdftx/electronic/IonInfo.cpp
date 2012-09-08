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
#include <cstdio>
#include <cmath>
#include <core/Units.h>

#define MIN_ION_DISTANCE 1e-10

IonInfo::IonInfo()
{	shouldPrintForceComponents = false;
}

void IonInfo::setup(const Everything &everything)
{	e = &everything;

	// Check atomic positions for problems
	checkPositions();

	 //Force output in same coordinate system as input forces
	if(forcesOutputCoords==ForcesCoordsPositions)
		forcesOutputCoords = coordsType==CoordsLattice ? ForcesCoordsLattice : ForcesCoordsCartesian;


	logPrintf("\n---------- Setting up pseudopotentials ----------\n");

	//Determine maximum G extents for local and non-local pseudopotentials:
	GmaxNL = sqrt(2.0*e->cntrl.Ecut);
	GmaxLoc = 0.0;
	vector3<int> c;
	for(c[0]=-1; c[0]<=1; c[0]+=2) for(c[1]=-1; c[1]<=1; c[1]+=2) for(c[2]=-1; c[2]<=1; c[2]+=2)
	{	vector3<> f; for(int k=0; k<3; k++) f[k] = c[k]*(e->gInfo.S[k]/2);
		double G = sqrt(e->gInfo.GGT.metric_length_squared(f));
		if(G>GmaxLoc) GmaxLoc=G;
	}

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
	logPrintf("Width of ionic core gaussian charges set to %lg\n", ionWidth);
	
	// Call the species setup routines
	int nAtomsTot=0;
	for(auto sp: species)
	{	nAtomsTot += sp->atpos.size();
		sp->setup(*e);
	}
	if(!nAtomsTot) logPrintf("Warning: no atoms in the calculation.\n");
	
	if(ionWidth && (e->eVars.fluidType != FluidNone) &&
		(e->eVars.fluidParams.ionicConcentration || e->eVars.fluidParams.hSIons.size()))
		logPrintf("\nCorrection to mu due to finite nuclear width = %lg\n", ionWidthMuCorrection());
}

void IonInfo::printPositions(FILE* fp) const
{	fprintf(fp, "# Ionic positions in %s coordinates:\n", coordsType==CoordsLattice ? "lattice" : "cartesian");
	for(auto sp: species) sp->print(fp);
}

// Check for overlapping atoms
void IonInfo::checkPositions() const
{
	bool err = false;
	double sizetest = 0;
	vector3<> vtest[2];

	for(auto sp: species)
		for(unsigned n=0; n < sp->atpos.size(); n++)
		{
			vtest[0] = sp->atpos[n];
			for(auto sp1: species)
				for (unsigned n1 = ((sp1==sp) ? (n+1) : 0); n1 < sp1->atpos.size(); n1++)
				{
					vtest[1] = vtest[0] - sp1->atpos[n1];
					sizetest = dot(vtest[1], vtest[1]);
					if (sizetest < MIN_ION_DISTANCE)
					{	err = true;
						logPrintf("Atom# %d of species %s and atom# %d of species %s coincide.\n",
							n, sp->name.c_str(), n1, sp1->name.c_str());
					}
				}
		}
	if(err) die("Coincident atoms found, please check lattice and atom positions.\n");
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
	
	//----------- update Vlocps, rhoIon, nCore and nChargeball --------------
	initZero(Vlocps, gInfo);
	initZero(rhoIon, gInfo);
	if(nChargeball) nChargeball->zero();
	DataGptr nCoreTilde, tauCoreTilde;
	for(auto sp: species) //collect contributions to the above from all species
		sp->updateLocal(Vlocps, rhoIon, nChargeball, nCoreTilde, tauCoreTilde);
	//Add long-range part to Vlocps and smoothen rhoIon:
	Vlocps += (*e->coulomb)(rhoIon);
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
	ener.Exc_core = nCore ? -e->exCorr(nCore, 0, false, &tauCore) : 0.0;
	
	//Ewald sum:
	ener.Eewald = ewaldAndGrad();
	
	//Pulay corrections:
	double dEtot_dnG = 0.0; //derivative of Etot w.r.t nG  (G-vectors/unit volume)
	for(auto sp: species)
		dEtot_dnG += sp->atpos.size() * sp->dE_dnG;
	ener.Epulay = dEtot_dnG * 
		( sqrt(2.0)*pow(e->cntrl.Ecut,1.5)/(3.0*M_PI*M_PI) //ideal nG
		-  e->basis[0].nbasis/e->gInfo.detR ); //actual nG
	
	//Update totals:
	ener.updateTotals();
}


double IonInfo::ionicEnergyAndGrad(IonicGradient& forces) const
{	const ElecInfo &eInfo = e->eInfo;
	const ElecVars &eVars = e->eVars;
	
	//---------- Ewald forces ---------
	IonicGradient forcesEwald; forcesEwald.init(*this);
	ewaldAndGrad(&forcesEwald);
	e->symm.symmetrize(forcesEwald);
	forces = forcesEwald;
	if(shouldPrintForceComponents)
		forcesEwald.print(*e, globalLog, "forceEwald");
	
	//---------- local part: Vlocps, chargeball, partial core etc.: --------------
	//compute the complex-conjugate gradient w.r.t the relevant densities/potentials:
	const DataGptr ccgrad_Vlocps = J(eVars.get_nTot()); //just the electron density for Vlocps
	const DataGptr ccgrad_nChargeball = eVars.V_cavity; //cavity potential for chargeballs
	DataGptr ccgrad_rhoIon = (*e->coulomb)(ccgrad_Vlocps); //long-range portion of Vlocps for rhoIon
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
		if(e->exCorr.needsKEdensity())
		{	DataRptr VtauAvg = (Vtau.size()==1) ? Vtau[0] : 0.5*(Vtau[0]+Vtau[1]);
			ccgrad_tauCore += J(VtauAvg - VtauCore);
		}
	}
	//Propagate those gradients to forces:
	IonicGradient forcesLoc; forcesLoc.init(*this);
	for(unsigned sp=0; sp<species.size(); sp++)
		forcesLoc[sp] = species[sp]->getLocalForces(ccgrad_Vlocps, ccgrad_rhoIon,
			ccgrad_nChargeball, ccgrad_nCore, ccgrad_tauCore);
	e->symm.symmetrize(forcesLoc);
	forces += forcesLoc;
	if(shouldPrintForceComponents)
		forcesLoc.print(*e, globalLog, "forceLoc");
	
	//--------- Forces due to pseudopotential density augmentation ---------
	IonicGradient forcesAug; forcesAug.init(*this);
	ColumnBundle nullHC; //a null column bundle (since we don't want elec gradient)
	for(int q=0; q<eInfo.nStates; q++)
		augmentDensityGrad(eVars.F[q], eVars.C[q], eVars.Vscloc[eInfo.qnums[q].index()],
			nullHC, &forcesAug, eVars.grad_CdagOC[q]);
	e->symm.symmetrize(forcesAug);
	forces += forcesAug;
	if(shouldPrintForceComponents)
		forcesAug.print(*e, globalLog, "forceAug");
	
	//---------------- non-local pseudopot contribution --------------------
	IonicGradient forcesNL; forcesNL.init(*this);
	for(int q=0; q<eInfo.nStates; q++)
		EnlAndGrad(eVars.F[q], eVars.C[q], nullHC, &forcesNL);
	e->symm.symmetrize(forcesNL);
	forces += forcesNL;
	if(shouldPrintForceComponents)
		forcesNL.print(*e, globalLog, "forceNL");
	
	return relevantFreeEnergy(*e);
}

double IonInfo::EnlAndGrad(const diagMatrix& Fq, const ColumnBundle& Cq, ColumnBundle& HCq, IonicGradient* forces) const
{	double Enlq = 0.0;
	for(unsigned sp=0; sp<species.size(); sp++)
		Enlq += species[sp]->EnlAndGrad(Fq, Cq, HCq, forces ? &(*forces)[sp] : 0);
	return Enlq;
}

void IonInfo::augmentOverlap(const ColumnBundle& Cq, ColumnBundle& OCq) const
{	for(auto sp: species)
		sp->augmentOverlap(Cq, OCq);
}

void IonInfo::augmentDensity(const diagMatrix& Fq, const ColumnBundle& Cq, DataRptr& n) const
{	for(auto sp: species)
		sp->augmentDensity(Fq, Cq, n);
}

void IonInfo::augmentDensityGrad(const diagMatrix& Fq, const ColumnBundle& Cq, const DataRptr& Vscloc,
	ColumnBundle& HCq, IonicGradient* forces, const matrix& gradCdagOCq) const
{	for(unsigned sp=0; sp<species.size(); sp++)
		species[sp]->augmentDensityGrad(Fq, Cq, Vscloc, HCq, forces ? &forces->at(sp) : 0, gradCdagOCq);
}


double IonInfo::ewaldAndGrad(IonicGradient* forces) const
{
	//Obtain the list of atomic positions and charges:
	std::vector<Coulomb::PointCharge> pointCharges;
	for(auto sp: species)
		for(const vector3<>& pos: sp->atpos)
			pointCharges.push_back({sp->Z, pos, vector3<>(0.,0.,0.)});
	//Compute Ewald sum and gradients
	double E = e->coulomb->energyAndGrad(pointCharges);
	//Store forces if requested:
	if(forces)
	{	auto pc = pointCharges.begin();
		for(unsigned sp=0; sp<species.size(); sp++)
			for(unsigned at=0; at<species[sp]->atpos.size(); at++)
				(*forces)[sp][at] = (pc++)->force;
	}
	return E;
}
