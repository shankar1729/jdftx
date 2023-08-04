/*
 * Copyright 2023 <copyright holder> <email>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <electronic/SpringConstant.h>
#include <electronic/IonInfo.h>
#include <electronic/Energies.h>
#include <electronic/ElecVars.h>
#include <electronic/PerturbationSolver.h>
#include <electronic/Everything.h>

SpringConstant::SpringConstant(Everything& e) : e(e), eVars(e.eVars), eInfo(e.eInfo), iInfo(e.iInfo), pInfo(e.vptInfo), ps(e) {}

void SpringConstant::setupModes() {
	vector3<> x(1,0,0), y(0,1,0), z(0,0,1);
	
	modes.clear();
	for (unsigned s = 0; s < iInfo.species.size(); s++) {
		auto sp = iInfo.species[s];
		for (unsigned i = 0; i < sp->atpos.size(); i++) {
			if (sp->perturbed[i]) {
				PerturbationInfo::Mode mx, my, mz;
				mx.sp = s; my.sp = s; mz.sp = s;
				mx.at = i; my.at = i; mz.at = i;
				mx.dirCartesian = x; my.dirCartesian = y; mz.dirCartesian = z;
				
				mx.dirLattice = e.gInfo.invR*mx.dirCartesian;
				my.dirLattice = e.gInfo.invR*my.dirCartesian;
				mz.dirLattice = e.gInfo.invR*mz.dirCartesian;
				
				modes.push_back(mx);
				modes.push_back(my);
				modes.push_back(mz);
			}
		}
	}
	
	//logPrintf("Listing degrees of freedom.\n");
	//for (unsigned i = 0; i < modes.size(); i++) {
	//	logPrintf("Species %d, atom %d\n", modes[i].sp, modes[i].at);
	//	vector3<> drCart = modes[i].dirCartesian;
	//	vector3<> drLatt = modes[i].dirLattice;
	//	logPrintf("Perturbation in cartesian coords: %g %g %g\n", drCart[0], drCart[1], drCart[2]);
	//	logPrintf("Perturbation in lattice coords: %g %g %g\n\n", drLatt[0], drLatt[1], drLatt[2]);
	//}
}

double SpringConstant::computeMatrixElement(int a, int b) {
	assert(!pInfo.incommensurate);
	assert(a == lastAtomPerturbed);
	
	
	PerturbationInfo::Mode modeA = modes[a];
	PerturbationInfo::Mode modeB = modes[b];
	
	Energies Eminusminus, Eminusplus, Eplusminus, Eplusplus;
	
	//double h = iInfo.species[modeA.sp]->coreRadius * 1e-7;
	//double k = iInfo.species[modeB.sp]->coreRadius * 1e-7;
	double h = 1e-3;
	double k = 1e-3;
	
	assert(h > 0 && k > 0);
	
	getPerturbedEnergy(Eminusminus, modeA, modeB, -h/2.0, -k/2.0);
	getPerturbedEnergy(Eminusplus, modeA, modeB, -h/2.0, k/2.0);
	getPerturbedEnergy(Eplusminus, modeA, modeB, h/2.0, -k/2.0);
	getPerturbedEnergy(Eplusplus, modeA, modeB, h/2.0, k/2.0);
	
	double dEtautau = (getIonDependentEnergy(Eplusplus) - getIonDependentEnergy(Eplusminus) - getIonDependentEnergy(Eminusplus) + getIonDependentEnergy(Eminusminus)) / (h*k);
	
	pInfo.atomdisplacement = modeB;
	//ps.updateHC();
	ps.updateNonlocalDerivs();
	ps.calcdGradTau();
	double kab = dEtautau + dot(pInfo.dGradTau, pInfo.dC, &pInfo, &eInfo);
	
	return kab;
}


void SpringConstant::getPerturbedEnergy(Energies& ener, PerturbationInfo::Mode modeA, PerturbationInfo::Mode modeB, double deltaA, double deltaB) {
	auto spA = iInfo.species[modeA.sp];
	auto spB = iInfo.species[modeB.sp];
	
	vector3<> posA_unperturbed = spA->atpos[modeA.at];
	spA->atpos[modeA.at] = posA_unperturbed + deltaA*modeA.dirLattice;
	mpiWorld->bcastData(spA->atpos);
	
	vector3<> posB_unperturbed = spB->atpos[modeB.at];
	spB->atpos[modeB.at] = posB_unperturbed + deltaB*modeB.dirLattice;
	mpiWorld->bcastData(spB->atpos);
	
	spA->sync_atpos();
	spB->sync_atpos();
	
	iInfo.update(ener);
	
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		if (ps.ultrasoftDensityAugRequired) {
			Ctmp[q] = eVars.C[q];
			ps.orthonormalize(q, eVars.C[q]);
		}
		
		e.iInfo.project(eVars.C[q], eVars.VdagC[q]);
	}
	eVars.elecEnergyAndGrad(ener, 0, 0, true);
	
	//pInfo.updateExcorrCache(e.exCorr, e.gInfo, ps.addnXC(eVars.n)[0]);
	//ps.updateHC();
	//ps.updateNonlocalDerivs();
	spB->atpos[modeB.at] = posB_unperturbed;
	mpiWorld->bcastData(spB->atpos);
	spA->atpos[modeA.at] = posA_unperturbed;
	mpiWorld->bcastData(spA->atpos);
	
	spB->sync_atpos();
	spA->sync_atpos();
	
	iInfo.update(e.ener);
	for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
		if (ps.ultrasoftDensityAugRequired)
			eVars.C[q] = Ctmp[q];
		e.iInfo.project(eVars.C[q], eVars.VdagC[q]);
	}
	eVars.elecEnergyAndGrad(e.ener, 0, 0, true);
}

double SpringConstant::getIonDependentEnergy(const Energies& ener) {
	return double(ener.E);
	//return ener.E["Eewald"];
}

void SpringConstant::computeMatrix() {
	pInfo.atposPerturbationExists = true;
	int n = modes.size();
	kmatrix = matrix(n, n);
	kmatrix.zero();
	
	ps.ultrasoftDensityAugRequired = true;
	if (ps.ultrasoftDensityAugRequired)
		init(Ctmp, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
	
	for (unsigned a = 0; a < n; a++) {
		pInfo.atomdisplacement = modes[a];
		ps.minimize(e.vptParams);
		lastAtomPerturbed = a;
				
		for (unsigned b = 0; b < n; b++) {
			kmatrix.set(a, b, computeMatrixElement(a, b));
		}
	}
}
