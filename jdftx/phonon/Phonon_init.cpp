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

#include <phonon/Phonon.h>
#include <electronic/ElecMinimizer.h>
#include <commands/parser.h>
#include <core/Units.h>

Phonon::Phonon()
: dr(0.01), T(298*Kelvin), Fcut(1e-8)
{
}

//Return size of stabilizer group of a Cartesian displacement (given Cartesian symmetry rotations)
inline int nStabilizer(const vector3<>& n, const std::vector< matrix3<> >& symCart)
{	int nStab = 0;
	for(const matrix3<>& m: symCart)
		if((n - m * n).length_squared() < symmThresholdSq * n.length_squared())
			nStab++;
	return nStab;
}

void Phonon::setup(bool printDefaults)
{
	//Parse input to initialize unit cell:
	parse(input, e, printDefaults);
	logSuspend();
	parse(input, eSupTemplate); //silently create a copy by re-parsing input (Everything is not trivially copyable)
	logResume();
	
	//Ensure phonon command specified:
	if(!sup.length())
		die("phonon supercell must be specified using the phonon command.\n");
	if(!e.gInfo.S.length_squared())
		die("Manual fftbox setting required for phonon. If supercell grid\n"
			"initialization fails, specify slightly larger manual fftbox.\n");
	//Check kpoint and supercell compatibility:
	if(e.eInfo.qnums.size()>1 || e.eInfo.qnums[0].k.length_squared())
		die("phonon requires a Gamma-centered uniform kpoint mesh.\n");
	for(int j=0; j<3; j++)
	{	if(!sup[j] || e.eInfo.kfold[j] % sup[j])
		{	die("kpoint folding %d is not a multiple of supercell count %d for lattice direction %d.\n",
				e.eInfo.kfold[j], sup[j], j);
		}
		eSupTemplate.eInfo.kfold[j] = e.eInfo.kfold[j] / sup[j];
	}
	
	logPrintf("########### Unit cell calculation #############\n");
	SpeciesInfo::Constraint constraintFull;
	constraintFull.moveScale = 0;
	constraintFull.type = SpeciesInfo::Constraint::None;
	for(size_t sp=0; sp<e.iInfo.species.size(); sp++)
		e.iInfo.species[sp]->constraints.assign(e.iInfo.species[sp]->atpos.size(), constraintFull);
	e.setup();
	if(!e.coulombParams.supercell) e.updateSupercell(true); //force supercell generation

	nSpins = e.eInfo.spinType==SpinZ ? 2 : 1;
	nSpinor = e.eInfo.spinorLength();

	//Initialize state of unit cell:
	if(e.cntrl.dumpOnly)
	{	//Single energy calculation so that all dependent quantities have been initialized:
		logPrintf("\n----------- Energy evaluation at fixed state -------------\n"); logFlush();
		e.eVars.elecEnergyAndGrad(e.ener, 0, 0, true);
	}
	else elecFluidMinimize(e);
	logPrintf("# Energy components:\n"); e.ener.print(); logPrintf("\n");

	//Determine optimum number of bands for supercell calculation:
	nBandsOpt = 0;
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{	int nBands_q = std::upper_bound(e.eVars.F[q].begin(), e.eVars.F[q].end(), Fcut, std::greater<double>()) - e.eVars.F[q].begin();
		nBandsOpt = std::max(nBandsOpt, nBands_q);
	}
	mpiUtil->allReduce(nBandsOpt, MPIUtil::ReduceMax);
	logPrintf("Fcut=%lg reduced nBands from %d to %d per unit cell.\n", Fcut, e.eInfo.nBands, nBandsOpt);

	//Make unit cell state available on all processes 
	//(since MPI division of qSup and q are different and independent of the map)
	for(int q=0; q<e.eInfo.nStates; q++)
	{	//Allocate:
		if(!e.eInfo.isMine(q))
		{	e.eVars.C[q].init(e.eInfo.nBands, e.basis[q].nbasis * e.eInfo.spinorLength(), &e.basis[q], &e.eInfo.qnums[q]);
			e.eVars.F[q].resize(e.eInfo.nBands);
			if(e.eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux)
				e.eVars.B[q].init(e.eInfo.nBands, e.eInfo.nBands);
		}
		//Broadcast from owner:
		int qSrc = e.eInfo.whose(q);
		e.eVars.C[q].bcast(qSrc);
		e.eVars.F[q].bcast(qSrc);
		if(e.eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux)
			e.eVars.B[q].bcast(qSrc);
	}

	logPrintf("\n------- Configuring supercell and perturbation modes -------\n");
	
	//Grid:
	eSupTemplate.gInfo.S = Diag(sup) * e.gInfo.S; //ensure exact supercell
	eSupTemplate.gInfo.R = e.gInfo.R * Diag(sup);
	prodSup = sup[0] * sup[1] * sup[2];
	
	//Replicate atoms (and related properties):
	for(size_t sp=0; sp<e.iInfo.species.size(); sp++)
	{	const SpeciesInfo& spIn = *(e.iInfo.species[sp]);
		SpeciesInfo& spOut = *(eSupTemplate.iInfo.species[sp]);
		spOut.atpos.clear();
		spOut.initialMagneticMoments.clear();
		matrix3<> invSup = inv(Diag(vector3<>(sup)));
		vector3<int> iR;
		for(iR[0]=0; iR[0]<sup[0]; iR[0]++)
		for(iR[1]=0; iR[1]<sup[1]; iR[1]++)
		for(iR[2]=0; iR[2]<sup[2]; iR[2]++)
		{	for(vector3<> pos: spIn.atpos)
				spOut.atpos.push_back(invSup * (pos + iR));
			for(vector3<> M: spIn.initialMagneticMoments)
				spOut.initialMagneticMoments.push_back(M); //needed only to determine supercell symmetries
		}
		spOut.constraints.assign(spOut.atpos.size(), constraintFull);
	}
	
	//Supercell symmetries:
	eSupTemplate.symm.setup(eSupTemplate);
	const std::vector< matrix3<int> >& symSup = eSupTemplate.symm.getMatrices();
	symSupCart.clear();
	eSupTemplate.gInfo.invR = inv(eSupTemplate.gInfo.R);
	for(const matrix3<int>& m: symSup)
		symSupCart.push_back(eSupTemplate.gInfo.R * m * eSupTemplate.gInfo.invR);
	
	//Pick maximally symmetric orthogonal basis:
	logPrintf("\nFinding maximally-symmetric orthogonal basis for displacements:\n");
	std::vector< vector3<> > dirBasis;
	{	std::multimap<int, vector3<> > dirList; //directions indexed by their stabilizer group cardinality
		vector3<int> iR;
		for(iR[0]=0; iR[0]<=+1; iR[0]++)
		for(iR[1]=-1; iR[1]<=+1; iR[1]++)
		for(iR[2]=-1; iR[2]<=+1; iR[2]++)
			if(iR.length_squared())
			{	//Try low-order lattice vector linear combination:
				vector3<> n = eSupTemplate.gInfo.R * iR; n *= (1./n.length());
				dirList.insert(std::make_pair(nStabilizer(n, symSupCart), n));
				//Try low-order reciprocal lattice vector linear combination:
				n = iR * eSupTemplate.gInfo.invR; n *= (1./n.length());
				dirList.insert(std::make_pair(nStabilizer(n, symSupCart), n));
			}
		dirBasis.push_back(dirList.rbegin()->second);
		//Pick second driection orthogonal to first:
		std::multimap<int, vector3<> > dirList2;
		for(auto entry: dirList)
		{	vector3<> n = entry.second;
			n -= dot(n, dirBasis[0]) * dirBasis[0];
			if(n.length_squared() < symmThresholdSq) continue;
			n *= (1./n.length());
			dirList2.insert(std::make_pair(nStabilizer(n, symSupCart), n));
		}
		dirBasis.push_back(dirList2.rbegin()->second);
		dirBasis.push_back(cross(dirBasis[0], dirBasis[1])); //third direction constrained by orthogonality
	}
	for(const vector3<>& n: dirBasis)
		logPrintf(" [ %+lf %+lf %+lf ] |Stabilizer|: %d\n", n[0], n[1], n[2], nStabilizer(n,symSupCart));
	
	//List all modes:
	modes.clear();
	for(size_t sp=0; sp<e.iInfo.species.size(); sp++)
		for(size_t at=0; at<e.iInfo.species[sp]->atpos.size(); at++) //only need to move atoms in first unit cell
			for(int iDir=0; iDir<3; iDir++)
			{	Mode mode;
				mode.sp = sp;
				mode.at = at;
				mode.dir[iDir] = 1.;
				modes.push_back(mode);
			}

	//Find irreducible modes:
	perturbations.clear();
	for(unsigned sp=0; sp<e.iInfo.species.size(); sp++)
	{	int nAtoms = e.iInfo.species[sp]->atpos.size();
		int nPert = nAtoms * dirBasis.size();
		//generate all perturbations first:
		std::vector<Perturbation> pertSp(nPert); //perturbations of this species
		std::vector<matrix> proj(nPert); //projection operator into subspace spanned by star of current perturbation
		matrix projTot;
		const auto& atomMap = eSupTemplate.symm.getAtomMap()[sp];
		for(int iPert=0; iPert<nPert; iPert++)
		{	pertSp[iPert].sp = sp;
			pertSp[iPert].at = iPert / dirBasis.size();
			pertSp[iPert].dir = dirBasis[iPert % dirBasis.size()];
			pertSp[iPert].weight = 1./symSupCart.size();
			for(unsigned iSym=0; iSym<symSupCart.size(); iSym++)
			{	int at = atomMap[pertSp[iPert].at][iSym] % nAtoms; //map back to first cell
				vector3<> dir = symSupCart[iSym] * pertSp[iPert].dir;
				matrix nHat = zeroes(nPert,1);
				for(int iDir=0; iDir<3; iDir++)
					nHat.set(at*3+iDir,0, dir[iDir]);
				proj[iPert] += pertSp[iPert].weight * nHat * dagger(nHat);
			}
			projTot += proj[iPert];
		}
		assert(nrm2(projTot - eye(nPert)) < symmThreshold);
		//only select perturbations with distinct subspace projections:
		std::vector<bool> irred(nPert, true); //whether each perturbation is in irreducible set
		for(int iPert=0; iPert<nPert; iPert++)
		{	for(int jPert=0; jPert<iPert; jPert++)
				if(irred[jPert] && nrm2(proj[iPert]-proj[jPert])<symmThreshold)
				{	pertSp[jPert].weight += pertSp[iPert].weight; //send weight of current mode to its image in irreducible set
					irred[iPert] = false; //this mode will be accounted for upon symmetrization
					break;
				}
		}
		for(int iPert=0; iPert<nPert; iPert++)
			if(irred[iPert])
				perturbations.push_back(pertSp[iPert]);
	}
	logPrintf("\n%d perturbations of the unit cell reduced to %d under symmetries:\n", int(modes.size()), int(perturbations.size()));
	for(const Perturbation& pert: perturbations)
		logPrintf("%s %d  [ %+lf %+lf %+lf ] %lf\n", e.iInfo.species[pert.sp]->name.c_str(),
			pert.at, pert.dir[0], pert.dir[1], pert.dir[2], pert.weight*symSupCart.size());
}
