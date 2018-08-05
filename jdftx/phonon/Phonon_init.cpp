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
#include <electronic/ColumnBundleTransform.h>
#include <commands/parser.h>
#include <core/Units.h>

PhononEverything::PhononEverything(Phonon& phonon) : phonon(phonon)
{
}

Phonon::Phonon()
: dr(0.1), T(298*Kelvin), Fcut(1e-8), rSmooth(1.), iPerturbation(-1), collectPerturbations(false), saveHsub(true), e(*this), eSupTemplate(*this)
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

inline bool isUnitary(const matrix& U) { return nrm2(U*dagger(U) - eye(U.nCols())) < symmThreshold; }

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

	nSpins = e.eInfo.nSpins();
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
	nBandsOpt = (e.eInfo.fillingsUpdate == ElecInfo::FillingsHsub)
		? 1+int(ceil(e.eInfo.nElectrons/e.eInfo.qWeightSum)) //make sure extra band present for smearing
		: 0; //no smearing: tightest number of bands acceptable
	nBandsOpt = std::min(e.eInfo.nBands, nBandsOpt); //in-case nBands was borderline
	for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
	{	int nBands_q = std::upper_bound(e.eVars.F[q].begin(), e.eVars.F[q].end(), Fcut, std::greater<double>()) - e.eVars.F[q].begin();
		nBandsOpt = std::max(nBandsOpt, nBands_q);
	}
	mpiWorld->allReduce(nBandsOpt, MPIUtil::ReduceMax);
	logPrintf("Fcut=%lg reduced nBands from %d to %d per unit cell.\n", Fcut, e.eInfo.nBands, nBandsOpt);

	//Make unit cell state available on all processes 
	//(since MPI division of qSup and q are different and independent of the map)
	for(int q=0; q<e.eInfo.nStates; q++)
	{	//Allocate:
		if(!e.eInfo.isMine(q))
		{	e.eVars.C[q].init(e.eInfo.nBands, e.basis[q].nbasis * e.eInfo.spinorLength(), &e.basis[q], &e.eInfo.qnums[q]);
			e.eVars.F[q].resize(e.eInfo.nBands);
			e.eVars.Hsub_eigs[q].resize(e.eInfo.nBands);
			if(e.eInfo.fillingsUpdate==ElecInfo::FillingsHsub)
				e.eVars.Haux_eigs[q].resize(e.eInfo.nBands);
		}
		//Broadcast from owner:
		int qSrc = e.eInfo.whose(q);
		mpiWorld->bcastData(e.eVars.C[q], qSrc);
		mpiWorld->bcastData(e.eVars.F[q], qSrc);
		mpiWorld->bcastData(e.eVars.Hsub_eigs[q], qSrc);
		if(e.eInfo.fillingsUpdate==ElecInfo::FillingsHsub)
			mpiWorld->bcastData(e.eVars.Haux_eigs[q], qSrc);
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
	
	//Store ideal energy and forces for unperturbed supercell:
	E0 = prodSup * relevantFreeEnergy(e);
	e.iInfo.ionicEnergyAndGrad(e.iInfo.forces); //compute forces in lattice coordinates
	IonicGradient grad0unit = -e.gInfo.invRT * e.iInfo.forces; //gradient in cartesian coordinates (and negative of force)
	grad0.assign(e.iInfo.forces.size(), std::vector<vector3<>>());
	for(size_t iSp=0; iSp<grad0.size(); iSp++)
		for(int iCell=0; iCell<prodSup; iCell++) //repeat for each unit cell in supercell
			grad0[iSp].insert(grad0[iSp].end(), grad0unit[iSp].begin(), grad0unit[iSp].end());
	
	//Supercell symmetries:
	//--- handle manual symmetries:
	if(eSupTemplate.symm.sym.size())
	{	eSupTemplate.symm.sym.clear();
		matrix3<> DiagSup = Diag(vector3<>(sup)), invDiagSup = inv(DiagSup);
		for(const SpaceGroupOp& op: e.symm.sym)
		{	matrix3<> rotSupTemp = DiagSup * op.rot * invDiagSup;
			SpaceGroupOp opSup;
			bool isSymSup = true;
			for(int i=0; i<3; i++)
				for(int j=0; j<3; j++)
				{	opSup.rot(i,j) = int(round(rotSupTemp(i,j)));
					if(fabs(rotSupTemp(i,j) - opSup.rot(i,j)) > symmThreshold)
						isSymSup = false; //not a supercell symmetry
				}
			if(!isSymSup) continue;
			opSup.a = op.a * invDiagSup;
			eSupTemplate.symm.sym.push_back(opSup);
		}
	}
	eSupTemplate.symm.sup = sup; //restrict space group to translations within unit cell
	eSupTemplate.symm.setup(eSupTemplate);
	symSup = eSupTemplate.symm.getMatrices();
	symSupCart.clear();
	eSupTemplate.gInfo.invR = inv(eSupTemplate.gInfo.R);
	for(const SpaceGroupOp& op: symSup)
		symSupCart.push_back(eSupTemplate.gInfo.R * op.rot * eSupTemplate.gInfo.invR);
	
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
		std::vector<bool> isEven(nPert); //whether perturbation is even under symmetries of calculation
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
				if(at==pertSp[iPert].at && (dir+pertSp[iPert].dir).length()<symmThreshold)
					isEven[iPert] = true; //symmetries make positive and negative displacements equivalent
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
			{	if(isEven[iPert]) perturbations.push_back(pertSp[iPert]);
				else //split into two perturbations to get a central difference formula
				{	pertSp[iPert].weight *= 0.5;
					for(int iSign=0; iSign<2; iSign++)
					{	perturbations.push_back(pertSp[iPert]);
						pertSp[iPert].dir *= -1.;
					}
				}
			}
	}
	logPrintf("\n%d signed perturbations of the unit cell reduced to %d under symmetries:\n", 2*int(modes.size()), int(perturbations.size()));
	for(const Perturbation& pert: perturbations)
		logPrintf("%s %d  [ %+lf %+lf %+lf ] %lf\n", e.iInfo.species[pert.sp]->name.c_str(),
			pert.at, pert.dir[0], pert.dir[1], pert.dir[2], pert.weight*symSupCart.size());
	if(iPerturbation>=int(perturbations.size()))
		die("Specified iPerturbation %d in command phonon is invalid (since it is > %lu).\n", iPerturbation+1, perturbations.size());
	
	//Determine wavefunction unitary rotations:
	logPrintf("\nCalculating unitary rotations of unit cell states under symmetries:\n");
	stateRot.resize(nSpins);
	double unitarityErr = 0.;
	for(int iSpin=0; iSpin<nSpins; iSpin++)
	{	//Find states involved in the supercell Gamma-point:
		struct Kpoint : public Supercell::KmeshTransform
		{	vector3<> k; //also store k-point for convenience (KmeshTransform doesn't have it)
		};
		std::vector<Kpoint> kpoints; kpoints.reserve(prodSup);
		const Supercell& supercell = *(e.coulombParams.supercell);
		for(unsigned ik=0; ik<supercell.kmesh.size(); ik++)
		{	double kSupErr; round(matrix3<>(Diag(sup)) * supercell.kmesh[ik], &kSupErr);
			if(kSupErr < symmThreshold) //maps to Gamma point
			{	Kpoint kpoint;
				(Supercell::KmeshTransform&)kpoint = supercell.kmeshTransform[ik]; //copy base class
				kpoint.k = supercell.kmesh[ik];
				kpoint.iReduced += iSpin*(e.eInfo.nStates/nSpins); //point to source k-point with appropriate spin
				kpoints.push_back(kpoint);
			}
		}
		assert(int(kpoints.size()) == prodSup);
		//Initialize basis and qnum for these states:
		std::vector<QuantumNumber> qnums(prodSup);
		std::vector<Basis> basis(prodSup);
		logSuspend();
		for(int ik=0; ik<prodSup; ik++)
		{	qnums[ik].k = kpoints[ik].k;
			qnums[ik].spin = (nSpins==1 ? 0 : (iSpin==0 ? +1 : -1));
			qnums[ik].weight = 1./prodSup;
			basis[ik].setup(e.gInfo, e.iInfo, e.cntrl.Ecut, kpoints[ik].k);
		}
		logResume();
		//Get wavefunctions for all these k-points:
		#define whose_ik(ik) (((ik) * mpiWorld->nProcesses())/prodSup) //local MPI division
		std::vector<ColumnBundle> C(prodSup);
		std::vector<std::shared_ptr<ColumnBundleTransform::BasisWrapper> > basisWrapper(prodSup);
		auto sym = e.symm.getMatrices(); //unit cell symmetries
		for(int ik=0; ik<prodSup; ik++)
		{	C[ik].init(e.eInfo.nBands, basis[ik].nbasis*nSpinor, &basis[ik], &qnums[ik], isGpuEnabled());
			if(whose_ik(ik) == mpiWorld->iProcess())
			{	int q = kpoints[ik].iReduced;
				C[ik].zero();
				basisWrapper[ik] = std::make_shared<ColumnBundleTransform::BasisWrapper>(basis[ik]);
				ColumnBundleTransform(e.eInfo.qnums[q].k, e.basis[q], qnums[ik].k, *(basisWrapper[ik]),
					nSpinor, sym[kpoints[ik].iSym], kpoints[ik].invert).scatterAxpy(1., e.eVars.C[q], C[ik],0,1);
			}
		}
		for(int ik=0; ik<prodSup; ik++)
			mpiWorld->bcastData(C[ik], whose_ik(ik)); //make available on all processes
		//Determine max eigenvalue:
		int nBands = e.eInfo.nBands;
		double Emax = -INFINITY;
		for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
			Emax = std::max(Emax, e.eVars.Hsub_eigs[q].back());
		mpiWorld->allReduce(Emax, MPIUtil::MPIUtil::ReduceMax);
		double EmaxValid = +INFINITY;
		//Loop over supercell symmetry operations:
		PeriodicLookup<QuantumNumber> plook(qnums, e.gInfo.GGT);
		stateRot[iSpin].resize(symSupCart.size());
		for(size_t iSym=0; iSym<symSupCart.size(); iSym++)
		{	matrix3<> rotUnitTmp = e.gInfo.invR * symSupCart[iSym] * e.gInfo.R; //in unit cell lattice coordinates
			#define SymmErrMsg \
				"Supercell symmetries do not map unit cell k-point mesh onto itself.\n" \
				"This implies that the supercell is more symmetric than the unit cell!\n" \
				"Please check to make sure that you have used the minimal unit cell.\n\n"
			SpaceGroupOp symUnit;
			for(int j1=0; j1<3; j1++)
				for(int j2=0; j2<3; j2++)
				{	symUnit.rot(j1,j2) = round(rotUnitTmp(j1,j2));
					if(fabs(symUnit.rot(j1,j2) - rotUnitTmp(j1,j2)) > symmThreshold)
						die(SymmErrMsg)
				}
			symUnit.a = Diag(sup) * symSup[iSym].a;
			//Find image kpoints under rotation: (do this for all k-points so that all processes exit together if necessary)
			std::vector<int> ikRot(prodSup);
			for(int ik=0; ik<prodSup; ik++)
			{	size_t ikRotCur = plook.find(qnums[ik].k * symUnit.rot);
				if(ikRotCur==string::npos) die(SymmErrMsg)
				ikRot[ik] = ikRotCur;
			}
			#undef SymmErrMsg
			//Calculate unitary transformation matrix:
			stateRot[iSpin][iSym].init(prodSup, nBands);
			for(int ik=0; ik<prodSup; ik++)
				if(whose_ik(ikRot[ik]) == mpiWorld->iProcess()) //MPI division by target k-point
				{	ColumnBundle Crot = C[ikRot[ik]].similar();
					Crot.zero();
					ColumnBundleTransform(qnums[ik].k, basis[ik], qnums[ikRot[ik]].k, *(basisWrapper[ikRot[ik]]),
						nSpinor, symUnit, +1).scatterAxpy(1., C[ik], Crot,0,1);
					matrix Urot = Crot ^ O(C[ikRot[ik]]); //will be unitary if Crot is a strict unitary rotation of C[ikRot[ik]]
					//Check maximal subspace that is unitary: (remainder must be incomplete degenerate subspace)
					int nBandsValid = nBands;
					while(nBandsValid && !isUnitary(Urot(0,nBandsValid, 0,nBandsValid)))
						nBandsValid--;
					if(nBandsValid<nBands)
					{	//Update energy range of validity:
						EmaxValid = std::min(EmaxValid, e.eVars.Hsub_eigs[kpoints[ik].iReduced][nBandsValid]);
						//Make valid subspace exactly unitary:
						matrix UrotSub = Urot(0,nBandsValid, 0,nBandsValid);
						matrix UrotOverlap = dagger(UrotSub) * UrotSub;
						UrotSub = UrotSub * invsqrt(UrotOverlap); //make exactly unitary
						unitarityErr += std::pow(nrm2(UrotOverlap - eye(nBandsValid)), 2);
						//Zero out invalid subspace:
						Urot.zero();
						Urot.set(0,nBandsValid, 0,nBandsValid, UrotSub);
					}
					stateRot[iSpin][iSym].set(ik, ikRot[ik], Urot);
				}
			stateRot[iSpin][iSym].allReduce();
		}
		#undef whose_ik
		mpiWorld->allReduce(EmaxValid, MPIUtil::ReduceMin);
		if(nSpins>1) logPrintf("\tSpin %+d: ", iSpin==0 ? +1 : -1);  else logPrintf("\t");
		logPrintf("Matrix elements valid for ");
		if(std::isfinite(EmaxValid)) logPrintf("E < %+.6lf (Emax = %+.6lf) due to incomplete degenerate subspaces.\n", EmaxValid, Emax);
		else logPrintf("all available states (all degenerate subspaces are complete).\n");
	}
	mpiWorld->allReduce(unitarityErr, MPIUtil::ReduceSum);
	unitarityErr = sqrt(unitarityErr / (nSpins * prodSup * symSupCart.size()));
	logPrintf("\tRMS unitarity error in valid subspaces: %le\n", unitarityErr);
}
