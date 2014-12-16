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
#include <core/Units.h>
#include <core/LatticeUtils.h>
#include <electronic/ElecMinimizer.h>
#include <electronic/operators.h>

Phonon::Phonon()
: dr(0.01), T(298*Kelvin), Fcut(1e-8)
{
	eSup.cntrl.dragWavefunctions = false; //wavefunction-drag doesn't always play nice with setSupState (especially with relativity)
}

inline vector3<> getCoord(const QuantumNumber& qnum) { return qnum.k; } //for k-point mapping
inline bool spinEqual(const QuantumNumber& qnum1, const QuantumNumber& qnum2) { return qnum1.spin == qnum2.spin; } //for k-point mapping (in spin polarized mode)

//Generate quadrature for BZ integrals with singularity at Gamma point:
//--- add contribution to quadrature due to box of size scale centered at gCenter (in reciprocal lattice coords)
void addQuadratureBZ_box(std::vector< std::pair<vector3<>,double> >& quad, double scale, vector3<> gCenter)
{	//Weights and abscissae of the 7-point gauss quadrature:
	const int N = 3;
	static const double w[N+1] = { 0.129484966168869693270611432679082, 0.279705391489276667901467771423780, 0.381830050505118944950369775488975, 0.417959183673469387755102040816327 };
	static const double x[N+1] = { 0.949107912342758524526189684047851, 0.741531185599394439863864773280788, 0.405845151377397166906606412076961, 0.000000000000000000000000000000000 };
	double h = 0.5 * scale; 
	vector3<> g;
	for(int i0=-N; i0<=N; i0++)
	{	g[0] = gCenter[0] + h*x[N-abs(i0)]*(i0>0?1:-1);
		double w0 = (N ? h*w[N-abs(i0)] : 1.);
		for(int i1=-N; i1<=N; i1++)
		{	g[1] = gCenter[1] + h*x[N-abs(i1)]*(i1>0?1:-1);
			double w01 = w0 * (N ? h*w[N-abs(i1)] : 1.);
			for(int i2=-N; i2<=N; i2++)
			{	g[2] = gCenter[2] + h*x[N-abs(i2)]*(i2>0?1:-1);
				double w012 = w01 * (N ? h*w[N-abs(i2)] : 1.);
				quad.push_back(std::make_pair(g, w012));
			}
		}
	}
}
//--- add contribution to quadrature between bozes of size scale and scale/3
void addQuadratureBZ_scale(std::vector< std::pair<vector3<>,double> >& quad, double scale) 
{	double scaleBy3 = scale/3.;
	vector3<int> ig;
	for(ig[0]=-1; ig[0]<=1; ig[0]++)
	for(ig[1]=-1; ig[1]<=1; ig[1]++)
	for(ig[2]=-1; ig[2]<=1; ig[2]++)
		if(ig.length_squared()) //except the center box
			addQuadratureBZ_box(quad, scaleBy3, scaleBy3*ig);
}
std::vector< std::pair<vector3<>,double> > getQuadratureBZ()
{	std::vector< std::pair<vector3<>,double> > quad;
	for(double scale=1.; scale>1e-3; scale/=3.)
		addQuadratureBZ_scale(quad, scale);
	return quad;
}

void Phonon::setup()
{
	//Ensure phonon command specified:
	if(!sup.length())
		die("phonon supercell must be specified using the phonon command.\n");
	if(!e.gInfo.S.length_squared())
		logPrintf(
			"WARNING: manual fftbox setting recommended for phonon. If supercell\n"
			"     grid initialization fails below, specify larger manual fftbox.\n");
	//Check symmetries
	if(e.symm.mode != SymmetriesNone)
		die("phonon does not support symmetries in the input file / unit cell calculation.\n"
			"Symmetries of the supercell will be applied to the force matrix automatically.\n");
	//Check kpoint and supercell compatibility:
	if(e.eInfo.qnums.size()>1 || e.eInfo.qnums[0].k.length_squared())
		die("phonon requires a Gamma-centered uniform kpoint mesh.\n");
	for(int j=0; j<3; j++)
	{	if(!sup[j] || e.eInfo.kfold[j] % sup[j])
		{	die("kpoint folding %d is not a multiple of supercell count %d for lattice direction %d.\n",
				e.eInfo.kfold[j], sup[j], j);
		}
		eSup.eInfo.kfold[j] = e.eInfo.kfold[j] / sup[j];
	}
	
	logPrintf("########### Initializing JDFTx for unit cell #############\n");
	SpeciesInfo::Constraint constraintFull;
	constraintFull.moveScale = 0;
	constraintFull.type = SpeciesInfo::Constraint::None;
	for(size_t sp=0; sp<e.iInfo.species.size(); sp++)
		e.iInfo.species[sp]->constraints.assign(e.iInfo.species[sp]->atpos.size(), constraintFull);
	e.setup();
	
	//Initialize state of unit cell:
	logPrintf("########### Unperturbed unit cell calculation #############\n");
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
	logPrintf("Fcut=%lg reduced nBands from %d to %d per unit cell.\n\n", Fcut, e.eInfo.nBands, nBandsOpt);

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

	logPrintf("########### Initializing JDFTx for supercell #############\n");
	//Grid:
	eSup.gInfo.S = Diag(sup) * e.gInfo.S; //ensure exact supercell
	eSup.gInfo.R = e.gInfo.R * Diag(sup);
	int prodSup = sup[0] * sup[1] * sup[2];
	//Replicate atoms:
	for(size_t sp=0; sp<e.iInfo.species.size(); sp++)
	{	eSup.iInfo.species[sp]->atpos.clear();
		eSup.iInfo.species[sp]->initialMagneticMoments.clear();
		matrix3<> invSup = inv(Diag(vector3<>(sup)));
		vector3<int> iR;
		for(iR[0]=0; iR[0]<sup[0]; iR[0]++)
		for(iR[1]=0; iR[1]<sup[1]; iR[1]++)
		for(iR[2]=0; iR[2]<sup[2]; iR[2]++)
			for(vector3<> pos: e.iInfo.species[sp]->atpos)
				eSup.iInfo.species[sp]->atpos.push_back(invSup * (pos + iR));
		eSup.iInfo.species[sp]->constraints.assign(eSup.iInfo.species[sp]->atpos.size(), constraintFull);
		for(SpeciesInfo::PlusU& Uparams: eSup.iInfo.species[sp]->plusU)
			Uparams.Vext.resize(eSup.iInfo.species[sp]->atpos.size());
	}
	//Remove initial state settings (incompatible with supercell):
	eSup.eVars.wfnsFilename.clear();
	eSup.eVars.HauxFilename.clear();
	eSup.eVars.eigsFilename.clear();
	eSup.eVars.fluidInitialStateFilename.clear();
	eSup.eInfo.initialFillingsFilename.clear();
	eSup.scfParams.historyFilename.clear();
	//ElecInfo:
	eSup.eInfo.nBands = nBandsOpt * prodSup;
	//ElecVars:
	eSup.eVars.initLCAO = false;
	eSup.setup();
	
	//Map states:
	PeriodicLookup<QuantumNumber> plook(eSup.eInfo.qnums, eSup.gInfo.GGT);
	std::vector<int> nqPrev(eSup.eInfo.qnums.size(), 0);
	for(int q=0; q<e.eInfo.nStates; q++)
	{	const QuantumNumber& qnum = e.eInfo.qnums[q];
		vector3<> kSup = matrix3<>(Diag(sup)) * qnum.k; //qnum.k in supercell reciprocal lattice coords
		size_t qSup = plook.find(kSup, qnum, &(eSup.eInfo.qnums), spinEqual);
		assert(qSup != string::npos);
		vector3<> iGtmp = kSup - eSup.eInfo.qnums[qSup].k;
		//Add to stateMap:
		auto sme = std::make_shared<StateMapEntry>();
		sme->qSup = qSup;
		for(int j=0; j<3; j++)
			sme->iG[j] = round(iGtmp[j]);
		sme->nqPrev = nqPrev[qSup];
		nqPrev[qSup]++;
		stateMap.push_back(sme);
	}
	for(int nq: nqPrev) assert(nq == prodSup); //each supercell k-point must be mapped to prod(sup) unit cell kpoints
	
	//Map corresponding basis objects:
	for(int qSup=eSup.eInfo.qStart; qSup<eSup.eInfo.qStop; qSup++)
	{	//Create lookup table for supercell indices:
		vector3<int> iGbox; //see Basis::setup(), this will bound the values of iG in the basis
		for(int i=0; i<3; i++)
			iGbox[i] = 1 + int(sqrt(2*eSup.cntrl.Ecut) * eSup.gInfo.R.column(i).length() / (2*M_PI));
		vector3<int> pitch;
		pitch[2] = 1;
		pitch[1] = pitch[2] * (2*iGbox[2]+1);
		pitch[0] = pitch[1] * (2*iGbox[1]+1);
		std::vector<int> supIndex(pitch[0] * (2*iGbox[0]+1), -1);
		const Basis& basisSup = eSup.basis[qSup];
		for(size_t n=0; n<basisSup.nbasis; n++)
			supIndex[dot(pitch,basisSup.iGarr[n]+iGbox)] = n;
		//Initialize index map for each unit cell k-point
		for(int q=0; q<e.eInfo.nStates; q++) if(stateMap[q]->qSup == qSup)
		{	const Basis& basis = e.basis[q];
			std::vector<int> indexVec(basis.nbasis);
			for(size_t n=0; n<basis.nbasis; n++)
			{	vector3<int> iGsup = basis.iGarr[n];
				for(int j=0; j<3; j++) iGsup[j] *= sup[j]; //to supercell reciprocal lattice coords
				iGsup += stateMap[q]->iG; //offset due to Bloch phase
				indexVec[n] = supIndex[dot(pitch,iGsup+iGbox)]; //lookup from above table
			}
			assert(*std::min_element(indexVec.begin(), indexVec.end()) >= 0); //make sure all entries were found
			stateMap[q]->setIndex(indexVec);
		}
	}
	
	//Symmetries:
	logPrintf("########### Initializing supercell symmetries #############\n");
	symm.mode = SymmetriesAutomatic;
	symm.setup(eSup);
}

void Phonon::dump()
{
	//List of displacements:
	struct Mode { int sp; int at; int dir; }; //specie, atom and cartesian directions for each displacement
	std::vector<Mode> modes;
	for(size_t sp=0; sp<e.iInfo.species.size(); sp++)
		for(size_t at=0; at<e.iInfo.species[sp]->atpos.size(); at++) //only need to move atoms in first unit cell
			for(int dir=0; dir<3; dir++)
			{	Mode mode = { int(sp), int(at), dir };
				modes.push_back(mode);
			}
	
	//Initialize state of supercell:
	logPrintf("########### Unperturbed supercell calculation #############\n");
	setSupState();
	IonicMinimizer imin(eSup);
	IonicGradient grad0;
	imin.compute(&grad0); //compute initial forces and energy
	logPrintf("# Energy components:\n"); eSup.ener.print(); logPrintf("\n");
	int prodSup = sup[0]*sup[1]*sup[2];
	double E0 = relevantFreeEnergy(eSup);
	logPrintf("Supercell energy discrepancy: %lg / unit cell\n", E0/prodSup - relevantFreeEnergy(e));
	logPrintf("RMS force in initial configuration: %lg\n", sqrt(dot(grad0,grad0)/(modes.size()*prodSup)));
	//subspace Hamiltonian of supercell Gamma point:
	std::vector<matrix> Hsub0;
	setSupState(&Hsub0);
	
	//Displaced modes:
	std::vector<IonicGradient> dgrad(modes.size()); //change in gradient for each perturbation
	std::vector< std::vector<matrix> > dHsub(modes.size()); //chaneg in Hsub for each perturbation
	for(size_t iMode=0; iMode<modes.size(); iMode++)
	{	const Mode& mode = modes[iMode];
		logPrintf("\n########### Perturbed supercell calculation %d of %d #############\n", int(iMode+1), int(modes.size()));
		//Move one atom:
		IonicGradient dir; dir.init(eSup.iInfo); //contains zeroes
		dir[mode.sp][mode.at][mode.dir] = 1.;
		imin.step(dir, dr);
		//Calculate energy and forces:
		IonicGradient grad;
		imin.compute(&grad);
		dgrad[iMode] = (grad - grad0) * (1./dr);
		logPrintf("Energy change: %lg / unit cell\n", (relevantFreeEnergy(eSup) - E0)/prodSup);
		logPrintf("RMS force: %lg\n", sqrt(dot(grad,grad)/(modes.size()*prodSup)));
		//Subspace hamiltonian change:
		std::vector<matrix> Hsub;
		setSupState(&Hsub);
		dHsub[iMode].resize(Hsub.size());
		for(size_t s=0; s<Hsub.size(); s++)
			dHsub[iMode][s] = (1./dr) * (Hsub[s] - Hsub0[s]);
		//Restore atom position:
		bool dragWfns = false;
		std::swap(dragWfns, eSup.cntrl.dragWavefunctions); //disable wave function drag because state has already been restored to unperturbed version
		imin.step(dir, -dr);
		std::swap(dragWfns, eSup.cntrl.dragWavefunctions); //restore wave function drag flag
	}
	
	//Symmetrize force matrix:
	//--- framework for organizing forces into 3x3 matrices per atom pairs:
	struct ForceMatrixIndex
	{	int sp1, at1, sp2, at2; 
		//Comparison operator for using in map:
		bool operator<(const ForceMatrixIndex& fmi) const
		{	if(sp1 < fmi.sp1) return true;
			if(sp1 == fmi.sp1)
			{	if(at1 < fmi.at1) return true;
				if(at1 == fmi.at1)
				{	if(sp2 < fmi.sp2) return true;
					if(sp2 == fmi.sp2)
						return (at2 < fmi.at2);
				}
			}
			return false;
		}
	};
	struct ForceMatrixIndexWrap //use periodicity to make sure at1 always belongs to unit cell
	{	const Phonon& phonon;
		ForceMatrixIndexWrap(const Phonon& phonon) : phonon(phonon) {}
		void operator()(ForceMatrixIndex& fmi)
		{	int nAtoms1 = phonon.e.iInfo.species[fmi.sp1]->atpos.size();
			int nAtoms2 = phonon.e.iInfo.species[fmi.sp2]->atpos.size();
			int unit1 = fmi.at1 / nAtoms1; if(!unit1) return; //already in unit cell
			int unit2 = fmi.at2 / nAtoms2;
			vector3<int> iR1 = getCell(unit1);
			vector3<int> iR2 = getCell(unit2);
			vector3<int> iR2new = iR2 - iR1; //so that atom goes to unit cell
			for(int j=0; j<3; j++)
				iR2new[j] = positiveRemainder(iR2new[j], phonon.sup[j]);
			int unit2new = (iR2new[0]*phonon.sup[1] + iR2new[1])*phonon.sup[2] + iR2new[2];
			fmi.at1 += (0 - unit1) * nAtoms1;
			fmi.at2 += (unit2new - unit2) * nAtoms2;
		}
	private:
		vector3<int> getCell(int unit)
		{	vector3<int> iR;
			iR[2] = unit % phonon.sup[2]; unit /= phonon.sup[2];
			iR[1] = unit % phonon.sup[1]; unit /= phonon.sup[1];
			iR[0] = unit;
			return iR;
		}
	}
	forceMatrixIndexWrap(*this);
	//---  collect forces with symmetric combinations:
	auto sym = symm.getMatrices();
	auto atomMap = symm.getAtomMap();
	std::vector< matrix3<> > symCart(sym.size());
	for(size_t iSym=0; iSym<sym.size(); iSym++)
		symCart[iSym] = eSup.gInfo.R * sym[iSym] * inv(eSup.gInfo.R);;
	std::map< ForceMatrixIndex, matrix3<> > forceMatrix;
	size_t iMode=0;
	for(size_t sp1=0; sp1<e.iInfo.species.size(); sp1++)
	for(size_t at1=0; at1<e.iInfo.species[sp1]->atpos.size(); at1++)
	{	for(size_t sp2=0; sp2<eSup.iInfo.species.size(); sp2++)
		for(size_t at2=0; at2<eSup.iInfo.species[sp2]->atpos.size(); at2++)
		{	//Pickup force for this combination:
			matrix3<> F;
			for(int j=0; j<3; j++)
				F.set_col(j, dgrad[iMode+j][sp2][at2]);
			//Save with all possible rotations:
			for(size_t iSym=0; iSym<sym.size(); iSym++)
			{	ForceMatrixIndex fmi;
				fmi.sp1 = sp1;
				fmi.sp2 = sp2;
				fmi.at1 = atomMap[sp1][at1][iSym];
				fmi.at2 = atomMap[sp1][at2][iSym];
				forceMatrixIndexWrap(fmi);
				forceMatrix[fmi] += symCart[iSym] * F * (~symCart[iSym]);
			}
		}
		iMode += 3;
	}
	//--- reconstitute dgrad:
	iMode=0;
	double nSymInv = 1./sym.size();
	for(size_t sp1=0; sp1<e.iInfo.species.size(); sp1++)
	for(size_t at1=0; at1<e.iInfo.species[sp1]->atpos.size(); at1++)
	{	for(size_t sp2=0; sp2<eSup.iInfo.species.size(); sp2++)
		for(size_t at2=0; at2<eSup.iInfo.species[sp2]->atpos.size(); at2++)
		{	ForceMatrixIndex fmi;
			fmi.sp1 = sp1; fmi.sp2 = sp2;
			fmi.at1 = at1; fmi.at2 = at2;
			const matrix3<>& F = forceMatrix[fmi]; //no wrapping needed for direct loop
			for(int j=0; j<3; j++)
				dgrad[iMode+j][sp2][at2] = nSymInv * F.column(j);
		}
		iMode += 3;
	}
	
	//Construct frequency-squared matrix:
	//--- convert forces to frequency-squared (divide by mass symmetrically):
	std::vector<double> invsqrtM;
	for(auto sp: e.iInfo.species)
		invsqrtM.push_back(1./sqrt(sp->mass * amu));
	for(size_t iMode=0; iMode<modes.size(); iMode++)
	{	dgrad[iMode] *= invsqrtM[modes[iMode].sp]; //divide by sqrt(M) on the left
		for(size_t sp=0; sp<e.iInfo.species.size(); sp++)
			for(vector3<>& f: dgrad[iMode][sp])
				f *= invsqrtM[sp]; // divide by sqrt(M) on the right
	}
	//--- remap to cells:
	std::map<vector3<int>,double> cellMap = getCellMap(e.gInfo.R, eSup.gInfo.R, e.dump.getFilename("phononCellMap"));
	std::vector<matrix> omegaSq(cellMap.size());
	auto iter = cellMap.begin();
	for(size_t iCell=0; iCell<cellMap.size(); iCell++)
	{	//Get cell map entry:
		vector3<int> iR = iter->first;
		double weight = iter->second;
		iter++;
		//Find index of cell in supercell:
		for(int j=0; j<3; j++)
			iR[j] = positiveRemainder(iR[j], sup[j]);
		int cellIndex =  (iR[0]*sup[1] + iR[1])*sup[2] + iR[2]; //corresponding to the order of atom replication in setup()
		//Collect omegSq entries:
		omegaSq[iCell].init(modes.size(), modes.size());
		for(size_t iMode1=0; iMode1<modes.size(); iMode1++)
		{	const IonicGradient& F1 = dgrad[iMode1];
			for(size_t iMode2=0; iMode2<modes.size(); iMode2++)
			{	const Mode& mode2 = modes[iMode2];
				size_t cellOffsetSp = cellIndex * e.iInfo.species[mode2.sp]->atpos.size(); //offset into atoms of current cell for current species
				omegaSq[iCell].set(iMode1,iMode2, weight * F1[mode2.sp][mode2.at + cellOffsetSp][mode2.dir]);
			}
		}
	}
	//--- enforce hermiticity:
	size_t nSymmetrizedCells = 0;
	auto iter1 = cellMap.begin();
	for(size_t iCell1=0; iCell1<cellMap.size(); iCell1++)
	{	auto iter2 = cellMap.begin();
		for(size_t iCell2=0; iCell2<cellMap.size(); iCell2++)
		{	vector3<int> iRsum = iter1->first + iter2->first;
			if(!iRsum.length_squared() && iCell2>=iCell1) //loop over iR1 + iR2 == 0 pairs
			{	matrix M = 0.5*(omegaSq[iCell1] + dagger(omegaSq[iCell2]));
				omegaSq[iCell1] = M;
				omegaSq[iCell2] = dagger(M);
				nSymmetrizedCells += (iCell1==iCell2 ? 1 : 2);
			}
			iter2++;
		}
		iter1++;
	}
	assert(nSymmetrizedCells == cellMap.size());
	//--- enforce translational invariance:
	matrix omegaSqSum;
	for(const matrix& M: omegaSq)
		omegaSqSum += M;
	matrix Fmean; //3x3 force matrix for all atoms moving together at Gamma point (should be zero)
	int nAtoms = modes.size()/3;
	for(int at=0; at<nAtoms; at++)
		Fmean += (1./(nAtoms*prodSup)) * omegaSqSum(3*at,3*(at+1), 3*at,3*(at+1)) * e.iInfo.species[modes[3*at].sp]->mass;
	matrix omegaSqCorrection = zeroes(modes.size(), modes.size());
	for(int at=0; at<nAtoms; at++)
		omegaSqCorrection.set(3*at,3*(at+1), 3*at,3*(at+1), Fmean * (1./e.iInfo.species[modes[3*at].sp]->mass));
	iter = cellMap.begin();
	for(size_t iCell=0; iCell<cellMap.size(); iCell++)
	{	omegaSq[iCell] -= iter->second * omegaSqCorrection;
		iter++;
	}
	//--- write to file
	if(mpiUtil->isHead())
	{	string fname = e.dump.getFilename("phononOmegaSq");
		logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
		FILE* fp = fopen(fname.c_str(), "w");
		for(const matrix& M: omegaSq)
			M.write_real(fp); //M is explicitly real by construction above
		fclose(fp);
		logPrintf("done.\n"); logFlush();
		//Write description of modes:
		fname = e.dump.getFilename("phononBasis");
		logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
		fp = fopen(fname.c_str(), "w");
		fprintf(fp, "#species atom dx dy dz [bohrs]\n");
		for(const Mode& mode: modes)
		{	vector3<> r; r[mode.dir] = invsqrtM[mode.sp];
			fprintf(fp, "%s %d %lg %lg %lg\n", e.iInfo.species[mode.sp]->name.c_str(), mode.at, r[0], r[1], r[2]);
		}
		fclose(fp);
		logPrintf("done.\n"); logFlush();
	}
	
	//Output electron-phonon matrix elements:
	if(mpiUtil->isHead())
	{	int nSpins = Hsub0.size();
		const int& nBands = e.eInfo.nBands;
		for(int s=0; s<nSpins; s++)
		{	string spinSuffix = (nSpins==1 ? "" : (s==0 ? "Up" : "Dn"));
			string fname = e.dump.getFilename("phononHsub" + spinSuffix);
			logPrintf("Dumping '%s' ... ", fname.c_str()); logFlush();
			FILE* fp = fopen(fname.c_str(), "w");
			for(size_t iMode=0; iMode<modes.size(); iMode++)
				for(int ik1=0; ik1<prodSup; ik1++)
					for(int ik2=0; ik2<prodSup; ik2++)
					{	matrix HePh = invsqrtM[modes[iMode].sp] * //incorporate mass factor in eigen-displacement
							dHsub[iMode][s](ik1*nBands,(ik1+1)*nBands, ik2*nBands,(ik2+1)*nBands); //split into k-point pairs
						HePh.write(fp);
					}
			fclose(fp);
			logPrintf("done.\n"); logFlush();
		}
	}

	//Calculate free energy (properly handling singularities at Gamma point):
	std::vector< std::pair<vector3<>,double> > quad = getQuadratureBZ();
	int ikStart = (quad.size() * mpiUtil->iProcess()) / mpiUtil->nProcesses();
	int ikStop = (quad.size() * (mpiUtil->iProcess()+1)) / mpiUtil->nProcesses();
	double ZPE = 0., Evib = 0., Avib = 0.;
	for(int ik=ikStart; ik<ikStop; ik++)
	{	//Calculate phonon omegaSq at current k:
		matrix omegaSq_k;
		iter = cellMap.begin();
		for(size_t iCell=0; iCell<cellMap.size(); iCell++)
		{	omegaSq_k += cis(2*M_PI*dot(iter->first, quad[ik].first)) * omegaSq[iCell];
			iter++;
		}
		//Diagonalize:
		diagMatrix omegaSqEigs; matrix omegaSqEvecs;
		omegaSq_k = dagger_symmetrize(omegaSq_k);
		omegaSq_k.diagonalize(omegaSqEvecs, omegaSqEigs);
		//Collect contributions:
		double w = quad[ik].second; //integration weight of current k-point
		for(double omegaSqEig: omegaSqEigs)
		{	double omega = sqrt(omegaSqEig);
			double expMomegaByT = exp(-omega/T);
			ZPE += w*( 0.5*omega );
			Evib += w*( 0.5*omega + omega * expMomegaByT / (1.-expMomegaByT) );
			Avib += w*( 0.5*omega + T * log(1.-expMomegaByT) );
		}
	}
	mpiUtil->allReduce(ZPE, MPIUtil::ReduceSum);
	mpiUtil->allReduce(Evib, MPIUtil::ReduceSum);
	mpiUtil->allReduce(Avib, MPIUtil::ReduceSum);
	double TSvib = Evib - Avib;
	logPrintf("\nPhonon free energy components (per unit cell) at T = %lg K:\n", T/Kelvin);
	logPrintf("\tZPE:   %15.6lf\n", ZPE);
	logPrintf("\tEvib:  %15.6lf\n", Evib);
	logPrintf("\tTSvib: %15.6lf\n", TSvib);
	logPrintf("\tAvib:  %15.6lf\n", Avib);
	if(isnan(ZPE))
		logPrintf("\tWARNING: free energies are undefined due to imaginary phonon frequencies.\n");
	logPrintf("\n");
}

void Phonon::setSupState(std::vector<matrix>* Hsub)
{
	int prodSup = sup[0]*sup[1]*sup[2];
	int nBandsSup = e.eInfo.nBands * prodSup; //Note >= eSup.eInfo.nBands, depending on e.eInfo.nBands >= nBandsOpt
	
	//Zero wavefunctions and auxiliary Hamiltonia (since a scatter used below)
	for(int qSup=eSup.eInfo.qStart; qSup<eSup.eInfo.qStop; qSup++)
	{	ColumnBundle& Cq = eSup.eVars.C[qSup];
		eSup.eVars.Y[qSup].free(); //to save memory (will be regenerated below, after Hsub calculation)
		if(Cq.nCols() != nBandsSup)
			Cq.init(nBandsSup, Cq.colLength(), Cq.basis, Cq.qnum, isGpuEnabled());
		Cq.zero();
		if(e.eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux)
			eSup.eVars.B[qSup] = zeroes(nBandsSup,nBandsSup);
	}
	
	//Update supercell quantities:
	int nSpinor = e.eInfo.spinorLength();
	double scaleFac = 1./sqrt(prodSup); //to account for normalization
	for(int q=0; q<e.eInfo.nStates; q++) if(eSup.eInfo.isMine(stateMap[q]->qSup))
	{	int nBandsPrev = e.eInfo.nBands * stateMap[q]->nqPrev;
		int qSup = stateMap[q]->qSup;
		//Wavefunctions:
		const ColumnBundle& C = e.eVars.C[q];
		ColumnBundle& Csup = eSup.eVars.C[qSup];
		for(int b=0; b<e.eInfo.nBands; b++)
			for(int s=0; s<nSpinor; s++)
				callPref(eblas_scatter_zdaxpy)(stateMap[q]->nIndices, scaleFac, stateMap[q]->indexPref,
					C.dataPref() + C.index(b,s*C.basis->nbasis),
					Csup.dataPref() + Csup.index(b + nBandsPrev, s*Csup.basis->nbasis));
		//Fillings:
		const diagMatrix& F = e.eVars.F[q];
		diagMatrix& Fsup = eSup.eVars.F[qSup];
		Fsup.resize(nBandsSup);
		Fsup.set(nBandsPrev,nBandsPrev+e.eInfo.nBands, F);
		//Auxiliary Hamiltonian (if necessary):
		if(e.eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux)
		{	const matrix& B = e.eVars.B[q];
			matrix& Bsup = eSup.eVars.B[qSup];
			Bsup.set(nBandsPrev,nBandsPrev+e.eInfo.nBands, nBandsPrev,nBandsPrev+e.eInfo.nBands, B);
		}
	}
	
	//Compute gamma point Hamiltonian if requested:
	if(Hsub)
	{	int nSpins = e.eInfo.spinType==SpinZ ? 2 : 1;
		Hsub->resize(nSpins);
		for(int s=0; s<nSpins; s++)
		{	int qSup = s*(eSup.eInfo.nStates/nSpins); //Gamma point is always first in the list for each spin
			assert(eSup.eInfo.qnums[qSup].k.length_squared() == 0); //make sure that above is true
			if(eSup.eInfo.isMine(qSup))
			{	ColumnBundle HC; Energies ener;
				eSup.iInfo.project(eSup.eVars.C[qSup], eSup.eVars.VdagC[qSup]); //update wavefunction projections
				eSup.eVars.applyHamiltonian(qSup, eye(nBandsSup), HC, ener, true);
				(*Hsub)[s] = eSup.eVars.Hsub[qSup] * prodSup; //account for scaling of wavefunctions above
			}
			else (*Hsub)[s].init(nBandsSup, nBandsSup);
			(*Hsub)[s].bcast(eSup.eInfo.whose(qSup));
		}
	}
	
	//Discard extra bands:
	for(int qSup=eSup.eInfo.qStart; qSup<eSup.eInfo.qStop; qSup++)
	{	ColumnBundle& Cq = eSup.eVars.C[qSup];
		ColumnBundle& Yq = eSup.eVars.Y[qSup];
		diagMatrix& Fq = eSup.eVars.F[qSup];
		matrix& Bq = eSup.eVars.B[qSup];
		Yq = Cq.similar(eSup.eInfo.nBands);
		diagMatrix Ftmp(eSup.eInfo.nBands);
		matrix Btmp = zeroes(eSup.eInfo.nBands, eSup.eInfo.nBands);
		for(int nqPrev=0; nqPrev<prodSup; nqPrev++)
		{	int offsIn = nqPrev * e.eInfo.nBands;
			int offsOut = nqPrev * nBandsOpt;
			callPref(eblas_copy)(Yq.data()+Yq.index(offsOut,0), Cq.data()+Cq.index(offsIn,0), nBandsOpt*Yq.colLength());
			Ftmp.set(offsOut,offsOut+nBandsOpt, Fq(offsIn,offsIn+nBandsOpt));
			if(e.eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux)
				Btmp.set(offsOut,offsOut+nBandsOpt, offsOut,offsOut+nBandsOpt, Bq(offsIn,offsIn+nBandsOpt, offsIn,offsIn+nBandsOpt));
		}
		Cq = Yq;
		std::swap(Fq, Ftmp);
		if(e.eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux)
			std::swap(Bq, Btmp);
		eSup.iInfo.project(Cq, eSup.eVars.VdagC[qSup]); //update wave function projections
	}
	
	//Update entropy contributions:
	if(eSup.eInfo.fillingsUpdate != ElecInfo::ConstantFillings)
		eSup.eInfo.updateFillingsEnergies(eSup.eVars.F, eSup.ener);
	
	if(e.eInfo.fillingsUpdate==ElecInfo::FermiFillingsAux)
		eSup.eVars.HauxInitialized = true;
}

Phonon::StateMapEntry::StateMapEntry() : nIndices(0), indexPref(0), index(0)
#ifdef GPU_ENABLED
, indexGpu(0)
#endif
{

}

Phonon::StateMapEntry::~StateMapEntry()
{
	if(index) delete[] index;
	#ifdef GPU_ENABLED
	if(indexGpu) cudaFree(indexGpu);
	#endif
}

void Phonon::StateMapEntry::setIndex(const std::vector<int>& indexVec)
{	nIndices = indexVec.size();
	//CPU version
	if(index) delete[] index;
	index = new int[nIndices];
	eblas_copy(index, indexVec.data(), nIndices);
	indexPref = index;
	//GPU version
	#ifdef GPU_ENABLED
	if(indexGpu) cudaFree(indexGpu);
	cudaMalloc(&indexGpu, sizeof(int)*nIndices); gpuErrorCheck();
	cudaMemcpy(indexGpu, index, sizeof(int)*nIndices, cudaMemcpyHostToDevice); gpuErrorCheck();
	indexPref = indexGpu;
	#endif
}
