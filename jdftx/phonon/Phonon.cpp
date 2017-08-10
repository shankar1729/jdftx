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
#include <core/WignerSeitz.h>

void Phonon::dump()
{	//Zero force matrix and electron-phonon matrix elements:
	IonicGradient zeroForce; zeroForce.init(eSupTemplate.iInfo);
	dgrad.assign(modes.size(), zeroForce);
	dHsub.assign(modes.size(), std::vector<matrix>(nSpins));
	
	//Accumulate contributions to force matrix and electron-phonon matrix elements for each irreducible perturbation:
	unsigned iPertStart = (iPerturbation>=0) ? iPerturbation : 0;
	unsigned iPertStop  = (iPerturbation>=0) ? iPerturbation+1 : perturbations.size();
	std::vector<int> nStatesPert(perturbations.size());
	for(unsigned iPert=iPertStart; iPert<iPertStop; iPert++)
	{	logPrintf("########### Perturbed supercell calculation %u of %d #############\n", iPert+1, int(perturbations.size()));
		ostringstream oss; oss << "phonon." << iPert+1 << ".$@#!"; //placeholder for $VAR
		string fnamePattern = e.dump.getFilename(oss.str()); //(because dump variable name cannot contain $VAR)
		fnamePattern.replace(fnamePattern.find("$@#!"), 4, "$VAR"); //replace placeholder with $VAR
		processPerturbation(perturbations[iPert], fnamePattern);
		nStatesPert[iPert] = eSup->eInfo.nStates;
		logPrintf("\n"); logFlush();
	}
	if(dryRun)
	{	logPrintf("\nParameter summary for supercell calculations:\n");
		for(unsigned iPert=iPertStart; iPert<iPertStop; iPert++)
			logPrintf("\tPerturbation: %u  nStates: %d\n", iPert+1, nStatesPert[iPert]);
		logPrintf("Use option iPerturbation of command phonon to run each supercell calculation separately.\n");
		return;
	}
	if(iPerturbation>=0)
	{	logPrintf("Completed supercell calculation for iPerturbation %d.\n", iPerturbation+1);
		logPrintf("After completing all supercells, rerun with option collectPerturbations in command phonon.\n");
		return;
	}
	
	//Generate phonon cell map:
	std::vector<vector3<>> xAtoms;  //lattice coordinates of all atoms in order
	for(const auto& sp: e.iInfo.species)
		xAtoms.insert(xAtoms.end(), sp->atpos.begin(), sp->atpos.end());
	assert(3*xAtoms.size() == modes.size());
	std::map<vector3<int>,matrix> cellMap = getCellMap(e.gInfo.R, eSupTemplate.gInfo.R,
		e.coulombParams.isTruncated(), xAtoms, xAtoms, rSmooth, e.dump.getFilename("phononCellMap"));
	
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
	std::vector<matrix> omegaSq(cellMap.size());
	auto iter = cellMap.begin();
	for(size_t iCell=0; iCell<cellMap.size(); iCell++)
	{	//Get cell map entry:
		vector3<int> iR = iter->first;
		const matrix& weight = iter->second;
		iter++;
		//Find index of cell in supercell:
		for(int j=0; j<3; j++)
			iR[j] = positiveRemainder(iR[j], sup[j]);
		int cellIndex =  (iR[0]*sup[1] + iR[1])*sup[2] + iR[2]; //corresponding to the order of atom replication in setup()
		//Collect omegaSq entries:
		omegaSq[iCell].init(modes.size(), modes.size());
		for(size_t iMode1=0; iMode1<modes.size(); iMode1++)
		{	const IonicGradient& F1 = dgrad[iMode1];
			for(size_t iMode2=0; iMode2<modes.size(); iMode2++)
			{	const Mode& mode2 = modes[iMode2];
				size_t cellOffsetSp = cellIndex * e.iInfo.species[mode2.sp]->atpos.size(); //offset into atoms of current cell for current species
				omegaSq[iCell].set(iMode1,iMode2, weight(iMode1/3,iMode2/3) * dot(F1[mode2.sp][mode2.at + cellOffsetSp], mode2.dir));
			}
		}
	}
	
	logPrintf("\nRefining force matrix:\n");
	forceMatrixDaggerSymmetrize(omegaSq, cellMap); //enforce hermiticity
	forceMatrixEnforceSumRule(omegaSq, cellMap, invsqrtM); //enforce translational invariance
	logPrintf("\n");
	
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
		{	vector3<> r = mode.dir * invsqrtM[mode.sp];
			fprintf(fp, "%s %d  %+lf %+lf %+lf\n", e.iInfo.species[mode.sp]->name.c_str(), mode.at, r[0], r[1], r[2]);
		}
		fclose(fp);
		logPrintf("done.\n"); logFlush();
	}
	
	//Output electron-phonon matrix elements:
	if(mpiUtil->isHead())
	{	const int& nBands = e.eInfo.nBands;
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
	std::vector< std::pair<vector3<>,double> > getQuadratureBZ(vector3<bool>); //implemented below
	std::vector< std::pair<vector3<>,double> > quad = getQuadratureBZ(e.coulombParams.isTruncated());
	int ikStart, ikStop;
	TaskDivision(quad.size(), mpiUtil).myRange(ikStart, ikStop);
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
	if(std::isnan(ZPE))
		logPrintf("\tWARNING: free energies are undefined due to imaginary phonon frequencies.\n");
	logPrintf("\n");
}

vector3<int> Phonon::getCell(int unit) const
{	vector3<int> cell;
	cell[2] = unit % sup[2]; unit /= sup[2];
	cell[1] = unit % sup[1]; unit /= sup[1];
	cell[0] = unit;
	return cell;
}

int Phonon::getUnit(const vector3<int>& cell) const
{	return
		( positiveRemainder(cell[0],sup[0]) * sup[1]
		+ positiveRemainder(cell[1],sup[1]) ) * sup[2]
		+ positiveRemainder(cell[2],sup[2]);
}


//Generate quadrature for BZ integrals with singularity at Gamma point:
//--- add contribution to quadrature due to box of size scale centered at gCenter (in reciprocal lattice coords)
void addQuadratureBZ_box(std::vector< std::pair<vector3<>,double> >& quad, double scale, vector3<> gCenter, vector3<int> dim)
{	//Weights and abscissae of the 7-point gauss quadrature:
	const int N = 3;
	static const double w[N+1] = { 0.129484966168869693270611432679082, 0.279705391489276667901467771423780, 0.381830050505118944950369775488975, 0.417959183673469387755102040816327 };
	static const double x[N+1] = { 0.949107912342758524526189684047851, 0.741531185599394439863864773280788, 0.405845151377397166906606412076961, 0.000000000000000000000000000000000 };
	double h = 0.5 * scale; 
	vector3<> g;
	vector3<int> i;
	#define LOOP(j) for(i[j]=-N*dim[j]; i[j]<=N*dim[j]; i[j]++)
	#define ABSCISSA(j) gCenter[j] + (dim[j] ? (h*x[N-abs(i[j])]*(i[j]>0?1:-1)) : 0)
	#define WEIGHT(j) (dim[j] ? h*w[N-abs(i[j])] : 1.)
	LOOP(0)
	{	g[0] = ABSCISSA(0);
		double w0 = WEIGHT(0);
		LOOP(1)
		{	g[1] = ABSCISSA(1);
			double w01 = w0 * WEIGHT(1);
			LOOP(2)
			{	g[2] = ABSCISSA(2);
				double w012 = w01 * WEIGHT(2);
				quad.push_back(std::make_pair(g, w012));
			}
		}
	}
	#undef LOOP
	#undef ABSCISSA
	#undef WEIGHT
}
//--- add contribution to quadrature between boxes of size scale and scale/3
//--- dim is 0 or 1 depending on whether that dimension should be sampled
void addQuadratureBZ_scale(std::vector< std::pair<vector3<>,double> >& quad, double scale, vector3<int> dim)
{	double scaleBy3 = scale/3.;
	vector3<int> ig;
	for(ig[0]=-dim[0]; ig[0]<=dim[0]; ig[0]++)
	for(ig[1]=-dim[1]; ig[1]<=dim[1]; ig[1]++)
	for(ig[2]=-dim[2]; ig[2]<=dim[2]; ig[2]++)
		if(ig.length_squared()) //except the center box
			addQuadratureBZ_box(quad, scaleBy3, scaleBy3*ig, dim);
}
std::vector< std::pair<vector3<>,double> > getQuadratureBZ(vector3<bool> isTruncated)
{	vector3<int> dim;
	for(int j=0; j<3; j++)
		dim[j] = (isTruncated[j] ? 0 : 1);
	std::vector< std::pair<vector3<>,double> > quad;
	for(double scale=1.; scale>1e-3; scale/=3.)
		addQuadratureBZ_scale(quad, scale, dim);
	return quad;
}

//-------------- class BlockRotationMatrix ---------------

void BlockRotationMatrix::init(int nBlocks, int blockSize)
{	this->nBlocks = nBlocks;
	this->blockSize = blockSize;
	colOffset.assign(nBlocks, -1);
	rots.assign(nBlocks, matrix());
}

void BlockRotationMatrix::set(int rowBlock, int colBlock, const matrix& rot)
{	assert(rowBlock >= 0 && rowBlock < nBlocks);
	assert(colBlock >= 0 && colBlock < nBlocks);
	assert(rot.nRows() == blockSize);
	assert(rot.nCols() == blockSize);
	colOffset[rowBlock] = colBlock;
	rots[rowBlock] = rot;
}

void BlockRotationMatrix::allReduce()
{	for(int rowBlock=0; rowBlock<nBlocks; rowBlock++)
	{	//Make sure exactly one process has rotation:
		bool myRot = colOffset[rowBlock]>=0;
		int haveRot = (myRot ? 1 : 0);
		mpiUtil->allReduce(haveRot, MPIUtil::ReduceSum);
		assert(haveRot == 1);
		//Reduce:
		mpiUtil->allReduce(colOffset[rowBlock], MPIUtil::ReduceMax);
		if(!myRot) rots[rowBlock] = zeroes(blockSize, blockSize);
		rots[rowBlock].allReduce(MPIUtil::ReduceSum);
	}
}

matrix BlockRotationMatrix::transform(const matrix& in) const
{	int matSize = blockSize * nBlocks;
	assert(in.nRows() == matSize);
	assert(in.nCols() == matSize);
	//First perform the left multiply by rot:
	matrix temp = zeroes(matSize, matSize);
	for(int rowBlock=0; rowBlock<nBlocks; rowBlock++)
	{	int colBlock = colOffset[rowBlock];
		temp.set(rowBlock*blockSize,(rowBlock+1)*blockSize, 0,matSize,
			rots[rowBlock] * in(colBlock*blockSize,(colBlock+1)*blockSize, 0,matSize));
	}
	//Next perform the right multiply by dagger(rot):
	matrix out = zeroes(matSize, matSize);
	for(int rowBlock=0; rowBlock<nBlocks; rowBlock++)
	{	int colBlock = colOffset[rowBlock];
		out.set(0,matSize, rowBlock*blockSize,(rowBlock+1)*blockSize,
			temp(0,matSize, colBlock*blockSize,(colBlock+1)*blockSize) * dagger(rots[rowBlock]));
	}
	return out;
}

//Enforce hermiticity of force matrix
void Phonon::forceMatrixDaggerSymmetrize(std::vector<matrix>& omegaSq, const std::map<vector3<int>,matrix>& cellMap) const
{	size_t nSymmetrizedCellsTot = 0;
	double hermErrNum = 0., hermErrDen = 0.;
	auto iter1 = cellMap.begin();
	for(size_t iCell1=0; iCell1<cellMap.size(); iCell1++)
	{	auto iter2 = cellMap.begin();
		for(size_t iCell2=0; iCell2<cellMap.size(); iCell2++)
		{	vector3<int> iRsum = iter1->first + iter2->first;
			if(!iRsum.length_squared() && iCell2>=iCell1) //loop over iR1 + iR2 == 0 pairs
			{	matrix M = 0.5*(omegaSq[iCell1] + dagger(omegaSq[iCell2]));
				matrix Merr = 0.5*(omegaSq[iCell1] - dagger(omegaSq[iCell2]));
				omegaSq[iCell1] = M;
				omegaSq[iCell2] = dagger(M);
				int nSymmetrizedCells = (iCell1==iCell2 ? 1 : 2);
				nSymmetrizedCellsTot += nSymmetrizedCells;
				hermErrNum += nSymmetrizedCells * std::pow(nrm2(Merr), 2);
				hermErrDen += nSymmetrizedCells * std::pow(nrm2(M), 2);
			}
			iter2++;
		}
		iter1++;
	}
	assert(nSymmetrizedCellsTot == cellMap.size());
	logPrintf("\tCorrected hermiticity relative error: %lg\n", sqrt(hermErrNum/hermErrDen));
}

//Enforce translational invariance sum rule on force matrix
void Phonon::forceMatrixEnforceSumRule(std::vector<matrix>& omegaSq, const std::map<vector3<int>,matrix>& cellMap, const std::vector<double>& invsqrtM) const
{	//Collect masses by mode:
	diagMatrix sqrtMmode, invsqrtMmode;
	for(const Mode& mode: modes)
	{	invsqrtMmode.push_back(invsqrtM[mode.sp]);
		sqrtMmode.push_back(1./invsqrtM[mode.sp]);
	}
	//Collect Gamma-point force matrix:
	matrix F0;
	for(const matrix& mat: omegaSq)
		F0 += sqrtMmode * mat * sqrtMmode;
	//Project out net force for each atom displacement:
	F0 = dagger_symmetrize(F0);
	matrix dF0 = zeroes(modes.size(), modes.size());
	int nAtoms = modes.size()/3;
	const complex* F0data = F0.data();
	complex* dF0data = dF0.data();
	for(int iAtom=0; iAtom<nAtoms; iAtom++)
	for(int iDir=0; iDir<3; iDir++)
	{	int iMode = 3*iAtom+iDir;
		for(int jAtom=0; jAtom<nAtoms; jAtom++)
		for(int jDir=0; jDir<3; jDir++)
		{	int jMode = 3*jAtom+jDir;
			dF0data[dF0.index(iMode,3*iAtom+jDir)] -= F0data[F0.index(iMode,jMode)];
		} //Note this makes dF0 block diagonal with block size 3
	}
	//Apply correction:
	auto iter = cellMap.begin();
	for(size_t iCell=0; iCell<cellMap.size(); iCell++)
	{	if(!iter->first.length_squared()) //apply correction to diagonal elements
			omegaSq[iCell] += (invsqrtMmode * dF0 * invsqrtMmode);
		iter++;
	}
	logPrintf("\tCorrected translational invariance relative error: %lg\n", nrm2(dF0)/nrm2(F0));
}
