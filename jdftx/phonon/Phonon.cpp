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
	
	//Process force matrix:
	//--- refine in reciprocal space
	dgradSymmetrize(dgrad);
	
	//--- generate phonon cell map:
	std::vector<vector3<>> xAtoms;  //lattice coordinates of all atoms in order
	for(const auto& sp: e.iInfo.species)
		xAtoms.insert(xAtoms.end(), sp->atpos.begin(), sp->atpos.end());
	assert(3*xAtoms.size() == modes.size());
	std::map<vector3<int>,matrix> cellMap = getCellMap(e.gInfo.R, eSupTemplate.gInfo.R,
		e.coulombParams.isTruncated(), xAtoms, xAtoms, rSmooth, e.dump.getFilename("phononCellMap"));
	
	//--- remap force matrix on cells (with duplicated ones on WS-supercell boundary):
	std::vector<matrix> F(cellMap.size());
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
		F[iCell].init(modes.size(), modes.size());
		for(size_t iMode1=0; iMode1<modes.size(); iMode1++)
		{	const IonicGradient& F1 = dgrad[iMode1];
			for(size_t iMode2=0; iMode2<modes.size(); iMode2++)
			{	const Mode& mode2 = modes[iMode2];
				size_t cellOffsetSp = cellIndex * e.iInfo.species[mode2.sp]->atpos.size(); //offset into atoms of current cell for current species
				F[iCell].set(iMode1,iMode2, weight(iMode1/3,iMode2/3) * dot(F1[mode2.sp][mode2.at + cellOffsetSp], mode2.dir));
			}
		}
	}
	
	//--- check force matrix
	logPrintf("\nFinalizing force matrix in real space:\n");
	forceMatrixGammaPrimeFix(F, cellMap, xAtoms); //fix first-derivative behavior near Gamma
	forceMatrixHermCheck(F, cellMap); //final hermiticity check
	forceMatrixSumRuleCheck(F, cellMap); //final translational invariance check
	logPrintf("\n");
	
	//--- collect species and mode masses:
	std::vector<double> invsqrtM; //by species
	for(auto sp: e.iInfo.species)
		invsqrtM.push_back(1./sqrt(sp->mass * amu));
	diagMatrix invsqrtMmode; //by mode
	for(const Mode& mode: modes)
		invsqrtMmode.push_back(invsqrtM[mode.sp]);
	
	//--- convert force matrix to omegaSq (divide symmetrically by sqrt(M)):
	std::vector<matrix> omegaSq;
	for(const matrix& Fi: F)
		omegaSq.push_back(invsqrtMmode * Fi * invsqrtMmode);
	
	//--- write to file
	if(mpiWorld->isHead())
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
		fprintf(fp, "#species atom dx[bohrs] dy[bohrs] dz[bohrs] M[amu]\n");
		for(const Mode& mode: modes)
		{	vector3<> r = mode.dir * invsqrtM[mode.sp];
			fprintf(fp, "%s %d  %+lf %+lf %+lf  %lf\n",
				e.iInfo.species[mode.sp]->name.c_str(), mode.at,
				r[0], r[1], r[2], e.iInfo.species[mode.sp]->mass);
		}
		fclose(fp);
		logPrintf("done.\n"); logFlush();
	}
	
	//Output electron-phonon matrix elements:
	if(mpiWorld->isHead() && saveHsub)
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
	TaskDivision(quad.size(), mpiWorld).myRange(ikStart, ikStop);
	double ZPE = 0., Evib = 0., Avib = 0.;
	double kSqMaxIm = 0.; //maximum |k|^2 with imaginary frequencies (if any)
	double omegaSqMin = 0.; int ikMin=-1; //most negative frequency squared and corresponding k index
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
		//Check for imaginary frequencies:
		if(omegaSqEigs[0] <= 0.) //eigenvalues in ascending order
		{	kSqMaxIm = std::max(kSqMaxIm, e.gInfo.GGT.metric_length_squared(quad[ik].first));
			if(omegaSqEigs[0] < omegaSqMin)
			{	omegaSqMin = omegaSqEigs[0];
				ikMin = ik;
			}
		}
		//Collect contributions:
		double w = quad[ik].second; //integration weight of current k-point
		for(double omegaSqEig: omegaSqEigs)
		{	if(omegaSqEig <= 0.) continue; //ignore imag frequency in free energy
			double omega = sqrt(omegaSqEig);
			double expMomegaByT = exp(-omega/T);
			ZPE += w*( 0.5*omega );
			Evib += w*( 0.5*omega + omega * expMomegaByT / (1.-expMomegaByT) );
			Avib += w*( 0.5*omega + T * log(1.-expMomegaByT) );
		}
	}
	mpiWorld->allReduce(ZPE, MPIUtil::ReduceSum);
	mpiWorld->allReduce(Evib, MPIUtil::ReduceSum);
	mpiWorld->allReduce(Avib, MPIUtil::ReduceSum);
	mpiWorld->allReduce(kSqMaxIm, MPIUtil::ReduceMax);
	mpiWorld->allReduce(omegaSqMin, ikMin, MPIUtil::ReduceMin);
	double TSvib = Evib - Avib;
	logPrintf("\nPhonon free energy components (per unit cell) at T = %lg K:\n", T/Kelvin);
	logPrintf("\tZPE:   %15.6lf\n", ZPE);
	logPrintf("\tEvib:  %15.6lf\n", Evib);
	logPrintf("\tTSvib: %15.6lf\n", TSvib);
	logPrintf("\tAvib:  %15.6lf\n", Avib);
	if(kSqMaxIm)
	{	logPrintf("\tWARNING: discarded imaginary frequencies at max |k| = %lg bohr^-1\n", sqrt(kSqMaxIm));
		const vector3<>& kMin = quad[ikMin].first;
		logPrintf("\tStrongest imaginary |omega|: %.6lf at k: [ %lg %lg %lg ]\n",
			sqrt(-omegaSqMin), kMin[0], kMin[1], kMin[2]);
	}
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
		mpiWorld->allReduce(haveRot, MPIUtil::ReduceSum);
		assert(haveRot == 1);
		//Reduce:
		mpiWorld->allReduce(colOffset[rowBlock], MPIUtil::ReduceMax);
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

//------------ force matrix symmetrization / check routines ---------------

//Enforce hermitian and translation invariance symmetry on dgrad (apply in reciprocal space)
void Phonon::dgradSymmetrize(std::vector<IonicGradient>& dgrad) const
{	logPrintf("Refining force matrix in reciprocal space:\n"); logFlush();
	int nModes= modes.size();
	int nModesSq = std::pow(nModes,2);
	//Collect into a nModesSq x nCells matrix:
	matrix F(nModesSq, prodSup);
	complex* Fdata = F.data();
	for(int iCell=0; iCell<prodSup; iCell++)
	{	for(int iMode=0; iMode<nModes; iMode++)
		{	int jMode = 0;
			for(const std::vector<vector3<>>& dgradSp: dgrad[iMode])
			{	int nAtomsSp = dgradSp.size()/prodSup;
				for(int jAtom=0; jAtom<nAtomsSp; jAtom++)
				{	const vector3<>& f = dgradSp[iCell*nAtomsSp+jAtom];
					for(int jDir=0; jDir<3; jDir++)
						Fdata[F.index(iMode*nModes+jMode++, iCell)] = f[jDir];
				}
			}
		}
	}
	//Calculate Fourier transform phase:
	matrix phase(prodSup, prodSup);
	matrix3<> invDiagSup = inv(Diag(vector3<>(sup)));
	complex* phaseData = phase.data();
	vector3<int> iR;
	for(iR[0]=0; iR[0]<sup[0]; iR[0]++)
	for(iR[1]=0; iR[1]<sup[1]; iR[1]++)
	for(iR[2]=0; iR[2]<sup[2]; iR[2]++)
	{	vector3<> xR_2pi = (2.*M_PI) * (invDiagSup * iR);
		vector3<int> ik;
		for(ik[0]=0; ik[0]<sup[0]; ik[0]++)
		for(ik[1]=0; ik[1]<sup[1]; ik[1]++)
		for(ik[2]=0; ik[2]<sup[2]; ik[2]++)
			*(phaseData++) = cis(dot(ik, xR_2pi));
	}
	//Forward transform:
	matrix Ftilde = F * phase;
	//Apply corrections per k:
	double Fnorm = 0., dFtransNorm = 0., dFhermNorm = 0.;
	for(int ik=0; ik<prodSup; ik++)
	{	matrix Fk = Ftilde(0,nModesSq, ik,ik+1);
		Fk.reshape(nModes, nModes);
		Fnorm += std::pow(nrm2(Fk),2);
		//Hermiticity correction:
		{	matrix FkNew = dagger_symmetrize(Fk);
			dFhermNorm += std::pow(nrm2(FkNew-Fk),2);
			Fk = FkNew;
		}
		//Translational invariance correction:
		if(ik==0) //Gamma point
		{	matrix mask = zeroes(nModes,3);
			int nAtoms = nModes/3;
			for(int iAtom=0; iAtom<nAtoms; iAtom++)
				for(int iDir=0; iDir<3; iDir++)
					mask.set(iAtom*3+iDir,iDir, 1.);
			matrix proj = (1./nAtoms) * (mask * dagger(mask));
			matrix dF0 = proj*Fk*proj - Fk*proj - proj*Fk; //double sum - row sum - column sum
			dFtransNorm = std::pow(nrm2(dF0),2);
			Fk += dF0;
		}
		//Store corrected version
		Fk.reshape(nModesSq, 1);
		Ftilde.set(0,nModesSq, ik,ik+1, Fk);
	}
	//Inverse transform:
	F = (Ftilde * dagger(phase)) * (1./prodSup);
	//Convert back to atoms:
	Fdata = F.data();
	for(int iCell=0; iCell<prodSup; iCell++)
	{	for(int iMode=0; iMode<nModes; iMode++)
		{	int jMode = 0;
			for(std::vector<vector3<>>& dgradSp: dgrad[iMode])
			{	int nAtomsSp = dgradSp.size()/prodSup;
				for(int jAtom=0; jAtom<nAtomsSp; jAtom++)
				{	vector3<>& f = dgradSp[iCell*nAtomsSp+jAtom];
					for(int jDir=0; jDir<3; jDir++)
						f[jDir] = Fdata[F.index(iMode*nModes+jMode++, iCell)].real();
				}
			}
		}
	}
	logPrintf("\tCorrected hermiticity relative error: %lg\n", sqrt(dFhermNorm/Fnorm));
	logPrintf("\tCorrected translational invariance relative error: %lg\n\n", sqrt(dFtransNorm/Fnorm));
}

//Fix first-derivative behavior near Gamma (to ensure finite sound velocities)
void Phonon::forceMatrixGammaPrimeFix(std::vector<matrix>& F, const std::map<vector3<int>,matrix>& cellMap, const std::vector<vector3<>>& xAtoms) const
{	double Fnorm=0., dFnorm=0.;
	//Collect the GammaPrime error and various moments of the weights:
	int nAtoms = modes.size()/3;
	std::vector<std::vector<vector3<>>> Err(3, std::vector<vector3<>>(3)); //Gamma prime errors
	std::vector<std::vector<vector3<>>> V(nAtoms, std::vector<vector3<>>(nAtoms)); //Sum(w.r)
	matrix3<> rho; //Sum(w.r.r)
	auto iter = cellMap.begin();
	for(size_t iCell=0; iCell<cellMap.size(); iCell++,iter++)
	{	const matrix& weight = iter->second;
		for(int iAtom=0; iAtom<nAtoms; iAtom++)
		for(int jAtom=0; jAtom<nAtoms; jAtom++)
		{	vector3<> dr = e.gInfo.R * (iter->first + xAtoms[jAtom] - xAtoms[iAtom]); //vector distance between atoms
			//Collect Gamma prime error:
			for(int iDir=0; iDir<3; iDir++)
			for(int jDir=0; jDir<3; jDir++)
				Err[iDir][jDir] += dr * F[iCell](3*iAtom+iDir,3*jAtom+jDir).real();
			//Collect weight moments:
			double w = weight(iAtom,jAtom).real();
			V[iAtom][jAtom] += w * dr;
			rho += w * outer(dr, dr);
		}
	}
	//Pre-account for Gamma-prime error due to subsequent translation-invariance correction:
	std::vector<vector3<>> Vr(nAtoms), Vc(nAtoms); //row and column means of V
	vector3<> Vrc; //overall mean of V
	double nAtomsInv = 1./nAtoms;
	for(int iAtom=0; iAtom<nAtoms; iAtom++)
		for(int jAtom=0; jAtom<nAtoms; jAtom++)
		{	vector3<> wV = nAtomsInv * V[iAtom][jAtom];
			Vr[iAtom] += wV;
			Vc[jAtom] += wV;
			Vrc += nAtomsInv * wV;
		}
	double prodSupInv = 1./prodSup;
	matrix3<> T; //translation invariance correction to rho
	for(int iAtom=0; iAtom<nAtoms; iAtom++)
	for(int jAtom=0; jAtom<nAtoms; jAtom++)
	{	vector3<> Vcur = V[iAtom][jAtom];
		vector3<> Ucur = prodSupInv*(Vrc - Vr[iAtom] - Vc[jAtom]);
		T += outer(Vcur, Ucur);
	}
	//Invert weight sums to calculate correction tensor:
	matrix3<> alphaInv = rho + T;
	vector3<bool> isTruncated = e.coulombParams.isTruncated();
	for(int iDir=0; iDir<3; iDir++)
		if(isTruncated[iDir])
			alphaInv(iDir,iDir) = 1.; //handle singularity in truncated directions
	matrix3<> alpha = inv(alphaInv);
	std::vector<std::vector<vector3<>>> g(3, std::vector<vector3<>>(3)); //correction tensor
	for(int iDir=0; iDir<3; iDir++)
	for(int jDir=0; jDir<3; jDir++)
		g[iDir][jDir] = -(alpha * Err[iDir][jDir]);
	//Apply correction:
	matrix dFsum;
	iter = cellMap.begin();
	for(size_t iCell=0; iCell<cellMap.size(); iCell++,iter++)
	{	const matrix& weight = iter->second;
		matrix dF(modes.size(), modes.size());
		for(int iAtom=0; iAtom<nAtoms; iAtom++)
		for(int jAtom=0; jAtom<nAtoms; jAtom++)
		{	vector3<> dr = e.gInfo.R * (iter->first + xAtoms[jAtom] - xAtoms[iAtom]); //vector distance between atoms
			double w = weight(iAtom,jAtom).real();
			for(int iDir=0; iDir<3; iDir++)
			for(int jDir=0; jDir<3; jDir++)
				dF.set(3*iAtom+iDir,3*jAtom+jDir, w * dot(dr, g[iDir][jDir]));
		}
		Fnorm += std::pow(nrm2(F[iCell]),2);
		dFnorm += std::pow(nrm2(dF),2);
		F[iCell] += dF;
		dFsum += dF;
	}
	//Restore exact translational invariance:
	matrix mask = zeroes(modes.size(),3);
	for(int iAtom=0; iAtom<nAtoms; iAtom++)
		for(int iDir=0; iDir<3; iDir++)
			mask.set(iAtom*3+iDir,iDir, 1.);
	matrix proj = (1./nAtoms) * (mask * dagger(mask));
	matrix dF0 = (1./prodSup) * (proj*dFsum*proj - dFsum*proj - proj*dFsum); //double sum - row sum - column sum
	iter = cellMap.begin();
	for(size_t iCell=0; iCell<cellMap.size(); iCell++,iter++)
	{	matrix dFcur(modes.size(), modes.size());
		const matrix& weight = iter->second;
		int iMode = 0;
		for(int iAtom=0; iAtom<nAtoms; iAtom++)
		for(int iDir=0; iDir<3; iDir++)
		{	int jMode = 0;
			for(int jAtom=0; jAtom<nAtoms; jAtom++)
			for(int jDir=0; jDir<3; jDir++)
			{	dFcur.set(iMode, jMode, dF0(iMode,jMode) * weight(iAtom,jAtom));
				jMode++;
			}
			iMode++;
		}
		F[iCell] += dFcur;
		dFnorm += std::pow(nrm2(dFcur),2);
	}
	logPrintf("\tCorrected gamma-point derivative relative error: %lg\n", sqrt(dFnorm/Fnorm));
}

//Check hermiticity of force matrix
void Phonon::forceMatrixHermCheck(const std::vector<matrix>& F, const std::map<vector3<int>,matrix>& cellMap) const
{	size_t nSymmetrizedCellsTot = 0;
	double hermErrNum = 0., hermErrDen = 0.;
	auto iter1 = cellMap.begin();
	for(size_t iCell1=0; iCell1<cellMap.size(); iCell1++)
	{	auto iter2 = cellMap.begin();
		for(size_t iCell2=0; iCell2<cellMap.size(); iCell2++)
		{	vector3<int> iRsum = iter1->first + iter2->first;
			if(!iRsum.length_squared() && iCell2>=iCell1) //loop over iR1 + iR2 == 0 pairs
			{	matrix M = 0.5*(F[iCell1] + dagger(F[iCell2]));
				matrix Merr = 0.5*(F[iCell1] - dagger(F[iCell2]));
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
	logPrintf("\tHermiticity relative error: %lg\n", sqrt(hermErrNum/hermErrDen));
}

//Check translational invariance sum rule of force matrix
void Phonon::forceMatrixSumRuleCheck(const std::vector<matrix>& F, const std::map<vector3<int>,matrix>& cellMap) const
{	//Collect Gamma-point force matrix:
	matrix F0; double Fnorm=0.;
	for(const matrix& Fi: F)
	{	F0 += Fi;
		Fnorm += std::pow(nrm2(Fi),2);
	}
	//Calculate correction matrix:
	int nModes = modes.size();
	matrix mask = zeroes(nModes,3);
	int nAtoms = nModes/3;
	for(int iAtom=0; iAtom<nAtoms; iAtom++)
		for(int iDir=0; iDir<3; iDir++)
			mask.set(iAtom*3+iDir,iDir, 1.);
	matrix proj = (1./nAtoms) * (mask * dagger(mask));
	matrix dF0 = (1./prodSup) * (proj*F0*proj - F0*proj - proj*F0); //double sum - row sum - column sum
	logPrintf("\tTranslational invariance relative error: %lg\n", nrm2(dF0)/sqrt(Fnorm/prodSup));
}
