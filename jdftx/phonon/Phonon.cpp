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

void Phonon::dump()
{	//Zero force matrix and electron-phonon matrix elements:
	IonicGradient zeroForce; zeroForce.init(eSupTemplate.iInfo);
	dgrad.assign(modes.size(), zeroForce);
	dHsub.assign(modes.size(), std::vector<matrix>(nSpins));
	
	//Accumulate contributions to force matrix and electron-phonon matrix elements for each irreducible perturbation:
	for(unsigned iPert=0; iPert<perturbations.size(); iPert++)
	{	logPrintf("########### Perturbed supercell calculation %u of %d #############\n", iPert+1, int(perturbations.size()));
		processPerturbation(perturbations[iPert]);
		logPrintf("\n"); logFlush();
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
	std::map<vector3<int>,double> cellMap = getCellMap(e.gInfo.R, eSupTemplate.gInfo.R, e.dump.getFilename("phononCellMap"));
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
		//Collect omegaSq entries:
		omegaSq[iCell].init(modes.size(), modes.size());
		for(size_t iMode1=0; iMode1<modes.size(); iMode1++)
		{	const IonicGradient& F1 = dgrad[iMode1];
			for(size_t iMode2=0; iMode2<modes.size(); iMode2++)
			{	const Mode& mode2 = modes[iMode2];
				size_t cellOffsetSp = cellIndex * e.iInfo.species[mode2.sp]->atpos.size(); //offset into atoms of current cell for current species
				omegaSq[iCell].set(iMode1,iMode2, weight * dot(F1[mode2.sp][mode2.at + cellOffsetSp], mode2.dir));
			}
		}
	}
	//--- enforce hermiticity:
	size_t nSymmetrizedCellsTot = 0;
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
	logPrintf("Corrected force-matrix hermiticity relative error: %lg\n", sqrt(hermErrNum/hermErrDen));
	//--- enforce translational invariance:
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
		//Project out each overall Cartesian displacement:
		int nAtoms = modes.size()/3;
		matrix proj(nAtoms, nAtoms); //set to ones(nAtoms)/nAtoms
		for(int iAtom=0; iAtom<nAtoms; iAtom++)
			for(int jAtom=0; jAtom<nAtoms; jAtom++)
				proj.set(iAtom, jAtom, 1./nAtoms);
		matrix dF0 = zeroes(modes.size(), modes.size());
		for(int iDir=0; iDir<3; iDir++)
			for(int jDir=0; jDir<3; jDir++)
			{	matrix F0sub = F0(iDir,3,modes.size(), jDir,3,modes.size()); //Gamma forces for current direction pair
				matrix dF0sub = proj*F0sub*proj - proj*F0sub - F0sub*proj; //correct net force, while preserving hermiticity
				dF0.set(iDir,3,modes.size(), jDir,3,modes.size(), dF0sub);
			}
		matrix domegaSqMean = (1./prodSup) * (invsqrtMmode * dF0 * invsqrtMmode); //correction to mean omegaSq
		//Apply correction:
		iter = cellMap.begin();
		for(size_t iCell=0; iCell<cellMap.size(); iCell++)
		{	omegaSq[iCell] += iter->second * domegaSqMean;
			iter++;
		}
		double nrm2omegaSqMean = hermErrDen / sqrt(cellMap.size());
		logPrintf("Corrected force-matrix translational invariance relative error: %lg\n", nrm2(domegaSqMean)/nrm2omegaSqMean);
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
	std::vector< std::pair<vector3<>,double> > getQuadratureBZ(void); //implemented below
	std::vector< std::pair<vector3<>,double> > quad = getQuadratureBZ();
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
